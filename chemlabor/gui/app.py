import threading
import queue
import math
import random
import json
import webbrowser

from collections import deque

import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog

from chemlabor.chemie import rdkit_bruecke as rb
from chemlabor.chemie import neuheit
from chemlabor.chemie import zufallsgenerator as idk
from chemlabor import labor

from chemlabor.labor import ChemLabor, formel
from chemlabor.gui.formatierung import formatiere_probleme
from chemlabor.ki.gedaechtnis import KIGedaechtnis, molekuel_signatur
from chemlabor.pfade import AUTOSAVE_PFAD


# ==========================================================
# LaborGUI
# ==========================================================
#
# GUI fuer das kuenstliche Chemie-Labor.
#
# Optimierungen:
#
# - Labor laeuft in Hintergrundthread
# - Tkinter wird nur im Hauptthread veraendert
# - Queue Kommunikation
# - Ergebnis Ringbuffer verhindert RAM Explosion
#
# ==========================================================


class LaborGUI:


    STATISTIK_FELDER = [
        "Experimente gesamt",
        "Versuche gesamt",
        "Durchschnittspunkte",
        "Minimum",
        "Maximum",
        "Anteil stabil (≥75)",
        "Anzahl Beobachtung",
        "Anzahl Automatik",
        "Verhältnis-Filter"
    ]


    ERGEBNIS_SPALTEN = (
        "nr",
        "phase",
        "formel",
        "atome",
        "bindungen",
        "punkte",
        "urteil",
        "rdkit"
    )


    SPALTEN_BESCHRIFTUNG = {

        "nr": "Nr",
        "phase": "Phase",
        "formel": "Formel",
        "atome": "Atome",
        "bindungen": "Bindungen",
        "punkte": "Punkte",
        "urteil": "Urteil",
        "rdkit": "RDKit"

    }


    SPALTEN_BREITE = {

        "nr": 50,
        "phase": 90,
        "formel": 100,
        "atome": 60,
        "bindungen": 70,
        "punkte": 70,
        "urteil": 150,
        "rdkit": 60

    }



    SPALTEN_SCHLUESSEL = {

        "nr":
            lambda e: e[0],

        "phase":
            lambda e: e[1],

        "formel":
            lambda e: formel(e[2]),

        "atome":
            lambda e: len(e[2]["atoms"]),

        "bindungen":
            lambda e: len(e[2]["bonds"]),

        "punkte":
            lambda e: e[3]["bewertung"]["punkte"],

        "urteil":
            lambda e: e[3]["bewertung"]["urteil"],

        # Unabhaengiges Validitaetssignal aus rdkit_bruecke.py (siehe
        # labor.ChemLabor.experiment()) - ob RDKits eigene, auf echter
        # Chemie basierende Sanitize-Pruefung diese Struktur akzeptieren
        # wuerde, unabhaengig vom Urteil unseres eigenen Kritikers.
        "rdkit":
            lambda e: "✓" if e[3].get("rdkit_gueltig") else "✗"

    }


    # Farbe je Element fuer die Lewis-Ansicht - rein kosmetisch, macht
    # groessere Strukturen aber deutlich leichter lesbar als einheitliche
    # weisse Kreise.
    ELEMENT_FARBEN = {
        "H": "#e8e8e8", "C": "#4a4a4a", "N": "#3060c8", "O": "#c83030",
        "F": "#3fae3f", "Cl": "#2f9e5f", "P": "#e08820", "S": "#d8c020",
        "Si": "#a8763f"
    }


    # Dunkles Farbschema fuer die gesamte Oberflaeche (siehe
    # _konfiguriere_style) - an einer Stelle gesammelt, damit ein
    # spaeteres helles Theme nur diese Werte austauschen muesste.
    FARBEN = {
        "bg": "#1e1e1e",
        "panel": "#252526",
        "panel_hell": "#2d2d30",
        "rand": "#3f3f46",
        "text": "#e0e0e0",
        "text_gedimmt": "#9a9a9a",
        "akzent": "#4fc3f7",
        "akzent_dunkel": "#2f9ee0",
        "gefahr": "#7a3030",
        "gefahr_hell": "#a23f3f",
        "auswahl": "#3a5f8a",
    }


    def __init__(self, root):


        self.root = root

        self.root.title(
            "Künstliches Chemie-Labor"
        )

        self.root.geometry(
            "1050x680"
        )


        # --------------------------------------------------
        # Thread sichere Queues
        # --------------------------------------------------

        self.ergebnis_queue = queue.Queue()

        self.fortschritt_queue = queue.Queue()

        # Fuer die asynchrone Lewis-Layout-Berechnung (siehe
        # _zeichne_lewis / _poll_lewis_queue) - getrennt von den
        # Labor-Queues, da unabhaengig vom Hintergrund-Lauf des Labors
        # benutzt wird (Klick auf eine Ergebniszeile).
        self.lewis_queue = queue.Queue()

        self._lewis_auftrag_id = 0



        # --------------------------------------------------
        # Labor Status
        # --------------------------------------------------

        self.labor = None

        self.laeuft = False



        # --------------------------------------------------
        # WICHTIG:
        #
        # Nur begrenzte Ergebnisse speichern.
        #
        # Millionen Experimente erzeugen keine
        # Millionen GUI Objekte mehr.
        #
        # --------------------------------------------------

        self.MAX_GUI_ERGEBNISSE = 5000


        self.alle_ergebnisse_lokal = deque(
            maxlen=self.MAX_GUI_ERGEBNISSE
        )

        # Wie viele Zeilen die Treeview beim Filter "Alle" hoechstens
        # gleichzeitig rendert (siehe _gefilterte_ergebnisse) - trennt
        # "wie viel wird gespeichert" (MAX_GUI_ERGEBNISSE) von "wie
        # viel wird auf einmal ins Treeview-Widget eingefuegt".
        self.MAX_ANGEZEIGTE_ERGEBNISSE = 500



        # "Zuletzt gerenderter Datenstand" je Tab - siehe
        # _render_wenn_noetig(). -1 sorgt dafuer, dass der allererste
        # Tick mit vorhandenen Daten garantiert einmal rendert (die
        # tatsaechliche Ergebnisanzahl ist nie negativ).
        self._statistik_render_stand = -1
        self._lernkurve_render_stand = -1
        self._tabelle_render_stand = -1


        self._sortier_spalte = None

        self._sortier_absteigend = False


        # Cache fuer molekuel_signatur() (siehe _gefilterte_ergebnisse):
        # ohne das wuerde bei aktiver "Wiederholungen ausblenden"-
        # Checkbox JEDES Tabellen-Update die Signatur fuer ALLE
        # gepufferten Ergebnisse (bis zu 5000, siehe
        # labor.ChemLabor.alle_ergebnisse) neu berechnen - eine der
        # Hauptursachen fuer spuerbares Laggen. Schluessel: die
        # Experiment-Nummer (eindeutig, aendert sich nie fuer ein
        # bereits gemeldetes Ergebnis).
        self._signatur_cache = {}

        # [patch_fixes] Monoton steigender Ergebniszaehler: len() der
        # deque bleibt ab MAX_GUI_ERGEBNISSE konstant und liess die
        # Anzeige danach einfrieren (siehe _render_wenn_noetig).
        self._ergebnis_zaehler = 0
        self._angezeigte_nummer = None



        self.aktuelles_molekuel = None

        # Von der Platte geladenes KIGedaechtnis, das beim naechsten
        # Start anstelle eines frischen verwendet wird (siehe
        # _gedaechtnis_laden_klick und _start_klick).
        self.geladenes_gedaechtnis = None

        # --------------------------------------------------
        # Neuheits-Filter (Tab "Entdeckungen", siehe neuheit.py)
        # --------------------------------------------------
        # Die Pruefung laeuft in einem Hintergrundthread (PubChem-
        # Anfragen dauern: ~0,25 s pro Kandidat), Rueckmeldungen
        # kommen ueber neuheit_queue in den Hauptthread.
        self.neuheit_queue = queue.Queue()
        self.neuheit_ergebnisse = []
        self._neuheit_laeuft = False
        self._neuheit_abbruch = threading.Event()

        self._konfiguriere_style()

        self._baue_oberflaeche()


        self._poll_queues()

        self._poll_lewis_queue()

        self._poll_neuheit_queue()





    # ======================================================
    # Oberfläche
    # ======================================================


    # ======================================================
    # Style / Farbschema
    # ======================================================
    #
    # Ein dunkles Theme fuer die gesamte Oberflaeche. ttk-Widgets
    # werden zentral ueber ttk.Style() eingefaerbt (Theme-Basis "clam",
    # da es - anders als die systemeigenen Themes - Hintergrund-/
    # Vordergrundfarben auf allen Widgets zuverlaessig annimmt). Reine
    # tk-Widgets (Canvas, ScrolledText) kennen kein ttk.Style und
    # bekommen ihre Farben stattdessen direkt beim Erzeugen ueber
    # self.FARBEN.
    # ======================================================

    def _konfiguriere_style(self):

        f = self.FARBEN

        self.root.configure(bg=f["bg"])

        style = ttk.Style(self.root)
        style.theme_use("clam")

        style.configure(
            ".",
            background=f["bg"],
            foreground=f["text"],
            fieldbackground=f["panel_hell"],
            bordercolor=f["rand"],
            font=("Segoe UI", 10)
        )

        style.configure("TFrame", background=f["bg"])

        style.configure("TLabel", background=f["bg"], foreground=f["text"])

        style.configure(
            "TLabelframe",
            background=f["bg"],
            foreground=f["text"],
            bordercolor=f["rand"]
        )

        style.configure(
            "TLabelframe.Label",
            background=f["bg"],
            foreground=f["akzent"],
            font=("Segoe UI", 10, "bold")
        )

        style.configure(
            "TButton",
            background=f["panel_hell"],
            foreground=f["text"],
            bordercolor=f["rand"],
            padding=6
        )

        style.map(
            "TButton",
            background=[("active", f["auswahl"]), ("disabled", f["panel"])],
            foreground=[("disabled", f["text_gedimmt"])]
        )

        # Hervorgehobene Buttons fuer die wichtigsten Aktionen (Start)
        # bzw. eine potenziell folgenreiche Aktion (Stoppen) - macht auf
        # einen Blick klar, welcher Button "los geht's" und welcher
        # "haelt an" bedeutet, statt dass beide gleich aussehen.
        style.configure(
            "Accent.TButton",
            background=f["akzent_dunkel"],
            foreground="#ffffff",
            bordercolor=f["akzent_dunkel"],
            padding=6,
            font=("Segoe UI", 10, "bold")
        )

        style.map(
            "Accent.TButton",
            background=[("active", f["akzent"]), ("disabled", f["panel"])]
        )

        style.configure(
            "Danger.TButton",
            background=f["gefahr"],
            foreground="#ffffff",
            bordercolor=f["gefahr"],
            padding=6
        )

        style.map(
            "Danger.TButton",
            background=[("active", f["gefahr_hell"]), ("disabled", f["panel"])]
        )

        style.configure(
            "TEntry",
            fieldbackground=f["panel_hell"],
            foreground=f["text"],
            insertcolor=f["text"],
            bordercolor=f["rand"]
        )

        style.configure("TCheckbutton", background=f["bg"], foreground=f["text"])

        style.map(
            "TCheckbutton",
            background=[("active", f["bg"])],
            foreground=[("disabled", f["text_gedimmt"])]
        )

        style.configure(
            "TNotebook",
            background=f["bg"],
            bordercolor=f["rand"]
        )

        style.configure(
            "TNotebook.Tab",
            background=f["panel"],
            foreground=f["text_gedimmt"],
            padding=(14, 7)
        )

        style.map(
            "TNotebook.Tab",
            background=[("selected", f["panel_hell"])],
            foreground=[("selected", f["akzent"])]
        )

        style.configure(
            "Treeview",
            background=f["panel_hell"],
            fieldbackground=f["panel_hell"],
            foreground=f["text"],
            bordercolor=f["rand"],
            borderwidth=0,
            rowheight=24
        )

        style.configure(
            "Treeview.Heading",
            background=f["panel"],
            foreground=f["akzent"],
            bordercolor=f["rand"],
            font=("Segoe UI", 9, "bold")
        )

        style.map(
            "Treeview",
            background=[("selected", f["auswahl"])],
            foreground=[("selected", "#ffffff")]
        )

        for orientierung in ("Vertical", "Horizontal"):

            style.configure(
                f"{orientierung}.TScrollbar",
                background=f["panel"],
                troughcolor=f["bg"],
                bordercolor=f["rand"],
                arrowcolor=f["text"]
            )


    def _baue_oberflaeche(self):


        f = self.FARBEN

        kopf = ttk.Frame(
            self.root,
            padding=(10, 10, 10, 4)
        )

        kopf.pack(
            side="top",
            fill="x"
        )


        # --------------------------------------------------
        # Gruppe: Generator-Einstellungen
        # --------------------------------------------------

        generator_frame = ttk.LabelFrame(
            kopf,
            text="Generator-Einstellungen",
            padding=(12, 8)
        )

        generator_frame.pack(
            side="left",
            fill="y",
            padx=(0, 8)
        )


        felder = [
            ("Beobachtungen:", "beobachtungen_var", str(labor.BEOBACHTUNGEN), 8),
            ("Automatik:", "automatik_var", str(labor.AUTOMATIK_EXPERIMENTE), 8),
            ("Kandidaten/Molekül:", "kandidaten_var", str(labor.KANDIDATEN_PRO_MOLEKUEL), 6),
            ("Max. Atome/Molekül:", "max_atome_var", str(idk.MOLEKUEL_MAX_ATOME), 6),
            ("Verhältnis Neg:Pos:", "verhaeltnis_var", str(labor.ZIEL_VERHAELTNIS_NEG_ZU_POS), 6),
        ]

        for spalte, (beschriftung, var_name, start_wert, breite) in enumerate(felder):

            ttk.Label(
                generator_frame,
                text=beschriftung
            ).grid(
                row=0,
                column=spalte,
                sticky="w",
                padx=(0 if spalte == 0 else 10, 4),
                pady=(0, 2)
            )

            var = tk.StringVar(value=start_wert)

            setattr(self, var_name, var)

            entry = ttk.Entry(
                generator_frame,
                textvariable=var,
                width=breite
            )

            entry.grid(
                row=1,
                column=spalte,
                sticky="w",
                padx=(0 if spalte == 0 else 10, 4)
            )

            if var_name == "automatik_var":
                self.automatik_entry = entry


        # Endlos-Modus: Automatikphase laeuft, bis "Stoppen" gedrueckt
        # wird, statt nach der obigen Anzahl von selbst zu enden.
        # Deaktiviert das Automatik-Zahlenfeld waehrenddessen, damit
        # nicht der Eindruck entsteht, die Zahl wuerde noch etwas
        # bewirken.

        self.endlos_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(
            generator_frame,
            text="Endlos-Automatik",
            variable=self.endlos_var,
            command=self._endlos_umschalten
        ).grid(
            row=1,
            column=len(felder),
            sticky="w",
            padx=(14, 0)
        )


        # Automatik-Modus: "Zufall" ist der bisherige Ansatz
        # (KIGenerator wuerfelt komplette Kandidaten und waehlt per
        # Best-of-N, siehe ki_gedaechtnis.KIGenerator). "Konstruktion"
        # baut jedes Molekuel schrittweise ueber KIKonstruktor
        # (konstruktor.py) auf und prueft jeden Bindungsschritt sofort
        # gegen RDKit (siehe ChemLabor.neues_molekuel_konstruktion()).

        ttk.Label(
            generator_frame,
            text="Automatik-Modus:"
        ).grid(
            row=0,
            column=len(felder) + 1,
            sticky="w",
            padx=(14, 4)
        )

        self.automatik_modus_var = tk.StringVar(value="Zufall")

        ttk.Combobox(
            generator_frame,
            textvariable=self.automatik_modus_var,
            values=["Zufall", "Konstruktion (RDKit)"],
            state="readonly",
            width=18
        ).grid(
            row=1,
            column=len(felder) + 1,
            sticky="w",
            padx=(14, 0)
        )


        # --------------------------------------------------
        # Gruppe: Lauf & Speicherung
        # --------------------------------------------------

        lauf_frame = ttk.LabelFrame(
            kopf,
            text="Lauf & Speicherung",
            padding=(12, 8)
        )

        lauf_frame.pack(
            side="left",
            fill="y"
        )


        self.start_button = ttk.Button(
            lauf_frame,
            text="▶ Start",
            style="Accent.TButton",
            command=self._start_klick
        )

        self.start_button.grid(
            row=0,
            column=0,
            sticky="w",
            padx=(0, 6),
            pady=(0, 6)
        )


        self.stop_button = ttk.Button(
            lauf_frame,
            text="■ Stoppen",
            style="Danger.TButton",
            command=self._stop_klick,
            state="disabled"
        )

        self.stop_button.grid(
            row=0,
            column=1,
            sticky="w",
            padx=(0, 6),
            pady=(0, 6)
        )


        # Autosave: sichert das Gedaechtnis automatisch alle
        # AUTOSAVE_INTERVALL Experimente in eine feste Datei - v.a. fuer
        # den Endlos-Modus gedacht, wo ohne das bei einem Absturz die
        # komplette Trainingszeit verloren waere.

        self.autosave_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(
            lauf_frame,
            text="Automatisch sichern",
            variable=self.autosave_var
        ).grid(
            row=0,
            column=2,
            sticky="w",
            padx=(6, 0),
            pady=(0, 6)
        )


        # Persistenz: gelerntes Wissen (KIGedaechtnis) auf Platte
        # sichern/wiederherstellen - ohne das ist bei jedem Neustart
        # der GUI die komplette Trainingszeit verloren.

        ttk.Button(
            lauf_frame,
            text="Gedächtnis speichern",
            command=self._gedaechtnis_speichern_klick
        ).grid(
            row=1,
            column=0,
            columnspan=2,
            sticky="we",
            padx=(0, 6)
        )

        ttk.Button(
            lauf_frame,
            text="Gedächtnis laden",
            command=self._gedaechtnis_laden_klick
        ).grid(
            row=1,
            column=2,
            sticky="we"
        )


        ttk.Button(
            lauf_frame,
            text="Gedächtnisse zusammenführen…",
            command=self._gedaechtnis_zusammenfuehren_klick
        ).grid(
            row=2,
            column=0,
            columnspan=3,
            sticky="we",
            pady=(6, 0)
        )


        # --------------------------------------------------
        # Statuszeile
        # --------------------------------------------------

        status_frame = ttk.Frame(
            self.root,
            padding=(10, 0, 10, 8)
        )

        status_frame.pack(
            side="top",
            fill="x"
        )

        self.status_label = ttk.Label(
            status_frame,
            text="Bereit.",
            foreground=f["akzent"]
        )

        self.status_label.pack(
            side="left"
        )


        self.notebook = ttk.Notebook(
            self.root
        )


        self.notebook.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=(0, 10)
        )

        self._baue_log_tab()

        self._baue_statistik_tab()

        self._baue_ergebnisse_tab()

        self._baue_lewis_tab()

        self._baue_lernkurve_tab()

        self._baue_entdeckungen_tab()


        # Kein Tab-Wechsel-Bind mehr noetig: _render_wenn_noetig() (Teil
        # von _poll_queues, laeuft alle 150ms fuer immer, unabhaengig
        # vom Lauf-Status) prueft bei jedem Tick ohnehin, welcher Tab
        # gerade sichtbar ist - ein Wechsel wird so binnen 150ms von
        # selbst erkannt, ganz ohne Event-Sonderfall.

    # ======================================================
    # Lewis Tab
    # ======================================================

    def _baue_lewis_tab(self):

        frame = ttk.Frame(
            self.notebook
        )

        self.notebook.add(
            frame,
            text="Lewis"
        )

        lewis_leiste = ttk.Frame(
            frame
        )

        lewis_leiste.pack(
            fill="x",
            padx=10,
            pady=5
        )

        self.lewis_status = ttk.Label(
            lewis_leiste,
            text="Noch kein Molekül ausgewählt."
        )

        self.lewis_status.pack(
            side="left"
        )

        ttk.Button(
            lewis_leiste,
            text="Bindungen als JSON anzeigen",
            command=self._zeige_json
        ).pack(
            side="right"
        )

        self.lewis_canvas = tk.Canvas(
            frame,
            width=800,
            height=500,
            bg=self.FARBEN["panel_hell"],
            highlightthickness=0
        )

        self.lewis_canvas.pack(
            fill="both",
            expand=True
        )
    # ======================================================
    # Log Tab
    # ======================================================


    def _baue_log_tab(self):


        frame = ttk.Frame(
            self.notebook
        )


        self.notebook.add(
            frame,
            text="Log"
        )


        self.log_text = scrolledtext.ScrolledText(
            frame,
            wrap="word",
            state="disabled",
            bg=self.FARBEN["panel_hell"],
            fg=self.FARBEN["text"],
            insertbackground=self.FARBEN["text"],
            borderwidth=0,
            highlightthickness=0
        )


        self.log_text.pack(
            fill="both",
            expand=True
        )





    # ======================================================
    # Statistik Tab Platzhalter
    # (kommt Teil 2)
    # ======================================================


    def _baue_statistik_tab(self):

        frame = ttk.Frame(
            self.notebook
        )

        self.notebook.add(
            frame,
            text="Statistik"
        )

        # Referenz fuer die Sichtbarkeitspruefung in
        # _aktualisiere_statistik() - Canvas-Neuzeichnen nur, wenn
        # dieser Tab gerade tatsaechlich zu sehen ist.
        self._statistik_tab_frame = frame


        self.stat_labels = {}


        for feld in self.STATISTIK_FELDER:

            label = ttk.Label(
                frame,
                text=feld + ": –"
            )


            label.pack(
                anchor="w"
            )


            self.stat_labels[feld] = label



        self.histogramm_canvas = tk.Canvas(
            frame,
            width=520,
            height=200,
            bg=self.FARBEN["panel_hell"],
            highlightthickness=0
        )


        self.histogramm_canvas.pack(
            pady=10
        )


    # ======================================================
    # Lernkurve Tab
    # ======================================================
    #
    # Zeigt KIGedaechtnis.lernkurve_daten() als einfaches Liniendiagramm:
    # Durchschnittspunktzahl je Intervall (siehe
    # KIGedaechtnis._lernkurve_abschliessen) ueber die Zahl der
    # Beobachtungen aufgetragen - macht sichtbar, ob (und wie schnell)
    # sich die erzeugten Molekuele im Lauf der Zeit verbessern, statt
    # nur den aktuellen Gesamtdurchschnitt zu zeigen.
    # ======================================================


    def _baue_lernkurve_tab(self):

        frame = ttk.Frame(
            self.notebook
        )

        self.notebook.add(
            frame,
            text="Lernkurve"
        )

        # Eigene Sichtbarkeitspruefung, da die Lernkurve ein SEPARATER
        # Tab ist (nicht Teil von "Statistik") - siehe
        # _aktualisiere_statistik().
        self._lernkurve_tab_frame = frame

        self.lernkurve_status = ttk.Label(
            frame,
            text="Noch keine Daten - Intervall wird erst nach den "
                 "ersten Beobachtungen gefüllt."
        )

        self.lernkurve_status.pack(
            anchor="w",
            padx=10,
            pady=(10, 0)
        )

        self.lernkurve_canvas = tk.Canvas(
            frame,
            width=800,
            height=400,
            bg=self.FARBEN["panel_hell"],
            highlightthickness=0
        )

        self.lernkurve_canvas.pack(
            fill="both",
            expand=True,
            padx=10,
            pady=10
        )


    def _zeichne_lernkurve(self):

        if self.labor is None:
            return

        punkte = self.labor.gedaechtnis.lernkurve_daten()

        canvas = self.lernkurve_canvas

        canvas.delete("all")

        if len(punkte) < 2:

            self.lernkurve_status.config(
                text=f"Noch zu wenig Daten für eine Kurve "
                     f"({len(punkte)} Intervall(e) abgeschlossen)."
            )

            return

        _differenz, trend_text = self.labor.gedaechtnis.lerntrend()

        self.lernkurve_status.config(
            text=f"{len(punkte)} Intervalle (je "
                 f"{self.labor.gedaechtnis.lernkurve_intervall} "
                 f"Beobachtungen) - {trend_text}"
        )

        breite = int(canvas["width"])
        hoehe = int(canvas["height"])

        rand_links = 45
        rand_unten = 25
        rand_oben = 15
        rand_rechts = 15

        werte = [durchschnitt for _stand, durchschnitt in punkte]

        y_min = min(0, min(werte))
        y_max = max(100, max(werte))

        # Sicherheitsabstand, falls alle Werte identisch sind (sonst
        # Division durch 0 bei der Skalierung unten).
        if y_max == y_min:
            y_max = y_min + 1

        x_min = punkte[0][0]
        x_max = punkte[-1][0]

        if x_max == x_min:
            x_max = x_min + 1

        plot_breite = breite - rand_links - rand_rechts
        plot_hoehe = hoehe - rand_oben - rand_unten

        def x_pos(stand):
            return rand_links + (stand - x_min) / (x_max - x_min) * plot_breite

        def y_pos(wert):
            return (
                rand_oben
                + (1 - (wert - y_min) / (y_max - y_min)) * plot_hoehe
            )

        # Achsen
        canvas.create_line(
            rand_links, rand_oben,
            rand_links, hoehe - rand_unten,
            fill=self.FARBEN["text_gedimmt"]
        )

        canvas.create_line(
            rand_links, hoehe - rand_unten,
            breite - rand_rechts, hoehe - rand_unten,
            fill=self.FARBEN["text_gedimmt"]
        )

        # Ein paar horizontale Hilfslinien mit Beschriftung.
        for anteil in (0, 0.25, 0.5, 0.75, 1.0):

            wert = y_min + anteil * (y_max - y_min)
            y = y_pos(wert)

            canvas.create_line(
                rand_links, y,
                breite - rand_rechts, y,
                fill=self.FARBEN["rand"]
            )

            canvas.create_text(
                rand_links - 8, y,
                text=f"{wert:.0f}",
                anchor="e",
                font=("", 8),
                fill=self.FARBEN["text_gedimmt"]
            )

        # Kurve selbst.
        koordinaten = []

        for stand, durchschnitt in punkte:
            koordinaten.append(x_pos(stand))
            koordinaten.append(y_pos(durchschnitt))

        canvas.create_line(
            *koordinaten,
            fill=self.FARBEN["akzent"],
            width=2
        )

        for stand, durchschnitt in punkte:

            x = x_pos(stand)
            y = y_pos(durchschnitt)

            canvas.create_oval(
                x - 2, y - 2,
                x + 2, y + 2,
                fill=self.FARBEN["akzent"],
                outline=""
            )

        canvas.create_text(
            breite - rand_rechts, hoehe - rand_unten + 12,
            text=f"Beobachtungen: {x_min} – {x_max}",
            anchor="e",
            font=("", 8),
            fill=self.FARBEN["text_gedimmt"]
        )


    # ======================================================
    # Ergebnisse Tab
    # ======================================================


    def _baue_ergebnisse_tab(self):


        frame = ttk.Frame(
            self.notebook
        )


        self.notebook.add(
            frame,
            text="Ergebnisse"
        )

        # Referenz fuer den zentralen Render-Tick (siehe
        # _render_wenn_noetig()).
        self._ergebnisse_tab_frame = frame



        toolbar = ttk.Frame(
            frame,
            padding=6
        )


        toolbar.pack(
            side="top",
            fill="x"
        )



        ttk.Label(
            toolbar,
            text="Filter:"
        ).pack(
            side="left"
        )



        self.filter_var = tk.StringVar(
            value="Alle"
        )


        filter_box = ttk.Combobox(
            toolbar,
            textvariable=self.filter_var,
            state="readonly",
            width=20,
            values=[
                "Alle",
                "Top 10",
                "Top 50",
                "Nur stabil (≥75)",
                "Nur instabil (<50)"
            ]
        )


        filter_box.pack(
            side="left",
            padx=(4,16)
        )


        filter_box.bind(
            "<<ComboboxSelected>>",
            lambda e:
            self._aktualisiere_ergebnisse_tabelle()
        )





        ttk.Label(
            toolbar,
            text="Phase:"
        ).pack(
            side="left"
        )



        self.phase_filter_var = tk.StringVar(
            value="Alle"
        )


        phase_box = ttk.Combobox(
            toolbar,
            textvariable=self.phase_filter_var,
            state="readonly",
            width=14,
            values=[
                "Alle",
                "Beobachtung",
                "Automatik"
            ]
        )


        phase_box.pack(
            side="left",
            padx=(4,16)
        )


        phase_box.bind(
            "<<ComboboxSelected>>",
            lambda e:
            self._aktualisiere_ergebnisse_tabelle()
        )




        # Blendet strukturell wiederholte Molekuele aus (gleiche
        # molekuel_signatur(), siehe ki_gedaechtnis.py) - zeigt je
        # Struktur nur den zuerst aufgetretenen Eintrag. Betrifft nur
        # die Anzeige hier; Statistik/Gedaechtnis lernen weiterhin aus
        # jedem Vorkommen.
        self.duplikate_ausblenden_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(
            toolbar,
            text="Wiederholungen ausblenden",
            variable=self.duplikate_ausblenden_var,
            command=self._aktualisiere_ergebnisse_tabelle
        ).pack(
            side="left",
            padx=(0, 16)
        )




        ttk.Button(
            toolbar,
            text="Aktualisieren",
            command=self._aktualisiere_ergebnisse_tabelle
        ).pack(
            side="left"
        )


        self.ergebnisse_info_label = ttk.Label(
            toolbar,
            text=""
        )

        self.ergebnisse_info_label.pack(
            side="left",
            padx=(12,0)
        )





        inhalt = ttk.Frame(
            frame
        )


        inhalt.pack(
            fill="both",
            expand=True
        )




        self.tabelle = ttk.Treeview(
            inhalt,
            columns=self.ERGEBNIS_SPALTEN,
            show="headings",
            selectmode="browse"
        )



        for spalte in self.ERGEBNIS_SPALTEN:


            self.tabelle.heading(
                spalte,
                text=self.SPALTEN_BESCHRIFTUNG[spalte],
                command=lambda s=spalte:
                    self._sortiere_tabelle(s)
            )



            self.tabelle.column(
                spalte,
                width=self.SPALTEN_BREITE[spalte],
                anchor="center"
            )




        scrollbar = ttk.Scrollbar(
            inhalt,
            orient="vertical",
            command=self.tabelle.yview
        )


        self.tabelle.configure(
            yscrollcommand=scrollbar.set
        )



        self.tabelle.pack(
            side="left",
            fill="both",
            expand=True
        )


        scrollbar.pack(
            side="left",
            fill="y"
        )





        self.detail_text = scrolledtext.ScrolledText(
            inhalt,
            wrap="word",
            width=45,
            state="disabled",
            bg=self.FARBEN["panel_hell"],
            fg=self.FARBEN["text"],
            insertbackground=self.FARBEN["text"],
            borderwidth=0,
            highlightthickness=0
        )


        self.detail_text.pack(
            side="left",
            fill="both"
        )



        self.tabelle.bind(
            "<<TreeviewSelect>>",
            self._zeige_detail
        )





    # ======================================================
    # Entdeckungen Tab (Neuheits-Filter)
    # ======================================================
    #
    # Prueft die gesammelten Ergebnisse gegen PubChem (siehe
    # neuheit.py): "unbekannt" heisst, das Molekuel ist dort NICHT
    # verzeichnet - ein Kandidat fuer etwas Neues. Hinweis zur Spalte
    # "Stereo": steht dort "offen", hat das Molekuel unfestgelegte
    # Stereozentren; ein "unbekannt" kann dann daran liegen, dass
    # PubChem nur konkrete Stereoisomere fuehrt - von Hand gegenpruefen.
    # ======================================================


    NEUHEIT_SPALTEN = (
        "status", "punkte", "formel", "molmasse", "smiles", "cid", "stereo"
    )

    NEUHEIT_BESCHRIFTUNG = {
        "status": "Status",
        "punkte": "Punkte",
        "formel": "Formel",
        "molmasse": "Molmasse",
        "smiles": "SMILES",
        "cid": "PubChem-CID",
        "stereo": "Stereo",
    }

    NEUHEIT_BREITE = {
        "status": 90,
        "punkte": 60,
        "formel": 100,
        "molmasse": 80,
        "smiles": 260,
        "cid": 90,
        "stereo": 60,
    }

    NEUHEIT_STATUS_TEXT = {
        "unbekannt": "unbekannt",
        "bekannt": "bekannt",
        "ungeprueft": "ungeprüft",
    }

    # Anzeige-Filter -> erlaubte Status (None = alle)
    NEUHEIT_FILTER = {
        "Alle": None,
        "Nur unbekannt": ("unbekannt",),
        "Nur bekannt": ("bekannt",),
        "Nur ungeprüft": ("ungeprueft",),
    }


    def _baue_entdeckungen_tab(self):

        frame = ttk.Frame(
            self.notebook
        )

        self.notebook.add(
            frame,
            text="Entdeckungen"
        )

        leiste = ttk.Frame(
            frame,
            padding=6
        )

        leiste.pack(
            side="top",
            fill="x"
        )

        self.neuheit_button = ttk.Button(
            leiste,
            text="Neuheit prüfen",
            style="Accent.TButton",
            command=self._neuheit_klick
        )

        self.neuheit_button.pack(
            side="left",
            padx=(0, 12)
        )

        ttk.Label(
            leiste,
            text="Min. Punkte:"
        ).pack(
            side="left"
        )

        self.neuheit_min_var = tk.StringVar(
            value=str(labor.STABIL_ANZEIGE_SCHWELLE)
        )

        ttk.Entry(
            leiste,
            textvariable=self.neuheit_min_var,
            width=5
        ).pack(
            side="left",
            padx=(4, 12)
        )

        ttk.Label(
            leiste,
            text="Max. Kandidaten:"
        ).pack(
            side="left"
        )

        # Leer = alle eindeutigen Kandidaten pruefen. Der Standard
        # bleibt bewusst begrenzt: jede neue Struktur kostet eine
        # PubChem-Anfrage (~0,25 s), 5000 Kandidaten waeren ~20 Minuten.
        self.neuheit_max_var = tk.StringVar(
            value="200"
        )

        ttk.Entry(
            leiste,
            textvariable=self.neuheit_max_var,
            width=6
        ).pack(
            side="left",
            padx=(4, 12)
        )

        self.neuheit_offline_var = tk.BooleanVar(
            value=False
        )

        ttk.Checkbutton(
            leiste,
            text="Nur Cache (offline)",
            variable=self.neuheit_offline_var
        ).pack(
            side="left",
            padx=(0, 12)
        )

        ttk.Label(
            leiste,
            text="Anzeige:"
        ).pack(
            side="left"
        )

        self.neuheit_filter_var = tk.StringVar(
            value="Alle"
        )

        neuheit_filter_box = ttk.Combobox(
            leiste,
            textvariable=self.neuheit_filter_var,
            state="readonly",
            width=14,
            values=list(self.NEUHEIT_FILTER)
        )

        neuheit_filter_box.pack(
            side="left",
            padx=(4, 12)
        )

        neuheit_filter_box.bind(
            "<<ComboboxSelected>>",
            lambda e: self._aktualisiere_entdeckungen_tabelle()
        )

        ttk.Button(
            leiste,
            text="Exportieren…",
            command=self._neuheit_export_klick
        ).pack(
            side="left"
        )

        self.neuheit_status = ttk.Label(
            frame,
            text="Noch nicht geprüft. Doppelklick auf eine Zeile mit "
                 "PubChem-CID öffnet den Eintrag im Browser.",
            padding=(8, 0, 8, 4)
        )

        self.neuheit_status.pack(
            side="top",
            fill="x"
        )

        inhalt = ttk.Frame(
            frame
        )

        inhalt.pack(
            fill="both",
            expand=True
        )

        self.neuheit_tabelle = ttk.Treeview(
            inhalt,
            columns=self.NEUHEIT_SPALTEN,
            show="headings",
            selectmode="browse"
        )

        for spalte in self.NEUHEIT_SPALTEN:

            self.neuheit_tabelle.heading(
                spalte,
                text=self.NEUHEIT_BESCHRIFTUNG[spalte]
            )

            self.neuheit_tabelle.column(
                spalte,
                width=self.NEUHEIT_BREITE[spalte],
                anchor="w" if spalte == "smiles" else "center"
            )

        # Unbekannte Kandidaten optisch hervorheben.
        self.neuheit_tabelle.tag_configure(
            "unbekannt",
            foreground=self.FARBEN["akzent"]
        )

        self.neuheit_tabelle.tag_configure(
            "ungeprueft",
            foreground=self.FARBEN["text_gedimmt"]
        )

        scrollbar = ttk.Scrollbar(
            inhalt,
            orient="vertical",
            command=self.neuheit_tabelle.yview
        )

        self.neuheit_tabelle.configure(
            yscrollcommand=scrollbar.set
        )

        self.neuheit_tabelle.pack(
            side="left",
            fill="both",
            expand=True
        )

        scrollbar.pack(
            side="left",
            fill="y"
        )

        self.neuheit_tabelle.bind(
            "<<TreeviewSelect>>",
            self._neuheit_auswahl
        )

        self.neuheit_tabelle.bind(
            "<Double-1>",
            self._neuheit_doppelklick
        )


    def _neuheit_klick(self):
        """Startet die Pruefung - oder bricht sie ab, falls sie schon
        laeuft (der Button wechselt waehrenddessen auf 'Abbrechen')."""

        if self._neuheit_laeuft:

            self._neuheit_abbruch.set()

            self.neuheit_button.config(state="disabled")

            self.neuheit_status.config(text="Breche ab...")

            return

        if not self.alle_ergebnisse_lokal:

            self.neuheit_status.config(
                text="Keine Ergebnisse vorhanden - erst einen Lauf starten."
            )

            return

        try:

            min_punkte = float(
                self.neuheit_min_var.get().strip().replace(",", ".")
            )

            max_text = self.neuheit_max_var.get().strip()

            limit = max(1, int(max_text)) if max_text else None

        except ValueError:

            self.neuheit_status.config(text="Ungültige Zahl.")

            return

        # Schnappschuss im Hauptthread - das Labor laeuft evtl. noch
        # weiter und fuellt alle_ergebnisse_lokal.
        ergebnisse = list(self.alle_ergebnisse_lokal)

        offline = self.neuheit_offline_var.get()

        self._neuheit_abbruch = threading.Event()

        self._neuheit_laeuft = True

        self.neuheit_button.config(text="■ Abbrechen")

        self.neuheit_status.config(text="Starte Neuheitsfilter...")

        threading.Thread(
            target=self._neuheit_thread,
            args=(
                ergebnisse,
                min_punkte,
                limit,
                offline,
                self._neuheit_abbruch
            ),
            daemon=True
        ).start()


    def _neuheit_thread(
            self,
            ergebnisse,
            min_punkte,
            limit,
            offline,
            abbruch
    ):

        pruefer = None

        try:

            # Der Pruefer (und damit seine SQLite-Verbindung) wird
            # HIER im Thread erzeugt - SQLite-Verbindungen duerfen
            # standardmaessig nur in dem Thread benutzt werden, der
            # sie angelegt hat.
            pruefer = neuheit.NeuheitsPruefer(
                offline=offline
            )

            ergebnis = pruefer.pruefe_ergebnisse(
                ergebnisse,
                min_punkte=min_punkte,
                limit=limit,
                on_fortschritt=lambda text:
                    self.neuheit_queue.put(("fortschritt", text)),
                abbruch=abbruch
            )

            self.neuheit_queue.put(
                ("fertig", ergebnis, abbruch.is_set())
            )

        except Exception as fehler:

            self.neuheit_queue.put(
                ("fehler", f"{type(fehler).__name__}: {fehler}")
            )

        finally:

            if pruefer is not None:
                pruefer.schliessen()


    def _poll_neuheit_queue(self):
        """Laeuft alle 150ms im Hauptthread und verarbeitet die
        Rueckmeldungen des Neuheits-Threads."""

        try:

            while True:

                nachricht = self.neuheit_queue.get_nowait()

                art = nachricht[0]

                if art == "fortschritt":

                    self.neuheit_status.config(text=nachricht[1])

                    self._log_anhaengen(nachricht[1])

                elif art == "fertig":

                    self._neuheit_beenden()

                    self.neuheit_ergebnisse = nachricht[1]

                    self._aktualisiere_entdeckungen_tabelle()

                    self.neuheit_status.config(
                        text=self._neuheit_zusammenfassung(
                            abgebrochen=nachricht[2]
                        )
                    )

                elif art == "fehler":

                    self._neuheit_beenden()

                    self.neuheit_status.config(
                        text=f"Neuheitsfilter fehlgeschlagen: {nachricht[1]}"
                    )

        except queue.Empty:

            pass

        self.root.after(150, self._poll_neuheit_queue)


    def _neuheit_beenden(self):

        self._neuheit_laeuft = False

        self.neuheit_button.config(
            text="Neuheit prüfen",
            state="normal"
        )


    def _neuheit_zusammenfassung(self, abgebrochen=False):

        zaehler = {}

        stereo_offen = 0

        for e in self.neuheit_ergebnisse:

            zaehler[e["status"]] = zaehler.get(e["status"], 0) + 1

            if e["status"] == "unbekannt" and e["stereo_offen"]:
                stereo_offen += 1

        text = (
            f"{zaehler.get('unbekannt', 0)} unbekannt "
            f"({stereo_offen} davon mit offenem Stereo), "
            f"{zaehler.get('bekannt', 0)} bekannt, "
            f"{zaehler.get('ungeprueft', 0)} ungeprüft."
        )

        if abgebrochen:
            text += " (abgebrochen)"

        return text


    def _aktualisiere_entdeckungen_tabelle(self):

        erlaubt = self.NEUHEIT_FILTER[self.neuheit_filter_var.get()]

        kinder = self.neuheit_tabelle.get_children()

        if kinder:
            self.neuheit_tabelle.delete(*kinder)

        for index, e in enumerate(self.neuheit_ergebnisse):

            if erlaubt is not None and e["status"] not in erlaubt:
                continue

            self.neuheit_tabelle.insert(
                "",
                "end",
                # iid = Position in neuheit_ergebnisse, damit Auswahl
                # und Doppelklick den Eintrag wiederfinden.
                iid=str(index),
                values=(
                    self.NEUHEIT_STATUS_TEXT.get(e["status"], e["status"]),
                    e["punkte"],
                    e["formel"],
                    f"{e['molmasse']:.2f}",
                    e["smiles"],
                    e["cid"] if e["cid"] else "–",
                    "offen" if e["stereo_offen"] else "–",
                ),
                tags=(e["status"],)
            )


    def _neuheit_eintrag_zur_auswahl(self):

        auswahl = self.neuheit_tabelle.selection()

        if not auswahl:
            return None

        return self.neuheit_ergebnisse[int(auswahl[0])]


    def _neuheit_auswahl(self, event):
        """Klick auf eine Zeile: Molekuel fuer den Lewis-Tab vormerken
        und dort zeichnen (wie im Ergebnisse-Tab)."""

        e = self._neuheit_eintrag_zur_auswahl()

        if e is None:
            return

        self.aktuelles_molekuel = e["molekuel"]

        # Damit ein spaeterer Klick im Ergebnisse-Tab auf dieselbe
        # Nummer wieder neu zeichnet (siehe _zeige_detail).
        self._angezeigte_nummer = None

        self._zeichne_lewis(e["molekuel"])

        self.neuheit_status.config(
            text=f"{e['formel']}  {e['smiles']} - Struktur im Tab "
                 f"'Lewis'."
        )


    def _neuheit_doppelklick(self, event):

        e = self._neuheit_eintrag_zur_auswahl()

        if e is None or not e["cid"]:
            return

        webbrowser.open(
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{e['cid']}"
        )


    def _neuheit_export_klick(self):
        """Exportiert genau die Eintraege, die der Anzeige-Filter
        gerade zeigt, als JSON (inklusive internem Molekuel-Dict, damit
        sich Entdeckungen spaeter wieder laden/pruefen lassen)."""

        if not self.neuheit_ergebnisse:

            self.neuheit_status.config(
                text="Nichts zu exportieren - erst 'Neuheit prüfen'."
            )

            return

        erlaubt = self.NEUHEIT_FILTER[self.neuheit_filter_var.get()]

        pfad = filedialog.asksaveasfilename(
            title="Entdeckungen exportieren",
            defaultextension=".json",
            initialfile="entdeckungen.json",
            filetypes=[("JSON-Dateien", "*.json"), ("Alle Dateien", "*.*")]
        )

        if not pfad:
            return

        try:

            anzahl = neuheit.NeuheitsPruefer.exportiere(
                self.neuheit_ergebnisse,
                pfad,
                nur_status=erlaubt
            )

            self.neuheit_status.config(
                text=f"{anzahl} Einträge exportiert: {pfad}"
            )

        except OSError as fehler:

            self.neuheit_status.config(
                text=f"Export fehlgeschlagen: {fehler}"
            )


    # ======================================================
    # Endlos-Modus
    # ======================================================


    def _endlos_umschalten(self):

        if self.endlos_var.get():

            self.automatik_entry.config(state="disabled")

        else:

            self.automatik_entry.config(state="normal")


    def _stop_klick(self):

        if self.labor is not None:

            self.labor.stoppen()

        self.stop_button.config(state="disabled")

        self.status_label.config(text="Stoppe...")


    # ======================================================
    # Gedächtnis speichern/laden
    # ======================================================


    def _gedaechtnis_speichern_klick(self):

        if self.labor is None:

            self.status_label.config(
                text="Kein Gedächtnis vorhanden - erst einen Lauf starten."
            )

            return

        pfad = filedialog.asksaveasfilename(
            title="Gedächtnis speichern",
            defaultextension=".json",
            filetypes=[("JSON-Dateien", "*.json"), ("Alle Dateien", "*.*")]
        )

        if not pfad:
            return

        try:

            self.labor.speichere_gedaechtnis(pfad)

            self.status_label.config(
                text=f"Gedächtnis gespeichert: {pfad}"
            )

        except OSError as fehler:

            self.status_label.config(
                text=f"Speichern fehlgeschlagen: {fehler}"
            )


    def _gedaechtnis_laden_klick(self):

        if self.laeuft:

            self.status_label.config(
                text="Erst den laufenden Lauf stoppen, bevor ein "
                     "Gedächtnis geladen wird."
            )

            return

        pfad = filedialog.askopenfilename(
            title="Gedächtnis laden",
            filetypes=[("JSON-Dateien", "*.json"), ("Alle Dateien", "*.*")]
        )

        if not pfad:
            return

        try:

            self.geladenes_gedaechtnis = KIGedaechtnis.laden(pfad)

            self.status_label.config(
                text=f"Gedächtnis geladen ({self.geladenes_gedaechtnis.beobachtungen} "
                     f"Beobachtungen) - wird beim naechsten Start verwendet."
            )

        except (OSError, ValueError, KeyError) as fehler:

            self.status_label.config(
                text=f"Laden fehlgeschlagen: {fehler}"
            )


    def _gedaechtnis_zusammenfuehren_klick(self):
        """Fasst mehrere separat gespeicherte Gedaechtnisdateien zu
        einem gemeinsamen 'Kollektivgedaechtnis' zusammen (siehe
        KIGedaechtnis.zusammenfuehren()) und bietet an, das Ergebnis
        gleich als neue Datei zu sichern. Das Ergebnis wird zusaetzlich
        als naechstes zu ladendes Gedaechtnis vorgemerkt, genau wie
        bei "Gedächtnis laden"."""

        if self.laeuft:

            self.status_label.config(
                text="Erst den laufenden Lauf stoppen, bevor "
                     "Gedächtnisse zusammengeführt werden."
            )

            return

        pfade = filedialog.askopenfilenames(
            title="Gedächtnisse zum Zusammenführen auswählen (mind. 2)",
            filetypes=[("JSON-Dateien", "*.json"), ("Alle Dateien", "*.*")]
        )

        if not pfade or len(pfade) < 2:

            if pfade:
                self.status_label.config(
                    text="Bitte mindestens 2 Dateien auswählen."
                )

            return

        try:

            kollektiv = KIGedaechtnis.zusammenfuehren(list(pfade))

        except (OSError, ValueError, KeyError) as fehler:

            self.status_label.config(
                text=f"Zusammenführen fehlgeschlagen: {fehler}"
            )

            return

        self.geladenes_gedaechtnis = kollektiv

        speicherpfad = filedialog.asksaveasfilename(
            title="Zusammengeführtes Gedächtnis speichern (optional)",
            defaultextension=".json",
            filetypes=[("JSON-Dateien", "*.json"), ("Alle Dateien", "*.*")]
        )

        if speicherpfad:

            try:
                kollektiv.speichern(speicherpfad)
            except OSError as fehler:
                self.status_label.config(
                    text=f"Zusammengeführt, aber Speichern fehlgeschlagen: {fehler}"
                )
                return

        self.status_label.config(
            text=f"{len(pfade)} Dateien zusammengeführt "
                 f"({kollektiv.beobachtungen} Beobachtungen gesamt) - "
                 f"wird beim nächsten Start verwendet."
        )


    # ======================================================
    # Start
    # ======================================================


    def _start_klick(self):


        if self.laeuft:

            return



        try:


            labor.BEOBACHTUNGEN = max(
                1,
                int(
                    self.beobachtungen_var.get()
                )
            )


            labor.AUTOMATIK_EXPERIMENTE = max(
                1,
                int(
                    self.automatik_var.get()
                )
            )



            labor.KANDIDATEN_PRO_MOLEKUEL = max(
                1,
                int(
                    self.kandidaten_var.get()
                )
            )



            idk.MOLEKUEL_MAX_ATOME = max(
                2,
                int(
                    self.max_atome_var.get()
                )
            )



            # Leer oder "aus" (Gross-/Kleinschreibung egal) schaltet den
            # Verhaeltnis-Filter komplett ab (altes Verhalten, siehe
            # labor.ChemLabor._verhaeltnis_pruefen) - sonst eine
            # nicht-negative Zahl.
            verhaeltnis_text = self.verhaeltnis_var.get().strip()

            if verhaeltnis_text == "" or verhaeltnis_text.lower() == "aus":

                labor.ZIEL_VERHAELTNIS_NEG_ZU_POS = None

            else:

                labor.ZIEL_VERHAELTNIS_NEG_ZU_POS = max(
                    0,
                    int(verhaeltnis_text)
                )



            # FIX: vorher "100000 wenn BEOBACHTUNGEN > 1000 sonst 1" -
            # dieselbe kaputte Schwelle wie der labor.py-Modul-Default
            # (siehe dortiger Kommentar): traf bei jedem typischen Lauf
            # praktisch nie, das Log blieb waehrend der kompletten
            # Beobachtungsphase stumm. Jetzt proportional zur
            # tatsaechlich eingestellten Laufgroesse - ~20 gleichmaessig
            # verteilte Meldungen, egal wie gross der Lauf ist.
            labor.BEOBACHTUNG_LOG_INTERVALL = max(
                1,
                labor.BEOBACHTUNGEN // 20
            )



        except ValueError:


            self.status_label.config(
                text="Ungültige Zahl."
            )

            return





        self._log_leeren()

        self._tabelle_leeren()


        # neuer Lauf

        self.alle_ergebnisse_lokal.clear()
        self._ergebnis_zaehler = 0
        self._angezeigte_nummer = None


        self._statistik_render_stand = -1
        self._lernkurve_render_stand = -1
        self._tabelle_render_stand = -1




        self.laeuft = True

        self._endlos_lauf = self.endlos_var.get()


        self.start_button.config(
            state="disabled"
        )


        # FIX: stoppen() wirkt jetzt in beiden Phasen und unabhaengig
        # vom Endlos-Modus (siehe labor.py) - der Stop-Button ist daher
        # waehrend JEDES Laufs aktiv, nicht mehr nur im Endlos-Modus.
        self.stop_button.config(
            state="normal"
        )


        self.status_label.config(
            text="Läuft (Endlos)..." if self._endlos_lauf else "Läuft..."
        )





        self.labor = ChemLabor(
            on_ergebnis=self._on_ergebnis,
            on_fortschritt=self._on_fortschritt
        )

        self.labor.automatik_modus = (
            "konstruktion"
            if self.automatik_modus_var.get().startswith("Konstruktion")
            else "zufall"
        )


        # Ein zuvor geladenes Gedaechtnis (siehe _gedaechtnis_laden_klick)
        # ersetzt das frische, das ChemLabor.__init__ gerade angelegt
        # hat - MUSS vor beobachtungsphase()/automatikphase() passieren,
        # da ki_generator erst danach gebaut wird und dann automatisch
        # auf das (jetzt schon richtige) self.labor.gedaechtnis zeigt.
        if self.geladenes_gedaechtnis is not None:

            self.labor.gedaechtnis = self.geladenes_gedaechtnis

            self._on_fortschritt(
                f"Geladenes Gedächtnis übernommen "
                f"({self.geladenes_gedaechtnis.beobachtungen} frühere Beobachtungen)."
            )


        if self.autosave_var.get():

            self.labor.autosave_pfad = AUTOSAVE_PFAD

            self._on_fortschritt(
                f"Automatisches Sichern aktiv - alle "
                f"{self.labor.autosave_intervall} Experimente nach "
                f"'{self.labor.autosave_pfad}'."
            )


        thread = threading.Thread(
            target=self._labor_thread,
            daemon=True
        )


        thread.start()







    def _labor_thread(self):


        try:


            self.labor.starten(endlos=self._endlos_lauf)



        except Exception as fehler:


            self.fortschritt_queue.put(
                f"Fehler: {fehler}"
            )



        finally:


            self.fortschritt_queue.put(
                "__FERTIG__"
            )







    # ======================================================
    # Callback vom Labor
    # ======================================================


    def _on_ergebnis(
            self,
            nummer,
            molekuel,
            analyse,
            phase
    ):


        self.ergebnis_queue.put(
            (
                nummer,
                phase,
                molekuel,
                analyse
            )
        )





    def _on_fortschritt(
            self,
            text
    ):


        self.fortschritt_queue.put(
            text
        )






    # ======================================================
    # Queue Verarbeitung
    # Hauptthread
    # ======================================================


    def _poll_queues(self):



        while not self.fortschritt_queue.empty():


            text = self.fortschritt_queue.get_nowait()



            if text == "__FERTIG__":


                self.laeuft = False


                self.start_button.config(
                    state="normal"
                )


                self.stop_button.config(
                    state="disabled"
                )


                self.status_label.config(
                    text="Fertig."
                )



            else:


                self._log_anhaengen(
                    text
                )






        while not self.ergebnis_queue.empty():



            eintrag = self.ergebnis_queue.get_nowait()



            self.alle_ergebnisse_lokal.append(
                eintrag
            )
            self._ergebnis_zaehler += 1



        # Rendering ist ab hier komplett von der Datensammlung oben
        # entkoppelt (siehe _render_wenn_noetig) - laeuft bei JEDEM
        # Tick, nicht nur wenn neue Ergebnisse da sind, damit auch ein
        # reiner Tab-Wechsel ohne neue Daten prompt sichtbar wird.
        self._render_wenn_noetig()



        self.root.after(
            150,
            self._poll_queues
        )


    # ======================================================
    # Zentrales Rendering
    # ======================================================
    #
    # FIX (grundlegend): ersetzt die vorherige Kombination aus
    # Zaehler-Schwellen (_statistik_counter >= 25, "alle 50 neuen
    # Ergebnisse") und Sichtbarkeits-Gates an mehreren Stellen, die sich
    # gegenseitig aus dem Tritt bringen konnten - z.B. konnte waehrend
    # der gesamten Beobachtungsphase kein einziges Mal beides
    # gleichzeitig zutreffen ("Schwelle erreicht" UND "richtiger Tab
    # sichtbar"), sodass ein Tab bis zum Phasenwechsel leer blieb.
    #
    # Stattdessen: EIN Ort, EINE einfache Regel, bei JEDEM Tick (alle
    # 150ms) neu geprueft, unabhaengig von Lauf-Status oder Phase:
    # "Welcher Tab ist gerade sichtbar, und hat sich die Datenmenge seit
    # dessen letztem Rendern geaendert?" Kein Zaehler kann hier
    # "verpasst" werden, weil nichts mehr auf einen bestimmten
    # Schwellenwert wartet - der naechste Tick prueft einfach wieder.

    def _render_wenn_noetig(self):

        if not self.alle_ergebnisse_lokal:
            return

        anzahl = self._ergebnis_zaehler

        aktueller_tab = self.notebook.select()

        if aktueller_tab == str(self._statistik_tab_frame):

            if anzahl != self._statistik_render_stand:

                self._aktualisiere_statistik()

                self._statistik_render_stand = anzahl

        elif aktueller_tab == str(self._lernkurve_tab_frame):

            if anzahl != self._lernkurve_render_stand:

                self._zeichne_lernkurve()

                self._lernkurve_render_stand = anzahl

        elif aktueller_tab == str(self._ergebnisse_tab_frame):

            if anzahl != self._tabelle_render_stand:

                self._aktualisiere_ergebnisse_tabelle()

                self._tabelle_render_stand = anzahl


    # ======================================================
    # Statistik
    # ======================================================


    def _aktualisiere_statistik(self):


        ergebnisse = list(
            self.alle_ergebnisse_lokal
        )



        if not ergebnisse:

            return



        # FIX (Performance): vorher liefen hier SIEBEN separate volle
        # Durchlaeufe ueber "ergebnisse"/"punkte" (Listenbildung, je ein
        # sum() fuer beobachtung_n/automatik_n/stabil_n, dann noch
        # sum()/min()/max() fuer die Punkte) - bei vollem Ringpuffer
        # (bis zu 5000 Eintraege) und einem Aufruf alle 25 neuen
        # Ergebnisse lief das waehrend eines schnellen Automatiklaufs
        # staendig mit. Jetzt EIN Durchlauf, der alles gleichzeitig
        # mitzaehlt.

        punkte = []
        punkte_summe = 0.0
        punkte_min = None
        punkte_max = None
        beobachtung_n = 0
        automatik_n = 0
        stabil_n = 0

        for _nr, phase, _mol, analyse in ergebnisse:

            p = analyse["bewertung"]["punkte"]

            punkte.append(p)
            punkte_summe += p

            if punkte_min is None or p < punkte_min:
                punkte_min = p

            if punkte_max is None or p > punkte_max:
                punkte_max = p

            if phase == "Beobachtung":
                beobachtung_n += 1
            elif phase == "Automatik":
                automatik_n += 1

            if p >= 75:
                stabil_n += 1





        self.stat_labels[
            "Experimente gesamt"
        ].config(
            text=f"Experimente gesamt: {self._ergebnis_zaehler}" + (
                f" (Statistik ueber die letzten {len(ergebnisse)})"
                if self._ergebnis_zaehler > len(ergebnisse) else ""
            )
        )


        self.stat_labels[
            "Versuche gesamt"
        ].config(
            text=
            f"Versuche gesamt: {self.labor.versuche if self.labor else '-'}"
        )



        self.stat_labels[
            "Durchschnittspunkte"
        ].config(
            text=
            f"Durchschnittspunkte: {punkte_summe/len(punkte):.2f}"
        )



        self.stat_labels[
            "Minimum"
        ].config(
            text=f"Minimum: {punkte_min}"
        )



        self.stat_labels[
            "Maximum"
        ].config(
            text=f"Maximum: {punkte_max}"
        )



        self.stat_labels[
            "Anteil stabil (≥75)"
        ].config(
            text=
            f"Anteil stabil (≥75): "
            f"{stabil_n}/{len(punkte)} "
            f"({stabil_n/len(punkte)*100:.1f}%)"
        )



        self.stat_labels[
            "Anzahl Beobachtung"
        ].config(
            text=f"Anzahl Beobachtung: {beobachtung_n}"
        )



        self.stat_labels[
            "Anzahl Automatik"
        ].config(
            text=f"Anzahl Automatik: {automatik_n}"
        )



        if self.labor is not None and self.labor.ziel_verhaeltnis_neg_zu_pos is not None:

            filter_text = (
                f"Verhältnis-Filter (Ziel {self.labor.ziel_verhaeltnis_neg_zu_pos}:1): "
                f"{self.labor._positiv_gezaehlt} positiv, "
                f"{self.labor._negativ_gezaehlt} negativ, "
                f"{self.labor._negativ_verworfen} verworfen"
            )

        else:

            filter_text = "Verhältnis-Filter: aus"

        self.stat_labels[
            "Verhältnis-Filter"
        ].config(
            text=filter_text
        )



        # Zeichnet immer, wenn aufgerufen - WANN das sinnvoll ist,
        # entscheidet jetzt zentral _render_wenn_noetig() (sichtbarer
        # Tab + neue Daten), nicht mehr diese Funktion selbst. Das
        # ersetzt die vorherige, fehleranfaellige Kombination aus
        # Zaehler-Schwelle (_statistik_counter) und lokalem
        # Sichtbarkeits-Check hier drin, die dazu fuehren konnte, dass
        # waehrend einer ganzen Phase nie beides gleichzeitig zutraf und
        # der Tab leer blieb.
        self._zeichne_histogramm(
            punkte
        )







    # ======================================================
    # Histogramm
    # ======================================================


    def _zeichne_histogramm(
            self,
            punkte
    ):


        canvas = self.histogramm_canvas


        canvas.delete(
            "all"
        )



        breite = 520

        hoehe = 200



        buckets = [
            0
            for _ in range(10)
        ]



        for p in punkte:


            index = min(
                int(p)//10,
                9
            )


            buckets[index] += 1





        maximum = max(
            buckets
        )


        balken = breite / 10





        for i, anzahl in enumerate(buckets):


            if maximum:


                h = (
                    anzahl / maximum
                ) * (
                    hoehe - 30
                )


            else:

                h = 0



            x0 = i*balken + 3

            x1 = (i+1)*balken - 3


            y1 = hoehe - 20

            y0 = y1 - h



            canvas.create_rectangle(
                x0,
                y0,
                x1,
                y1,
                fill=self.FARBEN["akzent_dunkel"],
                outline=self.FARBEN["akzent"]
            )



            canvas.create_text(
                (x0+x1)/2,
                hoehe-10,
                text=str(i*10),
                font=("",7),
                fill=self.FARBEN["text_gedimmt"]
            )



            if anzahl:


                canvas.create_text(

                    (x0+x1)/2,

                    y0-8,

                    text=str(anzahl),

                    font=("",7),

                    fill=self.FARBEN["text"]

                )







    # ======================================================
    # Ergebnis Filter
    # ======================================================


    def _gefilterte_ergebnisse(self):


        ergebnisse = list(
            self.alle_ergebnisse_lokal
        )



        phase = self.phase_filter_var.get()



        if phase != "Alle":


            ergebnisse = [

                e

                for e in ergebnisse

                if e[1] == phase

            ]





        filter_wahl = self.filter_var.get()



        if filter_wahl == "Nur stabil (≥75)":


            ergebnisse = [

                e

                for e in ergebnisse

                if e[3]["bewertung"]["punkte"] >= 75

            ]



        elif filter_wahl == "Nur instabil (<50)":


            ergebnisse = [

                e

                for e in ergebnisse

                if e[3]["bewertung"]["punkte"] < 50

            ]




        # FIX (Performance): "Alle" filterte bisher gar nicht weiter -
        # bei vollem Ringpuffer (5000 Eintraege) wurden bei jedem
        # Tabellen-Update alle 5000 sortiert UND neu in die Treeview
        # eingefuegt (~50ms, alle paar Sekunden waehrend eines Laufs
        # spuerbar). Fuer "Alle" werden jetzt nur die zuletzt
        # entstandenen MAX_ANGEZEIGTE_ERGEBNISSE angezeigt - die
        # vollstaendigen Daten bleiben in alle_ergebnisse_lokal
        # erhalten (Statistik/Speichern nutzen weiterhin alles), nur
        # die Tabellendarstellung ist begrenzt.
        #
        # WICHTIG: diese Begrenzung passiert jetzt VOR der
        # Duplikat-Ausblendung (vorher danach) - im "Alle"-Fall macht
        # das den ohnehin teuersten Schritt (Signaturberechnung, siehe
        # unten) auf hoechstens 500 statt bis zu 5000 Eintraegen noetig.
        if filter_wahl == "Alle" and len(ergebnisse) > self.MAX_ANGEZEIGTE_ERGEBNISSE:

            ergebnisse = ergebnisse[-self.MAX_ANGEZEIGTE_ERGEBNISSE:]


        if self.duplikate_ausblenden_var.get():

            # Signaturen aus dem Cache statt jedes Mal neu berechnet
            # (siehe self._signatur_cache in __init__) - der Cache wird
            # hier gleich auf die gerade tatsaechlich betrachteten
            # Eintraege zurechtgestutzt, damit er nicht unbegrenzt
            # waechst, wenn alte Ergebnisse aus dem Ringpuffer fallen.
            neuer_cache = {}

            gesehen = set()
            eindeutig = []

            for e in ergebnisse:

                nummer = e[0]

                signatur = self._signatur_cache.get(nummer)

                if signatur is None:
                    signatur = molekuel_signatur(e[2])

                neuer_cache[nummer] = signatur

                if signatur in gesehen:
                    continue

                gesehen.add(signatur)
                eindeutig.append(e)

            self._signatur_cache = neuer_cache

            ergebnisse = eindeutig








        # FIX (Performance): "Alle" filterte bisher gar nicht weiter -
        # bei vollem Ringpuffer (5000 Eintraege) wurden bei jedem
        # Tabellen-Update alle 5000 sortiert UND neu in die Treeview
        # eingefuegt (~50ms, alle paar Sekunden waehrend eines Laufs
        # spuerbar). Fuer "Alle" werden jetzt nur die zuletzt
        # entstandenen MAX_ANGEZEIGTE_ERGEBNISSE angezeigt - die
        # vollstaendigen Daten bleiben in alle_ergebnisse_lokal
        # erhalten (Statistik/Speichern nutzen weiterhin alles), nur
        # die Tabellendarstellung ist begrenzt. (Diese Begrenzung
        # passiert bereits weiter oben, vor der Duplikat-Ausblendung -
        # hier nur noch die Sortierung.)

        ergebnisse.sort(

            key=lambda e:
            e[3]["bewertung"]["punkte"],

            reverse=True

        )





        if filter_wahl == "Top 10":

            ergebnisse = ergebnisse[:10]


        elif filter_wahl == "Top 50":

            ergebnisse = ergebnisse[:50]






        if self._sortier_spalte:


            key = self.SPALTEN_SCHLUESSEL[
                self._sortier_spalte
            ]


            ergebnisse.sort(

                key=key,

                reverse=self._sortier_absteigend

            )



        return ergebnisse







    # ======================================================
    # Tabelle aktualisieren
    # ======================================================


    def _aktualisiere_ergebnisse_tabelle(self):
        # [patch_fixes] Auswahl ueber das Neuaufbauen hinweg erhalten -
        # sonst geht die markierte Zeile bei jedem Refresh waehrend
        # eines Laufs verloren und die Detailansicht springt.

        auswahl = self.tabelle.selection()

        gemerkt = auswahl[0] if auswahl else None

        self._tabelle_leeren()

        gefiltert = self._gefilterte_ergebnisse()

        gesamt = len(self.alle_ergebnisse_lokal)

        if len(gefiltert) < gesamt and self.filter_var.get() == "Alle":

            self.ergebnisse_info_label.config(
                text=f"(zeigt die neuesten {len(gefiltert)} von {gesamt})"
            )

        else:

            self.ergebnisse_info_label.config(text="")

        for e in gefiltert:

            self.tabelle.insert(
                "",
                "end",
                iid=str(e[0]),
                # Werte generisch aus SPALTEN_SCHLUESSEL bauen (statt
                # eines fest verdrahteten Tupels), damit neue Spalten
                # (wie "rdkit") nicht aus dem Tritt geraten.
                values=tuple(
                    self.SPALTEN_SCHLUESSEL[spalte](e)
                    for spalte in self.ERGEBNIS_SPALTEN
                )
            )

        if gemerkt is not None and self.tabelle.exists(gemerkt):

            # loest <<TreeviewSelect>> aus; _zeige_detail() erkennt,
            # dass dieses Molekuel schon angezeigt wird, und macht nichts.
            self.tabelle.selection_set(gemerkt)









    def _sortiere_tabelle(
            self,
            spalte
    ):


        if self._sortier_spalte == spalte:


            self._sortier_absteigend = not self._sortier_absteigend



        else:


            self._sortier_spalte = spalte

            self._sortier_absteigend = False



        self._aktualisiere_ergebnisse_tabelle()







    # ======================================================
    # Detailansicht
    # ======================================================


    def _zeige_detail(
            self,
            event
    ):


        auswahl = self.tabelle.selection()



        if not auswahl:

            return



        nummer = int(
            auswahl[0]
        )

        # [patch_fixes] Zeile schon angezeigt (z.B. nach Tabellen-Refresh)?
        if nummer == self._angezeigte_nummer:
            return

        self._angezeigte_nummer = nummer



        eintrag = next(

            (
                e

                for e

                in self.alle_ergebnisse_lokal

                if e[0] == nummer

            ),

            None

        )



        if not eintrag:

            return



        _, phase, molekuel, analyse = eintrag

        # aktuelles Molekül für Lewis-Ansicht merken

        self.aktuelles_molekuel = molekuel


        # Lewis-Zeichnung aktualisieren

        self._zeichne_lewis(
            molekuel
        )

        text = [

            f"Experiment {nummer}",

            f"Phase: {phase}",

            f"Formel: {formel(molekuel)}",

            "",

            f"Punkte: {analyse['bewertung']['punkte']}",

            f"Urteil: {analyse['bewertung']['urteil']}",

            "",

        ]


        # RDKit-Ansicht: unabhaengig vom eigenen Kritiker/formel(), auf
        # echten chemischen Regeln basierend (siehe rdkit_bruecke.py).
        # "rdkit_gueltig" wurde bereits von ChemLabor.experiment() in
        # jede Analyse geschrieben - hier nur noch anzeigen, nicht neu
        # berechnen, damit ein aus einer geladenen/aelteren Datei
        # stammendes Ergebnis ohne dieses Feld nicht abstuerzt.

        rdkit_gueltig = analyse.get("rdkit_gueltig")

        text.append(f"RDKit-gültig: {'Ja' if rdkit_gueltig else 'Nein'}")

        if rdkit_gueltig:

            formel_rdkit = rb.summenformel(molekuel)
            masse = rb.molmasse(molekuel)
            signatur = rb.kanonische_signatur(molekuel)

            if formel_rdkit:
                text.append(f"RDKit-Summenformel: {formel_rdkit}")

            if masse is not None:
                text.append(f"Molmasse: {masse:.2f} g/mol")

            if signatur:
                text.append(f"Kanonischer SMILES: {signatur}")

        text.append("")

        text.append("Probleme:")



        # Ausgelagert in die modulweite Funktion formatiere_probleme()
        # (siehe oben) - dort sass der Crash-Bug, und dort ist er jetzt
        # auch ohne Tkinter/GUI-Kontext testbar (siehe test_labor_gui.py).
        text.extend(formatiere_probleme(analyse["probleme"]))



        text.append("")

        text.append(
            "Erklärung:"
        )



        for e in analyse["erklaerung"]:

            text.append(
                " - " + e
            )





        self.detail_text.config(
            state="normal"
        )


        self.detail_text.delete(
            "1.0",
            "end"
        )


        self.detail_text.insert(
            "1.0",
            "\n".join(text)
        )


        self.detail_text.config(
            state="disabled"
        )



    # ======================================================
    # Lewis Analyse
    # ======================================================


    def _pruefe_lewis(
            self,
            molekuel
    ):


        # Keine Bindungen vorhanden

        if not molekuel.get("bonds"):

            return {
                "moeglich": False,
                "grund":
                    "Keine Bindungsdaten vorhanden."
            }





        # Atome prüfen

        for atom in molekuel.get("atoms", []):


            valenz = atom.get(
                "valence",
                None
            )


            benutzt = atom.get(
                "used_valence",
                None
            )



            if valenz is None or benutzt is None:


                return {
                    "moeglich": False,
                    "grund":
                        "Valenzdaten fehlen."
                }




            if benutzt > valenz:


                return {
                    "moeglich": False,
                    "grund":
                        f"Ungültige Valenz bei "
                        f"{atom['element']}."
                }





        # Bindungen prüfen

        for bindung in molekuel["bonds"]:


            if (
                "atom1" not in bindung
                or
                "atom2" not in bindung
            ):


                return {
                    "moeglich": False,
                    "grund":
                        "Unvollständige Bindungsdaten."
                }




        return {
            "moeglich": True,
            "grund":
                "Darstellung möglich."
        }








    # ======================================================
    # Layout-Berechnung (Force-Directed / Fruchterman-Reingold)
    # ======================================================
    #
    # Ersetzt die alte "Atom 0 in die Mitte, Rest auf einen Kreis mit
    # festem Radius"-Platzierung. Bei mehr als ~15 Atomen ueberlappten
    # sich damit zwangslaeufig Kreise und Linien, weil der Kreisumfang
    # nicht mit der Atomanzahl mitwuchs UND die tatsaechliche
    # Bindungstopologie komplett ignoriert wurde (ein langes Kettenmolekuel
    # sah genauso aus wie ein Ring).
    #
    # Fruchterman-Reingold-Prinzip, aus dem Stand implementiert (nur
    # math/random, keine neue Abhaengigkeit):
    #   - JEDES Atompaar stoesst sich ab (wie gleichnamige Ladungen) ->
    #     verhindert Ueberlappung, verteilt die Struktur im Raum.
    #   - NUR gebundene Atompaare ziehen sich zusaetzlich an (wie eine
    #     Feder) -> haelt tatsaechlich verbundene Atome nah beieinander,
    #     Ketten werden zu Ketten, Ringe zu Ringen.
    #   - Eine mit der Zeit sinkende Schrittweite ("Abkuehlung") sorgt
    #     dafuer, dass sich das System beruhigt statt endlos zu oszillieren.
    # ======================================================


    def _berechne_layout(
            self,
            atome,
            bindungen,
            breite,
            hoehe
    ):

        if not atome:
            return {}

        ids = [atom["id"] for atom in atome]
        n = len(ids)

        rand = 50  # Sicherheitsabstand zum Canvas-Rand

        # Ideale Distanz zwischen zwei Knoten, hergeleitet aus der
        # verfuegbaren Flaeche (klassische Fruchterman-Reingold-Formel).
        flaeche = max(breite - 2 * rand, 50) * max(hoehe - 2 * rand, 50)
        k = math.sqrt(flaeche / n) * 0.85

        # Weniger Iterationen bei sehr grossen Molekuelen, damit die
        # Berechnung (O(n²) pro Iteration wegen der paarweisen Abstossung)
        # auch bei hoher "Max. Atome"-Einstellung nicht spuerbar haengt.
        iterationen = 250 if n <= 40 else max(60, int(9000 / n))

        zufall = random.Random(42)  # deterministisch: gleiches Molekuel -> gleiches Bild

        pos = {
            aid: [
                zufall.uniform(rand, breite - rand),
                zufall.uniform(rand, hoehe - rand)
            ]
            for aid in ids
        }

        kanten = [
            (b["atom1"], b["atom2"])
            for b in bindungen
            if b["atom1"] in pos and b["atom2"] in pos
        ]

        temperatur = (breite + hoehe) / 20

        for _ in range(iterationen):

            kraft = {aid: [0.0, 0.0] for aid in ids}

            # Abstossung: alle Atompaare (O(n²), fuer Molekuelgroessen
            # hier unproblematisch)
            for i in range(n):
                for j in range(i + 1, n):

                    a, b = ids[i], ids[j]
                    dx = pos[a][0] - pos[b][0]
                    dy = pos[a][1] - pos[b][1]
                    dist = math.hypot(dx, dy) or 0.01

                    f = (k * k) / dist
                    fx, fy = dx / dist * f, dy / dist * f

                    kraft[a][0] += fx
                    kraft[a][1] += fy
                    kraft[b][0] -= fx
                    kraft[b][1] -= fy

            # Anziehung nur entlang tatsaechlicher Bindungen
            for a, b in kanten:

                dx = pos[a][0] - pos[b][0]
                dy = pos[a][1] - pos[b][1]
                dist = math.hypot(dx, dy) or 0.01

                f = (dist * dist) / k
                fx, fy = dx / dist * f, dy / dist * f

                kraft[a][0] -= fx
                kraft[a][1] -= fy
                kraft[b][0] += fx
                kraft[b][1] += fy

            # Positionen um die Nettokraft verschieben, auf maximal
            # "temperatur" pro Schritt begrenzt, und im Canvas halten
            for aid in ids:

                fx, fy = kraft[aid]
                dist = math.hypot(fx, fy) or 0.01
                schritt = min(dist, temperatur)

                pos[aid][0] += fx / dist * schritt
                pos[aid][1] += fy / dist * schritt

                pos[aid][0] = min(breite - rand, max(rand, pos[aid][0]))
                pos[aid][1] = min(hoehe - rand, max(rand, pos[aid][1]))

            temperatur *= 0.94  # abkuehlen, sonst oszilliert es endlos

        return {aid: (x, y) for aid, (x, y) in pos.items()}


    # ======================================================
    # Lewis Zeichner
    # ======================================================


    def _zeichne_lewis(
            self,
            molekuel
    ):
        """Einstiegspunkt: prueft, ob eine Lewis-Darstellung moeglich
        ist, und schickt die eigentliche Layout-BERECHNUNG (das teure,
        O(n²)-pro-Iteration Force-Directed-Verfahren) in einen
        Hintergrundthread.

        FIX (Performance): vorher lief die komplette Berechnung
        synchron hier im GUI-Thread - bei 100 Atomen ca. 170ms, bei 400
        Atomen ueber 1,4 SEKUNDEN, waehrend derer Tkinter komplett
        eingefroren war (nicht mal das Fenster liess sich verschieben).
        Jetzt berechnet ein Hintergrundthread die Positionen, das
        eigentliche Zeichnen (schnell, reine Canvas-Aufrufe) passiert
        weiterhin im Hauptthread, sobald das Ergebnis ueber
        lewis_queue eintrifft (siehe _poll_lewis_queue).

        Eine auftrag_id verhindert, dass ein spaeter verworfener Klick
        (Nutzer waehlt schnell mehrere Zeilen nacheinander) noch mit
        Verzoegerung ueber ein aktuelleres Ergebnis gezeichnet wird -
        nur der Auftrag mit der zuletzt vergebenen ID zaehlt."""

        canvas = self.lewis_canvas

        canvas.delete("all")

        pruefung = self._pruefe_lewis(molekuel)

        if not pruefung["moeglich"]:

            self.lewis_status.config(
                text="Lewis nicht verfügbar: " + pruefung["grund"]
            )

            self._zeige_lewis_rohdaten(molekuel)

            return

        atome = molekuel["atoms"]
        bindungen = molekuel["bonds"]
        n = len(atome)

        self.lewis_status.config(
            text=f"Layout wird berechnet ({n} Atome, {len(bindungen)} Bindungen)..."
        )

        canvas.update_idletasks()
        breite = max(canvas.winfo_width(), 400)
        hoehe = max(canvas.winfo_height(), 300)

        self._lewis_auftrag_id += 1
        auftrag_id = self._lewis_auftrag_id

        def berechne_im_hintergrund():

            positionen = self._berechne_layout(atome, bindungen, breite, hoehe)

            self.lewis_queue.put((auftrag_id, molekuel, positionen))

        threading.Thread(target=berechne_im_hintergrund, daemon=True).start()


    def _poll_lewis_queue(self):
        """Laeuft alle 50ms im Hauptthread (wie _poll_queues fuer das
        Labor) und holt fertige Layout-Ergebnisse aus lewis_queue.
        Verwirft veraltete Auftraege (auftrag_id passt nicht mehr zur
        zuletzt angeforderten Berechnung), damit bei schnellem
        Durchklicken mehrerer Molekuele nicht nachtraeglich ein altes
        Layout ueber das neue gezeichnet wird."""

        try:

            while True:

                auftrag_id, molekuel, positionen = self.lewis_queue.get_nowait()

                if auftrag_id == self._lewis_auftrag_id:
                    self._zeichne_lewis_mit_positionen(molekuel, positionen)

        except queue.Empty:

            pass

        self.root.after(50, self._poll_lewis_queue)


    def _zeichne_lewis_mit_positionen(
            self,
            molekuel,
            positionen
    ):
        """Reines Zeichnen (Canvas-Operationen) mit BEREITS berechneten
        Positionen - schnell genug, um im Hauptthread zu laufen, ohne
        die GUI spuerbar zu blockieren (siehe _zeichne_lewis fuer die
        ausgelagerte, teure Berechnung)."""

        canvas = self.lewis_canvas

        canvas.delete("all")

        atome = molekuel["atoms"]

        bindungen = molekuel["bonds"]

        n = len(atome)

        self.lewis_status.config(
            text=
            f"Lewis-Struktur erzeugt ({n} Atome, {len(bindungen)} Bindungen)."
        )

        # Atomradius und Schrift schrumpfen mit wachsender Atomanzahl,
        # damit sich benachbarte Atome bei grossen Molekuelen nicht
        # gegenseitig verdecken.
        atom_radius = max(10, min(25, 260 // max(n, 1)))
        schrift_groesse = max(7, min(14, 140 // max(n, 1)))



        # --------------------------------------------------
        # Bindungen zeichnen
        # --------------------------------------------------


        for bindung in bindungen:


            a1 = bindung["atom1"]

            a2 = bindung["atom2"]



            if (
                a1 not in positionen
                or
                a2 not in positionen
            ):

                continue



            x1,y1 = positionen[a1]

            x2,y2 = positionen[a2]



            ordnung = bindung.get(
                "order",
                1
            )


            # Bindungslinie(n) senkrecht zur Verbindungsachse versetzen
            # (nicht mehr nur vertikal wie vorher) - bei schraeg
            # verlaufenden Bindungen sonst optisch falsch.
            dx, dy = x2 - x1, y2 - y1
            laenge = math.hypot(dx, dy) or 1
            nx, ny = -dy / laenge, dx / laenge  # Normalenvektor



            for offset in range(
                ordnung
            ):


                versatz = (
                    offset
                    -
                    (ordnung-1)/2
                ) * 5



                canvas.create_line(

                    x1 + nx * versatz,

                    y1 + ny * versatz,

                    x2 + nx * versatz,

                    y2 + ny * versatz,

                    width=2,

                    fill=self.FARBEN["text"]

                )







        # --------------------------------------------------
        # Atome zeichnen
        # --------------------------------------------------


        for atom in atome:


            x,y = positionen[
                atom["id"]
            ]

            farbe = self.ELEMENT_FARBEN.get(atom["element"], "#ffffff")
            textfarbe = "#ffffff" if farbe in ("#4a4a4a", "#3060c8", "#c83030", "#2f9e5f", "#a8763f") else "#000000"



            canvas.create_oval(

                x-atom_radius,

                y-atom_radius,

                x+atom_radius,

                y+atom_radius,

                fill=farbe,
                outline=self.FARBEN["rand"]

            )



            canvas.create_text(

                x,

                y,

                text=atom["element"],

                font=("", schrift_groesse, "bold"),
                fill=textfarbe

            )









    # ======================================================
    # Rohdaten Fallback
    # ======================================================


    # ======================================================
    # JSON-Ansicht
    # ======================================================
    #
    # Zeigt das aktuell ausgewaehlte Molekuel (Atome + Bindungen) als
    # rohes JSON in einem eigenen Fenster - dieselben Daten, die auch
    # an den Kritiker gehen, aber ungefiltert und maschinenlesbar,
    # statt nur als gezeichnete Lewis-Struktur oder als Fliesstext in
    # der Detailansicht.
    # ======================================================


    def _zeige_json(self):

        if self.aktuelles_molekuel is None:

            self.status_label.config(
                text="Kein Molekül ausgewählt - erst eine Zeile in "
                     "'Ergebnisse' anklicken."
            )

            return

        fenster = tk.Toplevel(self.root)

        fenster.title("Molekül als JSON")

        fenster.geometry("500x600")

        fenster.configure(bg=self.FARBEN["bg"])

        text_widget = scrolledtext.ScrolledText(
            fenster,
            wrap="word",
            bg=self.FARBEN["panel_hell"],
            fg=self.FARBEN["text"],
            insertbackground=self.FARBEN["text"],
            borderwidth=0,
            highlightthickness=0
        )

        text_widget.pack(
            fill="both",
            expand=True
        )

        text_widget.insert(
            "1.0",
            json.dumps(
                self.aktuelles_molekuel,
                indent=2,
                ensure_ascii=False
            )
        )

        text_widget.config(
            state="disabled"
        )


    def _zeige_lewis_rohdaten(
            self,
            molekuel
    ):


        canvas = self.lewis_canvas


        y = 40



        canvas.create_text(

            400,

            y,

            text="Rohdaten:",

            font=("",14,"bold"),

            fill=self.FARBEN["akzent"]

        )


        y += 40



        for atom in molekuel["atoms"]:


            text = (

                f"Atom {atom['id']}: "
                f"{atom['element']} "
                f"(Valenz {atom.get('valence','?')}, "
                f"belegt {atom.get('used_valence','?')})"

            )


            canvas.create_text(

                400,

                y,

                text=text,

                fill=self.FARBEN["text"]

            )


            y += 25



        y += 20



        for bindung in molekuel["bonds"]:


            text = (

                f"Bindung: "
                f"{bindung.get('atom1')} - "
                f"{bindung.get('atom2')} "
                f"(Ordnung "
                f"{bindung.get('order','?')})"

            )


            canvas.create_text(

                400,

                y,

                text=text,

                fill=self.FARBEN["text"]

            )


            y += 25



    # ======================================================
    # Hilfsfunktionen
    # ======================================================


    # Deckelt die Zeilenzahl im Log-Textfeld - ohne das waechst es bei
    # einem langen (v.a. Endlos-)Lauf unbegrenzt, und Tkinters
    # Text-Widget wird mit wachsendem Inhalt spuerbar traeger beim
    # Einfuegen/Scrollen. Aeltere Zeilen werden einfach abgeschnitten,
    # der vollstaendige Verlauf ist ohnehin nicht der Zweck dieses
    # Live-Logs (die Ergebnisse selbst bleiben unabhaengig davon in
    # alle_ergebnisse_lokal erhalten).
    LOG_MAX_ZEILEN = 2000

    def _log_anhaengen(
            self,
            text
    ):


        self.log_text.config(
            state="normal"
        )


        self.log_text.insert(
            "end",
            text+"\n"
        )


        # "end-1c" zeigt auf den letzten Zeichenindex; die Zeilenzahl
        # davor ist die Gesamtzahl der Zeilen im Widget.
        zeilenzahl = int(
            self.log_text.index("end-1c").split(".")[0]
        )

        if zeilenzahl > self.LOG_MAX_ZEILEN:

            self.log_text.delete(
                "1.0",
                f"{zeilenzahl - self.LOG_MAX_ZEILEN}.0"
            )


        self.log_text.see(
            "end"
        )


        self.log_text.config(
            state="disabled"
        )





    def _log_leeren(self):


        self.log_text.config(
            state="normal"
        )


        self.log_text.delete(
            "1.0",
            "end"
        )


        self.log_text.config(
            state="disabled"
        )







    def _tabelle_leeren(self):

        # FIX (Performance): vorher wurde jede Zeile einzeln geloescht
        # (eine Tcl-Bridge-Rueckfrage pro Zeile - bei 500 Zeilen 500
        # einzelne Aufrufe). Treeview.delete() akzeptiert mehrere
        # Item-IDs auf einmal - EIN Aufruf statt bis zu 500 macht die
        # komplette Tabellenaktualisierung spuerbar fluessiger, gerade
        # weil sie bei laufendem Labor regelmaessig ausgefuehrt wird.

        kinder = self.tabelle.get_children()

        if kinder:
            self.tabelle.delete(*kinder)








# ==========================================================
# Start
# ==========================================================


def main():


    root = tk.Tk()


    LaborGUI(
        root
    )


    root.mainloop()



if __name__ == "__main__":

    main()