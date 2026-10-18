# ==========================================================
# tab_ergebnisse.py
# ==========================================================
#
# Tab "Ergebnisse": Tabelle mit Filtern, rechts Lewis-Struktur und
# Detailtext. Die PubChem-Pruefung (frueher eigener Tab "Entdeckungen")
# ist jetzt eine Spalte dieser Tabelle: "Neuheit prüfen" fragt die
# gerade SICHTBAREN Zeilen ab, "Anzeige: Nur unbekannt" zeigt danach nur
# die Kandidaten.
#
# Die Filter-Comboboxen zeigen uebersetzten Text, gelesen wird aber ueber
# den Index (Combobox.current()) und die internen Schluessel aus
# modell.py - so haengt die Logik nicht von der Sprache ab.

import json
import queue
import threading
import tkinter as tk
import webbrowser
from tkinter import filedialog, scrolledtext, ttk

from chemlabor.chemie import neuheit
from chemlabor.gui.kritiker_texte import uebersetze_urteil
from chemlabor.gui.lewis import LewisAnsicht
from chemlabor.gui.modell import (
    FILTER,
    NEUHEIT_FILTER,
    PHASEN,
    STABIL_AB,
    baue_detailtext,
    phase_text,
    punkte_von,
)
from chemlabor.gui.theme import FARBEN
from chemlabor.labor import formel
from chemlabor.texte import t


SPALTEN = ("nr", "phase", "formel", "atome", "punkte", "urteil", "rdkit", "pubchem", "qc")

BREITE = {
    "nr": 50, "phase": 90, "formel": 110, "atome": 55,
    "punkte": 60, "urteil": 140, "rdkit": 55, "pubchem": 90, "qc": 90,
}

# So viele verschiedene Strukturen werden pro Klick bei PubChem abgefragt
# (jede neue Struktur kostet ~0,25 s wegen des PubChem-Limits).
NEUHEIT_MAX_STRUKTUREN = 300

# Quantenchemie (ANI-2x / MACE): so viele verschiedene Strukturen pro Klick, beste Kritiker-Punkte zuerst,
# und nur bis zu dieser Zahl schwerer Atome (sonst dauert die Optimierung zu lange)
QC_MAX_STRUKTUREN = 20
QC_MAX_SCHWERE_ATOME = 25


def beschriftung(spalte):
    """Spaltenueberschrift in der aktiven Sprache."""

    if spalte == "qc":
        return "ΔHf (QC)"

    return t(f"spalte.{spalte}")


class TabErgebnisse(ttk.Frame):

    def __init__(self, master, modell, log=None):

        super().__init__(master)

        self.modell = modell
        self.log = log or (lambda text: None)

        self.sortierung = None
        self._sichtbar = []
        self._unterdruecke_auswahl = False

        self.neuheit_queue = queue.Queue()
        self._neuheit_laeuft = False
        self._neuheit_abbruch = threading.Event()
        self._neuheit_faden = None
        self._poll_id = None

        self.qc_ergebnisse = {}
        self.qc_queue = queue.Queue()
        self._qc_laeuft = False
        self._qc_abbruch = threading.Event()
        self._qc_faden = None

        self.duplikate_var = tk.BooleanVar(value=False)
        self.offline_var = tk.BooleanVar(value=False)

        self._neuheit_schluessel = tuple(NEUHEIT_FILTER)

        self._baue_leisten()
        self._baue_inhalt()

        self._poll_neuheit()

    # ------------------------------------------------------
    # Aufbau
    # ------------------------------------------------------

    def _baue_leisten(self):

        oben = ttk.Frame(self, padding=(8, 8, 8, 2))
        oben.pack(fill="x")

        ttk.Label(oben, text=t("erg.filter")).pack(side="left")

        self.filter_box = ttk.Combobox(
            oben,
            values=[t(f"filter.{k}", ab=STABIL_AB) for k in FILTER],
            state="readonly", width=40
        )
        self.filter_box.current(0)
        self.filter_box.pack(side="left", padx=(4, 12))
        self.filter_box.bind("<<ComboboxSelected>>", lambda _e: self.aktualisiere())

        ttk.Label(oben, text=t("erg.phase")).pack(side="left")

        self.phase_box = ttk.Combobox(
            oben,
            values=[t(f"phase.{p}") for p in PHASEN],
            state="readonly", width=14
        )
        self.phase_box.current(0)
        self.phase_box.pack(side="left", padx=(4, 12))
        self.phase_box.bind("<<ComboboxSelected>>", lambda _e: self.aktualisiere())

        ttk.Checkbutton(
            oben, text=t("erg.duplikate"),
            variable=self.duplikate_var, command=self.aktualisiere
        ).pack(side="left", padx=(0, 12))

        self.anzahl_label = ttk.Label(oben, text="", style="Gedimmt.TLabel")
        self.anzahl_label.pack(side="right")

        unten = ttk.Frame(self, padding=(8, 2, 8, 4))
        unten.pack(fill="x")

        self.neuheit_button = ttk.Button(unten, text=t("erg.neuheit_button"), command=self._neuheit_klick)
        self.neuheit_button.pack(side="left", padx=(0, 8))

        ttk.Checkbutton(unten, text=t("erg.offline"), variable=self.offline_var).pack(side="left", padx=(0, 12))

        ttk.Label(unten, text=t("erg.pubchem_anzeige")).pack(side="left")

        self.neuheit_box = ttk.Combobox(
            unten,
            values=[t(f"neufilter.{k}") for k in self._neuheit_schluessel],
            state="readonly", width=16
        )
        self.neuheit_box.current(0)
        self.neuheit_box.pack(side="left", padx=(4, 12))
        self.neuheit_box.bind("<<ComboboxSelected>>", lambda _e: self.aktualisiere())

        ttk.Button(unten, text=t("erg.export_button"), command=self._export_klick).pack(side="left")

        self.qc_button = ttk.Button(unten, text="Quantenchemie (Top 20)", command=self._qc_klick)
        self.qc_button.pack(side="left", padx=(12, 0))

        self.status = ttk.Label(self, text="", style="Gedimmt.TLabel", padding=(8, 0, 8, 4))
        self.status.pack(fill="x")

    def _baue_inhalt(self):

        teiler = ttk.PanedWindow(self, orient="horizontal")
        teiler.pack(fill="both", expand=True)

        links = ttk.Frame(teiler)
        teiler.add(links, weight=3)

        self.tabelle = ttk.Treeview(links, columns=SPALTEN, show="headings", selectmode="browse")

        for spalte in SPALTEN:

            self.tabelle.heading(
                spalte,
                text=beschriftung(spalte),
                command=(lambda s=spalte: self._sortiere(s)) if spalte != "qc" else (lambda: None)
            )

            self.tabelle.column(
                spalte, width=BREITE[spalte],
                anchor="w" if spalte in ("phase", "formel", "urteil") else "center"
            )

        self.tabelle.tag_configure("unbekannt", foreground=FARBEN["akzent"])
        self.tabelle.tag_configure("ungeprueft", foreground=FARBEN["text_gedimmt"])

        scroll = ttk.Scrollbar(links, orient="vertical", command=self.tabelle.yview)
        self.tabelle.configure(yscrollcommand=scroll.set)

        self.tabelle.pack(side="left", fill="both", expand=True)
        scroll.pack(side="left", fill="y")

        self.tabelle.bind("<<TreeviewSelect>>", self._auswahl)
        self.tabelle.bind("<Double-1>", self._doppelklick)

        rechts = ttk.PanedWindow(teiler, orient="vertical")
        teiler.add(rechts, weight=2)

        # Rechts oben: Reiter "Lewis" und "3D" (siehe ansicht3d.py)
        self._ansichten = ttk.Notebook(rechts)
        rechts.add(self._ansichten, weight=3)

        self.lewis = LewisAnsicht(self._ansichten)
        self._ansichten.add(self.lewis, text="Lewis")

        self.ansicht3d = None
        self._3d_offen = None

        try:

            from chemlabor.gui.ansicht3d import Ansicht3D

            self.ansicht3d = Ansicht3D(
                self._ansichten,
                hintergrund=FARBEN["panel_hell"],
                vordergrund=FARBEN["text"]
            )

            self._ansichten.add(self.ansicht3d, text="3D")
            self._ansichten.bind("<<NotebookTabChanged>>", self._3d_tab_gewechselt)

        except ImportError as fehler:

            hinweis = ttk.Frame(self._ansichten)

            ttk.Label(
                hinweis,
                text=f"3D-Ansicht nicht verfügbar: {fehler}\n(Abhilfe: pip install matplotlib numpy)"
            ).pack(padx=10, pady=10)

            self._ansichten.add(hinweis, text="3D")

        self.detail = scrolledtext.ScrolledText(
            rechts, wrap="word", height=12, state="disabled",
            bg=FARBEN["panel_hell"], fg=FARBEN["text"],
            insertbackground=FARBEN["text"], borderwidth=0, highlightthickness=0
        )
        rechts.add(self.detail, weight=2)

    # ------------------------------------------------------
    # Filterauswahl (Index -> interner Schluessel)
    # ------------------------------------------------------

    @staticmethod
    def _wahl(box, schluessel):
        """Interner Schluessel zur aktuellen Auswahl der Combobox."""

        index = box.current()

        return schluessel[index] if index >= 0 else schluessel[0]

    # ------------------------------------------------------
    # Tabelle
    # ------------------------------------------------------

    def leeren(self):

        self._sichtbar = []

        kinder = self.tabelle.get_children()

        if kinder:
            self.tabelle.delete(*kinder)

        self.qc_ergebnisse.clear()
        self.lewis.leeren()
        self._setze_detail("")
        self.anzahl_label.config(text="")
        self.status.config(text="")

    def aktualisiere(self):

        ausgewaehlt = self.tabelle.selection()

        self._sichtbar = self.modell.filtere(
            filter_wahl=self._wahl(self.filter_box, FILTER),
            phase=self._wahl(self.phase_box, PHASEN),
            duplikate_ausblenden=self.duplikate_var.get(),
            neuheit_filter=self._wahl(self.neuheit_box, self._neuheit_schluessel),
            sortierung=self.sortierung,
        )

        kinder = self.tabelle.get_children()

        if kinder:
            self.tabelle.delete(*kinder)

        for e in self._sichtbar:

            nummer, phase, molekuel, analyse = e

            pruefung = self.modell.neuheit.get(nummer)

            self.tabelle.insert(
                "", "end", iid=str(nummer),
                values=(
                    nummer,
                    phase_text(phase),
                    formel(molekuel),
                    len(molekuel["atoms"]),
                    punkte_von(e),
                    uebersetze_urteil(analyse["bewertung"]["urteil"]),
                    t("rdkit.ja_zelle") if analyse.get("rdkit_gueltig") else t("rdkit.nein_zelle"),
                    self.modell.pubchem_text(nummer),
                    self._qc_text(nummer),
                ),
                tags=(pruefung["status"],) if pruefung else ()
            )

        self.anzahl_label.config(
            text=t("erg.anzahl", n=len(self._sichtbar), gesamt=self.modell.zaehler)
        )

        # Auswahl nach dem Neuaufbau wiederherstellen, ohne Lewis neu zu berechnen
        if ausgewaehlt and self.tabelle.exists(ausgewaehlt[0]):

            self._unterdruecke_auswahl = True

            self.tabelle.selection_set(ausgewaehlt[0])

            self.after_idle(self._auswahl_freigeben)

    def _auswahl_freigeben(self):

        self._unterdruecke_auswahl = False

    def _sortiere(self, spalte):

        if self.sortierung is not None and self.sortierung[0] == spalte:
            absteigend = not self.sortierung[1]     # gleiche Spalte: Richtung umkehren
        else:
            absteigend = True                        # neue Spalte: zuerst absteigend

        self.sortierung = (spalte, absteigend)

        for s in SPALTEN:
            pfeil = (" ▼" if absteigend else " ▲") if s == spalte else ""
            self.tabelle.heading(s, text=beschriftung(s) + pfeil)

        self.aktualisiere()

    def _eintrag_zur_auswahl(self):

        auswahl = self.tabelle.selection()

        if not auswahl:
            return None

        return self.modell.eintrag(int(auswahl[0]))

    def _auswahl(self, _ereignis):

        if self._unterdruecke_auswahl:
            return

        eintrag = self._eintrag_zur_auswahl()

        if eintrag is None:
            return

        self.lewis.zeige(eintrag[2])
        self._melde_3d(eintrag[2], eintrag[3])

        self._setze_detail(
            baue_detailtext(eintrag, self.modell.neuheit.get(eintrag[0])) + self._qc_detail(eintrag[0])
        )

    # ------------------------------------------------------
    # 3D-Ansicht: gezeichnet wird erst, wenn der Reiter sichtbar ist
    # ------------------------------------------------------

    def _melde_3d(self, molekuel, analyse):

        if self.ansicht3d is None:
            return

        self._3d_offen = (molekuel, analyse)

        if self._ansichten.select() == str(self.ansicht3d):
            self._zeichne_3d()

    def _3d_tab_gewechselt(self, _ereignis=None):

        if self.ansicht3d is not None and self._ansichten.select() == str(self.ansicht3d):
            self._zeichne_3d()

    def _zeichne_3d(self):

        if self._3d_offen is None:
            return

        molekuel, analyse = self._3d_offen
        self._3d_offen = None

        try:
            self.ansicht3d.zeige(molekuel, analyse)
        except Exception as fehler:
            self.status.config(text=f"3D-Ansicht fehlgeschlagen: {type(fehler).__name__}: {fehler}")

    def _doppelklick(self, _ereignis):

        eintrag = self._eintrag_zur_auswahl()

        if eintrag is None:
            return

        pruefung = self.modell.neuheit.get(eintrag[0])

        if pruefung and pruefung.get("cid"):
            webbrowser.open(f"https://pubchem.ncbi.nlm.nih.gov/compound/{pruefung['cid']}")

    def _setze_detail(self, text):

        self.detail.config(state="normal")
        self.detail.delete("1.0", "end")
        self.detail.insert("1.0", text)
        self.detail.config(state="disabled")

    # ------------------------------------------------------
    # PubChem-Pruefung
    # ------------------------------------------------------

    def _neuheit_klick(self):
        """Startet die Pruefung der sichtbaren Zeilen - oder bricht sie
        ab, falls sie schon laeuft."""

        if self._neuheit_laeuft:

            self._neuheit_abbruch.set()

            self.neuheit_button.config(state="disabled")
            self.status.config(text=t("erg.breche_ab"))

            return

        if not self._sichtbar:

            self.status.config(text=t("erg.keine_zeilen"))

            return

        self._neuheit_abbruch = threading.Event()
        self._neuheit_laeuft = True

        self.neuheit_button.config(text=t("erg.abbrechen"))
        self.status.config(text=t("erg.starte"))

        self._neuheit_faden = threading.Thread(
            target=self._neuheit_worker,
            args=(list(self._sichtbar), self.offline_var.get(), self._neuheit_abbruch),
            daemon=True
        )

        self._neuheit_faden.start()

    def _neuheit_worker(self, eintraege, offline, abbruch):

        pruefer = None

        try:

            # Der Pruefer (und seine SQLite-Verbindung) wird HIER im Thread
            # erzeugt - SQLite-Verbindungen gehoeren dem Thread, der sie anlegt.
            pruefer = neuheit.NeuheitsPruefer(offline=offline)

            ergebnisse = {}
            je_struktur = {}
            uebersprungen = 0

            beste_zuerst = sorted(eintraege, key=punkte_von, reverse=True)

            for i, e in enumerate(beste_zuerst, start=1):

                if abbruch.is_set():
                    break

                nummer, _phase, molekuel, analyse = e

                if not analyse.get("rdkit_gueltig"):
                    continue

                schluessel = neuheit.struktur_schluessel(molekuel)

                if schluessel is None:
                    continue

                key = schluessel["inchikey"]

                if key not in je_struktur:

                    if len(je_struktur) >= NEUHEIT_MAX_STRUKTUREN:
                        uebersprungen += 1
                        continue

                    je_struktur[key] = pruefer.pruefe_schluessel(schluessel)

                pruefung = dict(je_struktur[key])
                pruefung["nr"] = nummer
                ergebnisse[nummer] = pruefung

                if i % 20 == 0:
                    self.neuheit_queue.put(("fortschritt", t("erg.fortschritt", i=i, n=len(beste_zuerst))))

            self.neuheit_queue.put(("fertig", ergebnisse, abbruch.is_set(), uebersprungen))

        except Exception as fehler:

            self.neuheit_queue.put(("fehler", f"{type(fehler).__name__}: {fehler}"))

        finally:

            if pruefer is not None:
                pruefer.schliessen()

    def _poll_neuheit(self):

        try:

            while True:

                nachricht = self.neuheit_queue.get_nowait()

                art = nachricht[0]

                if art == "fortschritt":

                    self.status.config(text=nachricht[1])

                elif art == "fertig":

                    self._neuheit_beenden()

                    ergebnisse, abgebrochen, uebersprungen = nachricht[1], nachricht[2], nachricht[3]

                    self.modell.neuheit.update(ergebnisse)

                    self.aktualisiere()

                    text = self._zusammenfassung(ergebnisse)

                    if abgebrochen:
                        text += t("erg.abgebrochen")

                    if uebersprungen:
                        text += t("erg.limit", n=uebersprungen, max=NEUHEIT_MAX_STRUKTUREN)

                    self.status.config(text=text)
                    self.log(t("erg.log_pruefung", text=text))

                elif art == "fehler":

                    self._neuheit_beenden()

                    self.status.config(text=t("erg.pruefung_fehler", fehler=nachricht[1]))

        except queue.Empty:

            pass

        try:
            self._poll_qc()
        except Exception as fehler:
            self.status.config(text=f"Quantenchemie: {type(fehler).__name__}: {fehler}")

        self._poll_id = self.after(150, self._poll_neuheit)

    def _neuheit_beenden(self):

        self._neuheit_laeuft = False

        self.neuheit_button.config(text=t("erg.neuheit_button"), state="normal")

    @staticmethod
    def _zusammenfassung(ergebnisse):

        zaehler = {}
        stereo_offen = 0

        for p in ergebnisse.values():

            zaehler[p["status"]] = zaehler.get(p["status"], 0) + 1

            if p["status"] == "unbekannt" and p.get("stereo_offen"):
                stereo_offen += 1

        return t(
            "erg.zusammenfassung",
            u=zaehler.get("unbekannt", 0),
            s=stereo_offen,
            b=zaehler.get("bekannt", 0),
            g=zaehler.get("ungeprueft", 0)
        )

    # ------------------------------------------------------
    # Quantenchemie (chemlabor/chemie/quantenchemie.py)
    # ------------------------------------------------------

    def _qc_text(self, nummer):

        r = self.qc_ergebnisse.get(nummer)

        if r is None:
            return ""

        if not r["anwendbar"]:
            return "n/a"

        if not r["intakt"]:
            return "zerfällt"

        return f"{r['dhf_kcal']:+.0f}"

    def _qc_detail(self, nummer):

        r = self.qc_ergebnisse.get(nummer)

        if r is None:
            return ""

        if not r["anwendbar"]:
            return f"\n\nQuantenchemie: nicht anwendbar - {r['grund']}"

        text = (f"\n\nQuantenchemie ({r.get('modell')}): Bildungsenthalpie ca. {r['dhf_kcal']:+.0f} kcal/mol "
                f"(Fehler typisch ±{r.get('fehler_kcal') or '10-30'}), {r['dauer_s']} s.")

        if r["intakt"]:
            text += " Die Struktur bleibt nach der Optimierung erhalten."
        else:
            text += (f" ACHTUNG: nach der Optimierung brechen {len(r['gebrochen'])} Bindung(en) auf, "
                     f"{len(r['neu'])} neue entstehen - die gezeichnete Struktur ist keine stabile Geometrie.")

        return text

    def _qc_klick(self):
        """Rechnet die besten sichtbaren Strukturen (oder bricht den Lauf ab)."""

        if self._qc_laeuft:
            self._qc_abbruch.set()
            self.qc_button.config(state="disabled")
            return

        if not self._sichtbar:
            self.status.config(text=t("erg.keine_zeilen"))
            return

        try:
            from chemlabor.chemie import quantenchemie
        except ImportError as fehler:
            self.status.config(text=f"Quantenchemie nicht verfügbar: {fehler}")
            return

        ok, hinweis = quantenchemie.verfuegbar()

        if not ok:
            self.status.config(text=hinweis)
            return

        self._qc_abbruch = threading.Event()
        self._qc_laeuft = True
        self.qc_button.config(text="Abbrechen")
        self.status.config(text="Quantenchemie: starte ...")

        self._qc_faden = threading.Thread(
            target=self._qc_worker, args=(list(self._sichtbar), self._qc_abbruch), daemon=True
        )
        self._qc_faden.start()

    def _qc_worker(self, eintraege, abbruch):

        try:
            from chemlabor.chemie import quantenchemie

            qc = quantenchemie.Quantenchemie()

            gruppen = {}      # Struktur -> Nummern (gleiche Struktur nur einmal rechnen)
            reihenfolge = []

            for e in sorted(eintraege, key=punkte_von, reverse=True):

                nummer, _phase, molekuel, analyse = e

                if not analyse.get("rdkit_gueltig"):
                    continue

                schluessel = neuheit.struktur_schluessel(molekuel)
                key = schluessel["inchikey"] if schluessel else f"nr{nummer}"

                if key not in gruppen:

                    if len(reihenfolge) >= QC_MAX_STRUKTUREN:
                        continue

                    gruppen[key] = []
                    reihenfolge.append((key, molekuel))

                gruppen[key].append(nummer)

            zerfallen = 0

            for i, (key, molekuel) in enumerate(reihenfolge, start=1):

                if abbruch.is_set():
                    break

                schwer = sum(1 for a in molekuel["atoms"] if a["element"] != "H")

                if schwer > QC_MAX_SCHWERE_ATOME:
                    r = {"anwendbar": False, "grund": f"zu gross ({schwer} schwere Atome, Grenze {QC_MAX_SCHWERE_ATOME})"}
                else:
                    r = qc.berechne(molekuel)

                if r["anwendbar"] and not r["intakt"]:
                    zerfallen += 1

                self.qc_queue.put(("ergebnis", gruppen[key], r))
                self.qc_queue.put(("fortschritt", f"Quantenchemie: {i} von {len(reihenfolge)} ..."))

            self.qc_queue.put(("fertig", len(reihenfolge), zerfallen, abbruch.is_set()))

        except Exception as fehler:

            self.qc_queue.put(("fehler", f"{type(fehler).__name__}: {fehler}"))

    def _poll_qc(self):

        neu = False

        while True:

            try:
                nachricht = self.qc_queue.get_nowait()
            except queue.Empty:
                break

            art = nachricht[0]

            if art == "ergebnis":

                for nummer in nachricht[1]:
                    self.qc_ergebnisse[nummer] = nachricht[2]

                neu = True

            elif art == "fortschritt":

                self.status.config(text=nachricht[1])

            elif art in ("fertig", "fehler"):

                self._qc_laeuft = False
                self.qc_button.config(text="Quantenchemie (Top 20)", state="normal")

                if art == "fertig":
                    text = f"Quantenchemie: {nachricht[1]} Strukturen gerechnet, {nachricht[2]} zerfallen bei der Optimierung."
                    if nachricht[3]:
                        text += " (abgebrochen)"
                    self.status.config(text=text)
                    self.log(text)
                else:
                    self.status.config(text=f"Quantenchemie: {nachricht[1]}")

        if neu:
            self.aktualisiere()

    def beenden_und_warten(self, sekunden=3):
        """Bricht eine laufende PubChem-Pruefung ab und wartet kurz auf
        den Thread (vor dem Schliessen des Fensters)."""

        self._neuheit_abbruch.set()
        self._qc_abbruch.set()

        faden = self._neuheit_faden

        if faden is not None and faden.is_alive():
            faden.join(sekunden)

    def destroy(self):

        if self._poll_id is not None:
            self.after_cancel(self._poll_id)
            self._poll_id = None

        super().destroy()

    # ------------------------------------------------------
    # Export
    # ------------------------------------------------------

    def _export_klick(self):
        """Exportiert genau die Zeilen, die die Tabelle gerade zeigt."""

        if not self._sichtbar:

            self.status.config(text=t("erg.nichts_export"))

            return

        pfad = filedialog.asksaveasfilename(
            title=t("erg.export_titel"),
            defaultextension=".json",
            initialfile="ergebnisse.json",
            filetypes=[(t("datei.json"), "*.json"), (t("datei.alle"), "*.*")]
        )

        if not pfad:
            return

        try:

            anzahl = self.exportiere(pfad)

            self.status.config(text=t("erg.exportiert", n=anzahl, pfad=pfad))

        except OSError as fehler:

            self.status.config(text=t("erg.export_fehler", fehler=fehler))

    def exportiere(self, pfad):
        """Die Exportdatei bleibt sprachunabhaengig: interne Schluessel
        und die urspruenglichen Kritiker-Texte, keine Anzeigetexte."""

        daten = []

        for e in self._sichtbar:

            nummer, phase, molekuel, analyse = e

            eintrag = {
                "nr": nummer,
                "phase": phase,
                "punkte": punkte_von(e),
                "urteil": analyse["bewertung"]["urteil"],
                "formel": formel(molekuel),
                "rdkit_gueltig": bool(analyse.get("rdkit_gueltig")),
                "molekuel": molekuel,
            }

            pruefung = self.modell.neuheit.get(nummer)

            if pruefung:
                eintrag["pubchem"] = {
                    k: pruefung.get(k)
                    for k in ("status", "cid", "smiles", "inchikey", "stereo_offen")
                }

            daten.append(eintrag)

        with open(pfad, "w", encoding="utf-8") as datei:
            json.dump(daten, datei, indent=2, ensure_ascii=False)

        return len(daten)