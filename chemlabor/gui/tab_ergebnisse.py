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


SPALTEN = ("nr", "phase", "formel", "atome", "punkte", "urteil", "rdkit", "pubchem")

BREITE = {
    "nr": 50, "phase": 90, "formel": 110, "atome": 55,
    "punkte": 60, "urteil": 140, "rdkit": 55, "pubchem": 90,
}

# So viele verschiedene Strukturen werden pro Klick bei PubChem abgefragt
# (jede neue Struktur kostet ~0,25 s wegen des PubChem-Limits).
NEUHEIT_MAX_STRUKTUREN = 300


def beschriftung(spalte):
    """Spaltenueberschrift in der aktiven Sprache."""

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
            state="readonly", width=20
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
                command=lambda s=spalte: self._sortiere(s)
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

        self.lewis = LewisAnsicht(rechts)
        rechts.add(self.lewis, weight=3)

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

        self._setze_detail(baue_detailtext(eintrag, self.modell.neuheit.get(eintrag[0])))

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

    def beenden_und_warten(self, sekunden=3):
        """Bricht eine laufende PubChem-Pruefung ab und wartet kurz auf
        den Thread (vor dem Schliessen des Fensters)."""

        self._neuheit_abbruch.set()

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