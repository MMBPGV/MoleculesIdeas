# ==========================================================
# fenster.py
# ==========================================================
#
# Hauptfenster der neuen, schlankeren Oberflaeche:
#
#   Kopfzeile:  Start / Stopp / Einstellungen / Gedaechtnis-Menue
#   Tabs:       Live (Kennzahlen + Lernkurve), Ergebnisse (Tabelle,
#               Lewis, Detail, PubChem)
#   unten:      einklappbares Log
#
# Einstellungen stehen in einer LaborConfig (siehe config.py) und
# werden ChemLabor uebergeben - nichts wird mehr ueber Modulkonstanten
# gesetzt. Die bisherige GUI bleibt in app.py erhalten (start_gui_alt.py).
#
# Alle sichtbaren Texte kommen aus chemlabor/texte.py.

import queue
import threading
import time
import tkinter as tk
from dataclasses import replace
from tkinter import filedialog, scrolledtext, ttk

from chemlabor.config import LaborConfig
from chemlabor.gui.dialog_einstellungen import zeige_einstellungen
from chemlabor.gui.modell import ErgebnisModell, config_text
from chemlabor.gui.tab_ergebnisse import TabErgebnisse
from chemlabor.gui.tab_live import TabLive
from chemlabor.gui.theme import FARBEN, ergaenze_style, konfiguriere_style
from chemlabor.ki.gedaechtnis import KIGedaechtnis
from chemlabor.labor import ChemLabor
from chemlabor.texte import t


LOG_MAX_ZEILEN = 2000

# Die Ergebnis-Tabelle wird hoechstens so oft neu aufgebaut (Sekunden)
TABELLE_MIN_ABSTAND = 0.4

FERTIG = "__FERTIG__"


def _json_dateitypen():

    return [(t("datei.json"), "*.json"), (t("datei.alle"), "*.*")]


class Hauptfenster:

    def __init__(self, root, config=None):

        self.root = root

        root.title(t("fenster.titel"))
        root.geometry("1200x780")
        root.minsize(920, 600)

        konfiguriere_style(root)
        ergaenze_style(root)

        self.config = config if config is not None else LaborConfig()

        self.labor = None
        self.laeuft = False
        self.geladenes_gedaechtnis = None
        self._labor_faden = None

        self.modell = ErgebnisModell()

        self.ergebnis_queue = queue.Queue()
        self.fortschritt_queue = queue.Queue()

        self._stand_live = -1
        self._stand_tabelle = -1
        self._tabelle_zeit = 0.0

        self._baue_kopfzeile()
        self._baue_tabs()
        self._baue_log()

        self._aktualisiere_beschriftungen()

        root.protocol("WM_DELETE_WINDOW", self.schliessen)

        self._poll_queues()

    # ------------------------------------------------------
    # Aufbau
    # ------------------------------------------------------

    def _baue_kopfzeile(self):

        kopf = ttk.Frame(self.root, padding=(10, 8, 10, 0))
        kopf.pack(fill="x")

        self.start_button = ttk.Button(kopf, text=t("fenster.start"), style="Accent.TButton", command=self.start)
        self.start_button.pack(side="left", padx=(0, 6))

        self.stop_button = ttk.Button(kopf, text=t("fenster.stopp"), style="Danger.TButton", command=self.stopp, state="disabled")
        self.stop_button.pack(side="left", padx=(0, 16))

        self.einstellungen_button = ttk.Button(kopf, text=t("fenster.einstellungen"), command=self.einstellungen_oeffnen)
        self.einstellungen_button.pack(side="left", padx=(0, 6))

        self.gedaechtnis_button = ttk.Menubutton(kopf, text=t("fenster.gedaechtnis"))
        self.gedaechtnis_button.pack(side="left")

        menue = tk.Menu(
            self.gedaechtnis_button, tearoff=0,
            bg=FARBEN["panel_hell"], fg=FARBEN["text"],
            activebackground=FARBEN["auswahl"], activeforeground=FARBEN["text"]
        )

        menue.add_command(label=t("fenster.menue_speichern"), command=self.gedaechtnis_speichern)
        menue.add_command(label=t("fenster.menue_laden"), command=self.gedaechtnis_laden)
        menue.add_command(label=t("fenster.menue_zusammenfuehren"), command=self.gedaechtnis_zusammenfuehren)
        menue.add_separator()
        menue.add_command(label=t("fenster.menue_verwerfen"), command=self.gedaechtnis_verwerfen)

        self.gedaechtnis_button["menu"] = menue

        self.status_label = ttk.Label(kopf, text=t("fenster.bereit"))
        self.status_label.pack(side="right")

        info = ttk.Frame(self.root, padding=(12, 4, 10, 6))
        info.pack(fill="x")

        self.config_label = ttk.Label(info, text="", style="Gedimmt.TLabel")
        self.config_label.pack(side="left")

        self.gedaechtnis_label = ttk.Label(info, text="", style="Gedimmt.TLabel")
        self.gedaechtnis_label.pack(side="right")

    def _baue_tabs(self):

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=6, pady=(0, 4))

        self.tab_live = TabLive(self.notebook, self.modell)
        self.notebook.add(self.tab_live, text=t("fenster.tab_live"))

        self.tab_ergebnisse = TabErgebnisse(self.notebook, self.modell, log=self.log)
        self.notebook.add(self.tab_ergebnisse, text=t("fenster.tab_ergebnisse"))

        self.notebook.bind("<<NotebookTabChanged>>", lambda _e: self._render(erzwingen=True))

    def _baue_log(self):

        self.log_rahmen = ttk.Frame(self.root, padding=(8, 0, 8, 6))
        self.log_rahmen.pack(fill="x")

        leiste = ttk.Frame(self.log_rahmen)
        leiste.pack(fill="x")

        self.log_sichtbar = True

        self.log_toggle = ttk.Button(leiste, text=t("fenster.log_offen"), command=self._log_umschalten, width=8)
        self.log_toggle.pack(side="left")

        ttk.Button(leiste, text=t("fenster.log_leeren"), command=self._log_leeren, width=8).pack(side="left", padx=6)

        self.log_text = scrolledtext.ScrolledText(
            self.log_rahmen, wrap="word", height=6, state="disabled",
            bg=FARBEN["panel_hell"], fg=FARBEN["text"],
            insertbackground=FARBEN["text"], borderwidth=0, highlightthickness=0
        )

        self.log_text.pack(fill="x", pady=(4, 0))

    # ------------------------------------------------------
    # Beschriftungen
    # ------------------------------------------------------

    def _aktualisiere_beschriftungen(self):

        self.config_label.config(text=config_text(self.config))

        if self.geladenes_gedaechtnis is not None:
            text = t("fenster.gedaechtnis_geladen", n=self.geladenes_gedaechtnis.beobachtungen)
        else:
            text = t("fenster.gedaechtnis_neu")

        self.gedaechtnis_label.config(text=text)

    def _setze_status(self, text):

        self.status_label.config(text=text)

    # ------------------------------------------------------
    # Steuerung
    # ------------------------------------------------------

    def einstellungen_oeffnen(self):

        if self.laeuft:
            return

        neu = zeige_einstellungen(self.root, self.config)

        if neu is not None:

            self.config = neu

            self._aktualisiere_beschriftungen()

    def start(self):

        if self.laeuft:
            return

        try:
            lauf_config = replace(self.config).validiert()   # Momentaufnahme
        except ValueError as fehler:
            self._setze_status(str(fehler))
            return

        try:
            labor = ChemLabor(
                on_ergebnis=self._on_ergebnis,
                on_fortschritt=self._on_fortschritt,
                config=lauf_config
            )
        except Exception as fehler:
            self._setze_status(t("fenster.labor_fehler_start", fehler=fehler))
            return

        self._log_leeren()
        self.modell.leeren()
        self.tab_live.leeren()
        self.tab_ergebnisse.leeren()

        self._stand_live = -1
        self._stand_tabelle = -1

        self.labor = labor

        if self.geladenes_gedaechtnis is not None:

            labor.gedaechtnis = self.geladenes_gedaechtnis

            self.log(
                t("fenster.gedaechtnis_uebernommen", n=self.geladenes_gedaechtnis.beobachtungen)
            )

        if lauf_config.autosave_pfad:

            self.log(
                t(
                    "fenster.autosave_aktiv",
                    n=lauf_config.autosave_intervall,
                    pfad=lauf_config.autosave_pfad
                )
            )

        self.laeuft = True

        self.start_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.einstellungen_button.config(state="disabled")

        self._setze_status(
            t("fenster.laeuft_endlos") if lauf_config.endlos else t("fenster.laeuft")
        )

        self._labor_faden = threading.Thread(target=self._labor_thread, daemon=True)
        self._labor_faden.start()

    def _labor_thread(self):

        try:
            self.labor.starten()
        except Exception as fehler:
            self.fortschritt_queue.put(t("fenster.fehler_praefix", fehler=fehler))
        finally:
            self.fortschritt_queue.put(FERTIG)

    def stopp(self):

        if self.labor is not None:
            self.labor.stoppen()

        self.stop_button.config(state="disabled")
        self._setze_status(t("fenster.stoppe"))

    def schliessen(self):
        """Stoppt Labor und PubChem-Pruefung und wartet kurz auf deren
        Threads, bevor das Fenster zerstoert wird - sonst raeumen sie
        evtl. noch Tk-Variablen aus einem Hintergrundthread ab und
        erzeugen beim Beenden Fehlermeldungen in der Konsole."""

        if self.laeuft and self.labor is not None:
            self.labor.stoppen()

        self.tab_ergebnisse.beenden_und_warten()

        if self._labor_faden is not None and self._labor_faden.is_alive():
            self._labor_faden.join(3)

        self.root.destroy()

    # ------------------------------------------------------
    # Callbacks aus dem Labor-Thread (nur Queue, nie direkt Widgets)
    # ------------------------------------------------------

    def _on_ergebnis(self, nummer, molekuel, analyse, phase):

        self.ergebnis_queue.put((nummer, phase, molekuel, analyse))

    def _on_fortschritt(self, text):

        self.fortschritt_queue.put(text)

    def _poll_queues(self):

        fertig = False

        while True:

            try:
                text = self.fortschritt_queue.get_nowait()
            except queue.Empty:
                break

            if text == FERTIG:
                fertig = True
            else:
                self.log(text)

        while True:

            try:
                eintrag = self.ergebnis_queue.get_nowait()
            except queue.Empty:
                break

            self.modell.hinzufuegen(eintrag)

        if fertig:

            self.laeuft = False

            self.start_button.config(state="normal")
            self.stop_button.config(state="disabled")
            self.einstellungen_button.config(state="normal")

            self._setze_status(t("fenster.fertig"))

            self._render(erzwingen=True)

        else:

            self._render()

        self.root.after(150, self._poll_queues)

    def _render(self, erzwingen=False):
        """Zeichnet nur den gerade sichtbaren Tab - und nur, wenn es neue
        Ergebnisse gibt (die Tabelle zusaetzlich hoechstens alle
        TABELLE_MIN_ABSTAND Sekunden)."""

        if not self.modell.eintraege:
            return

        anzahl = self.modell.zaehler
        aktueller = self.notebook.select()

        if aktueller == str(self.tab_live):

            if erzwingen or anzahl != self._stand_live:
                self.tab_live.aktualisiere(self.labor)
                self._stand_live = anzahl

        elif aktueller == str(self.tab_ergebnisse):

            jetzt = time.monotonic()

            if erzwingen or (
                anzahl != self._stand_tabelle
                and jetzt - self._tabelle_zeit >= TABELLE_MIN_ABSTAND
            ):
                self.tab_ergebnisse.aktualisiere()
                self._stand_tabelle = anzahl
                self._tabelle_zeit = jetzt

    # ------------------------------------------------------
    # Gedaechtnis
    # ------------------------------------------------------

    def gedaechtnis_speichern(self):

        if self.labor is None:
            self._setze_status(t("fenster.kein_gedaechtnis"))
            return

        pfad = filedialog.asksaveasfilename(
            title=t("fenster.speichern_titel"), defaultextension=".json",
            filetypes=_json_dateitypen()
        )

        if not pfad:
            return

        try:
            self.labor.speichere_gedaechtnis(pfad)
            self._setze_status(t("fenster.gespeichert", pfad=pfad))
        except OSError as fehler:
            self._setze_status(t("fenster.speichern_fehler", fehler=fehler))

    def gedaechtnis_laden(self):

        if self.laeuft:
            self._setze_status(t("fenster.laden_erst_stoppen"))
            return

        pfad = filedialog.askopenfilename(
            title=t("fenster.laden_titel"),
            filetypes=_json_dateitypen()
        )

        if not pfad:
            return

        try:
            self.geladenes_gedaechtnis = KIGedaechtnis.laden(pfad)
        except (OSError, ValueError, KeyError) as fehler:
            self._setze_status(t("fenster.laden_fehler", fehler=fehler))
            return

        self._aktualisiere_beschriftungen()

        self._setze_status(t("fenster.geladen"))

    def gedaechtnis_zusammenfuehren(self):

        if self.laeuft:
            self._setze_status(t("fenster.merge_erst_stoppen"))
            return

        pfade = filedialog.askopenfilenames(
            title=t("fenster.merge_titel"),
            filetypes=_json_dateitypen()
        )

        if not pfade:
            return

        if len(pfade) < 2:
            self._setze_status(t("fenster.merge_zu_wenig"))
            return

        try:
            kollektiv = KIGedaechtnis.zusammenfuehren(list(pfade))
        except (OSError, ValueError, KeyError) as fehler:
            self._setze_status(t("fenster.merge_fehler", fehler=fehler))
            return

        self.geladenes_gedaechtnis = kollektiv

        self._aktualisiere_beschriftungen()

        speicherpfad = filedialog.asksaveasfilename(
            title=t("fenster.merge_speichern_titel"), defaultextension=".json",
            filetypes=_json_dateitypen()
        )

        if speicherpfad:

            try:
                kollektiv.speichern(speicherpfad)
            except OSError as fehler:
                self._setze_status(t("fenster.merge_speichern_fehler", fehler=fehler))
                return

        self._setze_status(
            t("fenster.merge_fertig", n=len(pfade), b=kollektiv.beobachtungen)
        )

    def gedaechtnis_verwerfen(self):

        if self.laeuft:
            self._setze_status(t("fenster.verwerfen_erst_stoppen"))
            return

        self.geladenes_gedaechtnis = None

        self._aktualisiere_beschriftungen()

        self._setze_status(t("fenster.verworfen"))

    # ------------------------------------------------------
    # Log
    # ------------------------------------------------------

    def log(self, text):

        self.log_text.config(state="normal")
        self.log_text.insert("end", text + "\n")

        zeilen = int(self.log_text.index("end-1c").split(".")[0])

        if zeilen > LOG_MAX_ZEILEN:
            self.log_text.delete("1.0", f"{zeilen - LOG_MAX_ZEILEN}.0")

        self.log_text.see("end")
        self.log_text.config(state="disabled")

    def _log_leeren(self):

        self.log_text.config(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.config(state="disabled")

    def _log_umschalten(self):

        if self.log_sichtbar:
            self.log_text.pack_forget()
            self.log_toggle.config(text=t("fenster.log_zu"))
        else:
            self.log_text.pack(fill="x", pady=(4, 0))
            self.log_toggle.config(text=t("fenster.log_offen"))

        self.log_sichtbar = not self.log_sichtbar


def main():

    root = tk.Tk()

    Hauptfenster(root)

    root.mainloop()


if __name__ == "__main__":
    main()