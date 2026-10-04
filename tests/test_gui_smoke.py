# ==========================================================
# test_gui_smoke.py
# ==========================================================
#
# Rauchtest der neuen Oberflaeche: baut das Fenster, startet einen
# kleinen Lauf mit dem ECHTEN Labor, schaltet durch die Tabs, waehlt eine
# Zeile, filtert, sortiert, exportiert, prueft offline gegen den
# PubChem-Cache und probiert den Einstellungsdialog.
#
# Braucht ein Display (laeuft auf deinem Rechner normal; ohne Display
# meldet der Test "uebersprungen").
#
# Aufruf:   python tests/test_gui_smoke.py

import json
import os
import pathlib
import sys
import tempfile
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import tkinter as tk

from chemlabor.config import LaborConfig
from chemlabor.gui.dialog_einstellungen import EinstellungenDialog, zeige_einstellungen
from chemlabor.gui.fenster import Hauptfenster


def pumpe(root, bedingung, sekunden=90):
    """Haelt die Oberflaeche am Laufen, bis die Bedingung erfuellt ist."""

    ende = time.time() + sekunden

    while time.time() < ende:

        root.update()

        if bedingung():
            return True

        time.sleep(0.02)

    return False


def ruhe(root, sekunden=0.5):

    pumpe(root, lambda: False, sekunden)


def main():

    try:
        root = tk.Tk()
    except tk.TclError as fehler:
        print(f"uebersprungen: kein Display ({fehler})")
        return 0

    fehler = []

    def pruefe(bedingung, text):

        print(("ok       " if bedingung else "FEHLER   ") + text)

        if not bedingung:
            fehler.append(text)

    klein = LaborConfig(
        beobachtungen=40,
        automatik_experimente=15,
        kandidaten_pro_molekuel=2,
        max_atome=6,
        ziel_verhaeltnis_neg_zu_pos=None,
    )

    fenster = Hauptfenster(root, config=klein)

    tabs = [fenster.notebook.tab(t, "text") for t in fenster.notebook.tabs()]

    pruefe(tabs == ["Live", "Ergebnisse"], f"zwei Tabs: {tabs}")

    # --- Lauf ----------------------------------------------------------

    fenster.start()

    pruefe(fenster.laeuft and str(fenster.einstellungen_button["state"]) == "disabled",
           "Start sperrt die Einstellungen")

    fertig = pumpe(root, lambda: not fenster.laeuft)

    ruhe(root, 0.4)

    pruefe(fertig, "Lauf endet von selbst")
    pruefe(fenster.modell.zaehler == 55, f"55 Ergebnisse (40 + 15): {fenster.modell.zaehler}")
    pruefe(fenster.status_label.cget("text") == "Fertig.", "Status 'Fertig.'")

    # --- Live-Tab ----------------------------------------------------------

    fenster.notebook.select(fenster.tab_live)
    ruhe(root, 0.3)

    pruefe(fenster.tab_live.werte["experimente"].cget("text").startswith("55"),
           "Live: Experimente = " + fenster.tab_live.werte["experimente"].cget("text"))

    # --- Ergebnis-Tab ------------------------------------------------------------

    fenster.notebook.select(fenster.tab_ergebnisse)
    ruhe(root, 0.3)

    tab = fenster.tab_ergebnisse

    zeilen = tab.tabelle.get_children()

    pruefe(len(zeilen) == 55, f"Tabelle zeigt 55 Zeilen: {len(zeilen)}")

    tab.tabelle.selection_set(zeilen[0])
    ruhe(root, 0.8)

    detail = tab.detail.get("1.0", "end")

    pruefe("Experiment" in detail and "Punkte" in detail, "Detailtext erscheint nach Klick")
    pruefe(len(tab.lewis.lewis_canvas.find_all()) > 0, "Lewis-Ansicht hat etwas gezeichnet")

    tab.filter_var.set("Top 10")
    tab.aktualisiere()

    pruefe(len(tab.tabelle.get_children()) == 10, "Filter 'Top 10' zeigt 10 Zeilen")

    tab.filter_var.set("Alle")
    tab.phase_var.set("Automatik")
    tab.aktualisiere()

    pruefe(len(tab.tabelle.get_children()) == 15, "Phase 'Automatik' zeigt 15 Zeilen")

    tab.phase_var.set("Alle")
    tab._sortiere("punkte")
    tab._sortiere("punkte")
    tab.aktualisiere()

    punkte = [int(tab.tabelle.set(z, "punkte")) for z in tab.tabelle.get_children()]

    pruefe(punkte == sorted(punkte), "Sortierung nach Punkten (zweiter Klick: aufsteigend)")

    with tempfile.TemporaryDirectory() as ordner:

        pfad = os.path.join(ordner, "ergebnisse.json")

        anzahl = tab.exportiere(pfad)

        with open(pfad, encoding="utf-8") as datei:
            daten = json.load(datei)

        pruefe(anzahl == len(daten) == 55 and "molekuel" in daten[0], "Export schreibt alle sichtbaren Zeilen")

    # --- PubChem (offline: nur Cache, kein Netz) ---------------------------------

    tab.offline_var.set(True)
    tab._neuheit_klick()

    pruefe(pumpe(root, lambda: not tab._neuheit_laeuft, 60), "PubChem-Pruefung (offline) endet")

    ruhe(root, 0.3)

    status = {p["status"] for p in fenster.modell.neuheit.values()}

    pruefe(bool(status) and status <= {"bekannt", "unbekannt", "ungeprueft"}, f"PubChem-Status gesetzt: {status}")

    tab.neuheit_filter_var.set("Nur ungeprüft")
    tab.aktualisiere()

    pruefe(len(tab.tabelle.get_children()) == sum(
        1 for p in fenster.modell.neuheit.values() if p["status"] == "ungeprueft"
    ), "PubChem-Anzeigefilter passt zur Zahl der Status")

    # --- Einstellungsdialog (ueber den echten Weg zeige_einstellungen) ---------------------

    def im_dialog(aktion):
        """Fuehrt 'aktion(dialog)' kurz nach dem Oeffnen des modalen Dialogs aus."""

        def lauf():
            dialog = next(w for w in root.winfo_children() if isinstance(w, EinstellungenDialog))
            aktion(dialog)

        root.after(200, lauf)

    def gueltig(dialog):
        dialog.beobachtungen.set("1.500")
        dialog.max_atome.set("9")
        dialog.filter_aktiv.set(True)
        dialog.verhaeltnis.set("4")
        dialog.modus.set("konstruktion")
        dialog._uebernehmen()

    im_dialog(gueltig)

    neu = zeige_einstellungen(root, klein)

    pruefe(
        neu is not None and neu.beobachtungen == 1500 and neu.max_atome == 9
        and neu.ziel_verhaeltnis_neg_zu_pos == 4 and neu.automatik_modus == "konstruktion",
        "Dialog uebernimmt gueltige Werte"
    )

    meldung = []

    def ungueltig(dialog):
        dialog.max_atome.set("1")
        dialog._uebernehmen()
        meldung.append(dialog.fehler.cget("text"))
        dialog.destroy()

    im_dialog(ungueltig)

    pruefe(zeige_einstellungen(root, klein) is None and meldung and meldung[0] != "",
           "Dialog lehnt ungueltigen Wert mit Meldung ab")

    # --- Stoppen im Endlosmodus ------------------------------------------------------------

    fenster.config = klein.kopie(endlos=True)

    fenster.start()

    pumpe(root, lambda: fenster.modell.zaehler > 45, 60)

    fenster.stopp()

    pruefe(pumpe(root, lambda: not fenster.laeuft, 30), "Stopp beendet den Endlos-Lauf")

    fenster.schliessen()

    print("\nAlle Tests bestanden." if not fehler else f"\n{len(fehler)} FEHLER")

    return 1 if fehler else 0


if __name__ == "__main__":
    sys.exit(main())
