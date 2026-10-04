# ==========================================================
# pfade.py
# ==========================================================
#
# EINZIGE Stelle, an der Dateipfade des Projekts festgelegt werden.
# Alles ist relativ zum Projektordner (nicht zum Arbeitsverzeichnis):
# das Programm laeuft damit egal, von wo es gestartet wird.

from pathlib import Path

PROJEKT_WURZEL = Path(__file__).resolve().parent.parent

DATEN = PROJEKT_WURZEL / "daten"
GEDAECHTNIS = DATEN / "gedaechtnis"

PERIODENSYSTEM = str(DATEN / "periodensystem.json")
AUTOSAVE_PFAD = str(GEDAECHTNIS / "autosave_gedaechtnis.json")
NEUHEIT_CACHE_PFAD = str(DATEN / "neuheit_cache.sqlite")

GEDAECHTNIS.mkdir(parents=True, exist_ok=True)
