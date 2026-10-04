# ==========================================================
# test_config.py
# ==========================================================
#
# Prueft, dass ChemLabor seine Einstellungen aus LaborConfig bezieht
# (und nicht mehr aus ueberschriebenen Modulkonstanten).
#
# Aufruf:   python tests/test_config.py      (oder: pytest)

import os
import pathlib
import sys
import tempfile
import threading
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from chemlabor import labor as labor_modul
from chemlabor.config import LaborConfig
from chemlabor.labor import ChemLabor, standard_config


KONSTANTEN = (
    "BEOBACHTUNGEN", "AUTOMATIK_EXPERIMENTE", "KANDIDATEN_PRO_MOLEKUEL",
    "ZIEL_VERHAELTNIS_NEG_ZU_POS", "VERSUCHE_FAKTOR", "STABIL_ANZEIGE_SCHWELLE",
)


def klein(**aenderungen):

    werte = dict(
        beobachtungen=60,
        automatik_experimente=15,
        kandidaten_pro_molekuel=2,
        max_atome=6,
        ziel_verhaeltnis_neg_zu_pos=None,
    )

    werte.update(aenderungen)

    return LaborConfig(**werte)


def lauf(config):

    ergebnisse = []

    labor = ChemLabor(
        on_ergebnis=lambda nr, mol, ana, phase: ergebnisse.append((nr, phase, mol, ana)),
        on_fortschritt=lambda text: None,
        config=config,
    )

    labor.starten()

    return labor, ergebnisse


def test_config_wirkt():

    config = klein()

    labor, ergebnisse = lauf(config)

    assert labor.config is config

    beobachtung = [e for e in ergebnisse if e[1] == "Beobachtung"]
    automatik = [e for e in ergebnisse if e[1] == "Automatik"]

    assert len(beobachtung) == config.beobachtungen, len(beobachtung)
    assert len(automatik) == config.automatik_experimente, len(automatik)

    groesste = max(len(e[2]["atoms"]) for e in ergebnisse)

    assert groesste <= config.max_atome, f"{groesste} Atome > max_atome={config.max_atome}"


def test_modulkonstanten_bleiben_unveraendert():

    vorher = {name: getattr(labor_modul, name) for name in KONSTANTEN}

    lauf(klein(beobachtungen=10, automatik_experimente=5, max_atome=4))

    nachher = {name: getattr(labor_modul, name) for name in KONSTANTEN}

    assert vorher == nachher


def test_standard_config_liest_modulkonstanten():

    alt = labor_modul.AUTOMATIK_EXPERIMENTE

    try:
        labor_modul.AUTOMATIK_EXPERIMENTE = alt + 3
        assert standard_config().automatik_experimente == alt + 3
    finally:
        labor_modul.AUTOMATIK_EXPERIMENTE = alt


def test_validierung():

    kaputte = (
        dict(beobachtungen=0),
        dict(automatik_experimente=-1),
        dict(max_atome=1),
        dict(kandidaten_pro_molekuel=0),
        dict(automatik_modus="xyz"),
        dict(ziel_verhaeltnis_neg_zu_pos=-1),
        dict(beobachtungen="12"),
    )

    for werte in kaputte:

        try:
            LaborConfig(**werte).validiert()
        except ValueError:
            continue

        raise AssertionError(f"{werte} haette abgelehnt werden muessen")

    assert LaborConfig(ziel_verhaeltnis_neg_zu_pos=None).validiert() is not None


def test_stoppen_im_endlosmodus():

    config = klein(endlos=True)

    phasen = []

    labor = ChemLabor(
        on_ergebnis=lambda nr, mol, ana, phase: phasen.append(phase),
        on_fortschritt=lambda text: None,
        config=config,
    )

    faden = threading.Thread(target=labor.starten, daemon=True)
    faden.start()

    ende = time.time() + 20

    while time.time() < ende and phasen.count("Automatik") < 5:
        time.sleep(0.05)

    labor.stoppen()
    faden.join(10)

    assert not faden.is_alive(), "Labor hat nach stoppen() nicht aufgehoert"
    assert phasen.count("Automatik") >= 5


def test_konstruktionsmodus():

    config = klein(automatik_modus="konstruktion")

    labor, ergebnisse = lauf(config)

    assert len(ergebnisse) >= config.beobachtungen


def test_autosave_aus_config():

    with tempfile.TemporaryDirectory() as ordner:

        pfad = os.path.join(ordner, "auto.json")

        config = klein(
            beobachtungen=30,
            automatik_experimente=10,
            autosave_pfad=pfad,
            autosave_intervall=10,
        )

        lauf(config)

        assert os.path.exists(pfad), "Autosave-Datei wurde nicht geschrieben"


if __name__ == "__main__":

    fehler = 0

    for name, funktion in sorted(globals().items()):

        if not name.startswith("test_"):
            continue

        try:
            funktion()
            print(f"ok       {name}")
        except Exception as problem:
            fehler += 1
            print(f"FEHLER   {name}: {type(problem).__name__}: {problem}")

    print("\nAlle Tests bestanden." if not fehler else f"\n{fehler} FEHLER")

    sys.exit(1 if fehler else 0)
