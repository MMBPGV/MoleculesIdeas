# ==========================================================
# test_gui_modell.py
# ==========================================================
#
# Prueft die Daten- und Filterlogik der Oberflaeche (gui/modell.py)
# OHNE Fenster - laeuft also auch ohne Display.
#
# Aufruf:   python tests/test_gui_modell.py      (oder: pytest)

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from chemlabor.config import LaborConfig
from chemlabor.gui.modell import (
    ErgebnisModell,
    baue_detailtext,
    config_text,
    punkte_von,
)


def molekuel(*elemente):
    """Kette aus den angegebenen Elementen mit Einfachbindungen."""

    valenz = {"H": 1, "C": 4, "O": 2, "N": 3, "Cl": 1}

    benutzt = [0] * len(elemente)
    bindungen = []

    for i in range(len(elemente) - 1):
        bindungen.append({"atom1": i, "atom2": i + 1, "order": 1})
        benutzt[i] += 1
        benutzt[i + 1] += 1

    return {
        "atoms": [
            {"id": i, "element": e, "valence": valenz[e], "used_valence": benutzt[i]}
            for i, e in enumerate(elemente)
        ],
        "bonds": bindungen,
    }


def analyse(punkte, rdkit=True):

    return {
        "bewertung": {"punkte": punkte, "urteil": "x"},
        "rdkit_gueltig": rdkit,
        "probleme": [{"modell": "Radikal", "atom": 0}],
        "erklaerung": ["Beispiel"],
    }


def beispiel_modell():

    modell = ErgebnisModell()

    daten = [
        (1, "Beobachtung", molekuel("C", "H"), analyse(95)),
        (2, "Beobachtung", molekuel("C", "H"), analyse(40)),      # gleiche Struktur wie 1
        (3, "Automatik", molekuel("O", "H"), analyse(80)),
        (4, "Automatik", molekuel("N", "H"), analyse(60)),
        (5, "Automatik", molekuel("C", "Cl"), analyse(30, rdkit=False)),
    ]

    for eintrag in daten:
        modell.hinzufuegen(eintrag)

    return modell


def test_zaehler_und_leeren():

    modell = beispiel_modell()

    assert modell.zaehler == 5

    modell.leeren()

    assert modell.zaehler == 0 and not modell.eintraege and not modell.neuheit


def test_filter_phase_und_stabilitaet():

    modell = beispiel_modell()

    assert [e[0] for e in modell.filtere(phase="Automatik")] == [3, 4, 5]
    assert {e[0] for e in modell.filtere(filter_wahl="Nur stabil (≥75)")} == {1, 3}
    assert {e[0] for e in modell.filtere(filter_wahl="Nur instabil (<50)")} == {2, 5}


def test_top_n_sortiert_nach_punkten():

    modell = beispiel_modell()

    assert [e[0] for e in modell.filtere(filter_wahl="Top 10")] == [1, 3, 4, 2, 5]


def test_wiederholungen_ausblenden():

    modell = beispiel_modell()

    nummern = {e[0] for e in modell.filtere(duplikate_ausblenden=True)}

    assert not ({1, 2} <= nummern), "Beide gleichen Strukturen sind noch da"
    assert 1 in nummern or 2 in nummern


def test_sortierung_nach_spalte():

    modell = beispiel_modell()

    auf = [e[0] for e in modell.filtere(sortierung=("nr", False))]
    ab = [e[0] for e in modell.filtere(sortierung=("nr", True))]

    assert auf == [1, 2, 3, 4, 5] and ab == [5, 4, 3, 2, 1]


def test_pubchem_filter_und_text():

    modell = beispiel_modell()

    modell.neuheit[3] = {"status": "unbekannt", "cid": None}
    modell.neuheit[4] = {"status": "bekannt", "cid": 123}

    assert modell.pubchem_text(3) == "unbekannt"
    assert modell.pubchem_text(4) == "bekannt"
    assert modell.pubchem_text(1) == "–"

    assert [e[0] for e in modell.filtere(neuheit_filter="Nur unbekannt")] == [3]
    assert [e[0] for e in modell.filtere(neuheit_filter="Nur bekannt")] == [4]

    nach_pubchem = [e[0] for e in modell.filtere(sortierung=("pubchem", False))]

    assert len(nach_pubchem) == 5


def test_ringpuffer_begrenzt_eintraege_nicht_den_zaehler():

    modell = ErgebnisModell(max_eintraege=3)

    for i in range(10):
        modell.hinzufuegen((i, "Automatik", molekuel("C", "H"), analyse(50)))

    assert len(modell.eintraege) == 3 and modell.zaehler == 10


def test_kennzahlen():

    modell = beispiel_modell()

    k = modell.kennzahlen()

    assert k["gesamt"] == 5
    assert k["minimum"] == 30 and k["maximum"] == 95
    assert k["stabil"] == 2
    assert k["beobachtung"] == 2 and k["automatik"] == 3
    assert abs(k["mittel"] - (95 + 40 + 80 + 60 + 30) / 5) < 1e-9

    assert ErgebnisModell().kennzahlen() is None


def test_detailtext_enthaelt_wichtiges():

    modell = beispiel_modell()

    text = baue_detailtext(modell.eintrag(3), {
        "status": "unbekannt", "cid": None, "smiles": "O", "inchikey": "KEY", "stereo_offen": True
    })

    for teil in ("Experiment 3", "Punkte: 80", "PubChem: unbekannt", "Stereo offen", "Probleme:", "Erklärung:"):
        assert teil in text, teil

    assert baue_detailtext(modell.eintrag(1)).count("PubChem") == 0


def test_config_text():

    text = config_text(LaborConfig(beobachtungen=5000, max_atome=8, ziel_verhaeltnis_neg_zu_pos=None, endlos=True))

    assert "5.000 Beobachtungen" in text
    assert "max. 8 Atome" in text
    assert "Filter aus" in text
    assert "endlos" in text


def test_punkte_von():

    assert punkte_von((1, "Automatik", molekuel("C", "H"), analyse(77))) == 77


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
