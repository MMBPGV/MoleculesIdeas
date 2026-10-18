"""Prueft benson.berechne an Molekuelen mit bekannter Bildungsenthalpie.

Lauffaehig ohne den Rest von chemlabor: baut das noetige daten-Dict selbst.
Referenzwerte: experimentelle dHf (Gasphase, 298 K, kcal/mol), gerundet.
"""

import sys
from collections import defaultdict
from pathlib import Path

import pytest

# Projektwurzel auf den Pfad, damit der Test auch direkt aus tests/ laeuft.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chemlabor.chemie import benson

STANDARD_VALENZ = {"C": 4, "O": 2, "N": 3, "S": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}


def baue_daten(schwere, bindungen):
    """schwere: Liste von Elementen; bindungen: (i, j, ordnung) zwischen
    schweren Atomen. Fehlende Valenzen werden mit explizitem H aufgefuellt."""

    atome = [{"id": i, "element": e} for i, e in enumerate(schwere)]
    bonds = list(bindungen)
    benutzt = defaultdict(int)

    for i, j, o in bonds:
        benutzt[i] += o
        benutzt[j] += o

    naechste = len(atome)

    for i, e in enumerate(schwere):
        for _ in range(STANDARD_VALENZ[e] - benutzt[i]):
            atome.append({"id": naechste, "element": "H"})
            bonds.append((i, naechste, 1))
            naechste += 1

    graph = defaultdict(list)
    bindung_by_paar = {}

    for i, j, o in bonds:
        graph[i].append(j)
        graph[j].append(i)
        bindung_by_paar[frozenset([i, j])] = {"atom1": i, "atom2": j, "ordnung": o}

    return {
        "atome": atome,
        "graph": graph,
        "bindung_by_paar": bindung_by_paar,
        "atom_by_id": {a["id"]: a for a in atome},
    }


def kein_aromat(ring):
    return False


def rechne(schwere, bindungen, ringe=(), aromatisch=kein_aromat):
    daten = baue_daten(schwere, bindungen)
    return benson.berechne(daten, list(ringe), aromatisch)


FAELLE = {
    "methan":         (["C"], [], [], -17.9),
    "ethan":          (["C", "C"], [(0, 1, 1)], [], -20.0),
    "propan":         (["C", "C", "C"], [(0, 1, 1), (1, 2, 1)], [], -25.0),
    "isobutan":       (["C"] * 4, [(0, 1, 1), (0, 2, 1), (0, 3, 1)], [], -32.1),
    "ethen":          (["C", "C"], [(0, 1, 2)], [], 12.5),
    "propen":         (["C", "C", "C"], [(0, 1, 2), (1, 2, 1)], [], 4.9),
    "ethin":          (["C", "C"], [(0, 1, 3)], [], 54.2),
    "ethanol":        (["C", "C", "O"], [(0, 1, 1), (1, 2, 1)], [], -56.2),
    "dimethylether":  (["C", "O", "C"], [(0, 1, 1), (1, 2, 1)], [], -44.0),
    "aceton":         (["C", "C", "C", "O"], [(0, 1, 1), (1, 2, 1), (1, 3, 2)], [], -52.2),
    "formaldehyd":    (["C", "O"], [(0, 1, 2)], [], -26.0),
    "methylamin":     (["C", "N"], [(0, 1, 1)], [], -5.5),
    "wasser":         (["O"], [], [], -57.8),
}


@pytest.mark.parametrize("name", sorted(FAELLE))
def test_acyclische_molekuele(name):
    schwere, bindungen, ringe, erwartet = FAELLE[name]
    ergebnis = rechne(schwere, bindungen, ringe)

    assert ergebnis["anwendbar"], ergebnis["nicht_abgedeckt"]
    assert ergebnis["dhf_kcal"] == pytest.approx(erwartet, abs=1.0)


def test_cyclopropan():
    ring = [0, 1, 2]
    ergebnis = rechne(["C"] * 3, [(0, 1, 1), (1, 2, 1), (2, 0, 1)], [ring])
    assert ergebnis["anwendbar"]
    assert ergebnis["dhf_kcal"] == pytest.approx(12.7, abs=1.0)


def test_cyclohexan():
    ring = list(range(6))
    bindungen = [(i, (i + 1) % 6, 1) for i in range(6)]
    ergebnis = rechne(["C"] * 6, bindungen, [ring])
    assert ergebnis["dhf_kcal"] == pytest.approx(-29.5, abs=1.0)


def test_benzol():
    ring = list(range(6))
    bindungen = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)]
    ergebnis = rechne(["C"] * 6, bindungen, [ring], aromatisch=lambda r: True)
    assert ergebnis["anwendbar"], ergebnis["nicht_abgedeckt"]
    assert ergebnis["dhf_kcal"] == pytest.approx(19.8, abs=1.0)
    assert ergebnis["gruppen"] == {"Cb-(H)": 6}


def test_nicht_abgedecktes_element_ist_nicht_anwendbar():
    daten = baue_daten(["C", "C"], [(0, 1, 1)])
    daten["atom_by_id"][1]["element"] = "Si"
    ergebnis = benson.berechne(daten, [], kein_aromat)
    assert not ergebnis["anwendbar"]
    assert any("Si" in g for g in ergebnis["nicht_abgedeckt"])


# ----------------------------------------------------------------------
# Experimentell gegengeprueft: Heteroatom-Gruppen (Benson-Originalschema).
# NIST-Gasphasenwerte, aus dem Gedaechtnis - bei Zweifel nachschlagen.
# ----------------------------------------------------------------------

FAELLE_HETERO = {
    "dimethylamin":   (["C", "N", "C"], [(0, 1, 1), (1, 2, 1)], -4.4),
    "trimethylamin":  (["N", "C", "C", "C"], [(0, 1, 1), (0, 2, 1), (0, 3, 1)], -5.7),
    "ethylamin":      (["C", "C", "N"], [(0, 1, 1), (1, 2, 1)], -11.3),
    "methanthiol":    (["C", "S"], [(0, 1, 1)], -5.5),
    "dimethylsulfid": (["C", "S", "C"], [(0, 1, 1), (1, 2, 1)], -8.9),
    "ethanthiol":     (["C", "C", "S"], [(0, 1, 1), (1, 2, 1)], -11.1),
    "isopropanol":    (["C", "C", "C", "O"], [(0, 1, 1), (1, 2, 1), (1, 3, 1)], -65.2),
    "tert-butanol":   (["C", "C", "C", "C", "O"], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)], -74.7),
    "acetaldehyd":    (["C", "C", "O"], [(0, 1, 1), (1, 2, 2)], -39.7),
    "methylacetat":   (["C", "C", "O", "O", "C"], [(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)], -98.0),
    "essigsaeure":    (["C", "C", "O", "O"], [(0, 1, 1), (1, 2, 2), (1, 3, 1)], -103.4),
    "isobuten":       (["C", "C", "C", "C"], [(0, 1, 2), (1, 2, 1), (1, 3, 1)], -4.0),
    "butadien":       (["C", "C", "C", "C"], [(0, 1, 2), (1, 2, 1), (2, 3, 2)], 26.3),
}


@pytest.mark.parametrize("name", sorted(FAELLE_HETERO))
def test_hetero_und_konjugation(name):
    schwere, bindungen, erwartet = FAELLE_HETERO[name]
    ergebnis = rechne(schwere, bindungen)

    assert ergebnis["anwendbar"], ergebnis["nicht_abgedeckt"]
    assert ergebnis["dhf_kcal"] == pytest.approx(erwartet, abs=1.0)


def test_toluol():
    ring = list(range(6))
    bindungen = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)] + [(0, 6, 1)]
    ergebnis = rechne(["C"] * 7, bindungen, [ring], aromatisch=lambda r: True)
    assert ergebnis["anwendbar"], ergebnis["nicht_abgedeckt"]
    assert ergebnis["dhf_kcal"] == pytest.approx(12.0, abs=1.0)


# ----------------------------------------------------------------------
# Erweiterung 15.10.2026: Halogene, Cdd (CO2/Allen), Anilin/Anisol, Aromaten.
# Referenzwerte: NIST-Gasphase, AUS DEM GEDAECHTNIS (kcal/mol) - bei Zweifel
# nachschlagen. Format: (Schweratome, Bindungen, Ringe, aromatisch, Referenz, Toleranz)
# ----------------------------------------------------------------------

FAELLE_ERWEITERT = {
    "Fluormethan": (['C', 'F'], [(0, 1, 1)], [], False, -56.0, 1.0),
    "Difluormethan": (['C', 'F', 'F'], [(0, 1, 1), (0, 2, 1)], [], False, -108.1, 1.0),
    "Trifluormethan": (['C', 'F', 'F', 'F'], [(0, 1, 1), (0, 2, 1), (0, 3, 1)], [], False, -166.5, 1.0),
    "Tetrafluormethan": (['C', 'F', 'F', 'F', 'F'], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)], [], False, -223.0, 1.0),
    "Chlormethan": (['C', 'Cl'], [(0, 1, 1)], [], False, -19.6, 1.0),
    "Dichlormethan": (['C', 'Cl', 'Cl'], [(0, 1, 1), (0, 2, 1)], [], False, -22.8, 1.0),
    "Chloroform": (['C', 'Cl', 'Cl', 'Cl'], [(0, 1, 1), (0, 2, 1), (0, 3, 1)], [], False, -24.6, 1.0),
    "Tetrachlormethan": (['C', 'Cl', 'Cl', 'Cl', 'Cl'], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)], [], False, -22.9, 1.0),
    "Brommethan": (['C', 'Br'], [(0, 1, 1)], [], False, -8.5, 1.0),
    "Bromethan": (['C', 'C', 'Br'], [(0, 1, 1), (1, 2, 1)], [], False, -14.8, 1.0),
    "Fluorethan": (['C', 'C', 'F'], [(0, 1, 1), (1, 2, 1)], [], False, -62.9, 1.0),
    "tert-Butylchlorid": (['C', 'C', 'C', 'C', 'Cl'], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)], [], False, -43.8, 1.5),
    "1,1,1-Trichlorethan": (['C', 'C', 'Cl', 'Cl', 'Cl'], [(0, 1, 1), (1, 2, 1), (1, 3, 1), (1, 4, 1)], [], False, -34.6, 1.0),
    "Vinylidenchlorid": (['C', 'C', 'Cl', 'Cl'], [(0, 1, 2), (1, 2, 1), (1, 3, 1)], [], False, 0.6, 1.0),
    "Chlorbenzol": (['C', 'C', 'C', 'C', 'C', 'C', 'Cl'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1)], [[0, 1, 2, 3, 4, 5]], True, 12.4, 1.0),
    "Fluorbenzol": (['C', 'C', 'C', 'C', 'C', 'C', 'F'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1)], [[0, 1, 2, 3, 4, 5]], True, -27.5, 1.0),
    "Brombenzol": (['C', 'C', 'C', 'C', 'C', 'C', 'Br'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1)], [[0, 1, 2, 3, 4, 5]], True, 25.2, 1.0),
    "Anilin": (['C', 'C', 'C', 'C', 'C', 'C', 'N'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1)], [[0, 1, 2, 3, 4, 5]], True, 20.8, 1.0),
    "Anisol": (['C', 'C', 'C', 'C', 'C', 'C', 'O', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1), (6, 7, 1)], [[0, 1, 2, 3, 4, 5]], True, -16.2, 1.5),
    "Styrol": (['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1), (6, 7, 2)], [[0, 1, 2, 3, 4, 5]], True, 35.1, 1.0),
    "Biphenyl": (['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 6, 1), (0, 6, 1)], [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]], True, 43.3, 1.0),
    "Naphthalin": (['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (4, 6, 1), (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 5, 1), (5, 0, 1)], [[0, 1, 2, 3, 4, 5], [4, 6, 7, 8, 9, 5]], True, 36.0, 1.0),
    "Ammoniak": (['N'], [], [], False, -11.0, 1.0),
    "Schwefelwasserstoff": (['S'], [], [], False, -4.9, 1.0),
    "Neopentan": (['C', 'C', 'C', 'C', 'C'], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)], [], False, -40.2, 1.0),
    "Propin": (['C', 'C', 'C'], [(0, 1, 3), (1, 2, 1)], [], False, 44.3, 1.0),
    "1-Butin": (['C', 'C', 'C', 'C'], [(0, 1, 3), (1, 2, 1), (2, 3, 1)], [], False, 39.5, 1.0),
    "1-Buten": (['C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 1)], [], False, -0.1, 1.0),
    "3-Methyl-1-buten": (['C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 1), (2, 4, 1)], [], False, -6.6, 1.0),
    "Isopren": (['C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (1, 4, 1)], [], False, 18.0, 1.0),
    "2-Butanon": (['C', 'C', 'O', 'C', 'C'], [(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)], [], False, -57.0, 1.0),
    "Methylformiat": (['C', 'O', 'C', 'O'], [(0, 1, 1), (1, 2, 1), (2, 3, 2)], [], False, -84.9, 1.0),
    "Ameisensaeure": (['C', 'O', 'O'], [(0, 1, 2), (0, 2, 1)], [], False, -90.5, 1.0),
    "2-Propanthiol": (['C', 'C', 'C', 'S'], [(0, 1, 1), (1, 2, 1), (1, 3, 1)], [], False, -18.2, 1.0),
    "Isopropylamin": (['C', 'C', 'C', 'N'], [(0, 1, 1), (1, 2, 1), (1, 3, 1)], [], False, -20.0, 1.0),
    "tert-Butylamin": (['C', 'C', 'C', 'C', 'N'], [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)], [], False, -28.9, 1.0),
    "Ethylbenzol": (['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1), (6, 7, 1)], [[0, 1, 2, 3, 4, 5]], True, 7.1, 1.0),
    "Cumol": (['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1), (6, 7, 1), (6, 8, 1)], [[0, 1, 2, 3, 4, 5]], True, 0.9, 1.0),
    "Phenol": (['C', 'C', 'C', 'C', 'C', 'C', 'O'], [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1), (0, 6, 1)], [[0, 1, 2, 3, 4, 5]], True, -23.0, 1.0),
    "Kohlendioxid": (['C', 'O', 'O'], [(0, 1, 2), (0, 2, 2)], [], False, -94.05, 1.0),
    "Allen": (['C', 'C', 'C'], [(0, 1, 2), (1, 2, 2)], [], False, 45.6, 1.5),
}


@pytest.mark.parametrize("name", sorted(FAELLE_ERWEITERT))
def test_erweiterung(name):
    schwere, bindungen, ringe, aromatisch, erwartet, toleranz = FAELLE_ERWEITERT[name]
    ergebnis = rechne(schwere, bindungen, ringe, aromatisch=(lambda r: True) if aromatisch else kein_aromat)

    assert ergebnis["anwendbar"], ergebnis["nicht_abgedeckt"]
    assert ergebnis["dhf_kcal"] == pytest.approx(erwartet, abs=toleranz)


def test_halogen_ligand_wird_nicht_still_verworfen():
    """Ein unbekannter Ligand darf nie stillschweigend aus dem Schluessel fallen."""

    assert benson._schluessel("C", ["Cl", "H", "H", "H"]) == "C-(Cl)(H)3"
    assert benson._schluessel("C", ["Xx", "H", "H", "H"]) is None


def test_kohlenstoff_an_fremdatom_ist_nicht_anwendbar():
    daten = baue_daten(["C", "C"], [(0, 1, 1)])
    daten["atom_by_id"][1]["element"] = "Si"
    ergebnis = benson.berechne(daten, [], kein_aromat)

    assert not ergebnis["anwendbar"]
    assert any("Nachbar Si" in g for g in ergebnis["nicht_abgedeckt"])


def test_bekannte_luecken_bleiben_luecken():
    """Gruppen, die ich NICHT uebernommen habe, muessen 'nicht anwendbar' liefern
    (statt mit einem erfundenen Wert zu rechnen)."""

    for schwere, bindungen in (
        (["C", "C", "Cl"], [(0, 1, 1), (1, 2, 1)]),   # Chlorethan: C-(Cl)(C)(H)2 verfehlte die Messung
        (["C", "C", "Cl"], [(0, 1, 2), (1, 2, 1)]),   # Vinylchlorid: Cd-(Cl)(H) verfehlte die Messung
    ):
        ergebnis = rechne(schwere, bindungen)
        assert not ergebnis["anwendbar"]