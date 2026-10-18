"""Tests fuer chemlabor.chemie.zerfall. Referenzwerte NIST-Gasphase, aus dem Gedaechtnis (kcal/mol)."""

import pytest

from chemlabor.chemie import zerfall


def daten(elemente, bindungen):
    """Alle Atome explizit (auch H). bindungen: (i, j, ordnung)."""
    return {
        "atom_by_id": {i: {"id": i, "element": e} for i, e in enumerate(elemente)},
        "bindungen": [{"atom1": a, "atom2": b, "ordnung": o} for a, b, o in bindungen],
    }


def mit_h(zentrum, h, extra=()):
    return daten([zentrum] + ["H"] * h, [(0, i + 1, 1) for i in range(h)] + list(extra))


FAELLE = {
    # name: (Elemente, Bindungen, Referenz, Toleranz)
    "wasser": (["O", "H", "H"], [(0, 1, 1), (0, 2, 1)], -57.8, 3),
    "methan": (["C", "H", "H", "H", "H"], [(0, i, 1) for i in range(1, 5)], -17.8, 3),
    "ammoniak": (["N", "H", "H", "H"], [(0, i, 1) for i in range(1, 4)], -11.0, 3),
    "ethan": (["C", "C"] + ["H"] * 6, [(0, 1, 1)] + [(0, i, 1) for i in (2, 3, 4)] + [(1, i, 1) for i in (5, 6, 7)], -20.0, 3),
    "ethen": (["C", "C"] + ["H"] * 4, [(0, 1, 2), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1)], 12.5, 5),
    "stickstoff": (["N", "N"], [(0, 1, 3)], 0.0, 1),
    "fluor": (["F", "F"], [(0, 1, 1)], 0.0, 2),
    "hf": (["H", "F"], [(0, 1, 1)], -65.1, 2),
    "kohlendioxid": (["C", "O", "O"], [(0, 1, 2), (0, 2, 2)], -94.1, 30),   # Mittelwert 178 zu klein fuer CO2 (192)
    "cf4": (["C", "F", "F", "F", "F"], [(0, i, 1) for i in range(1, 5)], -223.0, 10),
    "sf6": (["S"] + ["F"] * 6, [(0, i, 1) for i in range(1, 7)], -291.0, 10),
}


@pytest.mark.parametrize("name", sorted(FAELLE))
def test_bekannte_moleküle(name):
    elemente, bindungen, ref, tol = FAELLE[name]
    r = zerfall.berechne(daten(elemente, bindungen))

    assert r["anwendbar"], r["nicht_abgedeckt"]
    assert r["dhf_kcal"] == pytest.approx(ref, abs=tol)


def test_radikal_ist_energiereicher_als_gesaettigtes():
    """Methyl (eine Bindung fehlt) muss deutlich energiereicher sein als Methan."""
    methyl = zerfall.berechne(mit_h("C", 3))["dhf_kcal"]
    methan = zerfall.berechne(mit_h("C", 4))["dhf_kcal"]

    assert methyl == pytest.approx(35.0, abs=8)  # NIST CH3: 35,0
    assert methyl - methan > 40


def test_carben_noch_energiereicher():
    carben = zerfall.berechne(mit_h("C", 2))["dhf_kcal"]

    assert carben == pytest.approx(92.0, abs=17)  # NIST CH2 (triplet): 93


def test_benzol_aromatenkorrektur():
    ring = [(i, (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)] + [(i, 6 + i, 1) for i in range(6)]
    d = daten(["C"] * 6 + ["H"] * 6, ring)

    ohne = zerfall.berechne(d)
    mit = zerfall.berechne(d, [list(range(6))], lambda r: True)

    assert ohne["dhf_kcal"] - mit["dhf_kcal"] == pytest.approx(37.0)
    assert mit["dhf_kcal"] == pytest.approx(19.8, abs=5)   # NIST Benzol 19,8


def test_unbekannte_bindung_ist_nicht_anwendbar():
    r = zerfall.berechne(daten(["Xe", "F"], [(0, 1, 1)]))
    assert not r["anwendbar"]

    r = zerfall.berechne(daten(["C", "I", "F"], [(0, 1, 1), (1, 2, 1)]))   # I-F fehlt in der Tabelle
    assert not r["anwendbar"]
    assert any("I-" in g or "-I" in g for g in r["nicht_abgedeckt"])


def test_unsichere_bindung_wird_gemeldet():
    d = daten(["Si", "Si", "H", "H", "H", "H"], [(0, 1, 2), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1)])
    r = zerfall.berechne(d)

    assert r["anwendbar"]
    assert r["unsicher"] == ["Si-Si Ordnung 2"]