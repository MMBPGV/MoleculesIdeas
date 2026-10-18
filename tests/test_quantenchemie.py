"""Tests fuer chemlabor.chemie.quantenchemie (ANI-2x). Werden uebersprungen, wenn torch/torchani fehlen.
Referenz: NIST-Bildungsenthalpien (Gasphase, kcal/mol), aus dem Gedaechtnis."""

import pytest

pytest.importorskip("torch")
pytest.importorskip("torchani")
pytest.importorskip("scipy")
pytest.importorskip("rdkit")

from chemlabor.chemie import quantenchemie as q  # noqa: E402


@pytest.fixture(scope="module")
def qc():
    return q.Quantenchemie()


def mol(*symbole_und_bindungen):
    atome, bindungen = symbole_und_bindungen
    return {"atoms": [{"id": i, "element": e} for i, e in enumerate(atome)],
            "bonds": [{"atom1": a, "atom2": b, "order": o} for a, b, o in bindungen]}


@pytest.mark.parametrize("smiles,referenz", [
    ("CCO", -56.2), ("c1ccccc1", 19.8), ("CC(C)=O", -52.2), ("CCl", -19.6), ("C1CC1", 12.7),
])
def test_bildungsenthalpie_im_fehlerbereich(qc, smiles, referenz):
    r = qc.berechne(q.molekuel_aus_smiles(smiles))

    assert r["anwendbar"], r["grund"]
    assert r["dhf_kcal"] == pytest.approx(referenz, abs=8.0)
    assert r["intakt"] is True


def test_carben_ist_nicht_anwendbar(qc):
    r = qc.berechne(mol(["C", "H", "H"], [(0, 1, 1), (0, 2, 1)]))

    assert not r["anwendbar"]
    assert "offene Schale" in r["grund"]


def test_nitren_trotz_gerader_elektronenzahl_nicht_anwendbar(qc):
    r = qc.berechne(mol(["N", "H"], [(0, 1, 1)]))      # 8 Elektronen, aber Triplett

    assert not r["anwendbar"]


def test_sauerstoff_ist_triplett(qc):
    r = qc.berechne(mol(["O", "O"], [(0, 1, 2)]))

    assert not r["anwendbar"]
    assert "Triplett" in r["grund"]


def test_fremdes_element_wird_gemeldet_nicht_geraten(qc):
    r = qc.berechne(mol(["Si", "H", "H", "H", "H"], [(0, i, 1) for i in range(1, 5)]))

    assert not r["anwendbar"]
    assert "Si" in r["grund"]


def test_unsinnige_struktur_ohne_startgeometrie(qc):
    r = qc.berechne(mol(["C", "H", "H", "H", "H", "H"], [(0, i, 1) for i in range(1, 6)]))   # CH5

    assert not r["anwendbar"]


def test_kalibrierung_deckt_alle_ani_elemente(qc):
    assert set(q.KALIBRIERUNG) == set(q.ANI_ELEMENTE)
    assert set(q.NORMALE_VALENZ) == set(q.ANI_ELEMENTE)


def test_verfuegbar_meldet_installierte_pakete():
    ok, hinweis = q.verfuegbar()

    assert ok and hinweis == ""