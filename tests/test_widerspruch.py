from chemlabor.chemie.widerspruch import klassifiziere, lernpunkte
from chemlabor.gui.modell import ErgebnisModell


def eintrag(nummer, punkte, rdkit):
    molekuel = {"atoms": [{"id": 0, "element": "H", "valence": 1, "used_valence": 0}], "bonds": []}
    analyse = {"bewertung": {"punkte": punkte, "urteil": "x"}, "probleme": [], "erklaerung": [],
               "rdkit_gueltig": rdkit}
    return (nummer, "Beobachtung", molekuel, analyse)


def test_klassifiziere():
    assert klassifiziere(90, False) == "milde"
    assert klassifiziere(75, False) == "milde"
    assert klassifiziere(74, False) is None
    assert klassifiziere(30, True) == "streng"
    assert klassifiziere(50, True) is None
    assert klassifiziere(90, True) is None
    assert klassifiziere(30, False) is None
    assert klassifiziere(90, None) is None


def test_lernpunkte_wirkt_nur_bei_milde_und_ist_klein():
    assert lernpunkte(90, False, 3) == 87
    assert lernpunkte(90, True, 3) == 90        # RDKit belohnt nie
    assert lernpunkte(30, True, 3) == 30        # 'streng' aendert nichts
    assert lernpunkte(30, False, 3) == 30
    assert lernpunkte(90, None, 3) == 90
    assert lernpunkte(90, False, 0) == 90       # aus
    assert lernpunkte(90, False, 99) == 80      # harte Obergrenze 10


def test_modell_filter_und_kennzahlen():
    modell = ErgebnisModell()
    for e in (eintrag(1, 90, False), eintrag(2, 90, True), eintrag(3, 30, True),
              eintrag(4, 30, False), eintrag(5, 85, False)):
        modell.hinzufuegen(e)

    milde = modell.filtere(filter_wahl="w_milde")
    streng = modell.filtere(filter_wahl="w_streng")

    assert sorted(e[0] for e in milde) == [1, 5]
    assert [e[0] for e in streng] == [3]

    k = modell.kennzahlen()
    assert k["widerspruch_milde"] == 2
    assert k["widerspruch_streng"] == 1