"""Zerfallsenergie-Term: Bildungsenthalpie aus mittleren Bindungsenergien.

dHf_schaetz = Summe(dHf der freien Atome) - Summe(Bindungsenergien)   [kcal/mol]

Das ist die Energie, die frei wird (negativ) bzw. aufgewendet werden muss (positiv),
wenn man das Molekuel aus seinen Elementen in Standardzustand bildet. Ein Molekuel mit
hohem positivem Wert zerfaellt unter Energiegewinn in stabile Kleinmolekuele.

Unterschied zu benson.py: keine Gruppenzerlegung. Es zaehlen nur Atompaar + Bindungsordnung,
deshalb sind auch Radikale (fehlende Bindung), Carbene, Si/B/P und N-Ketten abgedeckt.
Preis: Genauigkeit nur +-10..20 kcal/mol; Ringspannung und Resonanz fehlen (das macht Benson).

Werte: mittlere Bindungsenthalpien (CRC/Huheey), aus dem Gedaechtnis, auf ganze kcal/mol.
UNSICHER = Bindungen, deren Wert ich nur grob kenne (v.a. pi-Bindungen schwerer Elemente).
"""

# dHf der freien Atome, kcal/mol
ATOM_HF = {"H": 52.1, "B": 133.8, "C": 171.3, "N": 113.0, "O": 59.6, "F": 19.0, "Si": 107.6,
           "P": 75.6, "S": 66.6, "Cl": 29.1, "Br": 26.7, "I": 25.5}

# (Element, Element, Ordnung) -> kcal/mol ; Elementpaare alphabetisch sortiert
_BE = {
    ("H", "H", 1): 104, ("C", "H", 1): 99, ("H", "N", 1): 93, ("H", "O", 1): 111, ("F", "H", 1): 136,
    ("Cl", "H", 1): 103, ("Br", "H", 1): 88, ("H", "I", 1): 71, ("H", "S", 1): 81, ("H", "Si", 1): 76,
    ("H", "P", 1): 77, ("B", "H", 1): 93,
    ("C", "C", 1): 83, ("C", "C", 2): 146, ("C", "C", 3): 200,
    ("C", "N", 1): 73, ("C", "N", 2): 147, ("C", "N", 3): 213,
    ("C", "O", 1): 86, ("C", "O", 2): 178,
    ("C", "F", 1): 116, ("C", "Cl", 1): 81, ("Br", "C", 1): 68, ("C", "I", 1): 51,
    ("C", "S", 1): 65, ("C", "S", 2): 128,
    ("C", "Si", 1): 76, ("C", "P", 1): 70, ("B", "C", 1): 89,
    ("N", "N", 1): 39, ("N", "N", 2): 100, ("N", "N", 3): 226,
    ("N", "O", 1): 48, ("N", "O", 2): 145, ("F", "N", 1): 66, ("Cl", "N", 1): 46,
    ("O", "O", 1): 35, ("O", "O", 2): 119, ("F", "O", 1): 45, ("Cl", "O", 1): 49, ("Br", "O", 1): 48,
    ("O", "S", 1): 87, ("O", "S", 2): 126,
    ("S", "S", 1): 54, ("S", "S", 2): 100,
    ("Si", "Si", 1): 54, ("O", "Si", 1): 111, ("N", "Si", 1): 76, ("F", "Si", 1): 135, ("Cl", "Si", 1): 91,
    ("P", "P", 1): 50, ("O", "P", 1): 84, ("O", "P", 2): 130, ("N", "P", 1): 70,
    ("B", "O", 1): 125, ("B", "N", 1): 105, ("B", "F", 1): 150, ("B", "Cl", 1): 105, ("B", "B", 1): 71,
    ("F", "F", 1): 37, ("Cl", "Cl", 1): 58, ("Br", "Br", 1): 46, ("I", "I", 1): 36, ("Cl", "F", 1): 61,
    # unsicher (pi-Bindungen schwerer Elemente, grobe Schaetzung)
    ("F", "S", 1): 78, ("Cl", "S", 1): 61, ("P", "S", 2): 90,
    ("Si", "Si", 2): 75, ("C", "Si", 2): 95, ("O", "Si", 2): 190, ("C", "P", 2): 110, ("P", "P", 2): 85,
    ("B", "O", 2): 190, ("B", "N", 2): 140, ("B", "C", 2): 130,
}

UNSICHER = {("F", "S", 1), ("Cl", "S", 1), ("P", "S", 2), ("Si", "Si", 2), ("C", "Si", 2), ("O", "Si", 2), ("C", "P", 2), ("P", "P", 2),
            ("B", "O", 2), ("B", "N", 2), ("B", "C", 2)}


def _paar(a, b):
    return (a, b) if a <= b else (b, a)


# Resonanzstabilisierung je aromatischem Ring, kcal/mol (Benzol 36, Furan/Pyrrol ~ 16-22, Thiophen ~ 29)
AROMAT_6 = -37.0
AROMAT_5 = -25.0


def berechne(daten, ringe=None, ist_aromatisch=None):
    """ringe/ist_aromatisch optional (wie bei benson.berechne): aromatische Ringe bekommen eine
    Resonanzkorrektur, weil die Kekule-Struktur sonst zu energiereich erscheint.

    Gibt {anwendbar, dhf_kcal, schwere_atome, bindungen_gezaehlt, nicht_abgedeckt, unsicher}."""

    atome = daten["atom_by_id"]
    ergebnis = {"anwendbar": False, "dhf_kcal": None, "schwere_atome": 0,
                "nicht_abgedeckt": [], "unsicher": []}

    ergebnis["schwere_atome"] = sum(1 for a in atome.values() if a["element"] != "H")

    fehlt = sorted({a["element"] for a in atome.values() if a["element"] not in ATOM_HF})
    if fehlt:
        ergebnis["nicht_abgedeckt"].append("Element ohne Atomwert: " + ", ".join(fehlt))
        return ergebnis

    summe = sum(ATOM_HF[a["element"]] for a in atome.values())

    for b in daten["bindungen"]:
        e1, e2 = atome[b["atom1"]]["element"], atome[b["atom2"]]["element"]
        schluessel = _paar(e1, e2) + (b["ordnung"],)
        wert = _BE.get(schluessel)
        if wert is None:
            ergebnis["nicht_abgedeckt"].append(f"Bindung {e1}-{e2} Ordnung {b['ordnung']} ohne Energiewert")
            continue
        if schluessel in UNSICHER:
            ergebnis["unsicher"].append(f"{e1}-{e2} Ordnung {b['ordnung']}")
        summe -= wert

    if ergebnis["nicht_abgedeckt"]:
        return ergebnis

    ergebnis["dhf_ohne_aromat"] = round(summe, 1)

    if ringe and ist_aromatisch:
        for ring in ringe:
            if len(ring) in (5, 6) and ist_aromatisch(ring):
                summe += AROMAT_6 if len(ring) == 6 else AROMAT_5

    ergebnis["anwendbar"] = True
    ergebnis["dhf_kcal"] = round(summe, 1)
    return ergebnis