import json
import random
from chemlabor.pfade import PERIODENSYSTEM


# =====================================
# Einstellungen
# =====================================

MOLEKUEL_MAX_ATOME = 15
ANZAHL_GENERIERUNGEN = 10000


# Grundlegende Valenzen
VALENZEN = {
    "H": 1,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    # P und S duerfen ueber ihre "normale" Valenz (siehe
    # kritiker.py _NORMALE_VALENZ: P=3, S=2) hinausgehen - das ist die
    # deklarierte MAXIMALvalenz, nicht die uebliche. Reale, sehr
    # stabile Molekuele wie PCl5 (P mit 5 Bindungen) oder SF6 (S mit 6
    # Bindungen, "erweitertes Oktett") waren vorher strukturell gar
    # nicht erzeugbar, weil kann_binden() jeden Versuch ueber 3 bzw. 2
    # Bindungen hinaus blockiert hat. Der Kritiker straft die
    # Ueberschreitung ueber sein Hypervalenz-Modell bereits moderat ab
    # (dieselbe Logik, die vorher fuer P/S nie greifen konnte, weil
    # unerreichbar) - diese Grenze hier entscheidet nur, WAS UEBERHAUPT
    # gebaut werden darf, nicht, wie es bewertet wird.
    "P": 5,
    "S": 6,
    "Si": 4,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "B": 3
}


# Erlaubte Bindungsordnungen
# Noch kein Sonderfallwissen!
# Bewusst grosszuegig gehalten (siehe generiere_molekuel): der
# Generator soll auch chemisch fragwuerdige Molekuele erzeugen
# koennen, damit die KI-Gedaechtniskomponente an echten positiven
# UND negativen Beispielen lernt - der Kritiker sortiert danach aus.
BINDUNGSARTEN = {
    "H": [1],
    "F": [1],
    "Cl": [1],
    "Br": [1],
    "I": [1],
    "O": [1, 2],
    "S": [1, 2],
    "N": [1, 2, 3],
    "C": [1, 2, 3],
    "Si": [1, 2],
    "P": [1, 2, 3],
    "B": [1]
}


EDELGASE = {
    "He",
    "Ne",
    "Ar",
    "Kr",
    "Xe",
    "Rn",
    "Og"
}


# =====================================
# Elemente laden
# =====================================

def lade_elemente(datei):

    with open(datei, "r", encoding="utf-8") as f:
        daten = json.load(f)


    elemente = []

    for element in daten["elements"]:

        symbol = element["symbol"]

        if symbol in VALENZEN:

            if symbol not in EDELGASE:
                elemente.append(symbol)


    return elemente



# =====================================
# Atom
# =====================================

class Atom:

    def __init__(self, atom_id, symbol):

        self.id = atom_id
        self.symbol = symbol

        self.valenz = VALENZEN[symbol]

        self.verbrauchte_valenz = 0


    def kann_binden(self, ordnung):

        return (
            self.verbrauchte_valenz + ordnung
            <= self.valenz
        )


    def binde(self, ordnung):

        self.verbrauchte_valenz += ordnung



# =====================================
# Bindung
# =====================================

class Bindung:

    def __init__(self, atom1, atom2, ordnung):

        self.atom1 = atom1
        self.atom2 = atom2
        self.ordnung = ordnung



# =====================================
# Molekül erzeugen
# =====================================

def generiere_molekuel(elemente, max_atome=None):

    anzahl = random.randint(
        2,
        (max_atome if max_atome is not None else MOLEKUEL_MAX_ATOME)
    )


    atome = []

    for i in range(anzahl):

        symbol = random.choice(elemente)

        atome.append(
            Atom(i, symbol)
        )


    bindungen = []


    # FIX: vorher wurde bei jedem Versuch die komplette bindungen-Liste
    # linear durchsucht, um ein bereits verbundenes Atompaar zu erkennen
    # (O(Bindungen) pro Versuch). Ein Set aus frozenset-Paaren macht die
    # Pruefung O(1) - bei MOLEKUEL_MAX_ATOME=15 kaum messbar, aber sobald
    # der Wert (z.B. ueber die GUI-Einstellung) deutlich hochgesetzt
    # wird, faellt der Unterschied ins Gewicht.
    bindung_paare = set()


    versuche = anzahl * 10


    for _ in range(versuche):

        atom1 = random.choice(atome)
        atom2 = random.choice(atome)


        if atom1 == atom2:
            continue


        # gleiche Verbindung vermeiden

        paar = frozenset([atom1.id, atom2.id])

        if paar in bindung_paare:
            continue



        # erlaubte Bindungsordnung auswählen

        moegliche = BINDUNGSARTEN[
            atom1.symbol
        ]


        ordnung = random.choice(
            moegliche
        )


        if (
            atom1.kann_binden(ordnung)
            and
            atom2.kann_binden(ordnung)
        ):

            atom1.binde(ordnung)
            atom2.binde(ordnung)


            bindungen.append(
                Bindung(
                    atom1,
                    atom2,
                    ordnung
                )
            )

            bindung_paare.add(paar)


    return atome, bindungen



# =====================================
# Export für Kritiker
# =====================================

def exportiere_molekuel(atome, bindungen):

    daten = {

        "atoms": [],

        "bonds": []

    }


    for atom in atome:

        daten["atoms"].append(
            {
                "id": atom.id,

                "element": atom.symbol,

                "valence": atom.valenz,

                "used_valence":
                    atom.verbrauchte_valenz
            }
        )


    for bindung in bindungen:

        daten["bonds"].append(
            {
                "atom1":
                    bindung.atom1.id,

                "atom2":
                    bindung.atom2.id,

                "order":
                    bindung.ordnung
            }
        )


    return daten



# =====================================
# Formel
# =====================================

def formel(atome):

    zaehler = {}


    for atom in atome:

        if atom.symbol not in zaehler:
            zaehler[atom.symbol] = 0

        zaehler[atom.symbol] += 1


    ergebnis = ""

    for element, menge in zaehler.items():

        ergebnis += element

        if menge > 1:
            ergebnis += str(menge)


    return ergebnis



# =====================================
# Anzeige
# =====================================

def anzeigen(atome, bindungen):

    print("Formel:")
    print(formel(atome))


    print("\nAtome:")

    for atom in atome:

        print(
            atom.id,
            atom.symbol
        )


    print("\nBindungen:")

    for b in bindungen:

        print(
            b.atom1.id,
            "--",
            b.atom2.id,
            "Ordnung:",
            b.ordnung
        )


    print("\nKritiker-Daten:")

    daten = exportiere_molekuel(
        atome,
        bindungen
    )

    print(
        json.dumps(
            daten,
            indent=4
        )
    )



# =====================================
# Hauptprogramm
# =====================================

def main():

    elemente = lade_elemente(
        PERIODENSYSTEM
    )


    print(
        "Verfügbare Elemente:",
        elemente
    )


    print("\nGenerierte Moleküle:\n")


    erzeugt = 0


    while erzeugt < ANZAHL_GENERIERUNGEN:


        atome, bindungen = generiere_molekuel(
            elemente
        )


        if len(bindungen) > 0:

            anzeigen(
                atome,
                bindungen
            )


            print(
                "===================="
            )


            erzeugt += 1



if __name__ == "__main__":
    main()