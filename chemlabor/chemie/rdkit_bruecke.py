# ==========================================================
# rdkit_bruecke.py
# ==========================================================
#
# Bruecke zwischen dem internen Molekuel-Format (dicts mit "atoms"/
# "bonds", wie sie idk.exportiere_molekuel() liefert) und RDKit.
#
# Warum diese Datei existiert: idk.py/kritiker.py pruefen Chemie ueber
# selbstgeschriebene Heuristiken (Valenztabellen, VSEPR-Naeherung,
# eine eigene, NICHT graph-isomorphie-feste Signatur - siehe
# ki_gedaechtnis.molekuel_signatur()). Das funktioniert fuer den
# urspruenglichen Zweck (Lernsignal fuer die KI), stoesst aber an
# Grenzen, sobald man wissen will, ob ein Molekuel nach ECHTEN
# chemischen Regeln ueberhaupt existieren KOENNTE, oder ob zwei
# Strukturen wirklich dasselbe Molekuel sind (Ringsymmetrien,
# Automorphismen - siehe die Kommentare in molekuel_signatur()).
#
# RDKit uebernimmt hier NICHT den Generator oder den Kritiker (die
# bleiben wie sie sind) - es ist eine zusaetzliche, unabhaengige
# Pruef-/Auswertungsschicht:
#   - ist_gueltig(): wuerde dieses Molekuel RDKits Sanitize-Regeln
#     (Valenz, Aromatizitaet, Kekulisierung) ueberstehen?
#   - kanonische_signatur(): echte kanonische SMILES statt der
#     Naeherungssignatur - erkennt auch symmetrische Faelle korrekt
#     als identisch/verschieden.
#   - summenformel()/molmasse(): aus RDKit statt der eigenen,
#     einfachen formel()-Funktion in labor.py.
#
# HINWEIS zur Fehlerbehandlung: RDKits Sanitize-Fehler sind
# C++-Exceptions (ueber Boost.Python durchgereicht), ihre genauen
# Python-Klassen unterscheiden sich je nach RDKit-Version merklich.
# Ein zu enges except hier wuerde bei einem RDKit-Update leise wieder
# unbehandelte Abstuerze produzieren - deshalb bewusst breites
# "except Exception", nur an dieser einen, klar abgegrenzten Stelle.
# ==========================================================


from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


BINDUNGSTYP = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
}


def _zu_roh_rdkit(molekuel):
    """Baut ein (noch nicht sanitisiertes) RWMol aus dem internen
    Format. SetNoImplicit(True) auf jedem Atom, weil unser Format
    Wasserstoff immer als EIGENES Atom mit eigener Bindung fuehrt
    (siehe idk.py) - RDKit soll also keine zusaetzlichen impliziten
    H-Atome dazuerfinden, sonst wuerden ohnehin schon volle Valenzen
    doppelt gezaehlt.

    Rueckgabe: (mol, id_zu_index) - id_zu_index wird von
    aufrufenden Funktionen (z.B. dem Konstruktor) gebraucht, um
    RDKit-Atomindizes auf unsere eigenen Atom-IDs zurueckzufuehren."""

    mol = Chem.RWMol()

    id_zu_index = {}

    for atom in molekuel["atoms"]:

        rd_atom = Chem.Atom(atom["element"])

        rd_atom.SetNoImplicit(True)

        index = mol.AddAtom(rd_atom)

        id_zu_index[atom["id"]] = index

    for bindung in molekuel["bonds"]:

        typ = BINDUNGSTYP.get(bindung["order"])

        if typ is None:
            raise ValueError(
                f"Unbekannte Bindungsordnung: {bindung['order']}"
            )

        mol.AddBond(
            id_zu_index[bindung["atom1"]],
            id_zu_index[bindung["atom2"]],
            typ
        )

    return mol, id_zu_index


def pruefe(molekuel):
    """Zentrale Funktion, auf der alle anderen hier aufbauen.
    Rueckgabe: (gueltig: bool, fehlertext: str|None, rd_mol: Mol|None).

    'gueltig' heisst: RDKit haette gegen dieses Molekuel nach seinen
    eigenen (deutlich strengeren, auf echter Chemie basierenden)
    Regeln nichts einzuwenden - unabhaengig davon, was Kritiker.py
    dazu sagt."""

    try:
        roh, id_zu_index = _zu_roh_rdkit(molekuel)
    except (ValueError, KeyError) as fehler:
        return False, str(fehler), None

    m = roh.GetMol()

    try:
        Chem.SanitizeMol(m)
    except Exception as fehler:  # siehe Modul-Docstring oben
        return False, str(fehler), None

    return True, None, m


def ist_gueltig(molekuel):
    """Kurzform von pruefe(), wenn nur das Ja/Nein interessiert."""

    gueltig, _fehler, _m = pruefe(molekuel)

    return gueltig


def kanonische_signatur(molekuel):
    """Echte kanonische SMILES als Signatur - im Unterschied zu
    ki_gedaechtnis.molekuel_signatur() strukturell/graphentheoretisch
    korrekt (erkennt z.B. symmetrische Ringe zuverlaessig als
    gleich/verschieden). Rueckgabe None, wenn das Molekuel RDKits
    Sanitize-Pruefung nicht besteht - dann ist eine kanonische Form im
    chemischen Sinne ohnehin nicht sinnvoll definiert."""

    gueltig, _fehler, m = pruefe(molekuel)

    if not gueltig:
        return None

    return Chem.MolToSmiles(m)


def summenformel(molekuel):
    """RDKit-Summenformel (Hill-Notation) - Ersatz/Ergaenzung zu
    labor.formel(), die nur simpel durchzaehlt. None, wenn ungueltig."""

    gueltig, _fehler, m = pruefe(molekuel)

    if not gueltig:
        return None

    return rdMolDescriptors.CalcMolFormula(m)


def molmasse(molekuel):
    """Exakte Molmasse in g/mol, oder None, wenn ungueltig."""

    gueltig, _fehler, m = pruefe(molekuel)

    if not gueltig:
        return None

    return rdMolDescriptors.CalcExactMolWt(m)