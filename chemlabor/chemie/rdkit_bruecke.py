
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

import json
from functools import lru_cache

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


BINDUNGSTYP = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
}


def _molekuel_cache_key(molekuel):
    """Normalisiert Molekuel-Daten in einen stabilen String, damit
    wiederholte RDKit-Pruefungen zwischen derselben Struktur mit
    unterschiedlicher interner Reihenfolge/Dict-Reihenfolge nicht
    erneut berechnet werden muessen.
    """

    atome = tuple(
        sorted(
            (
                atom.get("id"),
                atom.get("element"),
                atom.get("valence"),
                atom.get("used_valence"),
            )
            for atom in molekuel.get("atoms", [])
        )
    )

    bindungen = tuple(
        sorted(
            (
                min(b.get("atom1"), b.get("atom2")),
                max(b.get("atom1"), b.get("atom2")),
                int(b.get("order", 1)),
            )
            for b in molekuel.get("bonds", [])
        )
    )

    return json.dumps(
        {"atoms": atome, "bonds": bindungen},
        separators=(",", ":"),
        sort_keys=True,
    )


@lru_cache(maxsize=20000)
def _pruefe_cached(serialisiert):
    """Zentrale RDKit-Pruefung mit LRU-Caching. Der Cache laeuft auf
    einem normalisierten Molekuel-Key und vermeidet wiederholte
    Sanitize-Aufrufe fuer dieselbe Struktur.
    """

    daten = json.loads(serialisiert)

    molekuel = {
        "atoms": [
            {"id": i, "element": e, "valence": v, "used_valence": u}
            for i, e, v, u in daten["atoms"]
        ],
        "bonds": [
            {"atom1": a, "atom2": b, "order": o}
            for a, b, o in daten["bonds"]
        ],
    }

    return _pruefe_raw(molekuel)


def _pruefe_raw(molekuel):
    """Rohlogik zur RDKit-Pruefung ohne Cache-Schicht."""

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

    key = _molekuel_cache_key(molekuel)
    return _pruefe_cached(key)


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


# ==========================================================
# End of file
# ==========================================================
