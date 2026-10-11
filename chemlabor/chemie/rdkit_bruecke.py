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
#
# Oeffentliche Funktionen (von labor.py, modell.py, neuheit.py u.a.
# benutzt): pruefe, ist_gueltig, kanonische_signatur, summenformel,
# molmasse. Intern: _molekuel_cache_key, _pruefe_cached, _pruefe_raw,
# _zu_roh_rdkit (der Konstruktor ruft _zu_roh_rdkit direkt auf).
# ==========================================================

import json
from functools import lru_cache

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

# Ungueltige Molekuele sind hier der Normalfall; RDKit schreibt sonst zu
# jedem eine Zeile auf stderr (spuerbar langsam, v.a. in der Windows-Konsole).
RDLogger.DisableLog("rdApp.*")


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


def _aus_key(serialisiert):
    """Baut aus dem Cache-Key wieder ein Molekuel im internen Format. JSON
    macht aus den Tupeln Listen - daher werden die dicts hier neu gebaut."""

    daten = json.loads(serialisiert)

    return {
        "atoms": [
            {"id": i, "element": e, "valence": v, "used_valence": u}
            for i, e, v, u in daten["atoms"]
        ],
        "bonds": [
            {"atom1": a, "atom2": b, "order": o}
            for a, b, o in daten["bonds"]
        ],
    }


# WICHTIG (Speicher): in den Caches stehen nur kleine Werte (bool, Text,
# Zahl), NIE RDKit-Mol-Objekte. Frueher hielt der Cache bis zu 20000 Mols
# fest - bei langen Automatik-/Endlos-Laeufen waren das viele hundert MB.
CACHE_GROESSE = 20000


@lru_cache(maxsize=CACHE_GROESSE)
def _pruefe_cached(serialisiert):
    """(gueltig, fehlertext) - ohne Mol-Objekt, siehe oben."""

    gueltig, fehler, _m = _pruefe_raw(_aus_key(serialisiert))

    return gueltig, fehler


@lru_cache(maxsize=CACHE_GROESSE)
def _smiles_cached(serialisiert):

    gueltig, _fehler, m = _pruefe_raw(_aus_key(serialisiert))

    return Chem.MolToSmiles(m) if gueltig else None


@lru_cache(maxsize=CACHE_GROESSE)
def _formel_cached(serialisiert):

    gueltig, _fehler, m = _pruefe_raw(_aus_key(serialisiert))

    return rdMolDescriptors.CalcMolFormula(m) if gueltig else None


@lru_cache(maxsize=CACHE_GROESSE)
def _masse_cached(serialisiert):

    gueltig, _fehler, m = _pruefe_raw(_aus_key(serialisiert))

    return rdMolDescriptors.CalcExactMolWt(m) if gueltig else None


def cache_leeren():
    """Leert alle Caches (z.B. zwischen zwei Laeufen)."""

    for funktion in (_pruefe_cached, _smiles_cached, _formel_cached, _masse_cached):
        funktion.cache_clear()


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
    dazu sagt.

    Das Ja/Nein kommt aus dem Cache; das Mol-Objekt wird nur bei einem
    gueltigen Molekuel neu gebaut (und nicht im Cache gehalten)."""

    key = _molekuel_cache_key(molekuel)

    gueltig, fehler = _pruefe_cached(key)

    if not gueltig:
        return False, fehler, None

    return _pruefe_raw(_aus_key(key))


def ist_gueltig(molekuel):
    """Nur Ja/Nein - der schnelle Weg ohne Mol-Objekt."""

    return _pruefe_cached(_molekuel_cache_key(molekuel))[0]


def kanonische_signatur(molekuel):
    """Echte kanonische SMILES als Signatur - strukturell/graphentheoretisch
    korrekt (erkennt z.B. symmetrische Ringe zuverlaessig als gleich/
    verschieden). None, wenn das Molekuel RDKits Sanitize-Pruefung nicht
    besteht."""

    return _smiles_cached(_molekuel_cache_key(molekuel))


def summenformel(molekuel):
    """RDKit-Summenformel (Hill-Notation); None, wenn ungueltig."""

    return _formel_cached(_molekuel_cache_key(molekuel))


def molmasse(molekuel):
    """Exakte Molmasse in g/mol, oder None, wenn ungueltig."""

    return _masse_cached(_molekuel_cache_key(molekuel))


# ==========================================================
# End of file
# ==========================================================