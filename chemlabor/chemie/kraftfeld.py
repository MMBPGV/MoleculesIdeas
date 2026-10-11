# ==========================================================
# kraftfeld.py
# ==========================================================
#
# Geometrie-Signal aus dem MMFF94-Kraftfeld von RDKit: Laesst sich das
# Molekuel raeumlich ueberhaupt bauen, und wie stark muessen Bindungen und
# Winkel dafuer verbogen werden?
#
# Gemessen wird nur die "Verbiegungsenergie" (Bindungsdehnung + Winkel +
# Dehnung-Biegung + Ebenheit) je Schweratom in kcal/mol - NICHT die
# Gesamtenergie (die waere zwischen verschiedenen Molekuelen nicht
# vergleichbar).
#
# Was das Signal kann und was nicht (gemessen am Referenzset):
#   + erkennt stark verzerrte Geometrien (Cyclopropen, Tetrahedran,
#     Cumulen im Ring) und nicht bauebare Strukturen (Cyclopropin)
#   - sieht NICHT, dass eine Bindung energetisch schwach ist (O-O-, N-N-
#     Ketten): Kraftfelder haben harmonische Bindungen
#   - MMFF hat fuer bekannte Ringe Sonderparameter: Cyclopropan zaehlt als
#     "normal". Hier bleibt die Ringspannungs-Heuristik des Kritikers
#     wichtig.
#
# Konservativ gedacht: nur ein Abzug, nie ein Bonus; fehlen Parameter oder
# lehnt RDKit das Molekuel ab, passiert nichts (kein Signal = kein Abzug).
# Ergebnisse stehen im Cache nur als kleine Zahlen, nie als Mol-Objekt.

from functools import lru_cache

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdForceFieldHelpers as FFH

from chemlabor.chemie import rdkit_bruecke as rb

RDLogger.DisableLog("rdApp.*")

# Ab dieser Verbiegungsenergie je Schweratom beginnt der Abzug
SCHWELLE_PRO_ATOM = 3.0
# Punkte Abzug je kcal/mol je Schweratom ueber der Schwelle
FAKTOR = 2.0
MAX_ABZUG = 20
# Abzug, wenn sich das Molekuel nicht in 3D einbetten laesst
ABZUG_NICHT_BAUBAR = 15

_VERBIEGUNG = (
    "SetMMFFBondTerm", "SetMMFFAngleTerm",
    "SetMMFFStretchBendTerm", "SetMMFFOopTerm",
)
_ALLE = _VERBIEGUNG + (
    "SetMMFFTorsionTerm", "SetMMFFVdWTerm", "SetMMFFEleTerm",
)


@lru_cache(maxsize=20000)
def _messe_cached(key):
    """(status, verbiegung_pro_schweratom). status: 'ok', 'nicht_baubar',
    'keine_parameter' oder 'abgelehnt' (RDKit-Sanitize schlaegt fehl)."""

    gueltig, _fehler, mol = rb._pruefe_raw(rb._aus_key(key))

    if not gueltig:
        return "abgelehnt", 0.0

    schwer = mol.GetNumHeavyAtoms()

    # Zwei bis drei Atome sind immer baubar (und ohne Ring nie verspannt)
    if mol.GetNumAtoms() <= 3:
        return "ok", 0.0

    if not FFH.MMFFHasAllMoleculeParams(mol):
        return "keine_parameter", 0.0

    parameter = AllChem.ETKDGv3()
    parameter.randomSeed = 7
    parameter.maxIterations = 200

    konf = AllChem.EmbedMolecule(mol, parameter)

    if konf < 0:
        parameter.useRandomCoords = True
        konf = AllChem.EmbedMolecule(mol, parameter)

    if konf < 0:
        return "nicht_baubar", 0.0

    FFH.MMFFOptimizeMolecule(mol, maxIters=500)

    eigenschaften = FFH.MMFFGetMoleculeProperties(mol)

    if eigenschaften is None:
        return "keine_parameter", 0.0

    for setter in _ALLE:
        getattr(eigenschaften, setter)(setter in _VERBIEGUNG)

    feld = FFH.MMFFGetMoleculeForceField(mol, eigenschaften, confId=konf)

    energie = max(0.0, feld.CalcEnergy())

    return "ok", energie / max(1, schwer)


def messe(molekuel):
    """(status, verbiegung_pro_schweratom) - siehe _messe_cached."""

    return _messe_cached(rb._molekuel_cache_key(molekuel))


def abzug(molekuel):
    """(punkte_abzug >= 0, status, verbiegung) - der Abzug fuer den
    Kritiker. 0, solange kein belastbares Signal vorliegt."""

    status, verbiegung = messe(molekuel)

    if status == "nicht_baubar":
        return ABZUG_NICHT_BAUBAR, status, verbiegung

    if status == "ok" and verbiegung > SCHWELLE_PRO_ATOM:
        punkte = min(MAX_ABZUG, round((verbiegung - SCHWELLE_PRO_ATOM) * FAKTOR))
        return punkte, status, verbiegung

    return 0, status, verbiegung


def cache_leeren():

    _messe_cached.cache_clear()