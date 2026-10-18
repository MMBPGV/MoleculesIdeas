"""Benson-Gruppenadditivitaet: Bildungsenthalpie (dHf) aus Gruppenbeitraegen.

Eigenstaendiges Modul ohne Abhaengigkeit vom Kritiker. Eingabe ist das
"daten"-Dict des Kritikers (atome, graph, bindung_by_paar, atom_by_id).

Gedacht als ZUSAETZLICHER Term, nicht als Ersatz: Wo die Tabelle ein Atom
oder einen Ring nicht abdeckt (Si, P, B, exotische Ringe ...), liefert
berechne() anwendbar=False. Der Aufrufer sollte das als "keine Aussage"
behandeln, nicht als Strafe.

Einheit: kcal/mol. Naeherung: keine gauche-/cis-/ortho-Korrekturen.

Schluesselformat der Tabelle: "Zentrum-(Ligand)n(Ligand)m ...", Liganden in
der Reihenfolge von _LIGANDEN_REIHENFOLGE. Typen: C = sp3-C, Cd = C=C,
Ct = C#C, Cb = aromatisches C, CO = Carbonyl-C, Cdd = Allen-/CO2-Zentrum,
Od = doppelt gebundener O am Cdd-Zentrum.

Halogengruppen sind (wie bei Benson) mit generischem Kohlenstoff-Liganden
"C" abgelegt; findet sich fuer einen Schluessel mit Cd/Ct/Cb/CO-Liganden
kein spezifischer Eintrag, wird derselbe Schluessel mit "C" gesucht.
"""

import json
import os
from collections import Counter

# Projektlayout: <wurzel>/chemlabor/chemie/benson.py  und  <wurzel>/daten/benson.json
_TABELLE_PFAD = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "daten", "benson.json")
)

# Reihenfolge der Liganden im Schluessel (muss zur JSON-Tabelle passen).
_LIGANDEN_REIHENFOLGE = ["CO", "Od", "Cd", "Ct", "Cb", "O", "N", "S", "F", "Cl", "Br", "I", "C", "H"]

_HALOGENE = ("F", "Cl", "Br", "I")
_C_TYPEN = ("CO", "Cd", "Ct", "Cb")

_UNBEKANNT = "?"  # Typ-Marker fuer Atome, die Benson hier nicht abdeckt

_cache = {}


def lade_tabelle(pfad=None):
    """Laedt die Tabelle einmalig (oder neu, wenn ein Pfad uebergeben wird)."""

    pfad = pfad or _TABELLE_PFAD

    if pfad not in _cache:
        with open(pfad, encoding="utf-8") as f:
            roh = json.load(f)

        gruppen = {k: v["wert"] for k, v in roh["gruppen"].items() if not k.startswith("_")}
        ringe = {k: v["wert"] for k, v in roh["ringkorrekturen"].items() if not k.startswith("_")}
        _cache[pfad] = {"gruppen": gruppen, "ringe": ringe, "meta": roh.get("_meta", {})}

    return _cache[pfad]


def _ordnung(daten, a, b):
    bindung = daten["bindung_by_paar"].get(frozenset([a, b]))
    return bindung["ordnung"] if bindung is not None else 1


def _atom_typen(daten, aromatische_c, aromatisch_unbekannt):
    """Ordnet jedem Atom einen Benson-Typ zu: C, Cd, Ct, Cb, CO, Cdd oder das
    Elementsymbol. Nicht abdeckbare Faelle bekommen _UNBEKANNT."""

    typen = {}
    graph = daten["graph"]

    for atom in daten["atome"]:

        aid = atom["id"]
        element = atom["element"]

        if aid in aromatisch_unbekannt:
            typen[aid] = _UNBEKANNT
            continue

        if element != "C":
            typen[aid] = element
            continue

        if aid in aromatische_c:
            typen[aid] = "Cb"
            continue

        ordnungen = [(n, _ordnung(daten, aid, n)) for n in graph[aid]]
        mehrfach = [(n, o) for n, o in ordnungen if o >= 2]

        if not mehrfach:
            typen[aid] = "C"
        elif len(mehrfach) > 1:
            # Zwei Doppelbindungen am selben C (Allen, CO2, Keten ...)
            if len(mehrfach) == 2 and all(o == 2 for _, o in mehrfach):
                typen[aid] = "Cdd"
            else:
                typen[aid] = _UNBEKANNT
        else:
            partner, o = mehrfach[0]
            partner_el = daten["atom_by_id"][partner]["element"]

            if o == 3 and partner_el == "C":
                typen[aid] = "Ct"
            elif o == 2 and partner_el == "C":
                typen[aid] = "Cd"
            elif o == 2 and partner_el == "O":
                typen[aid] = "CO"
            else:
                typen[aid] = _UNBEKANNT

    return typen


def _schluessel(zentrum, ligandenliste):
    """Baut den Tabellenschluessel. Gibt None zurueck, wenn ein Ligand nicht
    in _LIGANDEN_REIHENFOLGE steht - frueher wurde so ein Ligand stillschweigend
    weggelassen, und das Atom konnte zufaellig auf einen anderen, gueltigen
    Schluessel fallen."""

    zaehler = Counter(ligandenliste)

    if not set(zaehler) <= set(_LIGANDEN_REIHENFOLGE):
        return None

    teile = []

    for label in _LIGANDEN_REIHENFOLGE:
        n = zaehler.get(label, 0)
        if n:
            teile.append(f"({label}){n if n > 1 else ''}")

    return f"{zentrum}-" + "".join(teile)


def _alternativ_schluessel(zentrum, ligandenliste):
    """Generischer Schluessel fuer Halogengruppen: Cd/Ct/Cb/CO-Liganden -> C.
    Nur sinnvoll, wenn ein Halogen im Spiel ist (sonst None)."""

    if not any(l in _HALOGENE for l in ligandenliste):
        return None

    generisch = ["C" if l in _C_TYPEN else l for l in ligandenliste]

    if generisch == list(ligandenliste):
        return None

    return _schluessel(zentrum, generisch)


def _gruppe_fuer_atom(atom, typen, daten):
    """Gibt (schluessel, fehlergrund, alternativ_schluessel).
    schluessel None und fehlergrund None = Atom bildet keine eigene Gruppe
    (H, Halogen, Carbonyl-O). fehlergrund != None = nicht abdeckbar."""

    aid = atom["id"]
    element = atom["element"]
    typ = typen[aid]
    graph = daten["graph"]

    if typ == _UNBEKANNT:
        return None, f"{element} {aid}: Bindungs-/Ringsituation nicht von Benson-Tabelle abgedeckt", None

    if element == "H":
        return None, None, None

    nachbarn = list(graph[aid])
    ordnungen = {n: _ordnung(daten, aid, n) for n in nachbarn}

    if element in _HALOGENE:

        if len(nachbarn) == 1 and ordnungen[nachbarn[0]] == 1:
            return None, None, None  # der Beitrag steckt in der Gruppe des Nachbarn

        return None, f"{element} {aid}: Halogen mit mehr als einer Bindung nicht abgedeckt", None

    if element == "O":

        if any(o >= 2 for o in ordnungen.values()):
            # Carbonyl-O / Cdd-O gehoert zur Gruppe des Kohlenstoffs; alles andere unbekannt
            if (len(nachbarn) == 1 and ordnungen[nachbarn[0]] == 2
                    and typen[nachbarn[0]] in ("CO", "Cdd")):
                return None, None, None
            return None, f"O {aid}: Mehrfachbindung ausserhalb Carbonyl/Cdd", None

        ligand_typen = [typen[n] for n in nachbarn]

    elif element in ("N", "S"):

        if any(o != 1 for o in ordnungen.values()):
            return None, f"{element} {aid}: Mehrfachbindung nicht abgedeckt", None

        ligand_typen = [typen[n] for n in nachbarn]

    elif element == "C":

        if typ == "Cb":
            cb_nachbarn = [n for n in nachbarn if typen[n] == "Cb"]
            andere = [n for n in nachbarn if typen[n] != "Cb"]
            ligand_typen = [typen[n] for n in andere]
            if len(cb_nachbarn) >= 3:
                ligand_typen.append("Cb")  # Ringverknuepfung (Biphenyl o.ae.)

        elif typ == "Cdd":
            # Liganden sind die beiden Doppelbindungspartner
            ligand_typen = []
            for n in nachbarn:
                partner_el = daten["atom_by_id"][n]["element"]
                ligand_typen.append("Od" if partner_el == "O" else typen[n])

        else:
            ausgeschlossen = set()
            if typ in ("Cd", "Ct", "CO"):
                ausgeschlossen = {n for n in nachbarn if ordnungen[n] >= 2}
            ligand_typen = [typen[n] for n in nachbarn if n not in ausgeschlossen]

    else:
        return None, f"Element {element} nicht in Benson-Tabelle", None

    if _UNBEKANNT in ligand_typen:
        return None, f"{element} {aid}: Nachbaratom nicht abgedeckt", None

    zentrum = typ if element == "C" else element
    schluessel = _schluessel(zentrum, ligand_typen)

    if schluessel is None:
        fremd = sorted(set(ligand_typen) - set(_LIGANDEN_REIHENFOLGE))
        return None, f"{element} {aid}: Nachbar {'/'.join(fremd)} nicht abgedeckt", None

    return schluessel, None, _alternativ_schluessel(zentrum, ligand_typen)


def _ringkorrektur(ring, daten, ist_aromatisch, tabelle):
    """(wert, fehlergrund). wert None bei Fehler."""

    n = len(ring)
    elemente = [daten["atom_by_id"][a]["element"] for a in ring]
    fremd = [e for e in elemente if e != "C"]

    if ist_aromatisch(ring):
        if n == 6 and not fremd:
            return 0.0, None
        return None, f"aromatischer Ring ({n}, {''.join(sorted(set(elemente)))}) nicht abgedeckt"

    if not fremd:
        prefix = "C"
    elif len(fremd) == 1:
        prefix = fremd[0]
    else:
        return None, f"Ring mit {len(fremd)} Heteroatomen nicht abgedeckt"

    doppel = 0
    for i in range(n):
        o = _ordnung(daten, ring[i], ring[(i + 1) % n])
        if o == 3:
            return None, "Dreifachbindung im Ring nicht abgedeckt"
        if o == 2:
            doppel += 1

    schluessel = f"{prefix}:{n}:{doppel}"
    wert = tabelle["ringe"].get(schluessel)

    if wert is None:
        return None, f"Ringkorrektur {schluessel} fehlt"

    return wert, None


def berechne(daten, ringe, ist_aromatisch, tabelle_pfad=None):
    """Berechnet dHf nach Benson.

    daten           Kritiker-daten (atome, graph, bindung_by_paar, atom_by_id)
    ringe           Liste von Ring-Pfaden (Atom-IDs), z.B. Kritiker._alle_ringe
    ist_aromatisch  Funktion ring -> bool
    """

    tabelle = lade_tabelle(tabelle_pfad)

    aromatische_c = set()
    aromatisch_unbekannt = set()

    for ring in ringe:
        if ist_aromatisch(ring):
            elemente = {daten["atom_by_id"][a]["element"] for a in ring}
            if len(ring) == 6 and elemente == {"C"}:
                aromatische_c.update(ring)
            else:
                aromatisch_unbekannt.update(ring)

    typen = _atom_typen(daten, aromatische_c, aromatisch_unbekannt)

    gruppen = Counter()
    nicht_abgedeckt = []
    dhf_gruppen = 0.0

    for atom in daten["atome"]:

        schluessel, fehler, alternativ = _gruppe_fuer_atom(atom, typen, daten)

        if fehler:
            nicht_abgedeckt.append(fehler)
            continue

        if schluessel is None:
            continue

        wert = tabelle["gruppen"].get(schluessel)

        if wert is None and alternativ is not None:
            wert = tabelle["gruppen"].get(alternativ)
            if wert is not None:
                schluessel = alternativ

        if wert is None:
            nicht_abgedeckt.append(f"Gruppe {schluessel} fehlt in der Tabelle")
            continue

        gruppen[schluessel] += 1
        dhf_gruppen += wert

    ring_korrektur = 0.0
    warnungen = []

    for ring in ringe:

        wert, fehler = _ringkorrektur(ring, daten, ist_aromatisch, tabelle)

        if fehler:
            nicht_abgedeckt.append(fehler)
        else:
            ring_korrektur += wert

    # Kondensierte, nicht aromatische Ringsysteme: Summe der Einzelringe ist nur eine Naeherung.
    nicht_arom = [set(r) for r in ringe if not ist_aromatisch(r)]
    for i in range(len(nicht_arom)):
        for j in range(i + 1, len(nicht_arom)):
            if nicht_arom[i] & nicht_arom[j]:
                warnungen.append("Kondensierte/verbrueckte Ringe: Ringkorrektur nur als Summe genaehert")
                break
        else:
            continue
        break

    schwere_atome = sum(1 for a in daten["atome"] if a["element"] != "H")

    return {
        "anwendbar": not nicht_abgedeckt and schwere_atome > 0,
        "dhf_kcal": round(dhf_gruppen + ring_korrektur, 2),
        "dhf_ohne_ring": round(dhf_gruppen, 2),
        "ring_korrektur": round(ring_korrektur, 2),
        "schwere_atome": schwere_atome,
        "gruppen": dict(gruppen),
        "nicht_abgedeckt": nicht_abgedeckt,
        "warnungen": warnungen,
    }