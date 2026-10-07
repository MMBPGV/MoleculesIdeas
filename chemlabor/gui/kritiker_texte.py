# ==========================================================
# kritiker_texte.py
# ==========================================================
#
# Der Kritiker (kritiker.py, kritiker_regeln.py) liefert deutsche Texte:
# Urteil, Erklaerung und die Eintraege unter "probleme". Diese Texte
# bleiben im Datenmodell unveraendert (Gedaechtnis, gespeicherte Dateien
# und Export haengen daran) und werden erst bei der ANZEIGE uebersetzt.
#
# Bei deutscher Oberflaeche geben alle Funktionen den Text unveraendert
# zurueck. Bei jeder anderen Sprache wird nachgeschlagen; was nicht
# erkannt wird, bleibt im Original stehen (lieber ein deutscher Satz als
# ein falscher englischer).
#
# WICHTIG: die Satzmuster unten spiegeln die f-Strings im Kritiker. Wer
# dort einen Satz aendert oder ein neues Modell einbaut, muss hier das
# passende Muster nachziehen - tests/test_kritiker_texte.py laesst den
# echten Kritiker laufen und schlaegt an, wenn ein Satz unuebersetzt
# bleibt.

import re

from chemlabor.texte import sprache


URTEILE = {
    "sehr stabil": "very stable",
    "wahrscheinlich stabil": "probably stable",
    "fragwürdig": "questionable",
    "instabil": "unstable",
}

MODELLNAMEN = {
    "Valenz": "Valence",
    "Konnektivität": "Connectivity",
    "Radikal": "Radical",
    "Offene Valenzen": "Open valences",
    "Bindungsordnung": "Bond order",
    "Halogene": "Halogens",
    "Struktur": "Structure",
    "Bindungsenergie": "Bond energy",
    "Formalladung": "Formal charge",
    "Ladungsbilanz": "Charge balance",
    "Ringspannung": "Ring strain",
    "Hypervalenz": "Hypervalence",
    "Kumulen im Ring": "Cumulene in ring",
    "Dreifachbindung im Ring": "Triple bond in ring",
    "Hypohalogenit": "Hypohalite",
    "Schwere Mehrfachbindung": "Heavy multiple bond",
}

PARAMETER = {
    "atom": "atom",
    "atom1": "atom1",
    "atom2": "atom2",
    "atome": "atoms",
    "problem": "problem",
    "ungepaarte_elektronen": "unpaired_electrons",
    "freie_valenzen": "free_valences",
    "abzug": "penalty",
    "durchschnittsgrad": "avg_degree",
    "bindung": "bond",
    "energie_kj_mol": "energy_kj_mol",
    "ladung": "charge",
    "gesamtladung": "total_charge",
    "ringgroesse": "ring_size",
    "benutzt": "used",
    "normale_valenz": "normal_valence",
}

GEOMETRIEN = {
    "linear": "linear",
    "trigonal-planar": "trigonal planar",
    "gewinkelt": "bent",
    "tetraedrisch": "tetrahedral",
    "trigonal-pyramidal": "trigonal pyramidal",
    "trigonal-bipyramidal": "trigonal bipyramidal",
    "wippenförmig (see-saw)": "see-saw",
    "T-förmig": "T-shaped",
    "oktaedrisch": "octahedral",
    "quadratisch-pyramidal": "square pyramidal",
    "quadratisch-planar": "square planar",
}


# (deutsches Muster, englische Vorlage, Gruppe -> Nachschlagetabelle)
_MUSTER = [

    # --- kritiker.py: Erklaerungen ---
    (r"(?P<el>\w+) besitzt mehr Bindungen als erlaubt\.",
     "{el} has more bonds than allowed.", None),

    (r"Das Molekül besteht aus mehreren unabhängigen Teilstrukturen\.",
     "The molecule consists of several independent substructures.", None),

    (r"(?P<el>\w+) besitzt ein ungepaartes Elektron \(Radikalcharakter\) - reaktiv und instabil\.",
     "{el} has an unpaired electron (radical character) - reactive and unstable.", None),

    (r"(?P<el>\w+) besitzt (?P<n>\d+) ungesättigte Bindungsplätze\.",
     "{el} has {n} unsaturated bonding sites.", None),

    (r"(?P<el>\w+) erreicht kein vollständiges Oktett\.",
     "{el} does not reach a complete octet.", None),

    (r"(?P<e1>[A-Z][a-z]?)-(?P<e2>[A-Z][a-z]?)-Bindung ist extrem polar \(ΔEN=(?P<diff>[\d.]+)\)\.",
     "{e1}-{e2} bond is extremely polar (ΔEN={diff}).", None),

    (r"(?P<e1>[A-Z][a-z]?)-(?P<e2>[A-Z][a-z]?)-Bindung ist mit (?P<energie>[\d.]+) kJ/mol ungewöhnlich schwach\.",
     "{e1}-{e2} bond is unusually weak at {energie} kJ/mol.", None),

    (r"(?P<el>\w+) besitzt eine ungewöhnliche Formalladung\.",
     "{el} has an unusual formal charge.", None),

    (r"Summe aller Formalladungen ist (?P<g>-?[\d.]+), nicht 0 - für ein als neutral angenommenes Molekül unrealistisch\.",
     "The sum of all formal charges is {g}, not 0 - unrealistic for a molecule assumed to be neutral.", None),

    (r"(?P<n>\d+) konjugierte\(s\) Doppelbindungspaar\(e\) gefunden - Mesomeriestabilisierung möglich\.",
     "{n} conjugated double-bond pair(s) found - mesomeric stabilization possible.", None),

    (r"(?P<el>\w+) ist (?P<geo>.+) \((?P<hyb>sp\w*), (?P<sigma>\d+) Bindungspartner, (?P<freie>\d+) freie\(s\) Elektronenpaar\(e\)\)\.",
     "{el} is {geo} ({hyb}, {sigma} bonding partner(s), {freie} lone pair(s)).", {"geo": GEOMETRIEN}),

    (r"Ring mit (?P<n>\d+) Atomen gefunden - Ringspannung wahrscheinlich\.",
     "Ring of {n} atoms found - ring strain likely.", None),

    (r"Ring mit (?P<n>\d+) Atomen gefunden - günstige, praktisch spannungsfreie Ringgröße\.",
     "Ring of {n} atoms found - favorable, practically strain-free ring size.", None),

    (r"Ring mit (?P<n>\d+) Atomen ist aromatisch \((?P<pi>\d+) π-Elektronen, erfüllt die 4n\+2-Regel\) - zusätzlich stabilisiert\.",
     "Ring of {n} atoms is aromatic ({pi} π electrons, satisfies the 4n+2 rule) - additionally stabilized.", None),

    (r"(?P<el>\w+) nutzt (?P<b>\d+) Bindungen - mehr als die normale Valenz (?P<nv>\d+) - hypervalent\.",
     "{el} uses {b} bonds - more than the normal valence {nv} - hypervalent.", None),

    # --- kritiker_regeln.py: Erklaerungen ---
    (r"(?P<el>\w+) trägt zwei Doppelbindungen im (?P<n>\d+)er-Ring - will linear sein, im kleinen Ring stark verspannt\.",
     "{el} carries two double bonds in the {n}-membered ring - wants to be linear, strongly strained in a small ring.", None),

    (r"Dreifachbindung im (?P<n>\d+)er-Ring - will linear sein, im kleinen Ring stark verspannt\.",
     "Triple bond in the {n}-membered ring - wants to be linear, strongly strained in a small ring.", None),

    (r"O-(?P<x>[A-Z][a-z]?)-Bindung \(Hypohalogenit\): schwach und stark oxidierend, in der Praxis meist explosiv oder zersetzlich\.",
     "O-{x} bond (hypohalite): weak and strongly oxidizing, in practice usually explosive or prone to decomposition.", None),

    (r"(?P<name>\S+)-Mehrfachbindung ohne Schutzgruppen - reaktiv, real nur mit sperrigen Substituenten isolierbar\.",
     "{name} multiple bond without protecting groups - reactive, in reality only isolable with bulky substituents.", None),

    # --- Werte unter "problem" in den Eintraegen von "probleme" ---
    (r"Valenz überschritten", "valence exceeded", None),
    (r"(?P<n>\d+) Fragmente", "{n} fragments", None),
    (r"H kann keine Mehrfachbindung besitzen", "H cannot have a multiple bond", None),
    (r"Halogen besitzt zu viele Bindungen", "halogen has too many bonds", None),
    (r"Sehr geringe Vernetzung", "very low connectivity", None),
]

_KOMPILIERT = [
    (re.compile(f"^{muster}$"), vorlage, tabellen)
    for muster, vorlage, tabellen in _MUSTER
]


def _deutsch():

    return sprache() == "de"


def uebersetze_urteil(urteil):

    if _deutsch():
        return urteil

    return URTEILE.get(urteil, urteil)


def uebersetze_modell(name):

    if _deutsch():
        return name

    return MODELLNAMEN.get(name, name)


def uebersetze_parameter(schluessel):

    if _deutsch():
        return schluessel

    return PARAMETER.get(schluessel, schluessel)


def uebersetze_text(text):
    """Uebersetzt einen Satz aus Kritiker.erklaerung bzw. einen Text-Wert
    aus Kritiker.probleme. Nicht erkannte Saetze bleiben unveraendert."""

    if _deutsch() or not isinstance(text, str):
        return text

    for muster, vorlage, tabellen in _KOMPILIERT:

        treffer = muster.match(text)

        if treffer is None:
            continue

        gruppen = treffer.groupdict()

        for name, tabelle in (tabellen or {}).items():
            gruppen[name] = tabelle.get(gruppen[name], gruppen[name])

        return vorlage.format(**gruppen)

    return text


def uebersetze_wert(wert):
    """Wert eines Eintrags in 'probleme': Texte uebersetzen, Zahlen und
    Listen unveraendert lassen."""

    return uebersetze_text(wert)