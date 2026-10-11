# ==========================================================
# widerspruch.py
# ==========================================================
#
# Widersprueche zwischen dem eigenen Kritiker (Heuristik) und RDKit
# (Sanitize-Regeln der bekannten Chemie). Reine Funktionen, ohne
# Abhaengigkeiten - gut testbar.
#
# Zwei Faelle:
#
#   "milde"  - Kritiker findet das Molekuel gut (>= hoch), RDKit lehnt
#              es ab. Entweder hat der Kritiker eine Luecke - oder es
#              ist ungewoehnliche Chemie, die RDKit nicht kennt.
#   "streng" - Kritiker findet es schlecht (< niedrig), RDKit akzeptiert
#              es. Meist ein Hinweis, dass der Kritiker zu streng ist.
#
# Wie stark das in die BEWERTUNG eingeht, ist bewusst minimal:
#
#   - Die angezeigte Punktzahl bleibt immer die des Kritikers.
#   - Nur das LERNSIGNAL fuer das Gedaechtnis wird angepasst, und nur
#     im Fall "milde" um wenige Punkte (Standard 3, hoechstens 10).
#   - Der Fall "streng" und ein von RDKit akzeptiertes Molekuel aendern
#     NIE etwas: RDKit kennt nur die bekannte Chemie - ein Bonus dafuer
#     wuerde die KI darauf trimmen, nur Bekanntes nachzubauen.
#
# Mit abzug=0 ist die Wirkung komplett aus; die Liste bleibt trotzdem
# als reine Anzeige nutzbar.

KRITIKER_HOCH = 75
KRITIKER_NIEDRIG = 50

# Harte Obergrenze fuer den Lernabzug (siehe config.MAX_RDKIT_LERNABZUG)
OBERGRENZE_LERNABZUG = 10


def klassifiziere(punkte, rdkit_gueltig, hoch=KRITIKER_HOCH, niedrig=KRITIKER_NIEDRIG):
    """'milde', 'streng' oder None (kein Widerspruch). rdkit_gueltig=None
    (unbekannt, z.B. aus einer aelteren Datei) ergibt nie einen
    Widerspruch."""

    if rdkit_gueltig is None:
        return None

    if punkte >= hoch and not rdkit_gueltig:
        return "milde"

    if punkte < niedrig and rdkit_gueltig:
        return "streng"

    return None


def lernpunkte(punkte, rdkit_gueltig, abzug, hoch=KRITIKER_HOCH):
    """Punktzahl fuer das Lernsignal. Aendert sich nur im Fall 'milde'
    (Kritiker >= hoch, RDKit lehnt ab) und dort um hoechstens
    OBERGRENZE_LERNABZUG Punkte; sonst Rueckgabe = punkte."""

    if abzug <= 0 or rdkit_gueltig is not False or punkte < hoch:
        return punkte

    return punkte - min(abzug, OBERGRENZE_LERNABZUG)