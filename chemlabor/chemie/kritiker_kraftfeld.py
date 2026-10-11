# ==========================================================
# kritiker_kraftfeld.py
# ==========================================================
#
# Kritiker mit zusaetzlichem Geometrie-Abzug aus dem Kraftfeld (siehe
# kraftfeld.py). Der normale Kritiker bleibt unveraendert; diese Klasse
# haengt nach der normalen Analyse hoechstens einen Abzug an und vermerkt
# ihn unter analyse["probleme"] (Modell "Kraftfeld").

from chemlabor.chemie import kraftfeld
from chemlabor.chemie.kritiker import Kritiker


def _urteil(punkte):
    """Dieselben Schwellen wie Kritiker.stabilitaets_test()."""

    if punkte >= 90:
        return "sehr stabil"
    if punkte >= 75:
        return "wahrscheinlich stabil"
    if punkte >= 50:
        return "fragwürdig"
    return "instabil"


class KritikerKraftfeld(Kritiker):

    def analysiere(self, molekuel):

        daten = super().analysiere(molekuel)

        abzug, status, verbiegung = kraftfeld.abzug(molekuel)

        if abzug:

            punkte = max(0, daten["bewertung"]["punkte"] - abzug)

            daten["bewertung"] = {"punkte": punkte, "urteil": _urteil(punkte)}

            daten["probleme"].append({
                "modell": "Kraftfeld",
                "status": status,
                "verbiegung_kcal_pro_atom": round(verbiegung, 1),
                "abzug": abzug,
            })

        return daten