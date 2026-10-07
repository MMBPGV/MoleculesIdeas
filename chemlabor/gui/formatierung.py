# ==========================================================
# formatierung.py
# ==========================================================
#
# Reine Formatierungs-/Darstellungslogik OHNE tkinter-Abhaengigkeit
# (damit sie sich ohne Display testen laesst).

from chemlabor.gui.kritiker_texte import (
    uebersetze_modell,
    uebersetze_parameter,
    uebersetze_wert,
)


def formatiere_probleme(probleme):
    """Wandelt die von Kritiker.probleme gelieferten Dicts (z.B.
    {"modell": "Radikal", "atom": 3, "ungepaarte_elektronen": 1}) in
    lesbare Zeilen fuer die Detailansicht um.

    Modellnamen, Parameternamen und Text-Werte werden bei nicht-
    deutscher Oberflaeche uebersetzt (siehe kritiker_texte.py); bei
    deutscher Oberflaeche ist die Ausgabe wie bisher."""

    zeilen = []

    for p in probleme:

        modell = p.get("modell", "?")

        rest = {k: v for k, v in p.items() if k != "modell"}

        details = ", ".join(
            f"{uebersetze_parameter(k)}={uebersetze_wert(v)}"
            for k, v in rest.items()
        )

        zeilen.append(f" - [{uebersetze_modell(modell)}] {details}")

    return zeilen