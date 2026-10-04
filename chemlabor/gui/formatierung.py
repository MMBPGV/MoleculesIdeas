# ==========================================================
# labor_format.py
# ==========================================================
#
# Reine Formatierungs-/Darstellungslogik fuer das Chemie-Labor, OHNE
# jede Tkinter-Abhaengigkeit.
#
# Frueher lag formatiere_probleme() zwar schon ausserhalb der
# LaborGUI-Klasse in labor_gui.py, aber im selben MODUL wie
# "import tkinter" - jeder Import dieser Funktion (z.B. aus einem
# Test) zog damit tkinter als harte Abhaengigkeit mit, obwohl die
# Funktion selbst keins braucht. In einer Umgebung ohne Display/
# Tkinter (z.B. viele CI-Systeme) liess sie sich dadurch gar nicht
# importieren, geschweige denn testen - genau die Faehigkeit, die
# durch das Auslagern eigentlich hergestellt werden sollte.
#
# labor_gui.py importiert diese Funktion jetzt von hier, statt sie
# selbst zu definieren.
# ==========================================================


def formatiere_probleme(probleme):
    """Wandelt die von Kritiker.probleme gelieferten Dicts (z.B.
    {"modell": "Radikal", "atom": 3, "ungepaarte_elektronen": 1}) in
    lesbare Zeilen fuer die Detailansicht um.

    Hier sass der urspruengliche Crash-Bug: analyse["probleme"]
    enthaelt Dicts, keine Strings - " - " + p warf einen TypeError bei
    jedem Molekuel mit mindestens einem Problem."""

    zeilen = []

    for p in probleme:

        modell = p.get("modell", "?")

        rest = {k: v for k, v in p.items() if k != "modell"}

        details = ", ".join(f"{k}={v}" for k, v in rest.items())

        zeilen.append(f" - [{modell}] {details}")

    return zeilen