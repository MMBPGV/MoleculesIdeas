# ==========================================================
# layout.py
# ==========================================================
#
# Kleine, tkinter-freie Hilfsfunktionen fuer die Oberflaeche (testbar
# ohne Display).
#
# Hintergrund: Auf Windows mit Anzeigeskalierung (125 %/150 %) oder auf
# kleinen Laptop-Bildschirmen ist ein fest verdrahtetes 1200x780-Fenster
# hoeher als der nutzbare Bildschirm - der untere Teil (Log mit den
# Buttons "Log" und "Leeren") liegt dann ausserhalb und wirkt, als sei
# die Oberflaeche nicht vollstaendig geladen.


def fenstergroesse(
        bildschirm_b,
        bildschirm_h,
        wunsch_b=1200,
        wunsch_h=780,
        rand_b=60,
        rand_h=110
):
    """Gewuenschte Fenstergroesse, aber nie groesser als der Bildschirm
    abzueglich Rand (Titelleiste, Taskleiste). Rueckgabe: (breite, hoehe)."""

    breite = max(640, min(wunsch_b, bildschirm_b - rand_b))
    hoehe = max(480, min(wunsch_h, bildschirm_h - rand_h))

    return breite, hoehe


def mindestgroesse(breite, hoehe):
    """Mindestgroesse des Fensters - nie groesser als die tatsaechliche
    Startgroesse, damit das Fenster auf kleinen Bildschirmen nicht schon
    beim Start unter seinem eigenen Minimum liegt."""

    return min(920, breite), min(560, hoehe)


def geometrie(bildschirm_b, bildschirm_h, breite, hoehe):
    """tk-Geometriestring, mittig auf dem Bildschirm (etwas nach oben
    versetzt wegen der Taskleiste)."""

    x = max(0, (bildschirm_b - breite) // 2)
    y = max(0, (bildschirm_h - hoehe) // 2 - 20)

    return f"{breite}x{hoehe}+{x}+{y}"


def adaptiver_abstand(letzte_dauer, mindest, faktor=10, maximum=3.0):
    """Mindestabstand bis zum naechsten Neuzeichnen: mindestens 'mindest'
    Sekunden, aber so gross, dass das Zeichnen hoechstens etwa 1/faktor
    der Zeit belegt. Auf einem langsamen Rechner (oder bei vielen Zeilen)
    waechst der Abstand von selbst, auf einem schnellen bleibt er klein."""

    return min(maximum, max(mindest, letzte_dauer * faktor))