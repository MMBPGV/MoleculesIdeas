# ==========================================================
# lernkurve.py
# ==========================================================
#
# Lernkurve als eigenstaendiges Widget. Die Zeichenmethode ist aus der
# bisherigen GUI uebernommen; einzige Aenderung: sie nutzt die TATSAECH-
# LICHE Canvas-Groesse statt der festen 800x400 (das Diagramm passt sich
# damit dem Fenster an und wird beim Vergroessern neu gezeichnet).

import tkinter as tk
from tkinter import ttk

from chemlabor.gui.theme import FARBEN


class LernkurveAnsicht(ttk.Frame):

    FARBEN = FARBEN

    def __init__(self, master):

        super().__init__(master)

        self.labor = None

        self.lernkurve_status = ttk.Label(
            self,
            text="Noch keine Daten - die Kurve füllt sich nach den ersten Beobachtungen.",
            style="Gedimmt.TLabel"
        )

        self.lernkurve_status.pack(anchor="w", padx=6, pady=(4, 0))

        self.lernkurve_canvas = tk.Canvas(
            self,
            width=600,
            height=260,
            bg=FARBEN["panel_hell"],
            highlightthickness=0
        )

        self.lernkurve_canvas.pack(fill="both", expand=True, padx=6, pady=6)

        self._neuzeichnen_id = None

        self.lernkurve_canvas.bind("<Configure>", self._groesse_geaendert)

    def zeichne(self, labor):

        self.labor = labor

        self._zeichne_lernkurve()

    def _groesse_geaendert(self, _ereignis):

        # Beim Ziehen am Fenster nicht bei jedem Pixel neu zeichnen.
        if self._neuzeichnen_id is not None:
            self.after_cancel(self._neuzeichnen_id)

        self._neuzeichnen_id = self.after(120, self._neuzeichnen)

    def _neuzeichnen(self):

        self._neuzeichnen_id = None

        if self.labor is not None:
            self._zeichne_lernkurve()

    # ------------------------------------------------------------------
    # Ab hier: aus der bisherigen GUI (nur Groessenermittlung angepasst)
    # ------------------------------------------------------------------

    def _zeichne_lernkurve(self):

        if self.labor is None:
            return

        punkte = self.labor.gedaechtnis.lernkurve_daten()

        canvas = self.lernkurve_canvas

        canvas.delete("all")

        if len(punkte) < 2:

            self.lernkurve_status.config(
                text=f"Noch zu wenig Daten für eine Kurve "
                     f"({len(punkte)} Intervall(e) abgeschlossen)."
            )

            return

        _differenz, trend_text = self.labor.gedaechtnis.lerntrend()

        self.lernkurve_status.config(
            text=f"{len(punkte)} Intervalle (je "
                 f"{self.labor.gedaechtnis.lernkurve_intervall} "
                 f"Beobachtungen) - {trend_text}"
        )

        breite = max(canvas.winfo_width(), 300)
        hoehe = max(canvas.winfo_height(), 200)

        rand_links = 45
        rand_unten = 25
        rand_oben = 15
        rand_rechts = 15

        werte = [durchschnitt for _stand, durchschnitt in punkte]

        y_min = min(0, min(werte))
        y_max = max(100, max(werte))

        # Sicherheitsabstand, falls alle Werte identisch sind (sonst
        # Division durch 0 bei der Skalierung unten).
        if y_max == y_min:
            y_max = y_min + 1

        x_min = punkte[0][0]
        x_max = punkte[-1][0]

        if x_max == x_min:
            x_max = x_min + 1

        plot_breite = breite - rand_links - rand_rechts
        plot_hoehe = hoehe - rand_oben - rand_unten

        def x_pos(stand):
            return rand_links + (stand - x_min) / (x_max - x_min) * plot_breite

        def y_pos(wert):
            return (
                rand_oben
                + (1 - (wert - y_min) / (y_max - y_min)) * plot_hoehe
            )

        # Achsen
        canvas.create_line(
            rand_links, rand_oben,
            rand_links, hoehe - rand_unten,
            fill=self.FARBEN["text_gedimmt"]
        )

        canvas.create_line(
            rand_links, hoehe - rand_unten,
            breite - rand_rechts, hoehe - rand_unten,
            fill=self.FARBEN["text_gedimmt"]
        )

        # Ein paar horizontale Hilfslinien mit Beschriftung.
        for anteil in (0, 0.25, 0.5, 0.75, 1.0):

            wert = y_min + anteil * (y_max - y_min)
            y = y_pos(wert)

            canvas.create_line(
                rand_links, y,
                breite - rand_rechts, y,
                fill=self.FARBEN["rand"]
            )

            canvas.create_text(
                rand_links - 8, y,
                text=f"{wert:.0f}",
                anchor="e",
                font=("", 8),
                fill=self.FARBEN["text_gedimmt"]
            )

        # Kurve selbst.
        koordinaten = []

        for stand, durchschnitt in punkte:
            koordinaten.append(x_pos(stand))
            koordinaten.append(y_pos(durchschnitt))

        canvas.create_line(
            *koordinaten,
            fill=self.FARBEN["akzent"],
            width=2
        )

        for stand, durchschnitt in punkte:

            x = x_pos(stand)
            y = y_pos(durchschnitt)

            canvas.create_oval(
                x - 2, y - 2,
                x + 2, y + 2,
                fill=self.FARBEN["akzent"],
                outline=""
            )

        canvas.create_text(
            breite - rand_rechts, hoehe - rand_unten + 12,
            text=f"Beobachtungen: {x_min} – {x_max}",
            anchor="e",
            font=("", 8),
            fill=self.FARBEN["text_gedimmt"]
        )
