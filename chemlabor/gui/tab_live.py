# ==========================================================
# tab_live.py
# ==========================================================
#
# Tab "Live": die wichtigsten Kennzahlen des laufenden Labors plus
# Lernkurve. Ersetzt die Tabs Statistik und Lernkurve der bisherigen GUI
# (das Histogramm entfaellt).

from tkinter import ttk

from chemlabor.gui.lernkurve import LernkurveAnsicht
from chemlabor.gui.modell import STABIL_AB, filter_text
from chemlabor.texte import t, zahl


# Nur die Schluessel - die Titel werden erst beim Bauen des Tabs
# uebersetzt (nicht beim Import, sonst waere die Sprache zu frueh fest).
KENNZAHLEN = ("experimente", "mittel", "stabil", "spanne", "phasen")


class TabLive(ttk.Frame):

    def __init__(self, master, modell):

        super().__init__(master)

        self.modell = modell

        raster = ttk.Frame(self, padding=(10, 10, 10, 0))
        raster.pack(fill="x")

        self.werte = {}

        for spalte, schluessel in enumerate(KENNZAHLEN):

            karte = ttk.Frame(raster, padding=(0, 0, 24, 0))
            karte.grid(row=0, column=spalte, sticky="w")

            titel = t(f"live.{schluessel}", ab=STABIL_AB)

            ttk.Label(karte, text=titel, style="Gedimmt.TLabel").pack(anchor="w")

            wert = ttk.Label(karte, text="–", style="Kennzahl.TLabel")
            wert.pack(anchor="w")

            self.werte[schluessel] = wert

        self.filter_label = ttk.Label(self, text="", style="Gedimmt.TLabel", padding=(10, 8, 10, 0))
        self.filter_label.pack(anchor="w")

        self.trend_label = ttk.Label(self, text="", padding=(10, 2, 10, 0))
        self.trend_label.pack(anchor="w")

        self.lernkurve = LernkurveAnsicht(self)
        self.lernkurve.pack(fill="both", expand=True, pady=(6, 0))

    def leeren(self):

        for wert in self.werte.values():
            wert.config(text="–")

        self.filter_label.config(text="")
        self.trend_label.config(text="")

        self.lernkurve.labor = None
        self.lernkurve.lernkurve_canvas.delete("all")

    def aktualisiere(self, labor):

        k = self.modell.kennzahlen()

        if k is None:
            return

        ausschnitt = k["ausgewertet"] < k["gesamt"]

        experimente = zahl(k["gesamt"])

        if labor is not None:
            experimente += "  " + t("live.versuche", n=zahl(labor.versuche))

        self.werte["experimente"].config(text=experimente)
        self.werte["mittel"].config(
            text=f"{k['mittel']:.1f}" + ("*" if ausschnitt else "")
        )
        self.werte["stabil"].config(
            text=t(
                "live.stabil_wert",
                n=k["stabil"],
                p=f"{k['stabil'] / k['ausgewertet'] * 100:.0f}"
            )
        )
        self.werte["spanne"].config(text=f"{k['minimum']} – {k['maximum']}")
        self.werte["phasen"].config(text=f"{k['beobachtung']} / {k['automatik']}")

        filter_zeile = t("live.filter", text=filter_text(labor))

        if ausschnitt:
            filter_zeile += "   " + t("live.ausschnitt", n=k["ausgewertet"])

        self.filter_label.config(text=filter_zeile)

        if labor is not None:

            _differenz, trend = labor.gedaechtnis.lerntrend()

            self.trend_label.config(text=t("live.trend", text=trend))

            self.lernkurve.zeichne(labor)