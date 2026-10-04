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


KENNZAHLEN = (
    ("experimente", "Experimente"),
    ("mittel", "Ø Punkte"),
    ("stabil", f"Stabil (≥{STABIL_AB})"),
    ("spanne", "Min – Max"),
    ("phasen", "Beobachtung / Automatik"),
)


class TabLive(ttk.Frame):

    def __init__(self, master, modell):

        super().__init__(master)

        self.modell = modell

        raster = ttk.Frame(self, padding=(10, 10, 10, 0))
        raster.pack(fill="x")

        self.werte = {}

        for spalte, (schluessel, titel) in enumerate(KENNZAHLEN):

            karte = ttk.Frame(raster, padding=(0, 0, 24, 0))
            karte.grid(row=0, column=spalte, sticky="w")

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

        experimente = f"{k['gesamt']:,}".replace(",", ".")

        if labor is not None:
            experimente += f"  ({labor.versuche:,} Versuche)".replace(",", ".")

        self.werte["experimente"].config(text=experimente)
        self.werte["mittel"].config(text=f"{k['mittel']:.1f}")
        self.werte["stabil"].config(
            text=f"{k['stabil']} ({k['stabil'] / k['ausgewertet'] * 100:.0f} %)"
        )
        self.werte["spanne"].config(text=f"{k['minimum']} – {k['maximum']}")
        self.werte["phasen"].config(text=f"{k['beobachtung']} / {k['automatik']}")

        if k["ausgewertet"] < k["gesamt"]:
            self.werte["mittel"].config(text=f"{k['mittel']:.1f}*")

        self.filter_label.config(
            text="Verhältnis-Filter: " + filter_text(labor)
            + ("   (* Durchschnitt über die letzten " + str(k["ausgewertet"]) + ")"
               if k["ausgewertet"] < k["gesamt"] else "")
        )

        if labor is not None:

            _differenz, trend = labor.gedaechtnis.lerntrend()

            self.trend_label.config(text="Lerntrend: " + trend)

            self.lernkurve.zeichne(labor)
