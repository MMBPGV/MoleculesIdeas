# ==========================================================
# dialog_einstellungen.py
# ==========================================================
#
# Dialog fuer alle Einstellungen eines Laufs (LaborConfig). Ersetzt die
# dauerhaft sichtbaren Eingabefelder der bisherigen GUI.

import gc
import tkinter as tk
from tkinter import ttk

from chemlabor.config import LaborConfig
from chemlabor.pfade import AUTOSAVE_PFAD


def _ganzzahl(text, name):

    sauber = text.strip().replace(" ", "").replace(".", "").replace("_", "")

    try:
        return int(sauber)
    except ValueError:
        raise ValueError(f"{name}: '{text.strip()}' ist keine ganze Zahl.") from None


class EinstellungenDialog(tk.Toplevel):
    """Ergebnis steht nach dem Schliessen in self.ergebnis: eine neue,
    validierte LaborConfig - oder None bei Abbruch."""

    def __init__(self, master, config):

        super().__init__(master)

        self.title("Einstellungen")
        self.transient(master)
        self.resizable(False, False)

        self.config_alt = config
        self.ergebnis = None

        self.beobachtungen = tk.StringVar()
        self.automatik = tk.StringVar()
        self.kandidaten = tk.StringVar()
        self.max_atome = tk.StringVar()
        self.filter_aktiv = tk.BooleanVar()
        self.verhaeltnis = tk.StringVar()
        self.modus = tk.StringVar()
        self.endlos = tk.BooleanVar()
        self.autosave = tk.BooleanVar()

        self._baue()
        self._fuelle(config)

        self.bind("<Return>", lambda _e: self._uebernehmen())
        self.bind("<Escape>", lambda _e: self.destroy())

        self.protocol("WM_DELETE_WINDOW", self.destroy)

    # ------------------------------------------------------
    # Aufbau
    # ------------------------------------------------------

    def _baue(self):

        rahmen = ttk.Frame(self, padding=14)
        rahmen.pack(fill="both", expand=True)

        zeile = 0

        def feld(text, variable, hinweis):

            nonlocal zeile

            ttk.Label(rahmen, text=text).grid(row=zeile, column=0, sticky="w", pady=3)

            eingabe = ttk.Entry(rahmen, textvariable=variable, width=12)
            eingabe.grid(row=zeile, column=1, sticky="w", padx=(10, 10))

            ttk.Label(rahmen, text=hinweis, style="Gedimmt.TLabel").grid(
                row=zeile, column=2, sticky="w"
            )

            zeile += 1

            return eingabe

        self.erstes_feld = feld(
            "Beobachtungen", self.beobachtungen,
            "Zufallsmoleküle, aus denen die KI zuerst lernt"
        )

        feld("Automatik-Experimente", self.automatik, "danach von der KI selbst erzeugt")
        feld("Kandidaten pro Molekül", self.kandidaten, "Best-of-N im KI-Generator")
        feld("Max. Atome", self.max_atome, "Obergrenze je Molekül")

        ttk.Label(rahmen, text="Verhältnis-Filter").grid(row=zeile, column=0, sticky="w", pady=3)

        filterzeile = ttk.Frame(rahmen)
        filterzeile.grid(row=zeile, column=1, columnspan=2, sticky="w", padx=(10, 0))

        ttk.Checkbutton(filterzeile, text="aktiv, höchstens", variable=self.filter_aktiv).pack(side="left")
        ttk.Entry(filterzeile, textvariable=self.verhaeltnis, width=4).pack(side="left", padx=6)
        ttk.Label(filterzeile, text="schlechte je gutes Molekül", style="Gedimmt.TLabel").pack(side="left")

        zeile += 1

        ttk.Label(rahmen, text="Automatik-Modus").grid(row=zeile, column=0, sticky="nw", pady=3)

        modus = ttk.Frame(rahmen)
        modus.grid(row=zeile, column=1, columnspan=2, sticky="w", padx=(10, 0))

        ttk.Radiobutton(modus, text="Zufall (Best-of-N)", value="zufall", variable=self.modus).pack(anchor="w")
        ttk.Radiobutton(modus, text="Konstruktion (Schritt für Schritt mit RDKit)", value="konstruktion", variable=self.modus).pack(anchor="w")

        zeile += 1

        ttk.Checkbutton(
            rahmen,
            text="Endlos: Automatikphase läuft, bis du auf Stopp drückst",
            variable=self.endlos
        ).grid(row=zeile, column=0, columnspan=3, sticky="w", pady=(8, 0))

        zeile += 1

        ttk.Checkbutton(
            rahmen,
            text="Gedächtnis automatisch sichern",
            variable=self.autosave
        ).grid(row=zeile, column=0, columnspan=3, sticky="w")

        zeile += 1

        ttk.Label(rahmen, text=f"nach {AUTOSAVE_PFAD}", style="Gedimmt.TLabel").grid(
            row=zeile, column=0, columnspan=3, sticky="w", padx=(24, 0)
        )

        zeile += 1

        self.fehler = ttk.Label(rahmen, text="", foreground="#ff8a80")
        self.fehler.grid(row=zeile, column=0, columnspan=3, sticky="w", pady=(10, 0))

        zeile += 1

        knoepfe = ttk.Frame(rahmen)
        knoepfe.grid(row=zeile, column=0, columnspan=3, sticky="e", pady=(12, 0))

        ttk.Button(knoepfe, text="Standardwerte", command=lambda: self._fuelle(LaborConfig())).pack(side="left", padx=(0, 24))
        ttk.Button(knoepfe, text="Abbrechen", command=self.destroy).pack(side="left", padx=(0, 8))
        ttk.Button(knoepfe, text="OK", style="Accent.TButton", command=self._uebernehmen).pack(side="left")

    # ------------------------------------------------------
    # Werte
    # ------------------------------------------------------

    def _fuelle(self, config):

        self.beobachtungen.set(str(config.beobachtungen))
        self.automatik.set(str(config.automatik_experimente))
        self.kandidaten.set(str(config.kandidaten_pro_molekuel))
        self.max_atome.set(str(config.max_atome))

        self.filter_aktiv.set(config.ziel_verhaeltnis_neg_zu_pos is not None)
        self.verhaeltnis.set(
            str(config.ziel_verhaeltnis_neg_zu_pos)
            if config.ziel_verhaeltnis_neg_zu_pos is not None else "3"
        )

        self.modus.set(config.automatik_modus)
        self.endlos.set(config.endlos)
        self.autosave.set(bool(config.autosave_pfad))

        self.fehler.config(text="")

    def _werte(self):

        return {
            "beobachtungen": _ganzzahl(self.beobachtungen.get(), "Beobachtungen"),
            "automatik_experimente": _ganzzahl(self.automatik.get(), "Automatik-Experimente"),
            "kandidaten_pro_molekuel": _ganzzahl(self.kandidaten.get(), "Kandidaten pro Molekül"),
            "max_atome": _ganzzahl(self.max_atome.get(), "Max. Atome"),
            "ziel_verhaeltnis_neg_zu_pos": (
                _ganzzahl(self.verhaeltnis.get(), "Verhältnis-Filter")
                if self.filter_aktiv.get() else None
            ),
            "automatik_modus": self.modus.get(),
            "endlos": bool(self.endlos.get()),
            "autosave_pfad": AUTOSAVE_PFAD if self.autosave.get() else None,
        }

    def _uebernehmen(self):

        try:
            neu = self.config_alt.kopie(**self._werte())
        except ValueError as fehler:
            self.fehler.config(text=str(fehler))
            return

        self.ergebnis = neu

        self.destroy()


def zeige_einstellungen(master, config):
    """Oeffnet den Dialog modal. Rueckgabe: neue LaborConfig oder None."""

    dialog = EinstellungenDialog(master, config)

    dialog.update_idletasks()

    # mittig ueber dem Hauptfenster
    x = master.winfo_rootx() + (master.winfo_width() - dialog.winfo_width()) // 2
    y = master.winfo_rooty() + 80

    dialog.geometry(f"+{max(x, 0)}+{max(y, 0)}")

    dialog.grab_set()
    dialog.erstes_feld.focus_set()

    master.wait_window(dialog)

    ergebnis = dialog.ergebnis

    del dialog

    # Die Tk-Variablen des Dialogs haengen in Referenzzyklen. Wuerden sie
    # erst spaeter vom Muellsammler entsorgt, geschaehe das evtl. in einem
    # Hintergrundthread (z.B. dem Labor) - tkinter meldet dann "main thread
    # is not in main loop". Deshalb hier, im Hauptthread, aufraeumen.
    gc.collect()

    return ergebnis
