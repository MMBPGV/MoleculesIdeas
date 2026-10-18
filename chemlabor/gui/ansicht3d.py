"""3D-Ansicht fuer Molekuele des Labors (matplotlib, keine weiteren Abhaengigkeiten ausser numpy).

Eingabe: das Labor-Format {"atoms": [{id, element, ...}], "bonds": [{atom1, atom2, order}]} und
optional die Kritiker-Analyse. Atome und Bindungen, die in analyse["probleme"] vorkommen, werden rot
eingefaerbt, je groesser der Abzug, desto kraeftiger.

Koordinaten: zuerst RDKit (ETKDG), das liefert echte Geometrie. Scheitert RDKit (und das tut es bei
genau den unsinnigen Molekuelen, die der Kritiker bewerten soll), springt ein eigenes Federmodell ein,
damit IMMER etwas angezeigt wird.

Einbau in Tkinter:   ansicht = Ansicht3D(parent); ansicht.pack(fill="both", expand=True)
                     ansicht.zeige(molekuel, analyse)
"""

import queue
import threading

import numpy as np

# Farben (CPK) und Radien (relativ, nur fuer die Darstellung)
FARBEN = {"H": "#d9d9d9", "C": "#6b6b6b", "N": "#3050f8", "O": "#e8231f", "F": "#5fd35f", "Cl": "#1ff01f",
          "Br": "#a62929", "I": "#940094", "S": "#e6c800", "P": "#ff8000", "Si": "#b09070", "B": "#ffb5b5"}
RADIEN = {"H": 0.35, "C": 0.7, "N": 0.65, "O": 0.62, "F": 0.55, "Cl": 0.95, "Br": 1.1, "I": 1.3,
          "S": 1.0, "P": 1.0, "Si": 1.1, "B": 0.8}
STANDARD_FARBE, STANDARD_RADIUS = "#c060c0", 0.8

# Ab so vielen Atomen kein RDKit-Embedding mehr (dauert zu lange fuers GUI, Schemadarstellung genuegt)
MAX_ATOME_RDKIT = 80


# ---------------------------------------------------------------- Koordinaten

def _rdkit_koordinaten(molekuel):
    """Echte 3D-Geometrie oder None. Kein Sanitize noetig/erlaubt: ungueltige Molekuele sollen nicht scheitern."""
    if len(molekuel["atoms"]) > MAX_ATOME_RDKIT:
        return None

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None

    try:
        rw = Chem.RWMol()
        index = {}
        for a in molekuel["atoms"]:
            index[a["id"]] = rw.AddAtom(Chem.Atom(a["element"]))
        typen = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
        for b in molekuel["bonds"]:
            rw.AddBond(index[b["atom1"]], index[b["atom2"]], typen.get(b["order"], Chem.BondType.SINGLE))
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        parameter = AllChem.ETKDGv3()
        parameter.randomSeed = 7
        parameter.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, parameter) != 0:
            return None
        pos = mol.GetConformer().GetPositions()
        if not np.all(np.isfinite(pos)):
            return None
        return {a["id"]: pos[index[a["id"]]] for a in molekuel["atoms"]}
    except Exception:
        return None


def _feder_koordinaten(molekuel, schritte=None, seed=7):
    """Notloesung: Bindungen als Federn (Laenge 1.5), alle Atome stossen sich ab."""
    rng = np.random.default_rng(seed)
    ids = [a["id"] for a in molekuel["atoms"]]
    n = len(ids)
    if schritte is None:
        schritte = 400 if n <= 60 else 150
    nr = {i: k for k, i in enumerate(ids)}
    pos = rng.normal(size=(n, 3))
    kanten = [(nr[b["atom1"]], nr[b["atom2"]]) for b in molekuel["bonds"]]

    for _ in range(schritte):
        d = pos[:, None, :] - pos[None, :, :]
        r2 = (d ** 2).sum(-1) + np.eye(n)
        kraft = (d / r2[..., None] ** 1.5).sum(1) * 0.25
        for i, j in kanten:
            v = pos[j] - pos[i]
            l = np.linalg.norm(v) + 1e-9
            f = (l - 1.5) * v / l * 0.2
            kraft[i] += f
            kraft[j] -= f
        pos += np.clip(kraft, -0.3, 0.3)

    return {i: pos[nr[i]] for i in ids}


def koordinaten(molekuel):
    """(Koordinaten, quelle) mit quelle 'rdkit' oder 'feder'."""
    k = _rdkit_koordinaten(molekuel)
    if k is not None:
        return k, "rdkit"
    return _feder_koordinaten(molekuel), "feder"


# ---------------------------------------------------------------- Probleme aus der Analyse

_ATOM_SCHLUESSEL = ("atom", "atom1", "atom2")
_LISTEN_SCHLUESSEL = ("atome", "ring", "ring_atome", "pfad")


def problem_gewichte(analyse):
    """({atom_id: abzug}, {frozenset(a, b): abzug}) aus analyse["probleme"]."""
    atome, bindungen = {}, {}
    if not analyse:
        return atome, bindungen

    for p in analyse.get("probleme", []):
        if not isinstance(p, dict):
            continue
        abzug = abs(p.get("abzug", 0) or 0)
        ids = [p[k] for k in _ATOM_SCHLUESSEL if isinstance(p.get(k), int)]
        for k in _LISTEN_SCHLUESSEL:
            if isinstance(p.get(k), (list, tuple)):
                ids += [x for x in p[k] if isinstance(x, int)]
        for i in ids:
            atome[i] = atome.get(i, 0) + abzug
        if isinstance(p.get("atom1"), int) and isinstance(p.get("atom2"), int):
            paar = frozenset([p["atom1"], p["atom2"]])
            bindungen[paar] = bindungen.get(paar, 0) + abzug

    return atome, bindungen


def _rot(abzug, maximum):
    """Hellrot (kleiner Abzug) bis Dunkelrot (grosser Abzug) als RGB-Tupel."""
    t = min(1.0, abzug / max(maximum, 1.0))
    return (1.0 - 0.35 * t, 0.55 * (1 - t), 0.55 * (1 - t))


# ---------------------------------------------------------------- Zeichnen

def zeichne(ax, molekuel, analyse=None, beschriftung=True, geometrie=None):
    """Zeichnet das Molekuel in eine matplotlib-3D-Achse. Gibt die Koordinatenquelle zurueck.
    geometrie: fertige Koordinaten {atom_id: xyz} (z.B. von quantenchemie), sonst RDKit/Federmodell."""
    ax.clear()
    if geometrie is not None:
        pos, quelle = {i: np.asarray(p, float) for i, p in geometrie.items()}, "ani"
    else:
        pos, quelle = koordinaten(molekuel)
    gew_atom, gew_bindung = problem_gewichte(analyse)
    maximum = max(list(gew_atom.values()) + list(gew_bindung.values()) + [1.0])

    el = {a["id"]: a["element"] for a in molekuel["atoms"]}

    for b in molekuel["bonds"]:
        p1, p2 = pos[b["atom1"]], pos[b["atom2"]]
        paar = frozenset([b["atom1"], b["atom2"]])
        farbe = _rot(gew_bindung[paar], maximum) if paar in gew_bindung else "#808080"
        breite = 1.5 + 2.0 * b["order"] if paar in gew_bindung else 1.5 + 1.5 * b["order"]
        # Mehrfachbindungen als parallele Linien, sonst sieht man sie nicht
        v = p2 - p1
        senkr = np.cross(v, [0.3, 0.5, 0.8])
        senkr = senkr / (np.linalg.norm(senkr) + 1e-9) * 0.09
        for k in range(b["order"]):
            versatz = senkr * (k - (b["order"] - 1) / 2)
            ax.plot(*zip(p1 + versatz, p2 + versatz), color=farbe, linewidth=breite / b["order"] + 0.8,
                    solid_capstyle="round")

    for i, p in pos.items():
        e = el[i]
        farbe = _rot(gew_atom[i], maximum) if i in gew_atom else FARBEN.get(e, STANDARD_FARBE)
        flaeche = 180 * (RADIEN.get(e, STANDARD_RADIUS) ** 2) * 3
        ax.scatter(*p, s=flaeche, c=[farbe], edgecolors="black" if i in gew_atom else "#202020",
                   linewidths=1.6 if i in gew_atom else 0.6, depthshade=True)
        if beschriftung and e != "H":
            t = ax.text(*p, e, color="white" if e in ("C", "N", "O", "Si") else "black",
                        ha="center", va="center", fontsize=8, fontweight="bold")
            t.set_zorder(1000)

    alle = np.array(list(pos.values()))
    mitte = alle.mean(0)
    halb = max(1.0, np.abs(alle - mitte).max()) * 1.15
    ax.set_xlim(mitte[0] - halb, mitte[0] + halb)
    ax.set_ylim(mitte[1] - halb, mitte[1] + halb)
    ax.set_zlim(mitte[2] - halb, mitte[2] + halb)
    ax.set_axis_off()
    try:
        ax.set_box_aspect((1, 1, 1))
    except AttributeError:
        pass
    return quelle


def _kopfzeile(molekuel, analyse, quelle):
    summe = {}
    for a in molekuel["atoms"]:
        summe[a["element"]] = summe.get(a["element"], 0) + 1
    formel = "".join(e + (str(n) if n > 1 else "") for e, n in sorted(summe.items()))
    text = formel
    if analyse and "bewertung" in analyse:
        text += f"   {analyse['bewertung'].get('punkte', '?')} Punkte"
    if quelle == "feder":
        text += "   (Schemadarstellung, keine echte Geometrie)"
    elif quelle == "ani":
        text += "   (Geometrie optimiert mit ANI-2x)"
    return text


# ---------------------------------------------------------------- Quantenchemie (optional)

_QC = None


def _qc_verfuegbar():
    from chemlabor.chemie import quantenchemie
    return quantenchemie.verfuegbar()


def _qc():
    global _QC
    if _QC is None:
        from chemlabor.chemie import quantenchemie
        _QC = quantenchemie.Quantenchemie()
    return _QC


def _ergebnis_text(r):
    fehler = r.get("fehler_kcal")
    fehlertext = f"Fehler typisch +-{fehler:.0f}, bei exotischen Molekuelen mehr" if fehler else \
        "Fehler nur grob abschaetzbar (+-10 bis 30)"
    zeile = (f"{r.get('modell') or 'ANI-2x'}: Bildungsenthalpie ca. {r['dhf_kcal']:+.0f} kcal/mol ({fehlertext}), "
             f"{r['dauer_s']} s. ")

    if r["intakt"]:
        zeile += "Struktur bleibt nach der Optimierung erhalten."
    else:
        teile = []
        if r["gebrochen"]:
            teile.append(f"{len(r['gebrochen'])} Bindung(en) brechen auseinander")
        if r["neu"]:
            teile.append(f"{len(r['neu'])} neue Bindung(en) entstehen")
        zeile += "Struktur ist KEIN stabiles Minimum: " + ", ".join(teile) + " (Hinweis, kein Beweis: ANI kennt Hypervalenz und Bindungsbruch nur begrenzt)."

    return zeile


# ---------------------------------------------------------------- Tkinter-Widget

class Ansicht3D:
    """Tkinter-Widget: drehbar mit der Maus, Rad zum Zoomen. Wird erst beim Erzeugen importiert,
    damit das Modul ohne Tk/Display testbar bleibt."""

    def __new__(cls, parent, hintergrund=None, vordergrund=None):
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        class _Rahmen(ttk.Frame):
            def __init__(self, parent):
                super().__init__(parent)
                self.figur = Figure(figsize=(5, 5), dpi=100, facecolor=hintergrund)
                self.ax = self.figur.add_subplot(111, projection="3d")
                if hintergrund:
                    self.ax.set_facecolor(hintergrund)

                kopf = ttk.Frame(self)
                kopf.pack(fill="x", padx=6, pady=(4, 0))

                self.titel = ttk.Label(kopf, text="Kein Molekuel gewaehlt")
                self.titel.pack(side="left")

                self.knopf = ttk.Button(kopf, text="Berechnen (ANI-2x)", command=self.berechne)
                self.knopf.pack(side="right")

                self.ergebnis = ttk.Label(self, text="", wraplength=520, justify="left")
                self.ergebnis.pack(fill="x", padx=6)

                self.leinwand = FigureCanvasTkAgg(self.figur, master=self)
                self.leinwand.get_tk_widget().pack(fill="both", expand=True)
                self.molekuel = None
                self.analyse = None
                self._laeuft = False
                self._queue = queue.Queue()

            def zeige(self, molekuel, analyse=None, geometrie=None, ergebnis_text=""):
                self.molekuel = molekuel
                self.analyse = analyse
                quelle = zeichne(self.ax, molekuel, analyse, geometrie=geometrie)
                self.titel.config(text=_kopfzeile(molekuel, analyse, quelle))
                self.ergebnis.config(text=ergebnis_text)
                self.leinwand.draw_idle()

            # ---- Quantenchemie (ANI-2x) im Hintergrundthread

            def berechne(self):
                if self.molekuel is None or self._laeuft:
                    return

                ok, hinweis = _qc_verfuegbar()
                if not ok:
                    self.ergebnis.config(text="Quantenchemie nicht verfuegbar: " + hinweis)
                    return

                self._laeuft = True
                self.knopf.config(state="disabled")
                self.ergebnis.config(text="Rechne ... (beim ersten Mal laedt ANI-2x, ca. 5-10 s)")
                molekuel, analyse = self.molekuel, self.analyse

                def arbeit():
                    try:
                        self._queue.put(("fertig", molekuel, analyse, _qc().berechne(molekuel)))
                    except Exception as fehler:
                        self._queue.put(("fehler", f"{type(fehler).__name__}: {fehler}"))

                threading.Thread(target=arbeit, daemon=True).start()
                self.after(200, self._pruefe)

            def _pruefe(self):
                try:
                    nachricht = self._queue.get_nowait()
                except queue.Empty:
                    self.after(200, self._pruefe)
                    return

                self._laeuft = False
                self.knopf.config(state="normal")

                if nachricht[0] == "fehler":
                    self.ergebnis.config(text="Berechnung fehlgeschlagen: " + nachricht[1])
                    return

                _art, molekuel, analyse, r = nachricht

                if molekuel is not self.molekuel:      # inzwischen anderes Molekuel gewaehlt
                    self.ergebnis.config(text="")
                    return

                if not r["anwendbar"]:
                    self.ergebnis.config(text="Nicht berechenbar: " + r["grund"])
                    return

                self.zeige(molekuel, analyse, geometrie=r["geometrie"], ergebnis_text=_ergebnis_text(r))

        return _Rahmen(parent)