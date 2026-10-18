"""Echte Rechnung statt Heuristik: Geometrieoptimierung und Bildungsenthalpie mit neuronalen Potenzialen.

Drei Rechenmodelle, automatisch nach den vorkommenden Elementen gewaehlt:
  ANI-2x    (torchani)  H, C, N, O, F, S, Cl            - wB97X/6-31G*, am genauesten fuer diese Elemente
  MACE-OFF23 (mace)     H, C, N, O, F, P, S, Cl, Br, I  - organische Chemie inkl. P, Br, I
  MACE-MP-0 (mace)      u.a. Si, B (+ H, C, N, O, F, Cl) - PBE-Festkoerperdaten, fuer Molekuele UNGENAUER
Nur neutrale Molekuele mit geschlossener Schale. Alles andere (Radikale, Carbene, Elemente ausserhalb der
Modelle) meldet das Modul als "nicht anwendbar" mit Grund, statt eine Zahl zu erfinden.

Installation (alles pip, Windows ok):
    pip install torch numpy scipy h5py tqdm huggingface_hub rdkit
    pip install torchani --no-deps       (--no-deps: sonst zieht torchani zarr, das unter Windows nicht baut)
    pip install mace-torch               (nur fuer P, Br, I, Si, B; sonst reicht ANI-2x)
Gewichte werden beim ersten Start einmalig heruntergeladen.

Bildungsenthalpie (kcal/mol, Gasphase, 298 K):  dHf = E_Modell * Umrechnung + Summe n_Element * c_Element
Die Konstanten (KALIBRIERUNG je Modell) sind an Molekuele mit bekannter Bildungsenthalpie (NIST, aus dem
Gedaechtnis) angepasst; der Leave-one-out-Fehler steht in KALIBRIERUNG_FEHLER_KCAL. E ist die rein elektronische
Energie (keine Nullpunktsenergie/Thermik) - deshalb die Kalibrierung. Aussagekraft fuer die Stabilitaet: nur im
Bereich einiger kcal/mol (ANI) bzw. deutlich groeber (MACE, besonders MACE-MP bei Si/B).
"""

import math
import os
import time

import numpy as np

os.environ.setdefault("TORCHANI_NO_WARN_EXTENSIONS", "1")

HARTREE_KCAL = 627.5095

ANI_ELEMENTE = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16, "Cl": 17}

# Kovalenzradien (Angstrom), nur zur Pruefung, ob die Bindungen nach der Optimierung noch da sind
RADIEN = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05, "Cl": 1.02,
          "P": 1.07, "Br": 1.20, "I": 1.39, "Si": 1.11, "B": 0.84}

ORDNUNGSZAHL = {"H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16, "Cl": 17,
                "Br": 35, "I": 53}
OFF_ELEMENTE = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}
MP_ELEMENTE = {"H", "C", "N", "O", "F", "Cl", "Si", "B"}
EV_KCAL = 23.0605

# erlaubte Gesamtvalenz je Element (Summe der Bindungsordnungen); Schwefel auch hypervalent
NORMALE_VALENZ = {"H": (1,), "C": (4,), "N": (3,), "O": (2,), "F": (1,), "Cl": (1,), "S": (2, 4, 6),
                  "P": (3, 5), "Br": (1,), "I": (1,), "Si": (4,), "B": (3,)}

BINDUNG_GEBROCHEN_FAKTOR = 1.35   # Abstand > Faktor * Radiensumme: Bindung ist weg
BINDUNG_NEU_FAKTOR = 1.15         # nicht gebundenes Paar naeher als Faktor * Radiensumme: neue Bindung

# Atomkonstanten (kcal/mol) je Modell; werden von tools/kalibrierung_mace.py bestimmt (siehe unten)
KALIBRIERUNG_ANI = {"H": 371.440, "C": 23914.593, "N": 34355.764, "O": 47146.624, "F": 62581.126,
                    "Cl": 288754.381, "S": 249856.882}
KALIBRIERUNG_OFF = {"H": 370.329, "C": 23939.351, "N": 34393.426, "O": 47193.704, "F": 62634.764,
                    "P": 214246.966, "S": 249916.041, "Cl": 288808.862, "Br": 1615305.405, "I": 186872.785}
KALIBRIERUNG_MP = {"H": 80.66, "C": 214.168, "N": 195.472, "O": 107.299, "F": 36.706, "Cl": 40.115,
                   "Si": 117.35, "B": 154.507}
KALIBRIERUNG = {"ANI-2x": KALIBRIERUNG_ANI, "MACE-OFF23": KALIBRIERUNG_OFF, "MACE-MP-0": KALIBRIERUNG_MP}
KALIBRIERUNG_FEHLER_KCAL = {"ANI-2x": 4.6, "MACE-OFF23": 4.0, "MACE-MP-0": 8.0}
# OFF: Leave-one-out 3.7 (n=54); MP: LOO 4.6 (n=43), aber Si nur an 5 und B an 4 Molekuelen - eher 8 ansetzen


def waehle_modell(elemente):
    """Name des Rechenmodells fuer eine Elementmenge oder None."""
    el = set(elemente)
    if el <= set(ANI_ELEMENTE):
        return "ANI-2x"
    if el <= OFF_ELEMENTE and all(e in KALIBRIERUNG_OFF for e in el):
        return "MACE-OFF23"
    if el <= MP_ELEMENTE and all(e in KALIBRIERUNG_MP for e in el):
        return "MACE-MP-0"
    return None


def verfuegbar(modell="ANI-2x"):
    """(True, "") oder (False, Hinweis zur Installation)."""
    try:
        import scipy  # noqa: F401
        import torch  # noqa: F401
        from rdkit import Chem  # noqa: F401
        if modell == "ANI-2x":
            import torchani  # noqa: F401
        else:
            import ase  # noqa: F401
            import mace  # noqa: F401
    except ImportError as fehler:
        if modell == "ANI-2x":
            return False, (f"{fehler}. Installation: pip install torch numpy scipy h5py tqdm huggingface_hub rdkit "
                           f"und danach pip install torchani --no-deps")
        return False, f"{fehler}. Fuer {modell}: pip install mace-torch"
    return True, ""


def _rdkit_start(molekuel):
    """Startkoordinaten {atom_id: array} per RDKit-ETKDG oder None."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        rw = Chem.RWMol()
        index = {}
        for a in molekuel["atoms"]:
            index[a["id"]] = rw.AddAtom(Chem.Atom(a["element"]))
        typen = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
        for b in molekuel["bonds"]:
            rw.AddBond(index[b["atom1"]], index[b["atom2"]], typen.get(b["order"], Chem.BondType.SINGLE))
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        p = AllChem.ETKDGv3()
        p.randomSeed = 7
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p) != 0:
            return None
        pos = mol.GetConformer().GetPositions()
        if not np.all(np.isfinite(pos)):
            return None
        return {a["id"]: pos[index[a["id"]]] for a in molekuel["atoms"]}
    except Exception:
        return None


def molekuel_aus_smiles(smiles):
    """Labor-Format mit expliziten H aus SMILES (fuer die Hydrid-Referenzen und Tests)."""
    from rdkit import Chem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return {
        "atoms": [{"id": a.GetIdx(), "element": a.GetSymbol()} for a in mol.GetAtoms()],
        "bonds": [{"atom1": b.GetBeginAtomIdx(), "atom2": b.GetEndAtomIdx(),
                   "order": int(round(b.GetBondTypeAsDouble()))} for b in mol.GetBonds()],
    }


class Quantenchemie:

    def __init__(self, max_iterationen=300):
        self.max_iterationen = max_iterationen
        self._modelle = {}

    # ------------------------------------------------------ Modelle

    def _lade(self, modell="ANI-2x"):
        if modell not in self._modelle:
            if modell == "ANI-2x":
                import torch
                import torchani

                torch.set_default_dtype(torch.float64)
                self._modelle[modell] = torchani.models.ANI2x(periodic_table_index=True).double()
            elif modell == "MACE-OFF23":
                from mace.calculators import mace_off
                self._modelle[modell] = mace_off(model="small", device="cpu", default_dtype="float64")
            elif modell == "MACE-MP-0":
                from mace.calculators import mace_mp
                self._modelle[modell] = mace_mp(model="small", device="cpu", default_dtype="float64")
            else:
                raise ValueError(modell)

        return self._modelle[modell]

    def energie_und_kraft(self, z, xyz, modell="ANI-2x"):
        """E in Hartree und dE/dx in Hartree/Angstrom."""
        m = self._lade(modell)

        if modell == "ANI-2x":
            import torch

            arten = torch.tensor([z], dtype=torch.long)
            koord = torch.tensor(np.asarray(xyz, dtype=float).reshape(1, -1, 3), requires_grad=True)
            e = m((arten, koord)).energies
            g, = torch.autograd.grad(e.sum(), koord)
            return float(e.detach()), g.detach().numpy().reshape(-1, 3)

        from ase import Atoms

        atoms = Atoms(numbers=list(z), positions=np.asarray(xyz, dtype=float).reshape(-1, 3))
        atoms.calc = m
        e_ev = float(atoms.get_potential_energy())
        f_ev = atoms.get_forces()
        # MACE rechnet in eV: Umrechnung auf Hartree, Gradient = -Kraft
        return e_ev * EV_KCAL / HARTREE_KCAL, -f_ev * EV_KCAL / HARTREE_KCAL

    # ------------------------------------------------------ Optimierung

    def optimiere(self, z, xyz, modell="ANI-2x"):
        """L-BFGS auf der Modellenergie. Gibt (Koordinaten, Energie Ha, Iterationen) zurueck."""
        from scipy.optimize import minimize

        n = len(z)

        e_start = [None]

        def f(x):
            pos = x.reshape(n, 3)
            # Schutz vor unphysischen Linesearch-Schritten (Atome fast uebereinander): MACE liefert dort
            # absurde Energien. Dann eine einfache Abstossung statt der Modellenergie verwenden.
            d = pos[:, None, :] - pos[None, :, :]
            r = np.sqrt((d ** 2).sum(-1)) + np.eye(n) * 10.0
            if e_start[0] is not None and r.min() < 0.55:
                w = np.clip(0.55 - r, 0.0, None)
                e = e_start[0] + 10.0 * float((w ** 2).sum()) + 1.0
                g = -40.0 * ((w / r)[:, :, None] * d).sum(axis=1)
                return e, g.reshape(-1)
            e, g = self.energie_und_kraft(z, pos, modell)
            if e_start[0] is None:
                e_start[0] = e
            return e, g.reshape(-1)

        r = minimize(f, np.asarray(xyz, dtype=float).reshape(-1), jac=True, method="L-BFGS-B",
                     options={"maxiter": self.max_iterationen, "gtol": 1e-4})

        return r.x.reshape(n, 3), float(r.fun), int(r.nit)

    # ------------------------------------------------------ Hauptfunktion

    def berechne(self, molekuel):
        """Siehe Modul-Docstring. Rueckgabe-Schluessel: anwendbar, grund, dhf_kcal, energie_ha, intakt,
        gebrochen, neu, geometrie, iterationen, dauer_s, schwere_atome."""

        t0 = time.time()
        atome = molekuel["atoms"]

        erg = {"anwendbar": False, "grund": None, "dhf_kcal": None, "energie_ha": None, "intakt": None,
               "gebrochen": [], "neu": [], "geometrie": None, "iterationen": 0, "dauer_s": 0.0,
               "modell": None, "fehler_kcal": None, "schwere_atome": sum(1 for a in atome if a["element"] != "H")}

        elemente = {a["element"] for a in atome}
        fremd = sorted(e for e in elemente if e not in ORDNUNGSZAHL)
        modell = waehle_modell(elemente)
        if fremd or modell is None:
            erg["grund"] = ("kein Rechenmodell fuer: " + ", ".join(fremd or sorted(elemente)) +
                            " (abgedeckt: H, C, N, O, F, P, S, Cl, Br, I, Si, B)")
            return erg
        ok, hinweis = verfuegbar(modell)
        if not ok:
            erg["grund"] = hinweis
            return erg
        erg["modell"] = modell

        z = [ORDNUNGSZAHL[a["element"]] for a in atome]

        # Offene Schale: ein Atom mit fehlender Bindung (Radikal, Carben, Nitren). Die Elektronenzahl allein
        # reicht nicht - CH2 und NH haben gerade Zahlen und sind trotzdem Tripletts.
        valenz = {a["id"]: 0 for a in atome}
        for b in molekuel["bonds"]:
            valenz[b["atom1"]] += b["order"]
            valenz[b["atom2"]] += b["order"]

        elem = {a["id"]: a["element"] for a in atome}
        defizit = sorted({elem[i] for i, v in valenz.items() if v not in NORMALE_VALENZ[elem[i]]})

        if defizit:
            erg["grund"] = ("offene Schale oder ungewoehnliche Valenz an " + ", ".join(defizit) +
                            " (Radikal/Carben/Nitren): die Modelle kennen nur geschlossene Schalen")
            return erg

        if len(atome) == 2 and {a["element"] for a in atome} == {"O"}:
            erg["grund"] = "O2 hat einen Triplett-Grundzustand: nur Singuletts moeglich"
            return erg

        start = _rdkit_start(molekuel)
        if start is None:
            erg["grund"] = "keine 3D-Startgeometrie (RDKit scheitert an der Struktur)"
            return erg

        ids = [a["id"] for a in atome]
        xyz, e, nit = self.optimiere(z, [start[i] for i in ids], modell)

        erg["energie_ha"] = e
        erg["iterationen"] = nit
        erg["geometrie"] = {i: xyz[k] for k, i in enumerate(ids)}

        # Sind die Bindungen nach der Optimierung noch da, und sind keine neuen entstanden?
        pos = {i: xyz[k] for k, i in enumerate(ids)}
        gebunden = {frozenset([b["atom1"], b["atom2"]]) for b in molekuel["bonds"]}

        for p in gebunden:
            a, b = tuple(p)
            d = float(np.linalg.norm(pos[a] - pos[b]))
            if d > BINDUNG_GEBROCHEN_FAKTOR * (RADIEN[elem[a]] + RADIEN[elem[b]]):
                erg["gebrochen"].append((a, b, round(d, 2)))

        for k, a in enumerate(ids):
            for b in ids[k + 1:]:
                if frozenset([a, b]) in gebunden:
                    continue
                d = float(np.linalg.norm(pos[a] - pos[b]))
                if d < BINDUNG_NEU_FAKTOR * (RADIEN[elem[a]] + RADIEN[elem[b]]):
                    erg["neu"].append((a, b, round(d, 2)))

        erg["intakt"] = not erg["gebrochen"] and not erg["neu"]

        # Bildungsenthalpie aus der kalibrierten ANI-Energie
        dhf = e * HARTREE_KCAL
        for a in atome:
            dhf += KALIBRIERUNG[modell][a["element"]]
        erg["fehler_kcal"] = KALIBRIERUNG_FEHLER_KCAL[modell]

        erg["dhf_kcal"] = round(dhf, 1)

        erg["anwendbar"] = True
        erg["dauer_s"] = round(time.time() - t0, 2)
        return erg