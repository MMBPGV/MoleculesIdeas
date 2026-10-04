# ==========================================================
# lewis.py
# ==========================================================
#
# Lewis-Ansicht als eigenstaendiges Widget. Pruefung, Force-Directed-
# Layout (im Hintergrundthread) und Zeichnung sind WOERTLICH aus der
# bisherigen GUI uebernommen (app.py); neu ist nur der Rahmen darum
# (Widget statt Tab) und zeige()/leeren().

import math
import queue
import random
import threading
import tkinter as tk
from tkinter import ttk

from chemlabor.gui.theme import ELEMENT_FARBEN, FARBEN


class LewisAnsicht(ttk.Frame):

    FARBEN = FARBEN
    ELEMENT_FARBEN = ELEMENT_FARBEN

    def __init__(self, master):

        super().__init__(master)

        self.root = self.winfo_toplevel()

        self.lewis_queue = queue.Queue()
        self._lewis_auftrag_id = 0

        self.lewis_status = ttk.Label(
            self,
            text="Molekül in der Tabelle auswählen.",
            style="Gedimmt.TLabel"
        )

        self.lewis_status.pack(anchor="w", padx=6, pady=(4, 0))

        self.lewis_canvas = tk.Canvas(
            self,
            width=500,
            height=320,
            bg=FARBEN["panel_hell"],
            highlightthickness=0
        )

        self.lewis_canvas.pack(fill="both", expand=True, padx=6, pady=6)

        self._poll_lewis_queue()

    def zeige(self, molekuel):

        self._zeichne_lewis(molekuel)

    def leeren(self):

        self._lewis_auftrag_id += 1

        self.lewis_canvas.delete("all")

        self.lewis_status.config(text="Molekül in der Tabelle auswählen.")

    # ------------------------------------------------------------------
    # Ab hier: unveraendert aus der bisherigen GUI
    # ------------------------------------------------------------------

    def _pruefe_lewis(
            self,
            molekuel
    ):


        # Keine Bindungen vorhanden

        if not molekuel.get("bonds"):

            return {
                "moeglich": False,
                "grund":
                    "Keine Bindungsdaten vorhanden."
            }





        # Atome prüfen

        for atom in molekuel.get("atoms", []):


            valenz = atom.get(
                "valence",
                None
            )


            benutzt = atom.get(
                "used_valence",
                None
            )



            if valenz is None or benutzt is None:


                return {
                    "moeglich": False,
                    "grund":
                        "Valenzdaten fehlen."
                }




            if benutzt > valenz:


                return {
                    "moeglich": False,
                    "grund":
                        f"Ungültige Valenz bei "
                        f"{atom['element']}."
                }





        # Bindungen prüfen

        for bindung in molekuel["bonds"]:


            if (
                "atom1" not in bindung
                or
                "atom2" not in bindung
            ):


                return {
                    "moeglich": False,
                    "grund":
                        "Unvollständige Bindungsdaten."
                }




        return {
            "moeglich": True,
            "grund":
                "Darstellung möglich."
        }

    def _berechne_layout(
            self,
            atome,
            bindungen,
            breite,
            hoehe
    ):

        if not atome:
            return {}

        ids = [atom["id"] for atom in atome]
        n = len(ids)

        rand = 50  # Sicherheitsabstand zum Canvas-Rand

        # Ideale Distanz zwischen zwei Knoten, hergeleitet aus der
        # verfuegbaren Flaeche (klassische Fruchterman-Reingold-Formel).
        flaeche = max(breite - 2 * rand, 50) * max(hoehe - 2 * rand, 50)
        k = math.sqrt(flaeche / n) * 0.85

        # Weniger Iterationen bei sehr grossen Molekuelen, damit die
        # Berechnung (O(n²) pro Iteration wegen der paarweisen Abstossung)
        # auch bei hoher "Max. Atome"-Einstellung nicht spuerbar haengt.
        iterationen = 250 if n <= 40 else max(60, int(9000 / n))

        zufall = random.Random(42)  # deterministisch: gleiches Molekuel -> gleiches Bild

        pos = {
            aid: [
                zufall.uniform(rand, breite - rand),
                zufall.uniform(rand, hoehe - rand)
            ]
            for aid in ids
        }

        kanten = [
            (b["atom1"], b["atom2"])
            for b in bindungen
            if b["atom1"] in pos and b["atom2"] in pos
        ]

        temperatur = (breite + hoehe) / 20

        for _ in range(iterationen):

            kraft = {aid: [0.0, 0.0] for aid in ids}

            # Abstossung: alle Atompaare (O(n²), fuer Molekuelgroessen
            # hier unproblematisch)
            for i in range(n):
                for j in range(i + 1, n):

                    a, b = ids[i], ids[j]
                    dx = pos[a][0] - pos[b][0]
                    dy = pos[a][1] - pos[b][1]
                    dist = math.hypot(dx, dy) or 0.01

                    f = (k * k) / dist
                    fx, fy = dx / dist * f, dy / dist * f

                    kraft[a][0] += fx
                    kraft[a][1] += fy
                    kraft[b][0] -= fx
                    kraft[b][1] -= fy

            # Anziehung nur entlang tatsaechlicher Bindungen
            for a, b in kanten:

                dx = pos[a][0] - pos[b][0]
                dy = pos[a][1] - pos[b][1]
                dist = math.hypot(dx, dy) or 0.01

                f = (dist * dist) / k
                fx, fy = dx / dist * f, dy / dist * f

                kraft[a][0] -= fx
                kraft[a][1] -= fy
                kraft[b][0] += fx
                kraft[b][1] += fy

            # Positionen um die Nettokraft verschieben, auf maximal
            # "temperatur" pro Schritt begrenzt, und im Canvas halten
            for aid in ids:

                fx, fy = kraft[aid]
                dist = math.hypot(fx, fy) or 0.01
                schritt = min(dist, temperatur)

                pos[aid][0] += fx / dist * schritt
                pos[aid][1] += fy / dist * schritt

                pos[aid][0] = min(breite - rand, max(rand, pos[aid][0]))
                pos[aid][1] = min(hoehe - rand, max(rand, pos[aid][1]))

            temperatur *= 0.94  # abkuehlen, sonst oszilliert es endlos

        return {aid: (x, y) for aid, (x, y) in pos.items()}

    def _zeichne_lewis(
            self,
            molekuel
    ):
        """Einstiegspunkt: prueft, ob eine Lewis-Darstellung moeglich
        ist, und schickt die eigentliche Layout-BERECHNUNG (das teure,
        O(n²)-pro-Iteration Force-Directed-Verfahren) in einen
        Hintergrundthread.

        FIX (Performance): vorher lief die komplette Berechnung
        synchron hier im GUI-Thread - bei 100 Atomen ca. 170ms, bei 400
        Atomen ueber 1,4 SEKUNDEN, waehrend derer Tkinter komplett
        eingefroren war (nicht mal das Fenster liess sich verschieben).
        Jetzt berechnet ein Hintergrundthread die Positionen, das
        eigentliche Zeichnen (schnell, reine Canvas-Aufrufe) passiert
        weiterhin im Hauptthread, sobald das Ergebnis ueber
        lewis_queue eintrifft (siehe _poll_lewis_queue).

        Eine auftrag_id verhindert, dass ein spaeter verworfener Klick
        (Nutzer waehlt schnell mehrere Zeilen nacheinander) noch mit
        Verzoegerung ueber ein aktuelleres Ergebnis gezeichnet wird -
        nur der Auftrag mit der zuletzt vergebenen ID zaehlt."""

        canvas = self.lewis_canvas

        canvas.delete("all")

        pruefung = self._pruefe_lewis(molekuel)

        if not pruefung["moeglich"]:

            self.lewis_status.config(
                text="Lewis nicht verfügbar: " + pruefung["grund"]
            )

            self._zeige_lewis_rohdaten(molekuel)

            return

        atome = molekuel["atoms"]
        bindungen = molekuel["bonds"]
        n = len(atome)

        self.lewis_status.config(
            text=f"Layout wird berechnet ({n} Atome, {len(bindungen)} Bindungen)..."
        )

        canvas.update_idletasks()
        breite = max(canvas.winfo_width(), 400)
        hoehe = max(canvas.winfo_height(), 300)

        self._lewis_auftrag_id += 1
        auftrag_id = self._lewis_auftrag_id

        def berechne_im_hintergrund():

            positionen = self._berechne_layout(atome, bindungen, breite, hoehe)

            self.lewis_queue.put((auftrag_id, molekuel, positionen))

        threading.Thread(target=berechne_im_hintergrund, daemon=True).start()

    def _poll_lewis_queue(self):
        """Laeuft alle 50ms im Hauptthread (wie _poll_queues fuer das
        Labor) und holt fertige Layout-Ergebnisse aus lewis_queue.
        Verwirft veraltete Auftraege (auftrag_id passt nicht mehr zur
        zuletzt angeforderten Berechnung), damit bei schnellem
        Durchklicken mehrerer Molekuele nicht nachtraeglich ein altes
        Layout ueber das neue gezeichnet wird."""

        try:

            while True:

                auftrag_id, molekuel, positionen = self.lewis_queue.get_nowait()

                if auftrag_id == self._lewis_auftrag_id:
                    self._zeichne_lewis_mit_positionen(molekuel, positionen)

        except queue.Empty:

            pass

        self.root.after(50, self._poll_lewis_queue)

    def _zeichne_lewis_mit_positionen(
            self,
            molekuel,
            positionen
    ):
        """Reines Zeichnen (Canvas-Operationen) mit BEREITS berechneten
        Positionen - schnell genug, um im Hauptthread zu laufen, ohne
        die GUI spuerbar zu blockieren (siehe _zeichne_lewis fuer die
        ausgelagerte, teure Berechnung)."""

        canvas = self.lewis_canvas

        canvas.delete("all")

        atome = molekuel["atoms"]

        bindungen = molekuel["bonds"]

        n = len(atome)

        self.lewis_status.config(
            text=
            f"Lewis-Struktur erzeugt ({n} Atome, {len(bindungen)} Bindungen)."
        )

        # Atomradius und Schrift schrumpfen mit wachsender Atomanzahl,
        # damit sich benachbarte Atome bei grossen Molekuelen nicht
        # gegenseitig verdecken.
        atom_radius = max(10, min(25, 260 // max(n, 1)))
        schrift_groesse = max(7, min(14, 140 // max(n, 1)))



        # --------------------------------------------------
        # Bindungen zeichnen
        # --------------------------------------------------


        for bindung in bindungen:


            a1 = bindung["atom1"]

            a2 = bindung["atom2"]



            if (
                a1 not in positionen
                or
                a2 not in positionen
            ):

                continue



            x1,y1 = positionen[a1]

            x2,y2 = positionen[a2]



            ordnung = bindung.get(
                "order",
                1
            )


            # Bindungslinie(n) senkrecht zur Verbindungsachse versetzen
            # (nicht mehr nur vertikal wie vorher) - bei schraeg
            # verlaufenden Bindungen sonst optisch falsch.
            dx, dy = x2 - x1, y2 - y1
            laenge = math.hypot(dx, dy) or 1
            nx, ny = -dy / laenge, dx / laenge  # Normalenvektor



            for offset in range(
                ordnung
            ):


                versatz = (
                    offset
                    -
                    (ordnung-1)/2
                ) * 5



                canvas.create_line(

                    x1 + nx * versatz,

                    y1 + ny * versatz,

                    x2 + nx * versatz,

                    y2 + ny * versatz,

                    width=2,

                    fill=self.FARBEN["text"]

                )







        # --------------------------------------------------
        # Atome zeichnen
        # --------------------------------------------------


        for atom in atome:


            x,y = positionen[
                atom["id"]
            ]

            farbe = self.ELEMENT_FARBEN.get(atom["element"], "#ffffff")
            textfarbe = "#ffffff" if farbe in ("#4a4a4a", "#3060c8", "#c83030", "#2f9e5f", "#a8763f") else "#000000"



            canvas.create_oval(

                x-atom_radius,

                y-atom_radius,

                x+atom_radius,

                y+atom_radius,

                fill=farbe,
                outline=self.FARBEN["rand"]

            )



            canvas.create_text(

                x,

                y,

                text=atom["element"],

                font=("", schrift_groesse, "bold"),
                fill=textfarbe

            )

    def _zeige_lewis_rohdaten(
            self,
            molekuel
    ):


        canvas = self.lewis_canvas


        y = 40



        canvas.create_text(

            400,

            y,

            text="Rohdaten:",

            font=("",14,"bold"),

            fill=self.FARBEN["akzent"]

        )


        y += 40



        for atom in molekuel["atoms"]:


            text = (

                f"Atom {atom['id']}: "
                f"{atom['element']} "
                f"(Valenz {atom.get('valence','?')}, "
                f"belegt {atom.get('used_valence','?')})"

            )


            canvas.create_text(

                400,

                y,

                text=text,

                fill=self.FARBEN["text"]

            )


            y += 25



        y += 20



        for bindung in molekuel["bonds"]:


            text = (

                f"Bindung: "
                f"{bindung.get('atom1')} - "
                f"{bindung.get('atom2')} "
                f"(Ordnung "
                f"{bindung.get('order','?')})"

            )


            canvas.create_text(

                400,

                y,

                text=text,

                fill=self.FARBEN["text"]

            )


            y += 25
