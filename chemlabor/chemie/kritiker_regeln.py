# ==========================================================
# kritiker_regeln.py
# ==========================================================
#
# Zusatzmodelle fuer den Kritiker (Mixin-Klasse): schliesst Luecken,
# die der Neuheits-Filter sichtbar gemacht hat. Ein Kandidat mit
# Si=C=P im 6er-Ring, B-O-Cl und Si=C/P=C-Doppelbindungen bekam im
# alten Kritiker 95 Punkte - zu Unrecht: nur Cl-O (-4) und die
# Ringgroesse (-5) wurden bemerkt.
#
# Aufbau: Kritiker erbt von PlausibilitaetsModelle (siehe
# patch_kritiker.py) und nutzt dessen drei modell_*-Methoden wie
# seine eigenen. Sie brauchen von Kritiker nur:
#   - self._alle_ringe(daten)       (SSSR-Naeherung, existiert bereits)
#   - daten["bindungen"], daten["atom_by_id"], daten["molekuel"]
# Wie modell_halogene/modell_hypervalenz tragen sie nur "probleme"
# und "erklaerung" ein, KEINE neuen "regeln" - was KIGedaechtnis aus
# den Regeln lernt (z.B. "Ringgroesse"), bleibt dadurch unveraendert.
#
# Die Abzuege sind bewusst grob geschaetzt (wie die Ringspannungs-
# tabelle im Kritiker), keine gemessenen Werte.
# ==========================================================


from collections import defaultdict


class PlausibilitaetsModelle:


    # ======================================================
    # Kumulene und Alkine in kleinen Ringen
    # ======================================================
    #
    # Ein Atom mit zwei Doppelbindungen (C=C=C, Si=C=P ...) will
    # linear (180 Grad) sein, eine Dreifachbindung ebenso. In kleinen
    # Ringen geht das nicht: schon 1,2-Cyclohexadien ist nur ein
    # kurzlebiges Zwischenprodukt, Cyclooctin gerade noch isolierbar,
    # ab 9 Ringatomen ist es unproblematisch. Offene Ketten (Allen,
    # CO2, Keten) und Gruppen wie S(=O)(=O) im Ring sind NICHT
    # betroffen: bestraft wird nur, wenn BEIDE Mehrfachbindungen
    # selbst Ringbindungen sind.

    _KUMULEN_RING_ABZUG = {3: 50, 4: 50, 5: 50, 6: 45, 7: 40, 8: 20}
    _ALKIN_RING_ABZUG = {3: 50, 4: 50, 5: 50, 6: 45, 7: 40, 8: 15}

    def modell_kumulene_im_ring(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        m = daten["molekuel"]

        if m["anzahl_bindungen"] - m["anzahl_atome"] + m["fragmente"] <= 0:
            return e

        doppel_partner = defaultdict(set)
        dreifach = set()

        for b in daten["bindungen"]:

            if b["ordnung"] == 2:
                doppel_partner[b["atom1"]].add(b["atom2"])
                doppel_partner[b["atom2"]].add(b["atom1"])

            elif b["ordnung"] == 3:
                dreifach.add(frozenset([b["atom1"], b["atom2"]]))

        zentren = {
            atom: partner
            for atom, partner in doppel_partner.items()
            if len(partner) >= 2
        }

        if not zentren and not dreifach:
            return e

        kumulen_ring = {}
        alkin_ring = {}

        for ring in self._alle_ringe(daten):

            n = len(ring)

            for i, atom in enumerate(ring):

                if atom in zentren:

                    ring_nachbarn = {ring[i - 1], ring[(i + 1) % n]}

                    if len(zentren[atom] & ring_nachbarn) == 2:
                        kumulen_ring[atom] = min(n, kumulen_ring.get(atom, n))

                paar = frozenset([atom, ring[(i + 1) % n]])

                if paar in dreifach:
                    alkin_ring[paar] = min(n, alkin_ring.get(paar, n))

        for atom, groesse in sorted(kumulen_ring.items()):

            abzug = self._KUMULEN_RING_ABZUG.get(groesse, 0)

            if not abzug:
                continue

            element = daten["atom_by_id"][atom]["element"]

            e["punkte"] -= abzug
            e["probleme"].append({
                "modell": "Kumulen im Ring",
                "atom": atom,
                "ringgroesse": groesse,
                "abzug": abzug
            })
            e["erklaerung"].append(
                f"{element} trägt zwei Doppelbindungen im {groesse}er-Ring - "
                f"will linear sein, im kleinen Ring stark verspannt."
            )

        for paar, groesse in sorted(
                alkin_ring.items(), key=lambda p: sorted(p[0])
        ):

            abzug = self._ALKIN_RING_ABZUG.get(groesse, 0)

            if not abzug:
                continue

            e["punkte"] -= abzug
            e["probleme"].append({
                "modell": "Dreifachbindung im Ring",
                "atome": sorted(paar),
                "ringgroesse": groesse,
                "abzug": abzug
            })
            e["erklaerung"].append(
                f"Dreifachbindung im {groesse}er-Ring - will linear sein, "
                f"im kleinen Ring stark verspannt."
            )

        return e


    # ======================================================
    # Hypohalogenite (R-O-Cl, R-O-Br, R-O-I)
    # ======================================================
    #
    # Die O-X-Bindung ist schwach, solche Verbindungen sind starke
    # Oxidationsmittel (viele explosiv). Der alte Kritiker kannte nur
    # die O-Cl-Energie (205 kJ/mol -> 4 Punkte Abzug) und gar nichts
    # fuer O-Br/O-I. Hypervalente Halogene (ClO4 usw.) gibt es im
    # Generator nicht - modell_halogene bestraft mehr als eine
    # Bindung am Halogen ohnehin.

    _HYPOHALOGENIT_ABZUG = 20

    def modell_hypohalogenite(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for b in daten["bindungen"]:

            elemente = set(b["elemente"])

            if len(elemente) != 2 or "O" not in elemente:
                continue

            halogen = (elemente - {"O"}).pop()

            if halogen not in ("Cl", "Br", "I") or b["ordnung"] != 1:
                continue

            e["punkte"] -= self._HYPOHALOGENIT_ABZUG
            e["probleme"].append({
                "modell": "Hypohalogenit",
                "atom1": b["atom1"],
                "atom2": b["atom2"],
                "bindung": f"O-{halogen}",
                "abzug": self._HYPOHALOGENIT_ABZUG
            })
            e["erklaerung"].append(
                f"O-{halogen}-Bindung (Hypohalogenit): schwach und stark "
                f"oxidierend, in der Praxis meist explosiv oder zersetzlich."
            )

        return e


    # ======================================================
    # Mehrfachbindungen schwerer Elemente (ungeschuetzt)
    # ======================================================
    #
    # Si=C, Si=Si, P=C, P=P ... gibt es real nur, wenn sperrige
    # Gruppen sie kinetisch schuetzen (Silene, Phosphaalkene,
    # Diphosphene). Der Generator kennt keine solchen Schutzgruppen,
    # also ist jede dieser Bindungen hier ungeschuetzt und reaktiv.
    # NICHT betroffen: P=O, P=N, P=S, S=O, S=N (stabil, haeufig).

    def modell_schwere_mehrfachbindungen(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for b in daten["bindungen"]:

            if b["ordnung"] < 2:
                continue

            elemente = frozenset(b["elemente"])

            if "Si" in elemente:
                abzug = 15
            elif elemente == frozenset(["P"]):      # P=P
                abzug = 15
            elif elemente == frozenset(["P", "C"]):
                abzug = 10
            else:
                continue

            name = "=".join(b["elemente"]) if b["ordnung"] == 2 else "≡".join(b["elemente"])

            e["punkte"] -= abzug
            e["probleme"].append({
                "modell": "Schwere Mehrfachbindung",
                "atom1": b["atom1"],
                "atom2": b["atom2"],
                "bindung": name,
                "abzug": abzug
            })
            e["erklaerung"].append(
                f"{name}-Mehrfachbindung ohne Schutzgruppen - reaktiv, "
                f"real nur mit sperrigen Substituenten isolierbar."
            )

        return e