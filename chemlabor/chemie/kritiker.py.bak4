import json
from collections import defaultdict, deque
from chemlabor.chemie.kritiker_regeln import PlausibilitaetsModelle

# ==========================================================
# Kritiker v1.2 (gefixt)
# Stabilitätslehrer
# ==========================================================

ELEMENT_INFO = {

    "H": {"gruppe": "Nichtmetall", "oktett": 2, "en": 2.20},
    "C": {"gruppe": "Nichtmetall", "oktett": 8, "en": 2.55},
    "N": {"gruppe": "Nichtmetall", "oktett": 8, "en": 3.04},
    "O": {"gruppe": "Nichtmetall", "oktett": 8, "en": 3.44},
    "F": {"gruppe": "Halogen", "oktett": 8, "en": 3.98},
    "Cl": {"gruppe": "Halogen", "oktett": 8, "en": 3.16},
    "Br": {"gruppe": "Halogen", "oktett": 8, "en": 2.96},
    "I": {"gruppe": "Halogen", "oktett": 8, "en": 2.66},
    "P": {"gruppe": "Nichtmetall", "oktett": 8, "en": 2.19},
    "S": {"gruppe": "Nichtmetall", "oktett": 8, "en": 2.58},
    "Si": {"gruppe": "Halbmetall", "oktett": 8, "en": 1.90},

    # Bor ist ein legitimes "Elektronenmangelelement" - BH3/BF3 sind
    # stabile, reale Molekuele mit nur 3 Bindungspaaren und KEINEM
    # freien Elektronenpaar (das leere p-Orbital ist genau ihre
    # charakteristische Eigenschaft, nicht ein Fehler). Statt einer
    # Sonderregel in jedem einzelnen Modell wird hier derselbe Trick
    # wie bei Wasserstoff genutzt: "oktett" ist das Ziel fuer die
    # VOLLE Schale dieses Elements, nicht zwingend 8. H bekommt ein
    # Duett (2), Bor bekommt ein Sextett (6) - dadurch berechnen
    # modell_geometrie (freie_paare) und modell_formalladung (N)
    # automatisch korrekt "keine freien Elektronenpaare bei voller
    # Valenz", ohne dass Bor irgendwo explizit abgefragt werden muss.
    "B": {"gruppe": "Halbmetall", "oktett": 6, "en": 2.04}

}


# ==========================================================
# Kritiker
# ==========================================================

class Kritiker(PlausibilitaetsModelle):

    def __init__(self):

        # Reihenfolge ist bewusst gewählt.
        # Erst billige Modelle, dann aufwendige.

        self.modelle = [

            self.modell_valenz,
            self.modell_fragmente,
            self.modell_offene_stellen,
            self.modell_bindungen,
            self.modell_halogene,
            self.modell_struktur,
            self.modell_oktett,
            self.modell_elektronegativitaet,
            self.modell_formalladung,
            self.modell_mesomerie,
            self.modell_geometrie,
            self.modell_ringe,
            self.modell_hypervalenz,
            self.modell_kumulene_im_ring,
            self.modell_hypohalogenite,
            self.modell_schwere_mehrfachbindungen

        ]


    # ======================================================
    # Hauptanalyse
    # ======================================================

    def analysiere(self, molekuel):

        daten = {

            "molekuel": {
                "anzahl_atome": len(molekuel["atoms"]),
                "anzahl_bindungen": len(molekuel["bonds"]),
                "fragmente": 0
            },

            "atome": [],
            "bindungen": [],
            "muster": [],
            "graph": defaultdict(list),
            "wachstum": [],
            "regeln": [],
            "probleme": [],
            "erklaerung": [],
            "bewertung": {},

            # Einmalig aufgebaute Lookups statt wiederholter linearer
            # Suchen. WICHTIG: diese Dicts muessen auch tatsaechlich
            # ueberall benutzt werden, wo vorher linear gesucht wurde
            # (analyse_bindungen -> atom_by_id, _ring_bindungsordnungen
            # -> bindung_by_paar) - sonst ist der Lookup nur totes Gewicht.
            "atom_by_id": {a["id"]: a for a in molekuel["atoms"]},

        }

        self.baue_graph(molekuel, daten)
        self.analyse_atome(molekuel, daten)
        self.analyse_bindungen(molekuel, daten)
        self.analyse_fragmente(molekuel, daten)
        self.analyse_wachstum(daten)

        # Bindungs-Lookup erst NACH analyse_bindungen moeglich (braucht
        # daten["bindungen"]), aber vor stabilitaets_test, da mehrere
        # Modelle (Ringe, Konjugation) das brauchen.
        daten["bindung_by_paar"] = {
            frozenset([b["atom1"], b["atom2"]]): b for b in daten["bindungen"]
        }

        self.stabilitaets_test(daten)

        return daten


    # ======================================================
    # Graph erzeugen
    # ======================================================

    def baue_graph(self, molekuel, daten):

        for b in molekuel["bonds"]:
            daten["graph"][b["atom1"]].append(b["atom2"])
            daten["graph"][b["atom2"]].append(b["atom1"])


    # ======================================================
    # Atome analysieren
    # ======================================================

    def analyse_atome(self, molekuel, daten):

        for atom in molekuel["atoms"]:

            info = ELEMENT_INFO.get(atom["element"], {})
            frei = atom["valence"] - atom["used_valence"]

            daten["atome"].append({
                "id": atom["id"],
                "element": atom["element"],
                "typ": info.get("gruppe", "Unbekannt"),
                "elektronegativitaet": info.get("en"),
                "oktett": info.get("oktett"),
                "valenz": {
                    "max": atom["valence"],
                    "benutzt": atom["used_valence"],
                    "frei": frei
                }
            })


    # ======================================================
    # Bindungen
    # ======================================================

    def analyse_bindungen(self, molekuel, daten):

        namen = {1: "Einfachbindung", 2: "Doppelbindung", 3: "Dreifachbindung"}

        for b in molekuel["bonds"]:

            # FIX: vorher self.find_atom(molekuel, ...) - lineare Suche
            # ueber die komplette Atomliste, pro Bindung zweimal. Der
            # atom_by_id-Lookup wurde zwar in analysiere() gebaut, aber
            # hier nie benutzt. Jetzt O(1) statt O(n).
            a1 = daten["atom_by_id"][b["atom1"]]
            a2 = daten["atom_by_id"][b["atom2"]]

            daten["bindungen"].append({
                "atom1": a1["id"],
                "atom2": a2["id"],
                "elemente": [a1["element"], a2["element"]],
                "ordnung": b["order"],
                "typ": namen.get(b["order"], "Unbekannt")
            })

            daten["muster"].append({
                "muster": f"{a1['element']}-{a2['element']}",
                "ordnung": b["order"]
            })


    # ======================================================
    # Fragmente
    # ======================================================

    def analyse_fragmente(self, molekuel, daten):

        besucht = set()
        fragmente = 0

        for atom in molekuel["atoms"]:

            if atom["id"] in besucht:
                continue

            fragmente += 1
            self.dfs(atom["id"], daten["graph"], besucht)

        daten["molekuel"]["fragmente"] = fragmente


    def dfs(self, start, graph, besucht):

        besucht.add(start)

        for nachbar in graph[start]:
            if nachbar not in besucht:
                self.dfs(nachbar, graph, besucht)


    # ======================================================
    # Wachstum
    # ======================================================

    def analyse_wachstum(self, daten):

        for atom in daten["atome"]:

            frei = atom["valenz"]["frei"]

            if frei > 0:
                daten["wachstum"].append({
                    "atom": atom["id"],
                    "element": atom["element"],
                    "freie_valenz": frei
                })


    # ======================================================
    # Modelle ausführen
    # ======================================================

    def stabilitaets_test(self, daten):

        punkte = 100

        for modell in self.modelle:

            ergebnis = modell(daten)

            punkte += ergebnis["punkte"]
            daten["regeln"].extend(ergebnis["regeln"])
            daten["probleme"].extend(ergebnis["probleme"])
            daten["erklaerung"].extend(ergebnis["erklaerung"])

        punkte = max(0, min(100, punkte))

        if punkte >= 90:
            urteil = "sehr stabil"
        elif punkte >= 75:
            urteil = "wahrscheinlich stabil"
        elif punkte >= 50:
            urteil = "fragwürdig"
        else:
            urteil = "instabil"

        daten["bewertung"] = {"punkte": punkte, "urteil": urteil}


    # ======================================================
    # Hilfsfunktionen
    # ======================================================

    # [patch_fixes] Elemente mit mehreren ueblichen Valenzzustaenden.
    # idk.VALENZEN fuehrt fuer P/S die MAXIMALvalenz (5 bzw. 6) - daran
    # gemessen galten H2S (S mit 2 Bindungen) oder PH3 als "4 bzw. 2
    # offene Valenzen" und wurden hart bestraft, waehrend SF6 keine
    # offene Valenz hatte. Gemessen wird jetzt gegen den naechsten
    # erlaubten Zustand (RDKit akzeptiert genau diese: S 2/4/6, P 3/5).
    _ERLAUBTE_VALENZ = {"S": (2, 4, 6), "P": (3, 5)}
    _VALENZELEKTRONEN_EXPANDIERT = {"S": 6, "P": 5}

    def _offene_valenz(self, atom):
        """Freie Valenz relativ zum naechsten erlaubten Valenzzustand
        (bei allen Elementen ohne Zustandsliste: wie bisher
        max - benutzt)."""

        frei = atom["valenz"]["frei"]
        zustaende = self._ERLAUBTE_VALENZ.get(atom["element"])

        if zustaende is None or frei <= 0:
            return frei

        benutzt = atom["valenz"]["benutzt"]

        for zustand in zustaende:
            if zustand >= benutzt:
                return zustand - benutzt

        return frei

    def _nichtbindende_elektronen(self, atom):
        """Elektronen in freien Paaren (N in FC = V - (N + B/2)).
        Fuer normale Elemente unveraendert oktett - 2*benutzt. Bei S/P
        mit erweitertem Oktett (z.B. S in SO2: 4 Bindungen, 1 freies
        Paar) waere diese Formel bei 0 abgeschnitten worden und haette
        eine falsche Formalladung (+2) samt Ladungsbilanz-Abzug und
        falscher Geometrie (linear statt gewinkelt) ergeben."""

        benutzt = atom["valenz"]["benutzt"]
        v = self._VALENZELEKTRONEN_EXPANDIERT.get(atom["element"])

        if v is not None:
            return 2 * (max(v - benutzt, 0) // 2)

        return max(atom["oktett"] - 2 * benutzt, 0)

    def find_atom(self, molekuel, nummer):
        # Nur noch als Fallback fuer Aufrufer ohne fertiges daten-Dict.
        # Innerhalb der eigenen Pipeline immer daten["atom_by_id"] nutzen.

        for atom in molekuel["atoms"]:
            if atom["id"] == nummer:
                return atom

        return None


    # ==========================================================
    # MODELL 1 — VALENZ
    # ==========================================================

    def modell_valenz(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for atom in daten["atome"]:

            frei = atom["valenz"]["frei"]

            if frei < 0:

                e["punkte"] -= 40
                e["probleme"].append({
                    "modell": "Valenz",
                    "atom": atom["id"],
                    "problem": "Valenz überschritten"
                })
                e["erklaerung"].append(
                    f"{atom['element']} besitzt mehr Bindungen als erlaubt."
                )

            else:

                e["regeln"].append({"regel": "Valenzgrenze", "status": "OK"})

        return e


    # ==========================================================
    # MODELL 2 — FRAGMENTE
    # ==========================================================

    def modell_fragmente(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        fragmente = daten["molekuel"]["fragmente"]

        if fragmente > 1:

            abzug = (fragmente - 1) * 40
            e["punkte"] -= abzug
            e["probleme"].append({
                "modell": "Konnektivität",
                "problem": f"{fragmente} Fragmente",
                "abzug": abzug
            })
            e["erklaerung"].append(
                "Das Molekül besteht aus mehreren unabhängigen Teilstrukturen."
            )

        else:

            e["regeln"].append({"regel": "Zusammenhängende Struktur", "status": "OK"})

        return e


    # ==========================================================
    # MODELL 3 — OFFENE VALENZEN
    # ==========================================================

    def modell_offene_stellen(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for atom in daten["atome"]:

            frei = self._offene_valenz(atom)

            if frei == 0:
                continue

            if frei % 2 == 1:

                # Ungerade fehlende Valenz = ein ungepaartes Elektron =
                # Radikal. Radikale sind per Definition reaktiv und
                # instabil, egal wie "sauber" der Rest des Moleküls
                # aussieht - das muss den Score deutlich einbrechen
                # lassen, nicht nur um ein paar Punkte nach unten korrigieren.
                abzug = 50 + (frei - 1) * 25

                e["punkte"] -= abzug
                e["probleme"].append({
                    "modell": "Radikal",
                    "atom": atom["id"],
                    "ungepaarte_elektronen": (frei + 1) // 2
                })
                e["erklaerung"].append(
                    f"{atom['element']} besitzt ein ungepaartes Elektron (Radikalcharakter) - reaktiv und instabil."
                )

            else:

                abzug = frei * 10
                e["punkte"] -= abzug
                e["probleme"].append({
                    "modell": "Offene Valenzen",
                    "atom": atom["id"],
                    "freie_valenzen": frei
                })
                e["erklaerung"].append(
                    f"{atom['element']} besitzt {frei} ungesättigte Bindungsplätze."
                )

        return e


    # ==========================================================
    # MODELL 4 — BINDUNGSORDNUNGEN
    # ==========================================================

    def modell_bindungen(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for bindung in daten["bindungen"]:

            # FIX: vorher wurde nur die Dreifachbindung mit H geprueft -
            # eine Doppelbindung mit H (H kann ueberhaupt keine
            # Mehrfachbindung eingehen, da nur 1s-Orbital) rutschte
            # unbestraft durch.
            if bindung["ordnung"] >= 2 and "H" in bindung["elemente"]:

                e["punkte"] -= 30
                e["probleme"].append({
                    "modell": "Bindungsordnung",
                    "atom1": bindung["atom1"],
                    "atom2": bindung["atom2"],
                    "problem": "H kann keine Mehrfachbindung besitzen"
                })

        e["regeln"].append({"regel": "Bindungsordnungen geprüft", "status": "OK"})

        return e


    # ==========================================================
    # MODELL 5 — HALOGENE
    # ==========================================================

    def modell_halogene(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        halogene = {"F", "Cl", "Br", "I"}

        for atom in daten["atome"]:

            if atom["element"] not in halogene:
                continue

            if atom["valenz"]["benutzt"] > 1:

                e["punkte"] -= 20
                e["probleme"].append({
                    "modell": "Halogene",
                    "atom": atom["id"],
                    "problem": "Halogen besitzt zu viele Bindungen"
                })

        return e


    # ==========================================================
    # MODELL 6 — EINFACHE STRUKTURHEURISTIK
    # ==========================================================

    def modell_struktur(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        m = daten["molekuel"]

        grade = [len(daten["graph"][atom["id"]]) for atom in daten["atome"]]
        durchschnittsgrad = round(sum(grade) / m["anzahl_atome"], 2) if m["anzahl_atome"] > 0 else 0

        e["regeln"].append({"regel": "Durchschnittsgrad", "wert": durchschnittsgrad})

        if m["anzahl_bindungen"] < m["anzahl_atome"] - 1:

            e["punkte"] -= 15
            e["probleme"].append({
                "modell": "Struktur",
                "problem": "Sehr geringe Vernetzung",
                "durchschnittsgrad": durchschnittsgrad
            })

        else:

            e["regeln"].append({"regel": "Vernetzung ausreichend", "status": "OK"})

        return e


    # ==========================================================
    # MODELL 7 — OKTETTREGEL
    # ==========================================================

    def modell_oktett(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for atom in daten["atome"]:

            if atom["oktett"] is None:
                continue

            benutzt = atom["valenz"]["benutzt"]
            maximum = atom["valenz"]["max"]

            if benutzt >= maximum:

                e["regeln"].append({"regel": "Oktett wahrscheinlich erfüllt", "atom": atom["id"]})

            else:

                fehlend = (maximum - benutzt) * 2
                e["regeln"].append({
                    "regel": "Oktett unvollständig",
                    "atom": atom["id"],
                    "fehlende_elektronen": fehlend
                })
                e["erklaerung"].append(f"{atom['element']} erreicht kein vollständiges Oktett.")

        return e


    # ==========================================================
    # MODELL — ELEKTRONEGATIVITÄT
    # ==========================================================

    _BINDUNGSENERGIEN = {
        (frozenset(["H", "H"]), 1): 436, (frozenset(["H", "C"]), 1): 413,
        (frozenset(["H", "N"]), 1): 391, (frozenset(["H", "O"]), 1): 467,
        (frozenset(["H", "F"]), 1): 565, (frozenset(["H", "Cl"]), 1): 431,
        (frozenset(["H", "S"]), 1): 347, (frozenset(["H", "P"]), 1): 322,
        (frozenset(["H", "Si"]), 1): 393,
        (frozenset(["C", "C"]), 1): 347, (frozenset(["C", "N"]), 1): 305,
        (frozenset(["C", "O"]), 1): 358, (frozenset(["C", "F"]), 1): 485,
        (frozenset(["C", "Cl"]), 1): 339, (frozenset(["C", "S"]), 1): 259,
        (frozenset(["C", "Si"]), 1): 318,
        (frozenset(["N", "N"]), 1): 163, (frozenset(["N", "O"]), 1): 201,
        (frozenset(["N", "F"]), 1): 283,
        (frozenset(["O", "O"]), 1): 146, (frozenset(["O", "F"]), 1): 190,
        (frozenset(["O", "Cl"]), 1): 205, (frozenset(["O", "S"]), 1): 265,
        (frozenset(["O", "Si"]), 1): 452,
        (frozenset(["F", "F"]), 1): 155, (frozenset(["F", "Cl"]), 1): 253,
        (frozenset(["F", "Si"]), 1): 565,
        (frozenset(["Cl", "Cl"]), 1): 243, (frozenset(["Cl", "Si"]), 1): 381,
        (frozenset(["S", "S"]), 1): 226,
        (frozenset(["P", "P"]), 1): 201,
        (frozenset(["Si", "Si"]), 1): 222,
        # Brom/Iod - nur die gut belegten Standardwerte, bewusst ohne
        # Ergaenzungen wie Br-N/Br-O/I-S (dieselbe Zurueckhaltung wie
        # beim Rest der Tabelle: kein Eintrag ist besser als eine
        # geschaetzte Zahl).
        (frozenset(["H", "Br"]), 1): 366, (frozenset(["H", "I"]), 1): 298,
        (frozenset(["C", "Br"]), 1): 276, (frozenset(["C", "I"]), 1): 240,
        (frozenset(["Br", "Br"]), 1): 193, (frozenset(["I", "I"]), 1): 151,
        (frozenset(["F", "Br"]), 1): 249, (frozenset(["F", "I"]), 1): 278,
        (frozenset(["Cl", "Br"]), 1): 216, (frozenset(["Cl", "I"]), 1): 208,
        (frozenset(["Br", "I"]), 1): 175,
        # Bor - nur die gut belegten Standardwerte.
        (frozenset(["B", "H"]), 1): 389, (frozenset(["B", "F"]), 1): 613,
        (frozenset(["B", "Cl"]), 1): 456, (frozenset(["B", "B"]), 1): 293,
        # Doppelbindungen
        (frozenset(["C", "C"]), 2): 614, (frozenset(["C", "N"]), 2): 615,
        (frozenset(["C", "O"]), 2): 799,
        (frozenset(["N", "N"]), 2): 418, (frozenset(["N", "O"]), 2): 607,
        (frozenset(["O", "O"]), 2): 498,
        # Dreifachbindungen
        (frozenset(["C", "C"]), 3): 839, (frozenset(["C", "N"]), 3): 891,
        (frozenset(["N", "N"]), 3): 945,
    }

    _SCHWACHE_BINDUNG_SCHWELLE = 250  # kJ/mol

    def modell_elektronegativitaet(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for bindung in daten["bindungen"]:

            e1 = bindung["elemente"][0]
            e2 = bindung["elemente"][1]

            # FIX: .get(...) statt direktem [...]-Zugriff - vermeidet
            # einen KeyError, falls ein Molekuel mal ein Element ohne
            # EN-Wert in ELEMENT_INFO enthaelt. Bindung wird dann einfach
            # uebersprungen statt die ganze Analyse abstuerzen zu lassen.
            en1 = ELEMENT_INFO.get(e1, {}).get("en")
            en2 = ELEMENT_INFO.get(e2, {}).get("en")

            if en1 is None or en2 is None:
                continue

            diff = abs(en1 - en2)

            if diff < 0.3:
                e["regeln"].append({"regel": "Unpolare Bindung", "bindung": f"{e1}-{e2}"})
            elif diff < 1.8:
                e["regeln"].append({"regel": "Polare kovalente Bindung", "bindung": f"{e1}-{e2}"})
            else:
                abzug = min(round((diff - 1.8) * 10), 15)
                e["punkte"] -= abzug
                e["regeln"].append({"regel": "Sehr starke Polarität", "bindung": f"{e1}-{e2}"})
                e["erklaerung"].append(
                    f"{e1}-{e2}-Bindung ist extrem polar (ΔEN={diff:.2f})."
                )

            energie = self._BINDUNGSENERGIEN.get((frozenset([e1, e2]), bindung["ordnung"]))

            if energie is not None and energie < self._SCHWACHE_BINDUNG_SCHWELLE:

                abzug = round((self._SCHWACHE_BINDUNG_SCHWELLE - energie) / 10)
                e["punkte"] -= abzug
                e["probleme"].append({
                    "modell": "Bindungsenergie",
                    "bindung": f"{e1}-{e2}",
                    "energie_kj_mol": energie,
                    "abzug": abzug
                })
                e["erklaerung"].append(
                    f"{e1}-{e2}-Bindung ist mit {energie} kJ/mol ungewöhnlich schwach."
                )


        return e


    # ==========================================================
    # MODELL 8 — FORMALLADUNGEN
    # ==========================================================

    def modell_formalladung(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        valenzelektronen = {"H": 1, "C": 4, "N": 5, "O": 6, "F": 7, "Cl": 7, "Br": 7, "I": 7, "P": 5, "S": 6, "Si": 4, "B": 3}
        gesamtladung = 0

        for atom in daten["atome"]:

            element = atom["element"]

            if element not in valenzelektronen:
                continue

            if atom["oktett"] is None:
                continue

            benutzt = atom["valenz"]["benutzt"]

            V = valenzelektronen[element]
            N = self._nichtbindende_elektronen(atom)
            B = 2 * benutzt

            ladung = V - (N + B / 2)
            ladung = int(ladung) if ladung == int(ladung) else ladung

            gesamtladung += ladung

            if abs(ladung) > 2:

                e["punkte"] -= 10
                e["probleme"].append({
                    "modell": "Formalladung",
                    "atom": atom["id"],
                    "ladung": ladung
                })
                e["erklaerung"].append(f"{element} besitzt eine ungewöhnliche Formalladung.")

        if gesamtladung != 0:

            abzug = min(abs(gesamtladung) * 8, 40)
            e["punkte"] -= abzug
            e["probleme"].append({
                "modell": "Ladungsbilanz",
                "gesamtladung": gesamtladung,
                "abzug": abzug
            })
            e["erklaerung"].append(
                f"Summe aller Formalladungen ist {gesamtladung}, nicht 0 - "
                f"für ein als neutral angenommenes Molekül unrealistisch."
            )

        return e


    # ==========================================================
    # MODELL 9 — MESOMERIE
    # ==========================================================

    def modell_mesomerie(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        doppelbindungen = [b for b in daten["bindungen"] if b["ordnung"] >= 2]

        if len(doppelbindungen) < 2:
            return e

        einfachbindungs_paare = {
            frozenset([b["atom1"], b["atom2"]])
            for b in daten["bindungen"] if b["ordnung"] == 1
        }

        konjugierte_paare = 0

        for i in range(len(doppelbindungen)):

            atome1 = {doppelbindungen[i]["atom1"], doppelbindungen[i]["atom2"]}

            for j in range(i + 1, len(doppelbindungen)):

                atome2 = {doppelbindungen[j]["atom1"], doppelbindungen[j]["atom2"]}

                for a1 in atome1:
                    for a2 in atome2:
                        if a1 != a2 and frozenset([a1, a2]) in einfachbindungs_paare:
                            konjugierte_paare += 1

        if konjugierte_paare > 0:

            e["punkte"] += 5
            e["regeln"].append({
                "regel": "Konjugation",
                "konjugierte_paare": konjugierte_paare
            })
            e["erklaerung"].append(
                f"{konjugierte_paare} konjugierte(s) Doppelbindungspaar(e) gefunden - Mesomeriestabilisierung möglich."
            )

        return e


    # ==========================================================
    # MODELL 10 — GEOMETRIE (VSEPR)
    # ==========================================================

    _VSEPR_FORM = {
        (2, 0): ("linear", "sp"),
        (3, 0): ("trigonal-planar", "sp2"),
        (3, 1): ("gewinkelt", "sp2"),
        (4, 0): ("tetraedrisch", "sp3"),
        (4, 1): ("trigonal-pyramidal", "sp3"),
        (4, 2): ("gewinkelt", "sp3"),
        (5, 0): ("trigonal-bipyramidal", "sp3d"),
        (5, 1): ("wippenförmig (see-saw)", "sp3d"),
        (5, 2): ("T-förmig", "sp3d"),
        (5, 3): ("linear", "sp3d"),
        (6, 0): ("oktaedrisch", "sp3d2"),
        (6, 1): ("quadratisch-pyramidal", "sp3d2"),
        (6, 2): ("quadratisch-planar", "sp3d2"),
    }

    def modell_geometrie(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for atom in daten["atome"]:

            if atom["oktett"] is None:
                continue

            sigma = len(daten["graph"][atom["id"]])

            if sigma < 2:
                continue

            benutzt = atom["valenz"]["benutzt"]
            freie_paare = self._nichtbindende_elektronen(atom) // 2
            domänen = sigma + freie_paare

            form = self._VSEPR_FORM.get((domänen, freie_paare))

            if form is None:
                continue

            geometrie, hybridisierung = form

            e["regeln"].append({
                "regel": "Geometrie",
                "atom": atom["id"],
                "form": geometrie,
                "hybridisierung": hybridisierung
            })
            e["erklaerung"].append(
                f"{atom['element']} ist {geometrie} ({hybridisierung}, "
                f"{sigma} Bindungspartner, {freie_paare} freie(s) Elektronenpaar(e))."
            )

        return e


    # ==========================================================
    # MODELL 11 — RINGSPANNUNG
    # ==========================================================

    _RINGSPANNUNG = {3: -35, 4: -20, 5: -5, 6: 5, 7: -5, 8: -10}
    _AROMATIZITAETS_BONUS = 20

    def _alle_ringe(self, daten):
        """[patch_fixes] Kleinster Satz unabhaengiger Ringe (SSSR-
        Naeherung) statt der Fundamentalzyklen eines BFS-Spannbaums.
        Letztere koennen bei kondensierten Ringsystemen (Naphthalin)
        einen zu grossen Aussenring statt des kleinsten liefern und
        passten nicht zur Ringgroesse, die KIGenerator beim Bauen
        ueber die kuerzeste Pfadlaenge bestimmt.

        Vorgehen: fuer jede Kante der kuerzeste Zyklus durch sie
        (BFS ohne diese Kante), dazu die alten Fundamentalzyklen als
        Reserve; dann nach Laenge sortiert nur die linear unabhaengigen
        (GF(2)-Elimination ueber Kantenvektoren) behalten, bis die
        Ringzahl (Zyklomatische Zahl) erreicht ist."""

        graph = daten["graph"]

        kanten = set()
        for u, nachbarn in graph.items():
            for v in nachbarn:
                if u != v:
                    kanten.add(frozenset((u, v)))

        anzahl_ringe = (
            len(kanten) - len(daten["atome"]) + daten["molekuel"]["fragmente"]
        )

        if anzahl_ringe <= 0:
            return []

        kanten_index = {
            kante: i
            for i, kante in enumerate(sorted(kanten, key=lambda k: tuple(sorted(k))))
        }

        def kuerzester_zyklus(u, v):

            eltern = {u: None}
            queue = deque([u])

            while queue:

                knoten = queue.popleft()

                for nachbar in graph[knoten]:

                    if knoten == u and nachbar == v:
                        continue

                    if nachbar in eltern:
                        continue

                    eltern[nachbar] = knoten

                    if nachbar == v:

                        pfad = []
                        k = v

                        while k is not None:
                            pfad.append(k)
                            k = eltern[k]

                        pfad.reverse()

                        return pfad

                    queue.append(nachbar)

            return None

        def vektor_von(ring):

            vektor = 0
            n = len(ring)

            for i in range(n):
                kante = frozenset((ring[i], ring[(i + 1) % n]))
                vektor |= 1 << kanten_index[kante]

            return vektor

        kandidaten = {}

        for kante in kanten:

            u, v = tuple(kante)
            ring = kuerzester_zyklus(u, v)

            if ring is not None and len(ring) >= 3:
                kandidaten.setdefault(vektor_von(ring), ring)

        for ring in self._fundamentalzyklen(daten):

            try:
                kandidaten.setdefault(vektor_von(ring), ring)
            except KeyError:
                continue

        basis = {}
        ringe = []

        for vektor, ring in sorted(
            kandidaten.items(), key=lambda paar: len(paar[1])
        ):

            rest = vektor

            while rest:

                hoechstes = rest.bit_length() - 1

                if hoechstes in basis:
                    rest ^= basis[hoechstes]
                else:
                    basis[hoechstes] = rest
                    ringe.append(ring)
                    break

            if len(ringe) >= anzahl_ringe:
                break

        return ringe

    def _fundamentalzyklen(self, daten):
        """Findet ALLE unabhängigen Ringe im Molekül (nicht nur den
        kleinsten): pro Fragment wird ein Spannbaum aufgebaut (BFS), und
        jede Kante, die NICHT im Spannbaum liegt, schließt zusammen mit
        dem Baumpfad zwischen ihren beiden Enden genau einen
        Fundamentalzyklus."""

        graph = daten["graph"]
        alle_ids = [atom["id"] for atom in daten["atome"]]

        besucht = set()
        eltern = {}
        baum_kanten = set()

        for start in alle_ids:

            if start in besucht:
                continue

            besucht.add(start)
            eltern[start] = None
            queue = deque([start])

            while queue:

                knoten = queue.popleft()

                for nachbar in graph[knoten]:

                    if nachbar not in besucht:
                        besucht.add(nachbar)
                        eltern[nachbar] = knoten
                        baum_kanten.add(frozenset([knoten, nachbar]))
                        queue.append(nachbar)

        alle_kanten = set()
        for u, nachbarn in graph.items():
            for v in nachbarn:
                alle_kanten.add(frozenset([u, v]))

        ringe = []

        for kante in alle_kanten:

            if kante in baum_kanten:
                continue

            u, v = tuple(kante)

            pfad_u = []
            knoten = u
            while knoten is not None:
                pfad_u.append(knoten)
                knoten = eltern[knoten]

            pfad_v = []
            knoten = v
            while knoten is not None:
                pfad_v.append(knoten)
                knoten = eltern[knoten]

            menge_u = set(pfad_u)
            lca = next(k for k in pfad_v if k in menge_u)

            ring = []
            for knoten in pfad_u:
                ring.append(knoten)
                if knoten == lca:
                    break

            # FIX: der Pfad von v zur lca (pfad_v) verlaeuft in
            # Kind->Eltern-Richtung (v, eltern[v], eltern[eltern[v]],
            # ..., lca) - haengt man ihn in dieser Reihenfolge an, springt
            # der Ring direkt von der lca zu einem Knoten, der u.U. gar
            # nicht mit ihr benachbart ist (er ist mit dem naechsten
            # Element benachbart, nicht mit der lca). Der Ast muss
            # UMGEDREHT angehaengt werden, damit auf die lca zuerst ihr
            # tatsaechlicher Baum-Nachbar folgt und der Ringpfad
            # durchgehend aus benachbarten Knoten besteht.
            pfad_v_bis_lca = []
            for knoten in pfad_v:
                if knoten == lca:
                    break
                pfad_v_bis_lca.append(knoten)

            ring.extend(reversed(pfad_v_bis_lca))

            ringe.append(ring)

        return ringe


    def _ring_bindungsordnungen(self, ring_pfad, daten):
        """Bindungsordnungen entlang des Rings, in Ring-Reihenfolge
        (inkl. der schliessenden Bindung vom letzten zum ersten Atom)."""

        n = len(ring_pfad)
        ordnungen = []

        for i in range(n):

            a, b = ring_pfad[i], ring_pfad[(i + 1) % n]

            # FIX: vorher lineare Suche ueber daten["bindungen"] pro
            # Ringkante - genau die O(n)-Suche, die der bindung_by_paar-
            # Lookup in analysiere() eigentlich ersetzen sollte, aber
            # hier nie benutzt wurde. Jetzt O(1) statt O(Bindungen).
            bindung = daten["bindung_by_paar"].get(frozenset([a, b]))

            if bindung is not None:
                ordnungen.append(bindung["ordnung"])

        return ordnungen


    def _ist_aromatisch(self, ring_pfad, daten):
        """Hueckel-Pruefung ueber zwei Arten von pi-Elektronen-Beitraegen:

        1. Ring-Doppelbindungen (klassisches Kekule-Muster, z.B. Benzol)
           - jede liefert 2 pi-Elektronen.
        2. Ein freies Elektronenpaar an einem Ringatom OHNE eigene
           Ringdoppelbindung (z.B. das N in Pyrrol, das O in Furan) -
           liefert ebenfalls 2 pi-Elektronen.

        Jedes Ringatom muss ueber eine Ringdoppelbindung ODER ein freies
        Elektronenpaar am pi-System teilnehmen. Keine zwei
        Ringdoppelbindungen duerfen dasselbe Atom teilen (kumuliert).

        [patch_fixes] Zusaetzlich: mindestens eine Ring-Doppelbindung
        und hoechstens EIN Atom, das nur ueber ein freies Paar
        teilnimmt. Vorher galt z.B. ein O3-Ring oder ein S5-Ring (lauter
        Atome mit freien Paaren, keine Doppelbindung) als aromatisch und
        bekam +20 - echte Heteroaromaten (Pyrrol, Furan, Thiophen,
        Imidazol, Oxazol, ...) haben genau einen solchen Spender."""

        ordnungen = self._ring_bindungsordnungen(ring_pfad, daten)

        if len(ordnungen) != len(ring_pfad) or not ordnungen:
            return False, 0

        if not set(ordnungen) <= {1, 2}:
            return False, 0

        n = len(ordnungen)

        for i in range(n):
            if ordnungen[i] == 2 and ordnungen[(i - 1) % n] == 2:
                return False, 0  # kumulierte Doppelbindung am gemeinsamen Atom

        pi_elektronen = 0
        abgedeckt = set()
        ring_doppelbindungen = 0
        paar_spender = 0

        for i in range(n):

            if ordnungen[i] == 2:

                pi_elektronen += 2
                ring_doppelbindungen += 1
                abgedeckt.add(ring_pfad[i])
                abgedeckt.add(ring_pfad[(i + 1) % n])

        for atom_id in ring_pfad:

            if atom_id in abgedeckt:
                continue

            atom = daten["atom_by_id"][atom_id]
            oktett = ELEMENT_INFO.get(atom["element"], {}).get("oktett")

            if oktett is None:
                return False, 0  # unbekanntes Element - kein verlaessliches Urteil moeglich

            frei_paare = (oktett - 2 * atom["used_valence"]) // 2

            if frei_paare < 1:
                return False, 0  # weder Ringdoppelbindung noch freies Elektronenpaar

            pi_elektronen += 2
            paar_spender += 1

        if ring_doppelbindungen == 0 or paar_spender > 1:
            return False, 0

        hueckel = (pi_elektronen - 2) % 4 == 0

        return hueckel, pi_elektronen


    def modell_ringe(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        m = daten["molekuel"]
        anzahl_ringe = m["anzahl_bindungen"] - m["anzahl_atome"] + m["fragmente"]

        if anzahl_ringe <= 0:
            return e

        for ring_pfad in self._alle_ringe(daten):

            groesse = len(ring_pfad)

            abzug = self._RINGSPANNUNG.get(groesse)

            if abzug is None:
                abzug = -min(10 + (groesse - 8) * 2, 30)

            e["punkte"] += abzug

            # Unconditionaler, rein informativer Eintrag - unabhaengig
            # davon, ob die Ringgroesse einen Abzug oder Bonus gab.
            # Vorher wurde eine GUTE Ringgroesse (z.B. 6) nur als Text
            # in "erklaerung" erwaehnt, nirgends strukturiert erfasst -
            # ein Aufrufer konnte also nicht lernen, WELCHE Ringgroessen
            # sich lohnen, sondern nur, welche bestraft wurden. Aendert
            # die Bewertung (e["punkte"]) nicht, macht das Ergebnis nur
            # zusaetzlich maschinenlesbar (siehe KIGedaechtnis.beobachte).
            e["regeln"].append({
                "regel": "Ringgroesse",
                "ringgroesse": groesse
            })

            if abzug < 0:
                e["probleme"].append({
                    "modell": "Ringspannung",
                    "ringgroesse": groesse,
                    "abzug": abzug
                })
                e["erklaerung"].append(
                    f"Ring mit {groesse} Atomen gefunden - Ringspannung wahrscheinlich."
                )
            else:
                e["erklaerung"].append(
                    f"Ring mit {groesse} Atomen gefunden - günstige, praktisch spannungsfreie Ringgröße."
                )

            aromatisch, pi_elektronen = self._ist_aromatisch(ring_pfad, daten)

            if aromatisch:
                e["punkte"] += self._AROMATIZITAETS_BONUS
                e["regeln"].append({
                    "regel": "Aromatizität",
                    "ringgroesse": groesse,
                    "pi_elektronen": pi_elektronen
                })
                e["erklaerung"].append(
                    f"Ring mit {groesse} Atomen ist aromatisch ({pi_elektronen} π-Elektronen, "
                    f"erfüllt die 4n+2-Regel) - zusätzlich stabilisiert."
                )

        return e


    # ==========================================================
    # MODELL 12 — HYPERVALENZ
    # ==========================================================

    _NORMALE_VALENZ = {"H": 1, "C": 4, "N": 3, "O": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1, "P": 3, "S": 2, "Si": 4, "B": 3}

    def modell_hypervalenz(self, daten):

        e = {"punkte": 0, "regeln": [], "probleme": [], "erklaerung": []}

        for atom in daten["atome"]:

            element = atom["element"]
            normale_valenz = self._NORMALE_VALENZ.get(element)

            if normale_valenz is None:
                continue

            benutzt = atom["valenz"]["benutzt"]

            if benutzt > normale_valenz:

                e["punkte"] -= 20
                e["probleme"].append({
                    "modell": "Hypervalenz",
                    "atom": atom["id"],
                    "benutzt": benutzt,
                    "normale_valenz": normale_valenz
                })
                e["erklaerung"].append(
                    f"{element} nutzt {benutzt} Bindungen - mehr als die normale Valenz "
                    f"{normale_valenz} - hypervalent."
                )

        return e