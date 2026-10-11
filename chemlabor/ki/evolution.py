# ==========================================================
# evolution.py
# ==========================================================
#
# Evolutionaerer Generator: statt jedes Molekuel neu zu wuerfeln,
# verbessert er eine Population bereits bewerteter Molekuele durch
# Mutation und Rekombination. Bewertet wird weiterhin NUR vom Kritiker
# (hier als Punktzahl ueber melde() hereingereicht) - RDKit spielt in
# der Auswahl keine Rolle, damit nicht nur bekannte Chemie entsteht.
#
# Der Zufallsgenerator bleibt wichtig: er sorgt dafuer, dass die
# Ergebnisse nicht monoton werden. Er wird deshalb auf zwei Wegen
# eingebunden:
#
#   1. Grundrate: ein kleiner Anteil jeder Generation (Standard 5 %) sind
#      "Einwanderer" aus den Quellen (reiner Zufall und/oder KIGenerator).
#   2. Sackgasse: stagniert das beste Ergebnis ueber `sackgasse_generationen`
#      Generationen oder sinkt die Vielfalt (verschiedene Summenformeln)
#      unter `vielfalt_min`, wird ein groesserer Teil der naechsten
#      Generation durch Einwanderer ersetzt (Standard 30 %, bei weiterem
#      Stillstand mehr). Nach doppelter Stillstandsdauer folgt ein
#      Teilneustart: nur die Eliten bleiben.
#
# Einwanderer sind meist viel schlechter als die Elite und wuerden sofort
# aussortiert. Sie haben deshalb `schutz_generationen` Generationen Schutz:
# sie ueberleben die Auswahl und duerfen als Eltern mutieren.
#
# Benutzung (siehe ChemLabor.neues_molekuel_evolution):
#
#     kandidat = evolution.naechstes()                 # Molekuel erzeugen
#     ... Kritiker bewertet ...
#     evolution.melde(kandidat, punkte)                # Ergebnis zurueckmelden
#
# Ein nicht gemeldeter Kandidat (z.B. Fehler bei der Bewertung) zaehlt
# einfach nicht mit.

import copy
import random
from collections import Counter, defaultdict, deque


# (Name, Gewicht) der Mutationsoperatoren
OPERATOREN = (
    ("tausch", 3),
    ("anhaengen", 3),
    ("entfernen", 2),
    ("ordnung", 2),
    ("ring_zu", 1),
    ("ring_auf", 1),
)

# So viel besser muss das Beste werden, damit es nicht als Stillstand zaehlt
VERBESSERUNG_MIN = 0.5


# ----------------------------------------------------------
# Hilfsfunktionen auf dem Molekuelformat
# ----------------------------------------------------------

def _frei(atom):
    return atom["valence"] - atom["used_valence"]


def _komponente(molekuel, start, ohne=None):
    """IDs aller Atome, die von 'start' ueber Bindungen erreichbar sind.
    'ohne' ist der Index einer Bindung, die dabei ignoriert wird."""

    graph = defaultdict(set)

    for i, b in enumerate(molekuel["bonds"]):

        if i == ohne:
            continue

        graph[b["atom1"]].add(b["atom2"])
        graph[b["atom2"]].add(b["atom1"])

    besucht = {start}
    schlange = deque([start])

    while schlange:

        knoten = schlange.popleft()

        for nachbar in graph[knoten]:
            if nachbar not in besucht:
                besucht.add(nachbar)
                schlange.append(nachbar)

    return besucht


def _zusammenhaengend(molekuel, ohne=None):

    if not molekuel["atoms"]:
        return True

    start = molekuel["atoms"][0]["id"]

    return len(_komponente(molekuel, start, ohne)) == len(molekuel["atoms"])


def _brutto(molekuel):
    """Summenformel als Tupel - grobes Mass fuer die Vielfalt."""

    return tuple(sorted(Counter(a["element"] for a in molekuel["atoms"]).items()))


def _atom(molekuel, atom_id):

    return next(a for a in molekuel["atoms"] if a["id"] == atom_id)


# ----------------------------------------------------------
# Evolution
# ----------------------------------------------------------

class Evolution:

    def __init__(
            self,
            elemente,
            valenzen,
            bindungsarten,
            quellen,
            signatur,
            max_atome=15,
            min_atome=2,
            population=30,
            eliten=4,
            mutationsrate=0.5,
            rekombinationsrate=0.3,
            zufall_grundrate=0.05,
            sackgasse_generationen=12,
            sackgasse_anteil=0.3,
            schutz_generationen=3,
            vielfalt_min=0.4,
            plateau=95,
            nische_max=2,
            max_alter=8,
            on_generation=None,
            rng=None
    ):
        """elemente: Liste von Symbolen. valenzen: Symbol -> Valenz.
        bindungsarten: Symbol -> erlaubte Bindungsordnungen.
        quellen: Liste (name, funktion, gewicht); jede Funktion liefert ein
        neues Molekuel (Zufallsgenerator, KIGenerator ...).
        signatur: Funktion Molekuel -> hashbare Struktur-Signatur.
        on_generation: optional, bekommt nach jeder Generation ein dict."""

        self.elemente = [e for e in elemente if e in valenzen]
        self.valenzen = valenzen
        self.bindungsarten = bindungsarten
        self.quellen = list(quellen)
        self.signatur = signatur

        self.max_atome = max_atome
        self.min_atome = min_atome
        self.population_groesse = population
        self.eliten = max(1, min(eliten, population - 1))
        self.mutationsrate = mutationsrate
        self.rekombinationsrate = rekombinationsrate
        self.zufall_grundrate = zufall_grundrate
        self.sackgasse_generationen = max(2, sackgasse_generationen)
        self.sackgasse_anteil = sackgasse_anteil
        self.schutz_generationen = schutz_generationen
        self.vielfalt_min = vielfalt_min
        self.on_generation = on_generation

        # Punkte ab `plateau` gelten als gleich gut (Deckel der Bewertung).
        # Unter gleich guten gewinnt die seltener besuchte Summenformel
        # (Nische); je Nische bleiben hoechstens `nische_max` Individuen.
        # Plateau-Individuen gehen nach `max_alter` Generationen in die
        # Ruhmeshalle und werden gesperrt (0 = nie).
        self.plateau = plateau
        self.nische_max = max(1, nische_max)
        self.max_alter = max_alter

        self.ruhmeshalle = deque(maxlen=1000)
        self._tabu = set()
        self._besucht = Counter()
        self._top_nischen = set()

        self.rng = rng if rng is not None else random.Random()

        self.population = []          # ueberlebende Eintraege (bewertet)
        self.generation = 0

        self._warteschlange = deque()
        self._ausstehend = []
        self._sigs = set()

        self._bester = -1.0
        self._stillstand = 0
        self._vielfalt = 1.0
        self._ereignis = None

    # ------------------------------------------------------
    # Oeffentliche Schnittstelle
    # ------------------------------------------------------

    def besiedle(self, paare):
        """Fuegt bereits bewertete Molekuele (molekuel, punkte) als
        Startpopulation hinzu (z.B. die Besten der Beobachtungsphase)."""

        bekannt = {e["sig"] for e in self.population}

        for molekuel, punkte in sorted(paare, key=lambda p: p[1], reverse=True):

            eintrag = self._eintrag(molekuel, "start")

            if eintrag["sig"] in bekannt:
                continue

            eintrag["punkte"] = float(punkte)
            bekannt.add(eintrag["sig"])
            self.population.append(eintrag)

            if len(self.population) >= self.population_groesse:
                break

        if self.population:
            self._bester = max(e["punkte"] for e in self.population)

        for e in self.population:
            nische = _brutto(e["mol"])
            self._besucht[nische] += 1
            if e["punkte"] >= self.plateau:
                self._top_nischen.add(nische)

    def naechstes(self):
        """Das naechste zu bewertende Molekuel."""

        while not self._warteschlange:
            self._abschliessen()
            self._neue_generation()

        eintrag = self._warteschlange.popleft()
        self._ausstehend.append(eintrag)

        return eintrag["mol"]

    def melde(self, molekuel, punkte):
        """Bewertung eines zuvor mit naechstes() gelieferten Molekuels."""

        for eintrag in self._ausstehend:
            if eintrag["mol"] is molekuel:
                eintrag["punkte"] = float(punkte)
                return

    def bestes(self):
        """(molekuel, punkte) des besten Individuums oder None."""

        bewertet = [
            e for e in self.population + self._ausstehend
            if e["punkte"] is not None
        ]

        if not bewertet:
            return None

        bester = max(bewertet, key=lambda e: e["punkte"])

        return bester["mol"], bester["punkte"]

    # ------------------------------------------------------
    # Generationen
    # ------------------------------------------------------

    def _eintrag(self, molekuel, herkunft, schutz=0):

        return {
            "mol": molekuel,
            "sig": self.signatur(molekuel),
            "punkte": None,
            "herkunft": herkunft,
            "schutz": schutz,
            "alter": 0,
        }

    def _abschliessen(self):
        """Wertet die gerade fertige Generation aus und waehlt die
        Ueberlebenden."""

        bewertet = [e for e in self._ausstehend if e["punkte"] is not None]
        herkuenfte = [e["herkunft"] for e in bewertet]

        self._ausstehend = []

        if not bewertet and not self.population:
            return

        # je Struktur nur das beste Individuum behalten
        je_sig = {}

        for e in self.population + bewertet:

            alt = je_sig.get(e["sig"])

            if alt is None or e["punkte"] > alt["punkte"]:
                je_sig[e["sig"]] = e

        for e in bewertet:
            self._besucht[_brutto(e["mol"])] += 1

        kandidaten = sorted(je_sig.values(), key=self._schluessel, reverse=True)

        # Ruhestand: Plateau-Individuen, die lange genug dabei waren
        if self.max_alter:

            aktiv = []

            for e in kandidaten:
                if e["alter"] >= self.max_alter and e["punkte"] >= self.plateau:
                    self._pensioniere(e)
                else:
                    aktiv.append(e)

            kandidaten = aktiv or kandidaten

        # Nischenobergrenze
        gezaehlt = Counter()
        verteilt = []

        for e in kandidaten:

            nische = _brutto(e["mol"])

            if gezaehlt[nische] >= self.nische_max:
                continue

            gezaehlt[nische] += 1
            verteilt.append(e)

        kandidaten = verteilt

        eliten = kandidaten[:self.eliten]
        rest = kandidaten[self.eliten:]

        geschuetzt = [e for e in rest if e["schutz"] > 0][:max(1, self.population_groesse // 4)]

        gewaehlt = eliten + geschuetzt

        for e in rest:

            if len(gewaehlt) >= self.population_groesse:
                break

            if e not in gewaehlt:
                gewaehlt.append(e)

        for e in gewaehlt:
            e["alter"] += 1
            e["schutz"] = max(0, e["schutz"] - 1)

        self.population = gewaehlt

        bestes = max(e["punkte"] for e in gewaehlt)

        neue_nischen = {
            _brutto(e["mol"]) for e in bewertet if e["punkte"] >= self.plateau
        } - self._top_nischen

        self._top_nischen |= neue_nischen

        # Fortschritt: besseres Ergebnis ODER eine neue Summenformel, die
        # das Plateau erreicht (sonst waere bei Punktzahl 100 jede
        # Generation ein "Stillstand").
        if bestes > self._bester + VERBESSERUNG_MIN or neue_nischen:
            self._bester = max(self._bester, bestes)
            self._stillstand = 0
        else:
            self._stillstand += 1

        self._vielfalt = len({_brutto(e["mol"]) for e in gewaehlt}) / len(gewaehlt)

        if self.on_generation is not None and bewertet:

            self.on_generation({
                "generation": self.generation,
                "bestes": bestes,
                "mittel": sum(e["punkte"] for e in gewaehlt) / len(gewaehlt),
                "vielfalt": self._vielfalt,
                "einwanderer": sum(1 for h in herkuenfte if h in self._quellnamen()),
                "population": len(gewaehlt),
                "nischen": len(self._top_nischen),
                "ruhmeshalle": len(self.ruhmeshalle),
                "ereignis": self._ereignis,
            })

    def _schluessel(self, eintrag):
        """Sortierschluessel der Auswahl: Punkte (am Plateau gedeckelt),
        dann die seltener besuchte Nische."""

        return (
            min(eintrag["punkte"], self.plateau),
            -self._besucht[_brutto(eintrag["mol"])],
        )

    def _pensioniere(self, eintrag):

        self.ruhmeshalle.append((eintrag["mol"], eintrag["punkte"]))

        if len(self._tabu) > 200000:
            self._tabu.clear()

        self._tabu.add(eintrag["sig"])

    def _quellnamen(self):

        return {name for name, _f, _g in self.quellen}

    def _neue_generation(self):
        """Fuellt die Warteschlange mit der naechsten Generation."""

        self.generation += 1
        self._ereignis = None

        anteil = self.zufall_grundrate

        if self.population:

            if self._stillstand >= 2 * self.sackgasse_generationen:

                self._ereignis = "neustart"
                self.population = self.population[:self.eliten]
                anteil = 1.0 - self.eliten / self.population_groesse
                self._stillstand = 0

            elif (
                self._stillstand >= self.sackgasse_generationen
                or (self.generation > 3 and self._vielfalt < self.vielfalt_min)
            ):

                self._ereignis = "sackgasse"

                mehr = max(0, self._stillstand - self.sackgasse_generationen)
                anteil = min(0.9, self.sackgasse_anteil + 0.2 * mehr)

        self._sigs = {e["sig"] for e in self.population}

        for _ in range(self.population_groesse):

            eintrag = None

            if not self.population or self.rng.random() < anteil:
                eintrag = self._einwanderer()
            else:
                eintrag = self._nachkomme() or self._einwanderer()

            if eintrag is None:
                continue

            self._sigs.add(eintrag["sig"])
            self._warteschlange.append(eintrag)

        if not self._warteschlange:
            raise RuntimeError("Evolution: keine Kandidaten erzeugbar (Quellen liefern nichts).")

    # ------------------------------------------------------
    # Einwanderer (Zufall / KIGenerator)
    # ------------------------------------------------------

    def _einwanderer(self):

        for _ in range(6):

            name, quelle, _gewicht = self.rng.choices(
                self.quellen, weights=[g for _n, _f, g in self.quellen], k=1
            )[0]

            try:
                molekuel = quelle()
            except Exception:
                continue

            if not molekuel or not molekuel.get("atoms"):
                continue

            eintrag = self._eintrag(molekuel, name, schutz=self.schutz_generationen)

            if eintrag["sig"] in self._sigs or eintrag["sig"] in self._tabu:
                continue

            return eintrag

        return None

    # ------------------------------------------------------
    # Nachkommen
    # ------------------------------------------------------

    def _elternteil(self):

        geschuetzt = [e for e in self.population if e["schutz"] > 0]

        if geschuetzt and self.rng.random() < 0.25:
            return self.rng.choice(geschuetzt)

        turnier = [self.rng.choice(self.population) for _ in range(3)]

        return max(turnier, key=self._schluessel)

    def _nachkomme(self):

        for _ in range(5):

            herkunft = "mutation"

            if len(self.population) >= 2 and self.rng.random() < self.rekombinationsrate:

                a = self._elternteil()
                b = self._elternteil()

                kind = self.rekombiniere(a["mol"], b["mol"])
                herkunft = "rekombination"

                if kind is not None and self.rng.random() < self.mutationsrate:
                    kind = self.mutiere(kind) or kind

            else:

                kind = self.mutiere(self._elternteil()["mol"])

            if kind is None:
                continue

            eintrag = self._eintrag(kind, herkunft)

            if eintrag["sig"] in self._sigs or eintrag["sig"] in self._tabu:
                continue

            return eintrag

        return None

    # ------------------------------------------------------
    # Mutation
    # ------------------------------------------------------

    def mutiere(self, molekuel):
        """Kopie mit EINER Mutation oder None, wenn keine moeglich war."""

        namen = [n for n, _g in OPERATOREN]
        gewichte = [g for _n, g in OPERATOREN]

        for _ in range(8):

            kopie = copy.deepcopy(molekuel)
            operator = self.rng.choices(namen, weights=gewichte, k=1)[0]

            if getattr(self, f"_mut_{operator}")(kopie):
                return kopie

        return None

    def _erlaubt(self, element):

        return self.bindungsarten.get(element, [1])

    def _mut_tausch(self, mol):

        atome = list(mol["atoms"])
        self.rng.shuffle(atome)

        for atom in atome:

            ordnungen = [
                b["order"] for b in mol["bonds"]
                if atom["id"] in (b["atom1"], b["atom2"])
            ]

            moeglich = [
                e for e in self.elemente
                if e != atom["element"]
                and self.valenzen[e] >= atom["used_valence"]
                and all(o in self._erlaubt(e) for o in ordnungen)
            ]

            if moeglich:
                neu = self.rng.choice(moeglich)
                atom["element"] = neu
                atom["valence"] = self.valenzen[neu]
                return True

        return False

    def _mut_anhaengen(self, mol):

        if len(mol["atoms"]) >= self.max_atome:
            return False

        offene = [a for a in mol["atoms"] if _frei(a) >= 1]

        if not offene:
            return False

        ziel = self.rng.choice(offene)
        element = self.rng.choice(self.elemente)

        erlaubt = [
            o for o in self._erlaubt(element)
            if o in self._erlaubt(ziel["element"])
            and o <= min(_frei(ziel), self.valenzen[element])
        ] or [1]

        ordnung = 1 if (1 in erlaubt and self.rng.random() < 0.7) else self.rng.choice(erlaubt)

        neue_id = max(a["id"] for a in mol["atoms"]) + 1

        mol["atoms"].append({
            "id": neue_id,
            "element": element,
            "valence": self.valenzen[element],
            "used_valence": ordnung,
        })
        mol["bonds"].append({"atom1": ziel["id"], "atom2": neue_id, "order": ordnung})
        ziel["used_valence"] += ordnung

        return True

    def _mut_entfernen(self, mol):

        if len(mol["atoms"]) <= self.min_atome:
            return False

        grad = Counter()

        for b in mol["bonds"]:
            grad[b["atom1"]] += 1
            grad[b["atom2"]] += 1

        blaetter = [a for a in mol["atoms"] if grad[a["id"]] == 1]

        if not blaetter:
            return False

        blatt = self.rng.choice(blaetter)

        bindung = next(
            b for b in mol["bonds"] if blatt["id"] in (b["atom1"], b["atom2"])
        )

        partner = bindung["atom2"] if bindung["atom1"] == blatt["id"] else bindung["atom1"]

        _atom(mol, partner)["used_valence"] -= bindung["order"]

        mol["bonds"].remove(bindung)
        mol["atoms"].remove(blatt)

        return True

    def _mut_ordnung(self, mol):

        bindungen = list(mol["bonds"])
        self.rng.shuffle(bindungen)

        for b in bindungen:

            a1 = _atom(mol, b["atom1"])
            a2 = _atom(mol, b["atom2"])

            erlaubt = set(self._erlaubt(a1["element"])) & set(self._erlaubt(a2["element"]))

            optionen = []

            for neu in (b["order"] - 1, b["order"] + 1):

                if neu not in erlaubt or neu < 1:
                    continue

                diff = neu - b["order"]

                if diff > 0 and (_frei(a1) < diff or _frei(a2) < diff):
                    continue

                optionen.append(neu)

            if optionen:
                neu = self.rng.choice(optionen)
                diff = neu - b["order"]
                a1["used_valence"] += diff
                a2["used_valence"] += diff
                b["order"] = neu
                return True

        return False

    def _mut_ring_zu(self, mol):

        bestehend = {frozenset((b["atom1"], b["atom2"])) for b in mol["bonds"]}
        offene = [a for a in mol["atoms"] if _frei(a) >= 1]

        paare = [
            (a, b) for i, a in enumerate(offene) for b in offene[i + 1:]
            if frozenset((a["id"], b["id"])) not in bestehend
        ]

        if not paare:
            return False

        a, b = self.rng.choice(paare)

        mol["bonds"].append({"atom1": a["id"], "atom2": b["id"], "order": 1})
        a["used_valence"] += 1
        b["used_valence"] += 1

        return True

    def _mut_ring_auf(self, mol):

        kandidaten = [
            i for i in range(len(mol["bonds"]))
            if _zusammenhaengend(mol, ohne=i)
        ]

        if not kandidaten:
            return False

        bindung = mol["bonds"].pop(self.rng.choice(kandidaten))

        _atom(mol, bindung["atom1"])["used_valence"] -= bindung["order"]
        _atom(mol, bindung["atom2"])["used_valence"] -= bindung["order"]

        return True

    # ------------------------------------------------------
    # Rekombination
    # ------------------------------------------------------

    def _bruecken(self, mol):
        """Indizes der Bindungen, deren Entfernen das Molekuel teilt."""

        return [
            i for i in range(len(mol["bonds"]))
            if not _zusammenhaengend(mol, ohne=i)
        ]

    def _haelfte(self, mol):
        """Zufaellige Haelfte entlang einer Bruecke:
        (atome, bindungen, schnittatom_id) mit bereits freigegebener
        Valenz am Schnittatom - oder None."""

        bruecken = self._bruecken(mol)

        if not bruecken:
            return None

        index = self.rng.choice(bruecken)
        schnitt = mol["bonds"][index]

        start = self.rng.choice((schnitt["atom1"], schnitt["atom2"]))
        ids = _komponente(mol, start, ohne=index)

        atome = [copy.deepcopy(a) for a in mol["atoms"] if a["id"] in ids]
        bindungen = [
            copy.deepcopy(b) for i, b in enumerate(mol["bonds"])
            if i != index and b["atom1"] in ids and b["atom2"] in ids
        ]

        next(a for a in atome if a["id"] == start)["used_valence"] -= schnitt["order"]

        return atome, bindungen, start

    def rekombiniere(self, mol_a, mol_b):
        """Kind aus je einer Haelfte beider Eltern, verbunden durch eine
        Einfachbindung - oder None (keine Bruecke, zu gross/klein)."""

        a = self._haelfte(mol_a)
        b = self._haelfte(mol_b)

        if a is None or b is None:
            return None

        atome_a, bindungen_a, schnitt_a = a
        atome_b, bindungen_b, schnitt_b = b

        if not (self.min_atome <= len(atome_a) + len(atome_b) <= self.max_atome):
            return None

        versatz = max(x["id"] for x in atome_a) + 1
        abbildung = {x["id"]: versatz + i for i, x in enumerate(atome_b)}

        for x in atome_b:
            x["id"] = abbildung[x["id"]]

        for y in bindungen_b:
            y["atom1"] = abbildung[y["atom1"]]
            y["atom2"] = abbildung[y["atom2"]]

        kind = {"atoms": atome_a + atome_b, "bonds": bindungen_a + bindungen_b}

        ende_b = abbildung[schnitt_b]

        _atom(kind, schnitt_a)["used_valence"] += 1
        _atom(kind, ende_b)["used_valence"] += 1
        kind["bonds"].append({"atom1": schnitt_a, "atom2": ende_b, "order": 1})

        return kind