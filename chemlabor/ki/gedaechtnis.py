import math
import random
import json
from collections import defaultdict, deque, Counter
from datetime import datetime, timezone

from chemlabor.chemie.zufallsgenerator import VALENZEN, BINDUNGSARTEN, MOLEKUEL_MAX_ATOME
from chemlabor.texte import t

# ==========================================================
# KIGedaechtnis
# ==========================================================
#
# Kein LLM. Ein einfaches Gedaechtnis, das waehrend der
# Beobachtungsphase (Generator -> Kritiker) mitzaehlt:
#
#   - welche Bindungen (Element1-Element2, Bindungsordnung)
#     in eher stabilen Molekuelen vorkommen
#   - welche used_valence-Werte pro Element in eher stabilen
#     Molekuelen typisch sind
#   - welche Elemente ueberhaupt in stabilen Molekuelen auftauchen
#   - welche Ringgroessen im Schnitt gut abschneiden
#
# Aus diesen Zaehlungen baut KIGenerator danach eigene Molekuele.
#
# ----------------------------------------------------------
# v2: Qualitaet statt Haeufigkeit
# ----------------------------------------------------------
#
# Vorher wurde pro Muster eine punkte-gewichtete SUMME akkumuliert,
# mit einer binaeren Schwelle (stabilitaets_schwelle=70). Jetzt wird
# pro Muster/Wert (Summe, Anzahl) getrennt gefuehrt und daraus ein
# Bayes-geglaetteter DURCHSCHNITT gebildet (_glaetten): wenig
# beobachtete Muster werden Richtung Gesamtdurchschnitt gezogen, damit
# ein einzelner Zufallstreffer nicht sofort als "bestes Muster" gilt.
#
# ----------------------------------------------------------
# v3: Ringgroessen-Lernen
# ----------------------------------------------------------
#
# Der Kritiker bewertet Ringspannung und Aromatizitaet laengst
# (siehe kritiker.py, modell_ringe) und meldet seit dem juengsten Fix
# fuer JEDEN gefundenen Ring einen strukturierten
# {"regel": "Ringgroesse", "ringgroesse": n}-Eintrag in
# analyse["regeln"] - unabhaengig davon, ob die Groesse bestraft oder
# belohnt wurde. KIGedaechtnis wertet das jetzt aus (ringgroessen_
# nutzung) und KIGenerator nutzt es beim Bauen aktiv: bevor eine neue
# Bindung gesetzt wird, die einen Ring schliessen wuerde, prueft
# _zusatzbindungen(), ob die resultierende Ringgroesse laut gelernter
# Erfahrung gut abschneidet - und akzeptiert oder verwirft den
# Ringschluss entsprechend stochastisch, statt blind zu verbinden und
# erst hinterher per Best-of-N/Kritiker zu merken, dass z.B. ein
# Dreiring entstanden ist.
#
# ----------------------------------------------------------
# v4: Softmax/Temperatur und Baseline-Abzug
# ----------------------------------------------------------
#
# Vorher gingen die geglaetteten Durchschnittspunkte (0..100) direkt
# als Gewichte in random.choices() bzw. als Ringakzeptanz ein - ein
# Element mit Ø70 wurde dadurch nur ~1,2x so oft gewaehlt wie eins mit
# Ø58, und ein 3-Ring mit Ø35 wurde noch zu 35% akzeptiert. Jetzt
# zaehlt die ABWEICHUNG vom Gesamtmittel (Baseline-Abzug), und sie
# wird ueber exp(abweichung / T) in ein Gewicht umgerechnet (Softmax).
# Die Temperatur T faellt waehrend des Lernens bis temperatur_min:
# anfangs erkunden, spaeter ausnutzen.
# ==========================================================


def molekuel_signatur(molekuel):
    """Naeherungsweise kanonische Signatur eines Molekuels - robust
    gegenueber vertauschten Atom-IDs (zwei strukturell identische
    Molekuele mit unterschiedlicher ID-Reihenfolge ergeben dieselbe
    Signatur), aber KEINE vollstaendige Graph-Isomorphie-Pruefung.
    Jedes Atom wird durch (Element, sortierte Liste seiner
    (Nachbarelement, Bindungsordnung)-Paare) beschrieben; die Menge
    aller Atom-Signaturen, sortiert, ergibt die Molekuel-Signatur.
    Fuer die seltenen Faelle echter Automorphismen (z.B. hochsymmetrische
    Ringe) kann das zwei tatsaechlich unterschiedliche Strukturen
    faelschlich gleichsetzen - fuer eine Diversitaets-Heuristik (nicht
    fuer die eigentliche Chemie-Bewertung) ist das ein akzeptabler
    Kompromiss gegenueber einer vollen Isomorphie-Pruefung."""

    atome_by_id = {a["id"]: a for a in molekuel["atoms"]}

    nachbarn = defaultdict(list)

    for b in molekuel["bonds"]:

        e1 = atome_by_id[b["atom1"]]["element"]
        e2 = atome_by_id[b["atom2"]]["element"]

        nachbarn[b["atom1"]].append((e2, b["order"]))
        nachbarn[b["atom2"]].append((e1, b["order"]))

    atom_signaturen = [
        (atom["element"], tuple(sorted(nachbarn.get(atom["id"], []))))
        for atom in molekuel["atoms"]
    ]

    return tuple(sorted(atom_signaturen))


def _signatur_zu_json(signatur):
    """molekuel_signatur() liefert verschachtelte Tupel - JSON kennt
    weder Tupel noch Tupel-als-dict-Schluessel. Wandelt eine Signatur
    in eine reine Listenstruktur um (speichern())."""

    return [
        [element, [[nachbar_el, ordnung] for nachbar_el, ordnung in nachbarn]]
        for element, nachbarn in signatur
    ]


def _signatur_von_json(daten):
    """Kehrt _signatur_zu_json() um (laden())."""

    return tuple(
        (element, tuple((nachbar_el, ordnung) for nachbar_el, ordnung in nachbarn))
        for element, nachbarn in daten
    )


class KIGedaechtnis:

    def __init__(
        self,
        glaettung=5,
        ziel_perzentil=0.75,
        wiederholungsstrafe_basis=0.5,
        max_signaturen=20000,
        lernkurve_intervall=1000,
        max_lernkurve_punkte=2000,
        temperatur=10.0,
        temperatur_min=3.0,
        temperatur_zerfall=0.9999
    ):

        # Bayes-Glaettungskonstante: Anzahl "virtueller" Beobachtungen
        # am Gesamtdurchschnitt, mit denen jedes Muster startet. Hoeher
        # = vorsichtiger gegenueber wenig beobachteten Mustern.
        self.glaettung = glaettung

        # Perzentil, das stabilitaets_schwelle_adaptiv() standardmaessig
        # zurueckgibt (0.75 = oberes Quartil).
        self.ziel_perzentil = ziel_perzentil

        # Basis der exponentiellen Wiederholungsstrafe (0 < basis < 1).
        # Ein Molekuel, das schon n-mal beobachtet wurde, bekommt bei
        # der Best-of-N-Kandidatenauswahl den Faktor basis**n auf seine
        # Punktzahl - kleiner basis = haerte Strafe. 0.5 halbiert den
        # effektiven Auswahlwert bei jeder weiteren Wiederholung.
        self.wiederholungsstrafe_basis = wiederholungsstrafe_basis

        # Obergrenze fuer molekuel_signaturen (siehe beobachte() und
        # _signaturen_begrenzen()) - ohne diese Grenze waechst der
        # Zaehler in einem langen Endlos-Lauf unbegrenzt mit jeder neu
        # erzeugten Struktur, selbst wenn kaum je etwas wiederholt
        # wird (50000 realistisch variierte Testbeobachtungen erzeugten
        # ~25000 Eintraege - linear, kein Plateau). Alle anderen
        # wachsenden Sammlungen im Projekt sind bereits begrenzt
        # (alle_ergebnisse: deque(maxlen=5000), GUI-Ringpuffer analog) -
        # hier fehlte das bisher.
        self.max_signaturen = max_signaturen

        # Nach wie vielen Beobachtungen jeweils ein Lernkurven-Punkt
        # abgelegt wird (siehe beobachte()/_lernkurve_aktualisieren).
        self.lernkurve_intervall = lernkurve_intervall

        # Softmax-Temperatur in Punkten: Ein Muster, das T Punkte ueber
        # dem Gesamtmittel liegt, wird e-mal (~2,7x) haeufiger gewaehlt
        # als eines genau im Mittel. Hoch = erkunden, niedrig =
        # ausnutzen. Faellt in beobachte() pro Beobachtung um den
        # Faktor temperatur_zerfall bis hinunter zu temperatur_min.
        self.temperatur = temperatur
        self.temperatur_min = temperatur_min
        self.temperatur_zerfall = temperatur_zerfall

        self.beobachtungen = 0

        # Lernkurve: (beobachtungen_stand, durchschnitt_punkte_im_intervall)
        # - ein Sample alle lernkurve_intervall Beobachtungen, NICHT der
        # laufende Gesamtdurchschnitt (der wuerde sich mit wachsendem n
        # kaum noch bewegen und jede spaete Verbesserung/Verschlechterung
        # verschleiern). deque mit maxlen begrenzt den Speicher analog zu
        # max_signaturen - bei sehr langen Endlos-Laeufen soll das nicht
        # unbegrenzt wachsen.
        self.lernkurve = deque(maxlen=max_lernkurve_punkte)

        # Lauf-Historie: ein Eintrag pro starten()-Aufruf (siehe
        # ChemLabor.starten()/lauf_beginnen()) mit Zeitpunkt, den zu
        # dem Zeitpunkt geltenden Generator-Einstellungen und dem
        # Beobachtungsstand, ab dem dieser Lauf begann. Macht aus dem
        # Gedaechtnis kein anonymes Zahlenkonvolut mehr, sondern
        # nachvollziehbar, aus welchen Durchlaeufen mit welchen
        # Variablen sich das aktuelle "Kollektivwissen" zusammensetzt -
        # v.a. relevant, wenn man ueber viele Sitzungen hinweg immer
        # wieder dasselbe gespeicherte Gedaechtnis laedt und mit
        # unterschiedlichen Einstellungen weitertrainiert (siehe
        # zusammenfuehren() fuer das Verschmelzen mehrerer SEPARAT
        # gespeicherter Gedaechtnisdateien).
        self.lauf_historie = deque(maxlen=500)

        # Akkumulatoren fuer das laufende, noch nicht abgeschlossene
        # Intervall.
        self._intervall_summe = 0.0
        self._intervall_anzahl = 0

        # (element1, element2, ordnung) -> [summe_punkte, anzahl]
        self.bindungsmuster = defaultdict(lambda: [0.0, 0])

        # element -> {used_valence -> [summe_punkte, anzahl]}
        self.valenz_nutzung = defaultdict(lambda: defaultdict(lambda: [0.0, 0]))

        # element -> {Anzahl_Bindungspartner -> [summe_punkte, anzahl]}
        # (Domaenen im VSEPR-Sinn: zaehlt NACHBARATOME, nicht
        # Bindungsordnung - eine Doppelbindung ist immer noch 1 Domaene.
        # Bewusst getrennt von valenz_nutzung: zwei Molekuele koennen
        # dieselbe used_valence haben, aber unterschiedlich viele
        # Bindungspartner - und genau DAS entscheidet ueber die Form.)
        self.domänen_nutzung = defaultdict(lambda: defaultdict(lambda: [0.0, 0]))

        # ringgroesse -> [summe_punkte, anzahl] - aus den vom Kritiker
        # gemeldeten "Ringgroesse"-Eintraegen in analyse["regeln"].
        self.ringgroessen_nutzung = defaultdict(lambda: [0.0, 0])

        # element -> Haeufigkeit (reine Vorkommenszaehlung, fuer
        # Reporting - NICHT die Auswahlgrundlage fuer den Generator,
        # dafuer gibt es _element_qualitaet).
        self.elemente_haeufigkeit = Counter()

        # element -> [summe_punkte, anzahl] - Qualitaetsgrundlage fuer
        # bevorzugte_elemente().
        self._element_qualitaet = defaultdict(lambda: [0.0, 0])

        # Laufende Summe aller Punkte (einmal pro Molekuel, nicht pro
        # Atom/Bindung) als Bayes-Prior fuer wenig beobachtete Muster.
        self._punkte_summe_global = 0.0

        # Grobes Histogramm (11 Buckets, 0-100 in 10er-Schritten) fuer
        # eine adaptive Stabilitaetsschwelle, ohne alle Einzelwerte
        # speichern zu muessen.
        self._histogramm = [0] * 11

        # signatur -> Anzahl bisher beobachteter Molekuele mit exakt
        # dieser Struktur (siehe molekuel_signatur()) - Grundlage fuer
        # wiederholungsstrafe().
        self.molekuel_signaturen = Counter()


    # ------------------------------------------------------
    # Ein Experiment (Molekuel + Kritiker-Analyse) einlernen
    # ------------------------------------------------------

    def beobachte(self, molekuel, analyse, vollstaendig=True):
        """vollstaendig=False: zaehlt das Experiment nur in die globalen
        Statistiken (Durchschnitt, Histogramm, Lernkurve, adaptive
        Schwelle) ein, lernt aber KEINE Bindungsmuster/Valenz-/Domaenen-/
        Ringgroessen-/Signaturstatistik daraus. Gedacht fuer den
        Verhaeltnis-Filter in ChemLabor (siehe labor.py,
        ZIEL_VERHAELTNIS_NEG_ZU_POS): ueberschuessige schlechte
        Molekuele sollen die Verteilung/den Trend weiter korrekt
        widerspiegeln, aber nicht die gelernten Muster mit denselben
        paar schlechten Bindungsideen zumuellen."""

        self.beobachtungen += 1

        # Temperatur abkuehlen: erst erkunden, dann ausnutzen.
        self.temperatur = max(
            self.temperatur_min,
            self.temperatur * self.temperatur_zerfall
        )

        punkte = analyse["bewertung"]["punkte"]

        self._punkte_summe_global += punkte
        self._histogramm[min(int(punkte) // 10, 10)] += 1

        self._intervall_summe += punkte
        self._intervall_anzahl += 1

        if self.beobachtungen % self.lernkurve_intervall == 0:
            self._lernkurve_abschliessen()

        if vollstaendig:
            self._lerne_muster(molekuel, analyse, punkte)


    def _lerne_muster(self, molekuel, analyse, punkte):
        """Der teure Teil von beobachte(): fuellt Bindungsmuster-,
        Valenz-, Domaenen-, Ringgroessen- und Signaturstatistik.
        Ausgelagert, damit ChemLabor ihn gezielt ueberspringen kann
        (siehe beobachte())."""

        atome_by_id = {a["id"]: a for a in molekuel["atoms"]}

        nachbarn_pro_atom = defaultdict(set)
        for bindung in molekuel["bonds"]:
            nachbarn_pro_atom[bindung["atom1"]].add(bindung["atom2"])
            nachbarn_pro_atom[bindung["atom2"]].add(bindung["atom1"])

        for atom in molekuel["atoms"]:

            element = atom["element"]

            self.elemente_haeufigkeit[element] += 1

            eq = self._element_qualitaet[element]
            eq[0] += punkte
            eq[1] += 1

            ev = self.valenz_nutzung[element][atom["used_valence"]]
            ev[0] += punkte
            ev[1] += 1

            domänen = len(nachbarn_pro_atom.get(atom["id"], ()))
            if domänen > 0:
                ed = self.domänen_nutzung[element][domänen]
                ed[0] += punkte
                ed[1] += 1

        for bindung in molekuel["bonds"]:

            e1 = atome_by_id[bindung["atom1"]]["element"]
            e2 = atome_by_id[bindung["atom2"]]["element"]

            muster = tuple(sorted([e1, e2])) + (bindung["order"],)

            eb = self.bindungsmuster[muster]
            eb[0] += punkte
            eb[1] += 1

        # Ringgroessen aus den strukturierten Kritiker-Meldungen
        # herausziehen (siehe kritiker.py modell_ringe: ein
        # "Ringgroesse"-Eintrag pro gefundenem Ring, unabhaengig von
        # Bestrafung/Bonus).
        for regel in analyse.get("regeln", []):

            if regel.get("regel") == "Ringgroesse":

                er = self.ringgroessen_nutzung[regel["ringgroesse"]]
                er[0] += punkte
                er[1] += 1

        # Signatur dieses Molekuels mitzaehlen - Grundlage fuer die
        # exponentielle Wiederholungsstrafe bei der naechsten
        # Best-of-N-Auswahl (siehe wiederholungsstrafe()).
        self.molekuel_signaturen[molekuel_signatur(molekuel)] += 1

        self._signaturen_begrenzen()


    # ------------------------------------------------------
    # Bayes-geglaetteter Durchschnitt
    # ------------------------------------------------------

    def _globaler_mittelwert(self):
        """Durchschnittspunktzahl ueber alle bisher beobachteten
        Molekuele - die Baseline, von der Abweichungen gemessen werden
        (50.0 als neutraler Startwert, solange nichts beobachtet
        wurde)."""

        return (
            self._punkte_summe_global / self.beobachtungen
            if self.beobachtungen else 50.0
        )


    def _glaetten(self, summe, anzahl):
        """Zieht wenig beobachtete Muster Richtung Gesamtdurchschnitt,
        gewichtet mit 'glaettung' virtuellen Pseudo-Beobachtungen -
        verhindert, dass ein Muster, das nur ein- oder zweimal zufaellig
        mit hoher Punktzahl auftauchte, sofort als bestes Muster gilt."""

        globaler_mittelwert = self._globaler_mittelwert()

        return (
            (summe + self.glaettung * globaler_mittelwert)
            / (anzahl + self.glaettung)
        )


    def gewicht(self, summe, anzahl):
        """Softmax-Gewicht eines Musters: exp((geglaetteter Schnitt -
        Gesamtmittel) / T). Muster im Gesamtmittel bekommen Gewicht 1,
        T Punkte darueber ~2,7, T Punkte darunter ~0,37. Der Exponent
        wird begrenzt, damit nichts ueberlaeuft oder auf exakt 0
        faellt."""

        abweichung = self._glaetten(summe, anzahl) - self._globaler_mittelwert()

        return math.exp(max(-20.0, min(20.0, abweichung / self.temperatur)))


    # ------------------------------------------------------
    # Gelernte Muster abfragen
    # ------------------------------------------------------

    def bevorzugte_valenz(self, element, standard):
        """Used-valence-Wert mit der hoechsten geglaetteten
        Durchschnittspunktzahl fuer ein Element, oder 'standard' falls
        noch nichts beobachtet wurde. Waehlt die BESTE beobachtete
        Auspraegung, nicht mehr die haeufigste."""

        verteilung = self.valenz_nutzung.get(element)

        if not verteilung:
            return standard

        return max(
            verteilung,
            key=lambda wert: self._glaetten(*verteilung[wert])
        )


    def bevorzugte_domänen(self, element, standard):
        """Analog zu bevorzugte_valenz(), aber fuer die Anzahl der
        Bindungspartner (Domaenen im VSEPR-Sinn) - entscheidet ueber
        die FORM (viele schwache vs. wenige starke Bindungen), nicht
        nur ueber die genutzte Gesamtvalenz."""

        verteilung = self.domänen_nutzung.get(element)

        if not verteilung:
            return standard

        return max(
            verteilung,
            key=lambda wert: self._glaetten(*verteilung[wert])
        )


    def bevorzugte_elemente(self):
        """Liste von (element, gewicht), absteigend nach Gewicht
        sortiert. Das Gewicht ist das Softmax-Gewicht (siehe gewicht())
        der Bayes-geglaetteten Durchschnittspunktzahl der Molekuele, in
        denen das Element vorkam - NICHT die reine Haeufigkeit. None,
        falls noch nichts beobachtet wurde."""

        if not self.elemente_haeufigkeit:
            return None

        ergebnisse = [
            (element, self.gewicht(*self._element_qualitaet[element]))
            for element in self.elemente_haeufigkeit
        ]

        ergebnisse.sort(key=lambda eintrag: eintrag[1], reverse=True)

        return ergebnisse


    def bevorzugte_bindungen(self):
        """Liste von ((element1, element2, ordnung), gewicht),
        absteigend nach Softmax-Gewicht (siehe gewicht()) der
        Bayes-geglaetteten Durchschnittspunktzahl sortiert."""

        ergebnisse = [
            (muster, self.gewicht(*werte))
            for muster, werte in self.bindungsmuster.items()
        ]

        ergebnisse.sort(key=lambda eintrag: eintrag[1], reverse=True)

        return ergebnisse


    def ring_akzeptanz(self, ringgroesse):
        """Wahrscheinlichkeit (0..1), einen Ring dieser Groesse beim
        Bauen zu akzeptieren - logistische Funktion der Abweichung der
        Bayes-geglaetteten Durchschnittspunktzahl dieser Ringgroesse
        vom Gesamtmittel: im Mittel 0.5, T Punkte darueber ~0.73, T
        Punkte darunter ~0.27. Ohne Beobachtungen fuer diese Groesse
        neutral (0.5), damit unbekannte Ringgroessen nicht von
        vornherein blockiert werden - die KI muss sie erst einmal
        ausprobieren, um daraus zu lernen. Die Grenzen 0.05/0.95 halten
        die Erkundung auch bei sehr schlechten/guten Groessen offen."""

        werte = self.ringgroessen_nutzung.get(ringgroesse)

        if not werte:
            return 0.5

        abweichung = self._glaetten(*werte) - self._globaler_mittelwert()

        akzeptanz = 1.0 / (1.0 + math.exp(-max(-50.0, min(50.0, abweichung / self.temperatur))))

        return max(0.05, min(0.95, akzeptanz))


    def _lernkurve_abschliessen(self):
        """Schliesst das aktuelle Intervall ab: haengt (Beobachtungs-
        stand, Durchschnittspunktzahl DIESES Intervalls) an die
        Lernkurve an und setzt die Akkumulatoren zurueck. Getrennt von
        beobachte() gehalten, damit laden() denselben Mechanismus beim
        Wiederherstellen eines unvollstaendigen Intervalls nutzen kann
        (siehe laden())."""

        if self._intervall_anzahl == 0:
            return

        durchschnitt = self._intervall_summe / self._intervall_anzahl

        self.lernkurve.append((self.beobachtungen, durchschnitt))

        self._intervall_summe = 0.0
        self._intervall_anzahl = 0


    def lernkurve_daten(self):
        """Liste von (beobachtungen_stand, durchschnitt_punkte) Punkten,
        aeltester zuerst - Rohdaten fuer eine GUI-Grafik (z.B. Punkte
        gegen Beobachtungen auftragen)."""

        return list(self.lernkurve)


    def lerntrend(self, letzte_n=5):
        """Vergleicht den Durchschnitt der letzten 'letzte_n' Lernkurven-
        Punkte mit dem Durchschnitt der 'letzte_n' davor - ein simples,
        robusteres Signal als nur den allerersten mit dem allerletzten
        Punkt zu vergleichen (der jeweils auf einem einzelnen Intervall
        beruht und damit verrauscht sein kann). Liefert (differenz,
        text) - differenz ist None, wenn noch nicht genug Punkte fuer
        einen Vergleich vorliegen."""

        punkte = self.lernkurve

        if len(punkte) < 2:
            return None, t("ged.trend_zu_wenig")

        n = min(letzte_n, len(punkte) // 2) or 1

        fruehe = [p for _stand, p in list(punkte)[:n]]
        spaete = [p for _stand, p in list(punkte)[-n:]]

        fruehe_avg = sum(fruehe) / len(fruehe)
        spaete_avg = sum(spaete) / len(spaete)

        differenz = spaete_avg - fruehe_avg

        if differenz > 1:
            text = t("ged.trend_besser", d=f"{differenz:.1f}")
        elif differenz < -1:
            text = t("ged.trend_schlechter", d=f"{differenz:.1f}")
        else:
            text = t("ged.trend_stabil")

        return differenz, text


    def lauf_beginnen(self, einstellungen):
        """Traegt einen neuen Eintrag in die Lauf-Historie ein - von
        ChemLabor.starten() beim Start jedes Laufs aufgerufen.
        'einstellungen' ist ein flaches dict der zu dem Zeitpunkt
        geltenden Generator-/Filter-Einstellungen (siehe
        ChemLabor.starten())."""

        self.lauf_historie.append({
            "zeitpunkt": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "beobachtungen_ab": self.beobachtungen,
            "einstellungen": dict(einstellungen),
        })


    def lauf_historie_text(self):
        """Kurze, lesbare Zusammenfassung der Lauf-Historie fuer
        zusammenfassung()/die GUI: wie viele Laeufe stecken in diesem
        Gedaechtnis, und welche Werte wurden fuer die einzelnen
        Einstellungen jeweils benutzt (als Menge, nicht als Liste -
        zehnmal derselbe Wert soll nicht zehnmal auftauchen)."""

        if not self.lauf_historie:
            return t("ged.kein_lauf")

        werte_je_einstellung = defaultdict(set)

        for eintrag in self.lauf_historie:
            for schluessel, wert in eintrag["einstellungen"].items():
                werte_je_einstellung[schluessel].add(wert)

        teile = []

        for schluessel, werte in werte_je_einstellung.items():

            if len(werte) == 1:
                teile.append(f"{schluessel}={next(iter(werte))}")
            else:
                teile.append(f"{schluessel}∈{sorted(werte, key=str)}")

        return t(
            "ged.laeufe",
            n=len(self.lauf_historie),
            details=", ".join(teile)
        )


    def _signaturen_begrenzen(self):
        """Haelt molekuel_signaturen unter max_signaturen, indem beim
        Ueberschreiten nur die haeufigsten Signaturen behalten werden.
        Selten/einmalig aufgetauchte Molekuele fallen dabei aus der
        Wiederholungserkennung - das ist fuer den Zweck (Diversitaets-
        Malus fuer WIEDERHOLTE Molekuele) unproblematisch: ein
        Molekuel, das laengst aus dem Zaehler gefallen ist, wurde
        ohnehin selten genug beobachtet, dass es bei einer erneuten
        Erzeugung keine oder kaum Strafe verdient haette."""

        if len(self.molekuel_signaturen) <= self.max_signaturen:
            return

        self.molekuel_signaturen = Counter(
            dict(self.molekuel_signaturen.most_common(self.max_signaturen // 2))
        )


    def wiederholungsstrafe(self, molekuel):
        """Exponentiell fallender Faktor (0 < faktor <= 1) fuer bereits
        wiederholt beobachtete Molekuele - basis**anzahl, wobei anzahl
        die Zahl der bisherigen Beobachtungen mit exakt dieser Struktur
        ist. Ein Molekuel, das schon 3x aufgetaucht ist, bekommt bei
        basis=0.5 nur noch 1/8 seiner Punktzahl fuer die Auswahl. Ein
        komplett neues Molekuel bekommt Faktor 1 (keine Strafe).

        WICHTIG: das aendert NICHT die tatsaechliche Kritiker-
        Bewertung und NICHT die Lernstatistik (bindungsmuster,
        valenz_nutzung etc.) - ein Muster bleibt chemisch genauso gut
        oder schlecht, unabhaengig davon, wie oft es schon vorkam. Die
        Strafe wirkt ausschliesslich auf die Best-of-N-Auswahl in
        KIGenerator.erschaffe_molekuel(), um zu verhindern, dass die KI
        sich auf eine einmal gefundene 'sichere Bank' (z.B. ein
        winziges, immer stabiles Molekuel) einschiesst, statt neue
        Strukturen zu erkunden."""

        anzahl = self.molekuel_signaturen.get(molekuel_signatur(molekuel), 0)

        return self.wiederholungsstrafe_basis ** anzahl


    def stabilitaets_schwelle_adaptiv(self, perzentil=None):
        """Schaetzt den Punktwert am gegebenen Perzentil aus dem groben
        10er-Bucket-Histogramm - deutlich billiger als alle Einzelwerte
        zu sortieren, und ersetzt eine geratene feste Zahl (vorher: 70)
        durch einen Wert, der sich an der tatsaechlich beobachteten
        Verteilung orientiert. Default-Perzentil ist ziel_perzentil
        (0.75 = oberes Quartil)."""

        if perzentil is None:
            perzentil = self.ziel_perzentil

        gesamt = sum(self._histogramm)

        if gesamt == 0:
            return 70  # Fallback, solange nichts beobachtet wurde

        ziel = gesamt * (1 - perzentil)
        kumuliert = 0

        for bucket_index in range(10, -1, -1):
            kumuliert += self._histogramm[bucket_index]
            if kumuliert >= ziel:
                return bucket_index * 10

        return 0


    def zusammenfassung(self):
        """Kurzer Textbericht darueber, was die KI gelernt hat."""

        zeilen = [t("ged.beobachtet", n=self.beobachtungen)]

        if self.beobachtungen:

            globaler_mittelwert = self._punkte_summe_global / self.beobachtungen

            zeilen.append(
                t("ged.durchschnitt", x=f"{globaler_mittelwert:.1f}")
            )

            schwelle = self.stabilitaets_schwelle_adaptiv()

            zeilen.append(
                t("ged.schwelle", p=int(self.ziel_perzentil * 100), s=schwelle)
            )

            zeilen.append(
                t("ged.temperatur", wert=f"{self.temperatur:.2f}", minimum=f"{self.temperatur_min:.2f}")
            )

        if self.elemente_haeufigkeit:

            # Anzeige in Punkten (geglaetteter Schnitt), NICHT das
            # Softmax-Gewicht - sonst stuenden dort Gewichte statt
            # Punktzahlen.
            top_elemente = ", ".join(
                f"{el} (Ø{self._glaetten(*self._element_qualitaet[el]):.0f})"
                for el, _gewicht in self.bevorzugte_elemente()[:5]
            )

            zeilen.append(t("ged.beste_elemente", liste=top_elemente))

        if self.bindungsmuster:

            top_bindungen = ", ".join(
                t("ged.bindung_item", e1=e1, e2=e2, o=o,
                  x=f"{self._glaetten(*self.bindungsmuster[(e1, e2, o)]):.0f}")
                for (e1, e2, o), _gewicht in self.bevorzugte_bindungen()[:5]
            )

            zeilen.append(t("ged.beste_bindungen", liste=top_bindungen))

        if self.ringgroessen_nutzung:

            ring_zeilen = sorted(
                self.ringgroessen_nutzung.items(),
                key=lambda eintrag: self._glaetten(*eintrag[1]),
                reverse=True
            )

            top_ringe = ", ".join(
                t("ged.ring_item", g=groesse, x=f"{self._glaetten(*werte):.0f}")
                for groesse, werte in ring_zeilen[:5]
            )

            zeilen.append(t("ged.beste_ringe", liste=top_ringe))

        if self.lernkurve:

            _differenz, trend_text = self.lerntrend()

            zeilen.append(
                t("ged.lernkurve", n=len(self.lernkurve), i=self.lernkurve_intervall, trend=trend_text)
            )

        zeilen.append(t("ged.lauf_historie", text=self.lauf_historie_text()))

        return "\n".join(zeilen)


    # ------------------------------------------------------
    # Speichern / Laden
    # ------------------------------------------------------
    #
    # KIGedaechtnis lebt bisher nur im Arbeitsspeicher einer Sitzung -
    # schliesst man die GUI, ist das gesamte gelernte Wissen
    # (Bindungsmuster, Ringgroessen, Elementgewichte, ...) weg, und der
    # naechste Lauf faengt bei null an. Fuer ein System, dessen Wert
    # ueber viele Beobachtungen waechst, verschenkt das bei jedem
    # Neustart die gesamte Trainingszeit - daher jetzt persistierbar.
    # ------------------------------------------------------

    def speichern(self, pfad):
        """Schreibt den gesamten gelernten Zustand als JSON nach
        'pfad'. Tupel-Schluessel (Bindungsmuster, Signaturen) werden
        dafuer in reine Listenstrukturen umgewandelt, da JSON weder
        Tupel noch zusammengesetzte dict-Schluessel kennt - laden()
        baut sie exakt zurueck."""

        daten = {
            "version": 1,

            "einstellungen": {
                "glaettung": self.glaettung,
                "ziel_perzentil": self.ziel_perzentil,
                "wiederholungsstrafe_basis": self.wiederholungsstrafe_basis,
                "max_signaturen": self.max_signaturen,
                "lernkurve_intervall": self.lernkurve_intervall,
                "max_lernkurve_punkte": self.lernkurve.maxlen,
                "temperatur_min": self.temperatur_min,
                "temperatur_zerfall": self.temperatur_zerfall,
            },

            # Aktuelle (bereits abgekuehlte) Temperatur - ohne sie
            # startet ein geladener Lauf wieder bei voller Erkundung.
            "temperatur": self.temperatur,

            "beobachtungen": self.beobachtungen,
            "punkte_summe_global": self._punkte_summe_global,
            "histogramm": self._histogramm,

            # Abgeschlossene Intervalle plus der Stand des noch offenen
            # Intervalls (summe/anzahl) - so geht bei einem Save mitten
            # in einem Intervall keine Beobachtung verloren, und laden()
            # kann exakt dort weitermachen.
            "lernkurve": list(self.lernkurve),
            "intervall_summe": self._intervall_summe,
            "intervall_anzahl": self._intervall_anzahl,

            "lauf_historie": list(self.lauf_historie),

            "bindungsmuster": [
                [list(muster), werte] for muster, werte in self.bindungsmuster.items()
            ],

            # dict-Schluessel werden von JSON beim Schreiben automatisch
            # zu Strings - beim Laden muessen used_valence/Domaenenzahl
            # deshalb explizit wieder zu int gemacht werden.
            "valenz_nutzung": {
                element: {str(wert): werte for wert, werte in verteilung.items()}
                for element, verteilung in self.valenz_nutzung.items()
            },

            "domaenen_nutzung": {
                element: {str(wert): werte for wert, werte in verteilung.items()}
                for element, verteilung in self.domänen_nutzung.items()
            },

            "ringgroessen_nutzung": {
                str(groesse): werte for groesse, werte in self.ringgroessen_nutzung.items()
            },

            "elemente_haeufigkeit": dict(self.elemente_haeufigkeit),
            "element_qualitaet": dict(self._element_qualitaet),

            "molekuel_signaturen": [
                [_signatur_zu_json(signatur), anzahl]
                for signatur, anzahl in self.molekuel_signaturen.items()
            ],
        }

        with open(pfad, "w", encoding="utf-8") as f:
            json.dump(daten, f, ensure_ascii=False)


    @classmethod
    def laden(cls, pfad):
        """Liest eine mit speichern() geschriebene Datei und liefert
        ein neu aufgebautes KIGedaechtnis mit demselben Zustand.
        Aeltere Dateien ohne Temperatur-Felder bleiben ladbar: sie
        starten mit den Standardwerten (volle Anfangstemperatur)."""

        with open(pfad, "r", encoding="utf-8") as f:
            daten = json.load(f)

        einstellungen = daten.get("einstellungen", {})

        gedaechtnis = cls(
            glaettung=einstellungen.get("glaettung", 5),
            ziel_perzentil=einstellungen.get("ziel_perzentil", 0.75),
            wiederholungsstrafe_basis=einstellungen.get("wiederholungsstrafe_basis", 0.5),
            max_signaturen=einstellungen.get("max_signaturen", 20000),
            lernkurve_intervall=einstellungen.get("lernkurve_intervall", 1000),
            max_lernkurve_punkte=einstellungen.get("max_lernkurve_punkte", 2000),
            temperatur_min=einstellungen.get("temperatur_min", 3.0),
            temperatur_zerfall=einstellungen.get("temperatur_zerfall", 0.9999),
        )

        gedaechtnis.temperatur = daten.get("temperatur", gedaechtnis.temperatur)

        gedaechtnis.beobachtungen = daten.get("beobachtungen", 0)
        gedaechtnis._punkte_summe_global = daten.get("punkte_summe_global", 0.0)
        gedaechtnis._histogramm = daten.get("histogramm", [0] * 11)

        for stand, durchschnitt in daten.get("lernkurve", []):
            gedaechtnis.lernkurve.append((stand, durchschnitt))

        gedaechtnis._intervall_summe = daten.get("intervall_summe", 0.0)
        gedaechtnis._intervall_anzahl = daten.get("intervall_anzahl", 0)

        for eintrag in daten.get("lauf_historie", []):
            gedaechtnis.lauf_historie.append(eintrag)

        for muster, werte in daten.get("bindungsmuster", []):
            gedaechtnis.bindungsmuster[tuple(muster)] = werte

        for element, verteilung in daten.get("valenz_nutzung", {}).items():
            for wert_str, werte in verteilung.items():
                gedaechtnis.valenz_nutzung[element][int(wert_str)] = werte

        for element, verteilung in daten.get("domaenen_nutzung", {}).items():
            for wert_str, werte in verteilung.items():
                gedaechtnis.domänen_nutzung[element][int(wert_str)] = werte

        for groesse_str, werte in daten.get("ringgroessen_nutzung", {}).items():
            gedaechtnis.ringgroessen_nutzung[int(groesse_str)] = werte

        gedaechtnis.elemente_haeufigkeit = Counter(daten.get("elemente_haeufigkeit", {}))

        for element, werte in daten.get("element_qualitaet", {}).items():
            gedaechtnis._element_qualitaet[element] = werte

        for signatur_json, anzahl in daten.get("molekuel_signaturen", []):
            gedaechtnis.molekuel_signaturen[_signatur_von_json(signatur_json)] = anzahl

        return gedaechtnis


    @classmethod
    def zusammenfuehren(cls, dateipfade):
        """Fasst mehrere SEPARAT gespeicherte KIGedaechtnis-Dateien zu
        einem gemeinsamen 'Kollektivgedaechtnis' zusammen - additiv
        ueber alle Zaehler/Summen (Bindungsmuster, Valenz-, Domaenen-,
        Ringgroessen-, Elementstatistik, Signaturzaehler, Histogramm),
        mit vereinigter Lauf-Historie (siehe lauf_historie_text()).
        Gedacht fuer den Fall, dass mehrere unabhaengige Laeufe -
        moeglicherweise mit unterschiedlichen Einstellungen, z.B.
        verschiedenen Max-Atome-Werten - NICHT ueber
        "Gedaechtnis laden -> weitertrainieren" verkettet wurden,
        sondern als eigenstaendige Dateien vorliegen.

        Die Lernkurve wird bewusst NICHT numerisch verschmolzen: die
        Zeitachsen der Quelldateien sind nicht vergleichbar
        ("Beobachtung 500" in Datei A und Datei B sind zwei voneinander
        unabhaengige Zeitpunkte, ein gemeinsamer Verlauf daraus waere
        irrefuehrend). Das Ergebnis startet mit einer leeren
        Lernkurve; die Lauf-Historie bleibt dagegen vollstaendig
        erhalten, denn Zeitpunkt+Einstellungen jedes Quell-Laufs sind
        weiterhin unabhaengig voneinander gueltig.

        Konfigurationswerte (glaettung, ziel_perzentil,
        wiederholungsstrafe_basis, max_signaturen,
        lernkurve_intervall, temperatur_min, temperatur_zerfall) und
        die aktuelle Temperatur werden von der ERSTEN Datei
        uebernommen - bei abweichenden Werten in den anderen Quellen
        aendert das nichts an deren bereits gelernten Zahlen, wirkt
        sich aber ab jetzt auf neu hinzukommende Beobachtungen aus."""

        if not dateipfade:
            raise ValueError(t("ged.keine_datei"))

        quellen = [cls.laden(pfad) for pfad in dateipfade]

        basis = quellen[0]

        ergebnis = cls(
            glaettung=basis.glaettung,
            ziel_perzentil=basis.ziel_perzentil,
            wiederholungsstrafe_basis=basis.wiederholungsstrafe_basis,
            max_signaturen=basis.max_signaturen,
            lernkurve_intervall=basis.lernkurve_intervall,
            max_lernkurve_punkte=basis.lernkurve.maxlen,
            temperatur=basis.temperatur,
            temperatur_min=basis.temperatur_min,
            temperatur_zerfall=basis.temperatur_zerfall,
        )

        for quelle in quellen:

            ergebnis.beobachtungen += quelle.beobachtungen
            ergebnis._punkte_summe_global += quelle._punkte_summe_global

            for i, anzahl in enumerate(quelle._histogramm):
                ergebnis._histogramm[i] += anzahl

            for element, anzahl in quelle.elemente_haeufigkeit.items():
                ergebnis.elemente_haeufigkeit[element] += anzahl

            for element, (summe, anzahl) in quelle._element_qualitaet.items():
                eq = ergebnis._element_qualitaet[element]
                eq[0] += summe
                eq[1] += anzahl

            for element, werte in quelle.valenz_nutzung.items():
                for valenz, (summe, anzahl) in werte.items():
                    ev = ergebnis.valenz_nutzung[element][valenz]
                    ev[0] += summe
                    ev[1] += anzahl

            for element, werte in quelle.domänen_nutzung.items():
                for domänen, (summe, anzahl) in werte.items():
                    ed = ergebnis.domänen_nutzung[element][domänen]
                    ed[0] += summe
                    ed[1] += anzahl

            for muster, (summe, anzahl) in quelle.bindungsmuster.items():
                eb = ergebnis.bindungsmuster[muster]
                eb[0] += summe
                eb[1] += anzahl

            for groesse, (summe, anzahl) in quelle.ringgroessen_nutzung.items():
                er = ergebnis.ringgroessen_nutzung[groesse]
                er[0] += summe
                er[1] += anzahl

            for signatur, anzahl in quelle.molekuel_signaturen.items():
                ergebnis.molekuel_signaturen[signatur] += anzahl

            ergebnis.lauf_historie.extend(quelle.lauf_historie)

        ergebnis._signaturen_begrenzen()

        return ergebnis


# ==========================================================
# KIGenerator
# ==========================================================
#
# Baut neue Molekuele direkt (ohne den Zufallsgenerator aus 'idk'),
# auf Basis der in KIGedaechtnis gelernten Muster.
#
# 'elemente' ist das, was lade_elemente() liefert: eine einfache
# Liste von Symbolen (z.B. ["H", "C", "N", ...]). Valenzen und
# erlaubte Bindungsordnungen pro Element kommen direkt aus idk.py
# (VALENZEN, BINDUNGSARTEN) - dieselbe Quelle, die auch
# generiere_molekuel() benutzt, also garantiert konsistent.
# ==========================================================


class KIGenerator:

    def __init__(
        self,
        elemente,
        gedaechtnis,
        min_atome=2,
        max_atome=MOLEKUEL_MAX_ATOME,
        kritiker=None,
        kandidaten_pro_molekuel=5
    ):

        self.elemente = elemente
        self.gedaechtnis = gedaechtnis
        self.min_atome = min_atome
        self.max_atome = max_atome
        self.kritiker = kritiker
        self.kandidaten_pro_molekuel = kandidaten_pro_molekuel


    # ------------------------------------------------------
    # Elementinformationen
    # ------------------------------------------------------

    def _elementpool(self):
        """'elemente' ist bereits die fertige Symbolliste. Reiner
        Fallback fuer den Fall, dass noch nichts beobachtet wurde -
        dann sind alle verfuegbaren Elemente gleich gewichtet."""
        return list(self.elemente)


    def _elementpool_gewichtet(self):
        """Liefert (elemente_liste, gewichte_liste) - getrennte Listen,
        wie random.choices() sie erwartet. Nutzt gelernte Qualitaets-
        gewichte, falls vorhanden, sonst Gleichverteilung ueber alle
        bekannten Elemente."""

        gelernt = self.gedaechtnis.bevorzugte_elemente()

        if gelernt:
            return [el for el, _ in gelernt], [g for _, g in gelernt]

        pool = self._elementpool()
        return pool, [1.0] * len(pool)


    def _valenz(self, element):
        return VALENZEN.get(element, 1)


    def _erlaubte_ordnungen(self, element1, element2):
        """Schnittmenge der fuer beide Elemente erlaubten
        Bindungsordnungen (z.B. Si darf laut BINDUNGSARTEN nie eine
        Dreifachbindung eingehen, C nie eine Vierfachbindung - das
        geht ueber reine Valenzkapazitaet hinaus)."""

        a = set(BINDUNGSARTEN.get(element1, [1]))
        b = set(BINDUNGSARTEN.get(element2, [1]))

        return (a & b) or {1}


    # ------------------------------------------------------
    # Molekuel erschaffen
    # ------------------------------------------------------

    def erschaffe_molekuel(self):
        """Baut Molekuele und waehlt das beste aus mehreren Kandidaten
        (Best-of-N), FALLS ein Kritiker uebergeben wurde - die KI
        optimiert damit aktiv auf hohe Scores, statt nur passiv aus
        gelernten Haeufigkeiten heraus einmalig zu bauen. Ohne Kritiker
        (kritiker=None) faellt das auf einen einzelnen Kandidaten
        zurueck - z.B. praktisch fuer isolierte Tests."""

        if self.kritiker is None or self.kandidaten_pro_molekuel <= 1:
            return self._einzelnes_molekuel()

        bestes = None
        beste_bewertung = -1

        for _ in range(self.kandidaten_pro_molekuel):

            kandidat = self._einzelnes_molekuel()
            punkte = self.kritiker.analysiere(kandidat)["bewertung"]["punkte"]

            # Diversitaets-Malus: die tatsaechliche Kritiker-Punktzahl
            # bleibt unveraendert (die wird spaeter korrekt an Labor/
            # KIGedaechtnis gemeldet) - aber FUER DIE AUSWAHL zaehlt ein
            # bereits mehrfach beobachtetes Molekuel exponentiell
            # weniger, damit Best-of-N nicht immer wieder beim selben
            # bekannten Kandidaten landet.
            bewertung = punkte * self.gedaechtnis.wiederholungsstrafe(kandidat)

            if bewertung > beste_bewertung:
                beste_bewertung = bewertung
                bestes = kandidat

        return bestes


    def _einzelnes_molekuel(self):

        elemente_liste, gewichte_liste = self._elementpool_gewichtet()

        anzahl = random.randint(self.min_atome, self.max_atome)

        atome = []

        for i in range(anzahl):

            # FIX: vorher random.choice(pool) - eine GLEICHVERTEILTE
            # Wahl aus einer nach Qualitaet sortierten Liste. Die
            # Sortierung hatte dadurch keinerlei Einfluss auf die
            # tatsaechliche Auswahl. random.choices(..., weights=...)
            # gewichtet jetzt tatsaechlich nach der gelernten Qualitaet.
            element = random.choices(elemente_liste, weights=gewichte_liste, k=1)[0]

            valenz = self._valenz(element)
            ziel = min(self.gedaechtnis.bevorzugte_valenz(element, valenz), valenz)
            ziel_domänen = min(self.gedaechtnis.bevorzugte_domänen(element, valenz), valenz)

            atome.append({
                "id": i,
                "element": element,
                "valence": valenz,
                "used_valence": 0,
                "_ziel": max(ziel, 1),            # internes Bauziel (Valenz), wird vor Export entfernt
                "_ziel_domänen": max(ziel_domänen, 1)  # internes Bauziel (Form), wird vor Export entfernt
            })

        bindungen = self._verbinde(atome)

        for atom in atome:
            del atom["_ziel"]
            del atom["_ziel_domänen"]

        return {"atoms": atome, "bonds": bindungen}


    def _frei(self, atom):
        return atom["_ziel"] - atom["used_valence"]


    def _domänen(self, atom_id, bindungen):
        """Anzahl DISTINKTER Bindungspartner eines Atoms (VSEPR-Domaenen),
        nicht die Summe der Bindungsordnungen."""

        nachbarn = set()

        for b in bindungen:
            if b["atom1"] == atom_id:
                nachbarn.add(b["atom2"])
            elif b["atom2"] == atom_id:
                nachbarn.add(b["atom1"])

        return len(nachbarn)


    def _pfadlaenge(self, start, ziel, graph):
        """BFS-Distanz (Anzahl Kanten) zwischen zwei Atomen im
        uebergebenen Adjazenzgraphen, oder None falls (noch) nicht
        verbunden. 'graph' ist ein dict id -> set(nachbar_ids), das der
        Aufrufer pflegt (siehe _zusatzbindungen) - so muss hier nichts
        aus der Bindungsliste neu aufgebaut werden."""

        if start == ziel:
            return 0

        besucht = {start}
        queue = deque([(start, 0)])

        while queue:

            knoten, tiefe = queue.popleft()

            for nachbar in graph[knoten]:

                if nachbar == ziel:
                    return tiefe + 1

                if nachbar not in besucht:
                    besucht.add(nachbar)
                    queue.append((nachbar, tiefe + 1))

        return None


    def _verbinde(self, atome):

        bindungen = []
        gelernte_muster = self.gedaechtnis.bevorzugte_bindungen()

        # 1. Zusammenhang sicherstellen: jedes neue Atom mit einer
        #    EINFACHbindung an ein bereits verbundenes Atom haengen.
        #    Bewusst keine _passende_ordnung() hier: eine gleich zu
        #    Beginn vergebene Doppel-/Dreifachbindung kann die komplette
        #    freie Valenz eines Atoms auf einen Schlag aufbrauchen und
        #    andere Atome von vornherein aussperren (z.B. O=O verbraucht
        #    sofort beide freien Valenzen von O - ein drittes Atom findet
        #    dann nirgends mehr Platz). Ordnung 1 ist fuer jedes Element
        #    in BINDUNGSARTEN immer erlaubt, also nie ein Problem hier.
        #    Hoehere Ordnungen kommen erst in _zusatzbindungen(), wenn
        #    das Grundgeruest schon zusammenhaengt.
        verbunden = [atome[0]]
        unverbunden = atome[1:]
        uebrig_ohne_fortschritt = 0

        while unverbunden and uebrig_ohne_fortschritt < len(unverbunden) + 1:

            neu = unverbunden.pop(0)
            ziel = self._passendes_atom(neu, verbunden, gelernte_muster)

            if ziel is None:
                unverbunden.append(neu)
                uebrig_ohne_fortschritt += 1
                continue

            uebrig_ohne_fortschritt = 0

            bindungen.append({"atom1": neu["id"], "atom2": ziel["id"], "order": 1})
            neu["used_valence"] += 1
            ziel["used_valence"] += 1
            verbunden.append(neu)

        # 2. Erst jetzt, wo das Grundgeruest steht, freie Valenzen nach
        #    gelernten Mustern auffuellen - inklusive Hochstufen der
        #    gerade gesetzten Einfachbindungen zu Doppel-/Dreifachbindungen
        #    und gelerntem Ringgroessen-Feedback beim Ringschluss.
        self._zusatzbindungen(atome, bindungen, gelernte_muster)

        return bindungen


    def _passendes_atom(self, neu, kandidaten, gelernte_muster):

        gewichtung = {}

        for atom in kandidaten:

            if self._frei(atom) <= 0:
                continue

            muster = tuple(sorted([neu["element"], atom["element"]]))
            gewicht = sum(
                g for (e1, e2, _o), g in gelernte_muster
                if (e1, e2) == muster
            )

            gewichtung[atom["id"]] = gewicht + 1  # +1: auch ungelernte Paare bleiben moeglich

        if not gewichtung:
            return None

        ids = list(gewichtung.keys())
        gewichte = list(gewichtung.values())
        gewaehlte_id = random.choices(ids, weights=gewichte, k=1)[0]

        return next(a for a in kandidaten if a["id"] == gewaehlte_id)


    def _passende_ordnung(self, a1, a2, gelernte_muster):

        muster = tuple(sorted([a1["element"], a2["element"]]))
        erlaubt = self._erlaubte_ordnungen(a1["element"], a2["element"])

        max_frei = max(min(self._frei(a1), self._frei(a2)), 1)
        zulaessig = {o for o in erlaubt if o <= max_frei}

        if not zulaessig:
            return 1  # kein Platz mehr fuer eine gelernte Ordnung, Notbehelf

        optionen = [
            (o, g) for (e1, e2, o), g in gelernte_muster
            if (e1, e2) == muster and o in zulaessig
        ]

        if not optionen:
            return min(zulaessig)

        ordnungen = [o for o, _g in optionen]
        gewichte = [g for _o, g in optionen]

        return random.choices(ordnungen, weights=gewichte, k=1)[0]


    def _bindung_finden(self, bindungen, id1, id2):
        """Existierende Bindung zwischen zwei Atom-IDs, falls vorhanden.
        Nur noch als Fallback fuer Aufrufer ohne fertigen Index - intern
        nutzt _zusatzbindungen() jetzt bindung_index (siehe dort)."""

        paar = {id1, id2}

        for b in bindungen:
            if {b["atom1"], b["atom2"]} == paar:
                return b

        return None


    def _zusatzbindungen(self, atome, bindungen, gelernte_muster, max_versuche=20):

        # FIX (Performance): vorher durchsuchte _bindung_finden() bei
        # jedem der bis zu max_versuche Aufrufe die komplette
        # bindungen-Liste linear (O(Bindungen) pro Versuch). Ein Index
        # macht die Existenzpruefung O(1).
        bindung_index = {
            frozenset([b["atom1"], b["atom2"]]): b for b in bindungen
        }

        # Adjazenzgraph fuer die Ringgroessen-Erkennung (_pfadlaenge) -
        # einmal aus den bisherigen Bindungen aufgebaut und danach nur
        # noch inkrementell erweitert, statt bei jedem Versuch neu aus
        # bindungen zusammengesetzt zu werden.
        graph = defaultdict(set)
        for b in bindungen:
            graph[b["atom1"]].add(b["atom2"])
            graph[b["atom2"]].add(b["atom1"])

        versuch = 0

        while versuch < max_versuche:

            versuch += 1
            offene = [a for a in atome if self._frei(a) > 0]

            if len(offene) < 2:
                break

            a1, a2 = random.sample(offene, 2)
            bestehend = bindung_index.get(frozenset([a1["id"], a2["id"]]))

            a1_will_neue_domäne = self._domänen(a1["id"], bindungen) < a1["_ziel_domänen"]
            a2_will_neue_domäne = self._domänen(a2["id"], bindungen) < a2["_ziel_domänen"]

            if bestehend is not None:

                # Keine zweite, parallele Bindung zwischen denselben zwei
                # Atomen (chemisch bedeutungslos) - stattdessen die
                # bestehende zur Doppel-/Dreifachbindung hochstufen, aber
                # NUR wenn keins der beiden Atome noch eine neue Domaene
                # (einen zusaetzlichen Bindungspartner) anstrebt. Sonst
                # wuerde das Hochstufen genau die Valenz verbrauchen, die
                # besser fuer einen neuen Nachbarn reserviert bliebe -
                # z.B. lieber vier Einfachbindungen (tetraedrisch) als
                # zwei Doppelbindungen (linear), wenn Ersteres gelernt
                # das haeufigere Muster fuer dieses Element ist.
                if a1_will_neue_domäne or a2_will_neue_domäne:
                    continue

                erlaubt = self._erlaubte_ordnungen(a1["element"], a2["element"])
                neue_ordnung = bestehend["order"] + 1

                if neue_ordnung in erlaubt and self._frei(a1) >= 1 and self._frei(a2) >= 1:
                    bestehend["order"] = neue_ordnung
                    a1["used_valence"] += 1
                    a2["used_valence"] += 1

                continue

            # NEU: Ringgroessen-Lernen. a1/a2 sind noch nicht direkt
            # verbunden - sind sie aber ueber bestehende Bindungen
            # bereits ERREICHBAR, wuerde diese neue Bindung einen Ring
            # schliessen. Pfadlaenge + 1 = Ringgroesse. Die gelernte
            # Erfahrung entscheidet stochastisch, ob der Ringschluss
            # akzeptiert wird - bekannte gute Groessen (z.B. 6) werden
            # bevorzugt, bekannte schlechte (z.B. 3) meist vermieden,
            # unbekannte Groessen bleiben halb-offen (0.5), damit die KI
            # sie ueberhaupt erst ausprobieren kann.
            pfadlaenge = self._pfadlaenge(a1["id"], a2["id"], graph)

            if pfadlaenge is not None:

                ringgroesse = pfadlaenge + 1

                if random.random() > self.gedaechtnis.ring_akzeptanz(ringgroesse):
                    continue

            ordnung = self._passende_ordnung(a1, a2, gelernte_muster)

            # Haben beide Atome ihr Domaenen-Ziel schon erreicht, keine
            # weitere hochwertige Domaene mehr aufmachen - eine zusaetzliche
            # Einfachbindung ist okay (besser als offene Valenz), aber
            # keine neue Doppel-/Dreifachbindung, die von der gelernten
            # Form wegfuehrt.
            if not a1_will_neue_domäne and not a2_will_neue_domäne:
                ordnung = 1

            if ordnung > self._frei(a1) or ordnung > self._frei(a2):
                continue

            neue_bindung = {"atom1": a1["id"], "atom2": a2["id"], "order": ordnung}

            bindungen.append(neue_bindung)
            bindung_index[frozenset([a1["id"], a2["id"]])] = neue_bindung

            graph[a1["id"]].add(a2["id"])
            graph[a2["id"]].add(a1["id"])

            a1["used_valence"] += ordnung
            a2["used_valence"] += ordnung