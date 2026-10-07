import json
import threading

from collections import deque, Counter

from chemlabor.chemie import rdkit_bruecke as rb
from chemlabor.chemie import zufallsgenerator as idk

from chemlabor.chemie.zufallsgenerator import (
    lade_elemente,
    generiere_molekuel,
    exportiere_molekuel
)

from chemlabor.chemie.kritiker import Kritiker

from chemlabor.ki.gedaechtnis import (
    KIGedaechtnis,
    KIGenerator
)

from chemlabor.ki.konstruktor import KIKonstruktor
from chemlabor.pfade import PERIODENSYSTEM
from chemlabor.config import LaborConfig
from chemlabor.texte import t



# ==========================================================
# Einstellungen
# ==========================================================


BEOBACHTUNGEN = 100000

AUTOMATIK_EXPERIMENTE = 100


# maximale Fehlversuche pro Phase

VERSUCHE_FAKTOR = 5



# Ausgabeintervall
#
# FIX: vorher "100000 wenn BEOBACHTUNGEN > 1000 sonst 1" - bei jedem
# halbwegs typischen Lauf (BEOBACHTUNGEN zwischen ~1000 und 100000, der
# GUI-Standardbereich) traf die 100000er-Schwelle also praktisch NIE,
# das Log blieb waehrend der gesamten Beobachtungsphase stumm und
# wirkte dadurch, als wuerde die Anzeige weit hinter dem tatsaechlichen
# Fortschritt herhaenken (in Wahrheit gab es einfach keine
# Zwischenmeldungen). Jetzt proportional zur tatsaechlichen Laufgroesse:
# ~20 gleichmaessig verteilte Meldungen ueber die ganze Phase, egal ob
# sie 100 oder 100000 Beobachtungen umfasst.

BEOBACHTUNG_LOG_INTERVALL = max(1, BEOBACHTUNGEN // 20)



# stabile Moleküle anzeigen

STABIL_ANZEIGE_SCHWELLE = 90



# Best-of-N

KANDIDATEN_PRO_MOLEKUEL = 5


# Verhaeltnis-Filter: hoechstens so viele NEGATIVE (Punkte unter der
# adaptiven Stabilitaetsschwelle) Molekuele je POSITIVEM werden noch
# vollstaendig gelernt/angezeigt - der Rest zaehlt nur noch in die
# globalen Statistiken (Durchschnitt, Histogramm, Lernkurve) mit, um
# deren Aussagekraft nicht zu verzerren, taucht aber weder in den
# Ergebnissen noch in der Musterstatistik (Bindungsmuster, Valenz-,
# Domaenen-, Ringgroessenlernen) auf. None = Filter aus (altes
# Verhalten, jedes Experiment zaehlt vollstaendig).
ZIEL_VERHAELTNIS_NEG_ZU_POS = 3

# Sockelbetrag: so viele negative Experimente laufen IMMER vollstaendig
# mit, bevor der obige Verhaeltnis-Filter ueberhaupt zu greifen beginnt
# - sonst wuerde ganz am Lauf-Anfang (kaum Positive, adaptive Schwelle
# noch auf ihrem Fallback-Wert) fast alles sofort verworfen werden,
# bevor ueberhaupt ein verlaesslicher Massstab existiert.
VERHAELTNIS_SOCKEL = 20


# wie viele Fehler pro Fehlertyp maximal einzeln geloggt werden,
# bevor nur noch mitgezaehlt (aber nicht mehr ausgegeben) wird

FEHLER_LOG_LIMIT = 3


# Log-Intervall fuer die Endlos-Automatikphase (dort gibt es kein
# "X/Y" wie in der gezaehlten Phase, da Y unbekannt ist)

AUTOMATIK_ENDLOS_LOG_INTERVALL = 200


# Autosave: alle wie vielen Experimenten das Gedaechtnis automatisch
# gesichert wird - nur wirksam, wenn ChemLabor.autosave_pfad gesetzt
# ist (Standard: None, also aus).

AUTOSAVE_INTERVALL = 5000






# ==========================================================
# Hilfsfunktionen
# ==========================================================


def formel(molekuel):
    """
    Erstellt Summenformel.
    Beispiel:
    H2O
    C6H12O6
    """

    zaehler = {}


    for atom in molekuel["atoms"]:

        element = atom["element"]

        zaehler[element] = (
            zaehler.get(element, 0)
            + 1
        )



    ergebnis = ""


    for element in sorted(zaehler):


        ergebnis += element


        if zaehler[element] > 1:

            ergebnis += str(
                zaehler[element]
            )


    return ergebnis


def standard_config():
    """LaborConfig aus den Modulkonstanten oben in dieser Datei (und
    idk.MOLEKUEL_MAX_ATOME). Wird benutzt, wenn ChemLabor ohne
    config=... erzeugt wird - der alte Weg, auf dem die bisherige GUI
    Einstellungen durch Ueberschreiben der Konstanten setzte."""

    return LaborConfig(
        beobachtungen=BEOBACHTUNGEN,
        automatik_experimente=AUTOMATIK_EXPERIMENTE,
        kandidaten_pro_molekuel=KANDIDATEN_PRO_MOLEKUEL,
        max_atome=idk.MOLEKUEL_MAX_ATOME,
        ziel_verhaeltnis_neg_zu_pos=ZIEL_VERHAELTNIS_NEG_ZU_POS,
        verhaeltnis_sockel=VERHAELTNIS_SOCKEL,
        versuche_faktor=VERSUCHE_FAKTOR,
        stabil_anzeige_schwelle=STABIL_ANZEIGE_SCHWELLE,
        fehler_log_limit=FEHLER_LOG_LIMIT,
        automatik_endlos_log_intervall=AUTOMATIK_ENDLOS_LOG_INTERVALL,
        autosave_intervall=AUTOSAVE_INTERVALL,
    )









# ==========================================================
# ChemLabor
# ==========================================================


class ChemLabor:



    def __init__(
            self,
            on_ergebnis=None,
            on_fortschritt=None,
            config=None
    ):



        # Alle Einstellungen eines Laufs stehen in self.config (siehe
        # config.py). Ohne Uebergabe werden sie aus den Modulkonstanten
        # oben in dieser Datei gelesen (alter Weg, den die bisherige GUI
        # noch benutzt).
        self.config = (
            config if config is not None else standard_config()
        ).validiert()

        self.elemente = lade_elemente(
            PERIODENSYSTEM
        )



        self.kritiker = Kritiker()



        self.gedaechtnis = KIGedaechtnis()



        self.ki_generator = None

        # KIKonstruktor (siehe konstruktor.py) - alternative
        # Automatikphase: baut Molekuele SCHRITT FUER SCHRITT mit
        # sofortiger RDKit-Pruefung/Backtracking, statt wie KIGenerator
        # ganze Kandidaten zu wuerfeln und erst hinterher per Best-of-N
        # zu vergleichen. self.automatik_modus waehlt zwischen beiden
        # (siehe automatikphase()).
        self.konstruktor = None
        self.automatik_modus = self.config.automatik_modus





        # Statistik

        self.experimente = 0

        self.versuche = 0



        # FIX: vorher wurden Fehler in beobachtungsphase()/
        # automatikphase() nur einzeln geloggt und dann komplett
        # vergessen - ein systematischer Bug (z.B. inkonsistente
        # Daten zwischen Generator und Kritiker) konnte sich so als
        # "ein paar rote Log-Zeilen" tarnen, ohne dass am Ende sichtbar
        # wurde, WIE OFT welcher Fehlertyp aufgetreten ist.

        self.fehler_counter = Counter()

        # Verhaeltnis-Filter (siehe ZIEL_VERHAELTNIS_NEG_ZU_POS/
        # VERHAELTNIS_SOCKEL und _verhaeltnis_pruefen()).
        self.ziel_verhaeltnis_neg_zu_pos = self.config.ziel_verhaeltnis_neg_zu_pos
        self._positiv_gezaehlt = 0
        self._negativ_gezaehlt = 0
        self._negativ_verworfen = 0



        # Signal fuer die Endlos-Automatikphase (siehe automatikphase()
        # und stoppen()) - threading.Event statt eines einfachen bool,
        # weil das Labor im Hintergrundthread laeuft und die GUI aus
        # dem Hauptthread heraus stoppen koennen muss.

        self.stop_event = threading.Event()

        # Autosave (siehe _pruefe_autosave/_melde_ergebnis) - aus per
        # Default, GUI/Konsole koennen autosave_pfad vor starten()
        # setzen, um es zu aktivieren.
        self.autosave_pfad = self.config.autosave_pfad
        self.autosave_intervall = self.config.autosave_intervall





        # ==================================================
        # Speicherbegrenzung
        #
        # vorher:
        #   Millionen Ergebnisse im RAM
        #
        # jetzt:
        #   maximal 5000
        #
        # ==================================================


        self.MAX_SPEICHER_ERGEBNISSE = 5000



        self.alle_ergebnisse = deque(
            maxlen=self.MAX_SPEICHER_ERGEBNISSE
        )






        # ==================================================
        # Callback-Puffer
        #
        # GUI bekommt nicht jedes Molekül einzeln,
        # sondern Pakete.
        #
        # ==================================================


        self.callback_puffer = []


        self.CALLBACK_BATCH = 25





        self.on_ergebnis = (

            on_ergebnis

            if on_ergebnis

            else self._konsolen_ausgabe

        )



        self.on_fortschritt = (

            on_fortschritt

            if on_fortschritt

            else lambda text:
            print(text)

        )







    # ------------------------------------------------------
    # Molekül-Erzeugung Zufall
    # ------------------------------------------------------


    def neues_molekuel_zufall(self):


        atome, bindungen = generiere_molekuel(self.elemente, max_atome=self.config.max_atome)


        return exportiere_molekuel(
            atome,
            bindungen
        )







    # ------------------------------------------------------
    # Molekül-Erzeugung KI
    # ------------------------------------------------------


    def neues_molekuel_ki(self):


        return self.ki_generator.erschaffe_molekuel()




    def neues_molekuel_konstruktion(self):
        """Wrapper um KIKonstruktor.baue_molekuel() mit demselben
        Null-Argument-Vertrag wie neues_molekuel_zufall()/
        neues_molekuel_ki() (siehe experiment()). Ein gescheiterter Bau
        wird als Exception gemeldet, statt still ein unvollstaendiges
        Molekuel zurueckzugeben - so greift dieselbe Fehlerbehandlung
        (_behandle_fehler(), einfach naechster Versuch) wie bei jedem
        anderen Erzeuger auch, ohne die Aufrufer (automatikphase()) zu
        verzweigen."""

        molekuel, erfolgreich = self.konstruktor.baue_molekuel()

        if not erfolgreich:
            raise RuntimeError(t("lab.konstruktion_gescheitert"))

        return molekuel







    # ------------------------------------------------------
    # Experiment
    # ------------------------------------------------------


    def experiment(
            self,
            erzeuge_molekuel
    ):


        molekuel = erzeuge_molekuel()



        analyse = self.kritiker.analysiere(
            molekuel
        )


        # Zusaetzliches, vom eigenen Kritiker UNABHAENGIGES
        # Validitaetssignal (siehe rdkit_bruecke.py): unser Kritiker
        # bewertet nach eigenen Heuristiken, RDKit prueft nach echten
        # chemischen Sanitize-Regeln. Beide koennen auseinanderlaufen -
        # gerade das macht diesen zweiten Blick nuetzlich.
        analyse["rdkit_gueltig"] = rb.ist_gueltig(molekuel)


        return molekuel, analyse






    # ------------------------------------------------------
    # Konsolen-Ausgabe
    # ------------------------------------------------------


    def _konsolen_ausgabe(
            self,
            nummer,
            molekuel,
            analyse,
            phase
    ):



        if (
            phase == "Automatik"
            and
            analyse["bewertung"]["punkte"]
            <
            self.config.stabil_anzeige_schwelle
        ):

            return





        print(
            f"\n[{phase}] Experiment {nummer}\n"
        )



        print("Molekül:")



        print(
            json.dumps(
                molekuel,
                indent=4,
                ensure_ascii=False
            )
        )



        print("\nKritiker:")



        print(
            json.dumps(
                analyse,
                indent=4,
                ensure_ascii=False
            )
        )



        print(
            "\n"
            +
            "-" * 50
        )




    # ------------------------------------------------------
    # Ergebnis melden
    #
    # Speicherung:
    # - Labor behält nur begrenzte Historie
    # - GUI bekommt Ergebnisse in Paketen
    #
    # ------------------------------------------------------


    def _melde_ergebnis(
            self,
            molekuel,
            analyse,
            phase
    ):
        # [patch_fixes] FIX: Der Batch-Flush stand versehentlich am Ende
        # von _pruefe_autosave() und wurde ohne Autosave nie erreicht -
        # die GUI bekam waehrend einer Phase nichts, und der
        # callback_puffer wuchs unbegrenzt (v.a. im Endlos-Modus).

        self.experimente += 1

        eintrag = (
            self.experimente,
            phase,
            molekuel,
            analyse
        )

        # begrenzter Speicher
        self.alle_ergebnisse.append(eintrag)

        self._pruefe_autosave()

        # Callback sammeln
        self.callback_puffer.append(eintrag)

        # erst ab Batch-Groesse senden
        if len(self.callback_puffer) >= self.CALLBACK_BATCH:
            self._puffer_leeren()



    def _pruefe_autosave(self):
        """Sichert das Gedaechtnis automatisch alle autosave_intervall
        Experimente, FALLS autosave_pfad gesetzt ist (Standard: None,
        also aus). Ein fehlgeschlagenes Autosave (z.B. Platte voll,
        Pfad nicht beschreibbar) bricht den Lauf nicht ab - nur eine
        Meldung, dann geht's weiter."""

        if self.autosave_pfad is None:
            return

        if self.experimente % self.autosave_intervall != 0:
            return

        try:

            self.speichere_gedaechtnis(self.autosave_pfad)

            self.on_fortschritt(
                t("lab.autosave_ok", n=self.experimente, pfad=self.autosave_pfad)
            )

        except OSError as fehler:

            self.on_fortschritt(
                t("lab.autosave_fehler", fehler=fehler)
            )



    # ------------------------------------------------------
    # Puffer leeren
    #
    # FIX: vorher wurde der callback_puffer nur bei Erreichen von
    # CALLBACK_BATCH ODER ganz am Ende von starten() geleert. Endete
    # z.B. die Beobachtungsphase mit 17 ungesendeten Eintraegen im
    # Puffer, blieben die bis zum Ende der Automatikphase liegen und
    # wurden dann zusammen mit den ersten Automatik-Ergebnissen in
    # einem Rutsch ausgeliefert - die GUI bekam Beobachtung/Automatik
    # dadurch nicht sauber in Echtzeit getrennt. Jetzt wird der Puffer
    # explizit am Ende jeder Phase geleert.
    # ------------------------------------------------------


    def _puffer_leeren(self):

        for nr, ph, mol, ana in self.callback_puffer:

            self.on_ergebnis(
                nr,
                mol,
                ana,
                ph
            )

        self.callback_puffer.clear()



    # ------------------------------------------------------
    # Fehler behandeln
    #
    # FIX: zaehlt Fehler nach Typ mit, statt sie nur zu loggen und zu
    # vergessen. Die ersten FEHLER_LOG_LIMIT Vorkommen eines Fehlertyps
    # werden weiterhin einzeln ausgegeben (fuer sofortige Sichtbarkeit),
    # danach nur noch stumm mitgezaehlt - sonst wuerde ein
    # systematischer Fehler bei 100.000 Durchlaeufen den Log fluten.
    # Die Gesamtzahl je Typ steht am Ende in der Abschlussmeldung.
    # ------------------------------------------------------


    def _behandle_fehler(self, fehler):

        typ = type(fehler).__name__

        self.fehler_counter[typ] += 1

        if self.fehler_counter[typ] <= self.config.fehler_log_limit:

            self.on_fortschritt(
                t("lab.fehler", typ=typ, fehler=fehler)
            )

        elif self.fehler_counter[typ] == self.config.fehler_log_limit + 1:

            self.on_fortschritt(
                t("lab.weitere_fehler", typ=typ)
            )






    # ------------------------------------------------------
    # Verhaeltnis-Filter
    #
    # Reduziert den "Muell" ueberschuessiger schlechter Molekuele:
    # positive Experimente (Punkte >= adaptiver Stabilitaetsschwelle)
    # laufen immer vollstaendig mit, negative nur noch bis zum
    # eingestellten Verhaeltnis - der Ueberschuss wird verworfen
    # (weder in die Ergebnisliste noch in die Musterstatistik), fliesst
    # aber weiterhin in die globalen Statistiken der KIGedaechtnis ein
    # (siehe KIGedaechtnis.beobachte(vollstaendig=)), damit Durchschnitt/
    # Histogramm/Lernkurve/adaptive Schwelle repraesentativ bleiben.
    #
    # Nebeneffekt: da verworfene Experimente weder die teure
    # Musterstatistik fuellen noch in alle_ergebnisse/GUI landen, wird
    # dadurch auch der gesamte Beobachtungs-/Automatiklauf schneller
    # und die Ergebnisliste bleibt uebersichtlicher.
    # ------------------------------------------------------

    def _verhaeltnis_pruefen(self, analyse):
        """Rueckgabe: (positiv, vollstaendig). 'positiv' sagt nur etwas
        ueber die Punktzahl relativ zur aktuellen Schwelle aus,
        'vollstaendig' ob dieses Experiment voll gelernt/gemeldet werden
        soll."""

        if self.ziel_verhaeltnis_neg_zu_pos is None:
            return True, True

        punkte = analyse["bewertung"]["punkte"]
        schwelle = self.gedaechtnis.stabilitaets_schwelle_adaptiv()

        positiv = punkte >= schwelle

        if positiv:
            self._positiv_gezaehlt += 1
            return True, True

        self._negativ_gezaehlt += 1

        erlaubt = (
            self._positiv_gezaehlt * self.ziel_verhaeltnis_neg_zu_pos
            + self.config.verhaeltnis_sockel
        )

        if self._negativ_gezaehlt <= erlaubt:
            return False, True

        self._negativ_verworfen += 1

        return False, False


    # ------------------------------------------------------
    # Phase 1:
    #
    # KI beobachtet Zufallsmoleküle
    #
    # ------------------------------------------------------


    def beobachtungsphase(self):



        self.on_fortschritt(

            t("lab.beobachtungsphase", n=self.config.beobachtungen)

        )




        beobachtet = 0


        versuche_in_phase = 0



        max_versuche_phase = (

            self.config.beobachtungen
            *
            self.config.versuche_faktor

        )







        while (


            beobachtet < self.config.beobachtungen


            and


            versuche_in_phase < max_versuche_phase


            and


            not self.stop_event.is_set()


        ):



            versuche_in_phase += 1


            self.versuche += 1






            try:



                molekuel, analyse = self.experiment(

                    self.neues_molekuel_zufall

                )



            except Exception as fehler:



                self._behandle_fehler(fehler)


                continue







            beobachtet += 1




            positiv, vollstaendig = self._verhaeltnis_pruefen(analyse)

            # KI bekommt Erfahrung

            self.gedaechtnis.beobachte(

                molekuel,

                analyse,

                vollstaendig=vollstaendig

            )





            if vollstaendig:

                self._melde_ergebnis(

                    molekuel,

                    analyse,

                    "Beobachtung"

                )







            # =================================================
            # Fortschritt begrenzen
            #
            # FIX: vorher wurde hier bei jedem Aufruf ein eigenes,
            # lokales "wenn > 1000 dann alle 5000 sonst alle 100"
            # berechnet - eine Kopie der Logik, die oben als Modul-
            # konstante BEOBACHTUNG_LOG_INTERVALL bereits existiert,
            # aber nirgends verwendet wurde (auch die GUI setzte sie
            # nur, ohne dass hier je jemand draufgeschaut hat). Jetzt
            # wird die tatsaechliche Konstante benutzt - die GUI kann
            # sie weiterhin vor dem Start ueberschreiben und hat damit
            # wirklich Wirkung.
            # =================================================



            if beobachtet % self.config.beobachtung_log_intervall == 0:



                self.on_fortschritt(

                    t("lab.beobachtung_fortschritt", i=beobachtet, n=self.config.beobachtungen,
                    p=analyse['bewertung']['punkte'])

                )









        self.on_fortschritt(

            self.gedaechtnis.zusammenfassung()

        )



        # Reste dieser Phase sofort ausliefern, statt sie bis zum
        # Ende der Automatikphase im Puffer liegen zu lassen.

        self._puffer_leeren()






        # KI Generator wird nach Beobachtung erzeugt


        self.ki_generator = KIGenerator(


            self.elemente,


            self.gedaechtnis,


            max_atome=self.config.max_atome,


            kritiker=self.kritiker,


            kandidaten_pro_molekuel=
            self.config.kandidaten_pro_molekuel


        )


        self.konstruktor = KIKonstruktor(
            self.elemente,
            self.gedaechtnis,
            max_atome=self.config.max_atome
        )




    # ------------------------------------------------------
    # Phase 2:
    #
    # KI erschafft eigene Moleküle
    #
    # ------------------------------------------------------


    def stoppen(self):
        """Signalisiert dem laufenden Labor, sich nach dem aktuellen
        Versuch zu beenden - wirkt in BEIDEN Phasen (Beobachtung und
        Automatik) und unabhaengig vom Endlos-Modus. Frueher wirkte das
        nur im Endlos-Modus; ein gezaehlter Lauf mit sehr hoher Zahl
        (z.B. 1.000.000 Beobachtungszyklen, wie in einem frueheren
        Trainingslauf tatsaechlich genutzt) liess sich davor gar nicht
        vorzeitig abbrechen."""

        self.stop_event.set()


    def speichere_gedaechtnis(self, pfad):
        """Duenner Wrapper um KIGedaechtnis.speichern() - Einstiegspunkt
        fuer sowohl Konsole als auch GUI, damit gelerntes Wissen einen
        Sitzungsende ueberlebt."""

        self.gedaechtnis.speichern(pfad)


    def lade_gedaechtnis(self, pfad):
        """Ersetzt self.gedaechtnis durch den gespeicherten Zustand aus
        'pfad'. Ein bereits existierender ki_generator (nach einer
        Beobachtungsphase) referenziert weiterhin das ALTE gedaechtnis-
        Objekt - lade_gedaechtnis() ist daher nur sinnvoll VOR
        beobachtungsphase()/automatikphase(), wenn ki_generator noch
        None ist. Ein Aufruf danach wuerde das Gedaechtnis zwar
        ersetzen, der bereits gebaute Generator wuerde davon aber
        nichts mitbekommen - deshalb hier bewusst kein stiller
        Auto-Fix, sondern ein Hinweis im Docstring."""

        self.gedaechtnis = KIGedaechtnis.laden(pfad)


    def automatikphase(self, endlos=False):
        """endlos=True: laeuft, bis stoppen() aufgerufen wird (siehe
        stop_event), statt nach AUTOMATIK_EXPERIMENTE Durchlaeufen von
        selbst zu enden. Alles andere (Fehlerbehandlung, Lernen,
        Ergebnismeldung) ist identisch zum gezaehlten Modus."""

        if endlos:

            self.on_fortschritt(t("lab.automatik_endlos"))

        else:

            self.on_fortschritt(
                t("lab.automatikphase", n=self.config.automatik_experimente)
            )


        automatisch = 0

        versuche_in_phase = 0

        max_versuche_phase = self.config.automatik_experimente * self.config.versuche_faktor  # nur im gezaehlten Modus relevant


        erzeuger = (
            self.neues_molekuel_konstruktion
            if self.automatik_modus == "konstruktion"
            else self.neues_molekuel_ki
        )


        while True:

            # FIX: stop_event wird jetzt in BEIDEN Zweigen geprueft -
            # vorher konnte ein gezaehlter Lauf (endlos=False) nur ueber
            # sein eigenes Limit enden, stoppen() hatte darauf keine
            # Wirkung.
            if self.stop_event.is_set():
                break

            if not endlos:

                if (
                    automatisch >= self.config.automatik_experimente
                    or versuche_in_phase >= max_versuche_phase
                ):
                    break


            versuche_in_phase += 1

            self.versuche += 1


            try:

                molekuel, analyse = self.experiment(
                    erzeuger
                )

            except Exception as fehler:

                self._behandle_fehler(fehler)

                continue


            automatisch += 1


            positiv, vollstaendig = self._verhaeltnis_pruefen(analyse)

            # Auch die KI lernt weiter
            self.gedaechtnis.beobachte(
                molekuel,
                analyse,
                vollstaendig=vollstaendig
            )


            if vollstaendig:
                self._melde_ergebnis(
                    molekuel,
                    analyse,
                    "Automatik"
                )


            if endlos:

                if automatisch % self.config.automatik_endlos_log_intervall == 0:

                    self.on_fortschritt(
                        t("lab.automatik_endlos_fortschritt", n=automatisch,
                        p=analyse['bewertung']['punkte'])
                    )

            else:

                # FIX: bisher gab es im NICHT-Endlos-Fall waehrend der
                # ganzen Automatikphase keine einzige Zwischenmeldung -
                # bei einem gross eingestellten AUTOMATIK_EXPERIMENTE
                # (z.B. mehrere Tausend) wirkte das Log genauso
                # "eingefroren" wie die Beobachtungsphase vor deren Fix
                # (siehe BEOBACHTUNG_LOG_INTERVALL). Dieselbe Idee: ~20
                # gleichmaessig verteilte Meldungen ueber die Phase.
                automatik_log_intervall = max(1, self.config.automatik_experimente // 20)

                if automatisch % automatik_log_intervall == 0:

                    self.on_fortschritt(
                        t("lab.automatik_fortschritt", i=automatisch, n=self.config.automatik_experimente,
                        p=analyse['bewertung']['punkte'])
                    )


        # Reste dieser Phase sofort ausliefern.

        self._puffer_leeren()







    # ------------------------------------------------------
    # Labor starten
    # ------------------------------------------------------


    def starten(self, endlos=None):
        """endlos=True: die Automatikphase laeuft, bis stoppen()
        aufgerufen wird, statt nach AUTOMATIK_EXPERIMENTE von selbst
        zu enden (siehe automatikphase())."""



        if endlos is None:
            endlos = self.config.endlos

        self.on_fortschritt(

            t("lab.titel")

        )


        # Traegt diesen Lauf in die Lauf-Historie des Gedaechtnisses
        # ein (siehe KIGedaechtnis.lauf_beginnen()) - macht sichtbar,
        # mit welchen Einstellungen die verschiedenen Abschnitte des
        # aktuell geladenen "Kollektivwissens" jeweils entstanden sind.
        self.gedaechtnis.lauf_beginnen({
            "beobachtungen_soll": self.config.beobachtungen,
            "automatik_soll": self.config.automatik_experimente,
            "endlos": endlos,
            "kandidaten_pro_molekuel": self.config.kandidaten_pro_molekuel,
            "max_atome": self.config.max_atome,
            "verhaeltnis_neg_zu_pos": self.ziel_verhaeltnis_neg_zu_pos,
            "automatik_modus": self.automatik_modus,
        })


        self.beobachtungsphase()



        self.automatikphase(endlos=endlos)





        # Etwaige Reste ausliefern (durch die Flushes am Ende jeder
        # Phase sollte der Puffer hier bereits leer sein - dieser
        # Aufruf bleibt als Sicherheitsnetz bestehen).


        self._puffer_leeren()



        # Fehler-Zusammenfassung, statt die Information stillschweigend
        # verloren gehen zu lassen.

        if self.fehler_counter:

            fehler_text = ", ".join(
                f"{typ}: {anzahl}"
                for typ, anzahl in self.fehler_counter.most_common()
            )

        else:

            fehler_text = t("lab.keine")



        self.on_fortschritt(

            t("lab.beendet", v=self.versuche, e=self.experimente, f=fehler_text)

        )

        if self.ziel_verhaeltnis_neg_zu_pos is not None:

            self.on_fortschritt(
                t("lab.verhaeltnis_ende", ziel=self.ziel_verhaeltnis_neg_zu_pos,
                pos=self._positiv_gezaehlt, neg=self._negativ_gezaehlt,
                verw=self._negativ_verworfen)
            )








# ==========================================================
# Start (Konsolenmodus)
# ==========================================================


if __name__ == "__main__":



    labor = ChemLabor()



    labor.starten()