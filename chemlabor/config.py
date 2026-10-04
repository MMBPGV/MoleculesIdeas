# ==========================================================
# config.py
# ==========================================================
#
# Alle Einstellungen eines Laufs an EINER Stelle. Wird ChemLabor
# uebergeben:
#
#     config = LaborConfig(beobachtungen=2000, max_atome=10)
#     ChemLabor(config=config).starten()
#
# Die Standardwerte hier gelten fuer neuen Code. Die Konstanten oben in
# labor.py sind nur noch Rueckwaertskompatibilitaet fuer Aufrufer, die
# sie noch ueberschreiben (standard_config() in labor.py).

from dataclasses import dataclass, replace
from typing import Optional


AUTOMATIK_MODI = ("zufall", "konstruktion")


def _ganzzahl(name, wert, minimum):

    if isinstance(wert, bool) or not isinstance(wert, int):
        raise ValueError(f"{name} muss eine ganze Zahl sein (ist: {wert!r}).")

    if wert < minimum:
        raise ValueError(f"{name} muss mindestens {minimum} sein (ist: {wert}).")


@dataclass
class LaborConfig:

    # Beobachtungsphase: so viele Zufallsmolekuele lernt die KI zuerst
    beobachtungen: int = 100000

    # Automatikphase: so viele Molekuele erzeugt die KI danach selbst
    automatik_experimente: int = 100

    # Best-of-N: so viele Kandidaten pro Molekuel (KIGenerator)
    kandidaten_pro_molekuel: int = 5

    # Obergrenze der Atomzahl je Molekuel (Zufall, KIGenerator, Konstruktor)
    max_atome: int = 15

    # Hoechstens so viele schlechte je gutes Molekuel werden voll gelernt;
    # None = Filter aus
    ziel_verhaeltnis_neg_zu_pos: Optional[int] = 3
    verhaeltnis_sockel: int = 20

    # Fehlversuche pro Phase = Anzahl * versuche_faktor
    versuche_faktor: int = 5

    # Ab dieser Punktzahl zeigt die Konsolenausgabe ein Automatik-Molekuel
    stabil_anzeige_schwelle: int = 90

    # Logging
    fehler_log_limit: int = 3
    automatik_endlos_log_intervall: int = 200

    # "zufall" = KIGenerator (Best-of-N), "konstruktion" = KIKonstruktor
    automatik_modus: str = "zufall"

    # True: Automatikphase laeuft, bis ChemLabor.stoppen() aufgerufen wird
    endlos: bool = False

    # Autosave des Gedaechtnisses: Pfad (None = aus) und Intervall
    autosave_pfad: Optional[str] = None
    autosave_intervall: int = 5000

    @property
    def beobachtung_log_intervall(self):
        """~20 gleichmaessig verteilte Meldungen ueber die Beobachtungsphase."""

        return max(1, self.beobachtungen // 20)

    def validiert(self):
        """Prueft alle Werte und gibt self zurueck (zum Verketten).
        Wirft ValueError mit verstaendlicher Meldung."""

        _ganzzahl("Beobachtungen", self.beobachtungen, 1)
        _ganzzahl("Automatik-Experimente", self.automatik_experimente, 1)
        _ganzzahl("Kandidaten pro Molekuel", self.kandidaten_pro_molekuel, 1)
        _ganzzahl("Max. Atome", self.max_atome, 2)
        _ganzzahl("Verhaeltnis-Sockel", self.verhaeltnis_sockel, 0)
        _ganzzahl("Versuchsfaktor", self.versuche_faktor, 1)
        _ganzzahl("Anzeigeschwelle", self.stabil_anzeige_schwelle, 0)
        _ganzzahl("Fehler-Log-Limit", self.fehler_log_limit, 0)
        _ganzzahl("Endlos-Log-Intervall", self.automatik_endlos_log_intervall, 1)
        _ganzzahl("Autosave-Intervall", self.autosave_intervall, 1)

        if self.ziel_verhaeltnis_neg_zu_pos is not None:
            _ganzzahl("Verhaeltnis Neg:Pos", self.ziel_verhaeltnis_neg_zu_pos, 0)

        if self.automatik_modus not in AUTOMATIK_MODI:
            raise ValueError(
                f"Automatik-Modus muss einer von {AUTOMATIK_MODI} sein "
                f"(ist: {self.automatik_modus!r})."
            )

        return self

    def kopie(self, **aenderungen):
        """Neue, validierte Config mit einzelnen geaenderten Feldern."""

        return replace(self, **aenderungen).validiert()
