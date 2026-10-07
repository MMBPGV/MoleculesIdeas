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
#
# Fehlermeldungen (ValueError) kommen aus chemlabor/texte.py und sind
# daher in der aktiven Sprache formuliert.

from dataclasses import dataclass, replace
from typing import Optional

from chemlabor.texte import t


AUTOMATIK_MODI = ("zufall", "konstruktion")


def _ganzzahl(feld, wert, minimum):
    """feld ist der Schluessel des Feldnamens in texte.py (ohne 'feld.')."""

    name = t(f"feld.{feld}")

    if isinstance(wert, bool) or not isinstance(wert, int):
        raise ValueError(t("config.ganzzahl", name=name, wert=repr(wert)))

    if wert < minimum:
        raise ValueError(
            t("config.minimum", name=name, minimum=minimum, wert=wert)
        )


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

        _ganzzahl("beobachtungen", self.beobachtungen, 1)
        _ganzzahl("automatik_experimente", self.automatik_experimente, 1)
        _ganzzahl("kandidaten", self.kandidaten_pro_molekuel, 1)
        _ganzzahl("max_atome", self.max_atome, 2)
        _ganzzahl("verhaeltnis_sockel", self.verhaeltnis_sockel, 0)
        _ganzzahl("versuche_faktor", self.versuche_faktor, 1)
        _ganzzahl("anzeigeschwelle", self.stabil_anzeige_schwelle, 0)
        _ganzzahl("fehler_log_limit", self.fehler_log_limit, 0)
        _ganzzahl("endlos_log", self.automatik_endlos_log_intervall, 1)
        _ganzzahl("autosave_intervall", self.autosave_intervall, 1)

        if self.ziel_verhaeltnis_neg_zu_pos is not None:
            _ganzzahl("verhaeltnis", self.ziel_verhaeltnis_neg_zu_pos, 0)

        if self.automatik_modus not in AUTOMATIK_MODI:
            raise ValueError(
                t(
                    "config.modus",
                    modi=AUTOMATIK_MODI,
                    wert=repr(self.automatik_modus)
                )
            )

        return self

    def kopie(self, **aenderungen):
        """Neue, validierte Config mit einzelnen geaenderten Feldern."""

        return replace(self, **aenderungen).validiert()