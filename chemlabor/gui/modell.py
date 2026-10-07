# ==========================================================
# modell.py
# ==========================================================
#
# Daten und Rechenlogik der Oberflaeche - OHNE tkinter, damit sie sich
# ohne Fenster testen laesst (tests/test_gui_modell.py).
#
# Ein Ergebnis ist wie im ganzen Projekt ein Tupel
#   (nummer, phase, molekuel, analyse)
# so wie ChemLabor es liefert.
#
# Filter werden ueber INTERNE Schluessel angesprochen ("alle", "top10",
# "stabil", ...), nicht ueber ihren Anzeigetext - der haengt von der
# Sprache ab (texte.py). Die alten deutschen Anzeigetexte werden von
# filtere() weiterhin akzeptiert.

from collections import deque

from chemlabor.chemie import rdkit_bruecke as rb
from chemlabor.gui.formatierung import formatiere_probleme
from chemlabor.gui.kritiker_texte import uebersetze_text, uebersetze_urteil
from chemlabor.ki.gedaechtnis import molekuel_signatur
from chemlabor.labor import formel
from chemlabor.texte import t, t_oder, zahl


# Phasen-Filter: "alle" oder ein Phasenwert, wie ChemLabor ihn liefert
PHASEN = ("alle", "Beobachtung", "Automatik")

FILTER = ("alle", "top10", "top50", "stabil", "instabil")

# Anzeige-Filter fuer die PubChem-Spalte -> erlaubte Status (None = alle)
NEUHEIT_FILTER = {
    "alle": None,
    "unbekannt": ("unbekannt",),
    "bekannt": ("bekannt",),
    "ungeprueft": ("ungeprueft",),
}

# Fruehere deutsche Anzeigetexte -> interner Schluessel
_ALTE_NAMEN = {
    "Alle": "alle",
    "Top 10": "top10",
    "Top 50": "top50",
    "Nur stabil (≥75)": "stabil",
    "Nur instabil (<50)": "instabil",
    "Nur unbekannt": "unbekannt",
    "Nur bekannt": "bekannt",
    "Nur ungeprüft": "ungeprueft",
}

# Mehr als so viele Zeilen zeigt die Tabelle bei Filter "alle" nicht
MAX_ZEILEN = 500

STABIL_AB = 75


def punkte_von(eintrag):

    return eintrag[3]["bewertung"]["punkte"]


def phase_text(phase):
    """Anzeigetext einer Phase ('Beobachtung' -> 'Observation')."""

    return t_oder(f"phase.{phase}", phase)


def status_text(status):
    """Anzeigetext eines PubChem-Status ('unbekannt' -> 'unknown')."""

    return t_oder(f"status.{status}", status)


class ErgebnisModell:

    def __init__(self, max_eintraege=5000):

        self.eintraege = deque(maxlen=max_eintraege)

        # alle Ergebnisse seit Start (die Liste oben ist ein Ringpuffer)
        self.zaehler = 0

        # ergebnis-nummer -> Pruefergebnis aus neuheit.py (status, cid, ...)
        self.neuheit = {}

        self._signaturen = {}

    def leeren(self):

        self.eintraege.clear()
        self.zaehler = 0
        self.neuheit.clear()
        self._signaturen.clear()

    def hinzufuegen(self, eintrag):

        self.eintraege.append(eintrag)
        self.zaehler += 1

    def eintrag(self, nummer):

        return next((e for e in self.eintraege if e[0] == nummer), None)

    # ------------------------------------------------------
    # Filtern und Sortieren
    # ------------------------------------------------------

    def sortierschluessel(self, spalte):

        schluessel = {
            "nr": lambda e: e[0],
            "phase": lambda e: e[1],
            "formel": lambda e: formel(e[2]),
            "atome": lambda e: len(e[2]["atoms"]),
            "punkte": punkte_von,
            "urteil": lambda e: e[3]["bewertung"]["urteil"],
            "rdkit": lambda e: bool(e[3].get("rdkit_gueltig")),
            "pubchem": lambda e: self.pubchem_text(e[0]),
        }

        return schluessel[spalte]

    def pubchem_text(self, nummer):

        pruefung = self.neuheit.get(nummer)

        if pruefung is None:
            return "–"

        return status_text(pruefung["status"])

    def filtere(
            self,
            filter_wahl="alle",
            phase="alle",
            duplikate_ausblenden=False,
            neuheit_filter="alle",
            max_zeilen=MAX_ZEILEN,
            sortierung=None
    ):
        """sortierung: None oder (spalte, absteigend)."""

        filter_wahl = _ALTE_NAMEN.get(filter_wahl, filter_wahl)
        phase = _ALTE_NAMEN.get(phase, phase)
        neuheit_filter = _ALTE_NAMEN.get(neuheit_filter, neuheit_filter)

        ergebnisse = list(self.eintraege)

        if phase != "alle":
            ergebnisse = [e for e in ergebnisse if e[1] == phase]

        if filter_wahl == "stabil":
            ergebnisse = [e for e in ergebnisse if punkte_von(e) >= STABIL_AB]
        elif filter_wahl == "instabil":
            ergebnisse = [e for e in ergebnisse if punkte_von(e) < 50]

        erlaubt = NEUHEIT_FILTER.get(neuheit_filter)

        if erlaubt is not None:
            ergebnisse = [
                e for e in ergebnisse
                if self.neuheit.get(e[0], {}).get("status") in erlaubt
            ]

        if filter_wahl == "alle" and len(ergebnisse) > max_zeilen:
            ergebnisse = ergebnisse[-max_zeilen:]

        if duplikate_ausblenden:

            neuer_cache = {}
            gesehen = set()
            eindeutig = []

            for e in ergebnisse:

                signatur = self._signaturen.get(e[0])

                if signatur is None:
                    signatur = molekuel_signatur(e[2])

                neuer_cache[e[0]] = signatur

                if signatur in gesehen:
                    continue

                gesehen.add(signatur)
                eindeutig.append(e)

            self._signaturen = neuer_cache
            ergebnisse = eindeutig

        ergebnisse.sort(key=punkte_von, reverse=True)

        if filter_wahl == "top10":
            ergebnisse = ergebnisse[:10]
        elif filter_wahl == "top50":
            ergebnisse = ergebnisse[:50]

        if sortierung:
            spalte, absteigend = sortierung
            ergebnisse.sort(key=self.sortierschluessel(spalte), reverse=absteigend)

        return ergebnisse

    # ------------------------------------------------------
    # Kennzahlen
    # ------------------------------------------------------

    def kennzahlen(self):
        """None, solange es keine Ergebnisse gibt."""

        eintraege = list(self.eintraege)

        if not eintraege:
            return None

        punkte = [punkte_von(e) for e in eintraege]

        return {
            "gesamt": self.zaehler,
            "ausgewertet": len(eintraege),
            "mittel": sum(punkte) / len(punkte),
            "minimum": min(punkte),
            "maximum": max(punkte),
            "stabil": sum(1 for p in punkte if p >= STABIL_AB),
            "beobachtung": sum(1 for e in eintraege if e[1] == "Beobachtung"),
            "automatik": sum(1 for e in eintraege if e[1] == "Automatik"),
        }


def filter_text(labor):
    """Eine Zeile zum Verhaeltnis-Filter des Labors."""

    if labor is None or labor.ziel_verhaeltnis_neg_zu_pos is None:
        return t("verh.aus")

    return t(
        "verh.stand",
        ziel=labor.ziel_verhaeltnis_neg_zu_pos,
        pos=labor._positiv_gezaehlt,
        neg=labor._negativ_gezaehlt,
        verw=labor._negativ_verworfen
    )


def config_text(config):
    """Kurzfassung der Einstellungen fuer die Kopfzeile."""

    teile = [
        t("cfg.beobachtungen", n=zahl(config.beobachtungen)),
        t("cfg.endlos") if config.endlos
        else t("cfg.automatik", n=zahl(config.automatik_experimente)),
        t("cfg.max_atome", n=config.max_atome),
        t("cfg.konstruktion") if config.automatik_modus == "konstruktion"
        else t("cfg.zufall"),
        t("cfg.filter_aus") if config.ziel_verhaeltnis_neg_zu_pos is None
        else t("cfg.filter", n=config.ziel_verhaeltnis_neg_zu_pos),
    ]

    if config.autosave_pfad:
        teile.append(t("cfg.autosave"))

    return "  ·  ".join(teile)


def baue_detailtext(eintrag, pruefung=None):
    """Text fuer die Detailansicht eines Ergebnisses. pruefung: Eintrag
    aus ErgebnisModell.neuheit oder None."""

    nummer, phase, molekuel, analyse = eintrag

    zeilen = [
        t("detail.experiment", nummer=nummer, phase=phase_text(phase)),
        t("detail.formel", formel=formel(molekuel)),
        "",
        t(
            "detail.punkte",
            punkte=analyse["bewertung"]["punkte"],
            urteil=uebersetze_urteil(analyse["bewertung"]["urteil"])
        ),
        "",
    ]

    rdkit_gueltig = analyse.get("rdkit_gueltig")

    zeilen.append(
        t(
            "detail.rdkit_gueltig",
            antwort=t("detail.ja") if rdkit_gueltig else t("detail.nein")
        )
    )

    if rdkit_gueltig:

        formel_rdkit = rb.summenformel(molekuel)
        masse = rb.molmasse(molekuel)
        signatur = rb.kanonische_signatur(molekuel)

        if formel_rdkit:
            zeilen.append(t("detail.rdkit_formel", formel=formel_rdkit))

        if masse is not None:
            zeilen.append(t("detail.molmasse", masse=f"{masse:.2f}"))

        if signatur:
            zeilen.append(t("detail.smiles", smiles=signatur))

    if pruefung is not None:

        zeilen.append("")
        zeilen.append(t("detail.pubchem", status=status_text(pruefung["status"])))

        if pruefung.get("cid"):
            zeilen.append(t("detail.cid", cid=pruefung["cid"]))

        if pruefung.get("smiles"):
            zeilen.append(t("detail.smiles_ohne_h", smiles=pruefung["smiles"]))

        if pruefung.get("inchikey"):
            zeilen.append(t("detail.inchikey", key=pruefung["inchikey"]))

        if pruefung["status"] == "unbekannt" and pruefung.get("stereo_offen"):
            zeilen.append(t("detail.stereo"))

        if pruefung.get("fehler"):
            zeilen.append(t("detail.abfrage_fehler", fehler=pruefung["fehler"]))

    zeilen.append("")
    zeilen.append(t("detail.probleme"))

    zeilen.extend(formatiere_probleme(analyse["probleme"]) or [t("detail.keine")])

    zeilen.append("")
    zeilen.append(t("detail.erklaerung"))

    for text in analyse["erklaerung"]:
        zeilen.append(" - " + uebersetze_text(text))

    return "\n".join(zeilen)