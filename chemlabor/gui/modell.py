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

from collections import deque

from chemlabor.chemie import rdkit_bruecke as rb
from chemlabor.gui.formatierung import formatiere_probleme
from chemlabor.ki.gedaechtnis import molekuel_signatur
from chemlabor.labor import formel


PHASEN = ("Alle", "Beobachtung", "Automatik")

FILTER = ("Alle", "Top 10", "Top 50", "Nur stabil (≥75)", "Nur instabil (<50)")

# Anzeige-Filter fuer die PubChem-Spalte -> erlaubte Status (None = alle)
NEUHEIT_FILTER = {
    "Alle": None,
    "Nur unbekannt": ("unbekannt",),
    "Nur bekannt": ("bekannt",),
    "Nur ungeprüft": ("ungeprueft",),
}

STATUS_TEXT = {
    "unbekannt": "unbekannt",
    "bekannt": "bekannt",
    "ungeprueft": "ungeprüft",
}

# Mehr als so viele Zeilen zeigt die Tabelle bei Filter "Alle" nicht
MAX_ZEILEN = 500

STABIL_AB = 75


def punkte_von(eintrag):

    return eintrag[3]["bewertung"]["punkte"]


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

        return STATUS_TEXT.get(pruefung["status"], pruefung["status"])

    def filtere(
            self,
            filter_wahl="Alle",
            phase="Alle",
            duplikate_ausblenden=False,
            neuheit_filter="Alle",
            max_zeilen=MAX_ZEILEN,
            sortierung=None
    ):
        """sortierung: None oder (spalte, absteigend)."""

        ergebnisse = list(self.eintraege)

        if phase != "Alle":
            ergebnisse = [e for e in ergebnisse if e[1] == phase]

        if filter_wahl == "Nur stabil (≥75)":
            ergebnisse = [e for e in ergebnisse if punkte_von(e) >= STABIL_AB]
        elif filter_wahl == "Nur instabil (<50)":
            ergebnisse = [e for e in ergebnisse if punkte_von(e) < 50]

        erlaubt = NEUHEIT_FILTER.get(neuheit_filter)

        if erlaubt is not None:
            ergebnisse = [
                e for e in ergebnisse
                if self.neuheit.get(e[0], {}).get("status") in erlaubt
            ]

        if filter_wahl == "Alle" and len(ergebnisse) > max_zeilen:
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

        if filter_wahl == "Top 10":
            ergebnisse = ergebnisse[:10]
        elif filter_wahl == "Top 50":
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
        return "aus"

    return (
        f"{labor.ziel_verhaeltnis_neg_zu_pos}:1 – "
        f"{labor._positiv_gezaehlt} positiv, {labor._negativ_gezaehlt} negativ, "
        f"{labor._negativ_verworfen} verworfen"
    )


def config_text(config):
    """Kurzfassung der Einstellungen fuer die Kopfzeile."""

    def zahl(n):
        return f"{n:,}".replace(",", ".")

    teile = [
        f"{zahl(config.beobachtungen)} Beobachtungen",
        "Automatik endlos" if config.endlos else f"{zahl(config.automatik_experimente)} Automatik",
        f"max. {config.max_atome} Atome",
        "Konstruktion" if config.automatik_modus == "konstruktion" else "Zufall (Best-of-N)",
        "Filter aus" if config.ziel_verhaeltnis_neg_zu_pos is None
        else f"Filter {config.ziel_verhaeltnis_neg_zu_pos}:1",
    ]

    if config.autosave_pfad:
        teile.append("Autosave")

    return "  ·  ".join(teile)


def baue_detailtext(eintrag, pruefung=None):
    """Text fuer die Detailansicht eines Ergebnisses. pruefung: Eintrag
    aus ErgebnisModell.neuheit oder None."""

    nummer, phase, molekuel, analyse = eintrag

    zeilen = [
        f"Experiment {nummer}  ({phase})",
        f"Formel: {formel(molekuel)}",
        "",
        f"Punkte: {analyse['bewertung']['punkte']}  –  {analyse['bewertung']['urteil']}",
        "",
    ]

    rdkit_gueltig = analyse.get("rdkit_gueltig")

    zeilen.append(f"RDKit-gültig: {'Ja' if rdkit_gueltig else 'Nein'}")

    if rdkit_gueltig:

        formel_rdkit = rb.summenformel(molekuel)
        masse = rb.molmasse(molekuel)
        signatur = rb.kanonische_signatur(molekuel)

        if formel_rdkit:
            zeilen.append(f"RDKit-Summenformel: {formel_rdkit}")

        if masse is not None:
            zeilen.append(f"Molmasse: {masse:.2f} g/mol")

        if signatur:
            zeilen.append(f"Kanonischer SMILES: {signatur}")

    if pruefung is not None:

        zeilen.append("")
        zeilen.append(f"PubChem: {STATUS_TEXT.get(pruefung['status'], pruefung['status'])}")

        if pruefung.get("cid"):
            zeilen.append(f"CID {pruefung['cid']} (Doppelklick auf die Zeile öffnet den Eintrag)")

        if pruefung.get("smiles"):
            zeilen.append(f"SMILES ohne H: {pruefung['smiles']}")

        if pruefung.get("inchikey"):
            zeilen.append(f"InChIKey: {pruefung['inchikey']}")

        if pruefung["status"] == "unbekannt" and pruefung.get("stereo_offen"):
            zeilen.append(
                "Achtung: Stereo offen – PubChem kennt evtl. nur konkrete "
                "Stereoisomere, 'unbekannt' bitte von Hand gegenprüfen."
            )

        if pruefung.get("fehler"):
            zeilen.append(f"Fehler bei der Abfrage: {pruefung['fehler']}")

    zeilen.append("")
    zeilen.append("Probleme:")

    zeilen.extend(formatiere_probleme(analyse["probleme"]) or [" - keine"])

    zeilen.append("")
    zeilen.append("Erklärung:")

    for text in analyse["erklaerung"]:
        zeilen.append(" - " + text)

    return "\n".join(zeilen)
