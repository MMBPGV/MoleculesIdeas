# ==========================================================
# neuheit.py
# ==========================================================
#
# Neuheits-Filter: prueft, ob ein vom Labor erzeugtes Molekuel
# bereits in PubChem bekannt ist.
#
# Ablauf pro Molekuel:
#   1. rdkit_bruecke.pruefe()  -> nur RDKit-gueltige Molekuele weiter
#   2. Wasserstoffe entfernen  -> unser Format fuehrt H als eigene
#      Atome (siehe rdkit_bruecke._zu_roh_rdkit), PubChem/InChI
#      erwarten die normale Schreibweise ohne explizite H-Atome
#   3. kanonische SMILES + InChIKey berechnen
#   4. Cache (SQLite) befragen, sonst PubChem PUG REST:
#        /compound/inchikey/<KEY>/cids/JSON
#      Treffer -> "bekannt" (mit CID), HTTP 404 -> "unbekannt"
#
# Status-Werte im Ergebnis:
#   "bekannt"     in PubChem gefunden (cid gesetzt)
#   "unbekannt"   in PubChem NICHT gefunden (siehe Stereo-Hinweis unten)
#   "ungeprueft"  offline und nicht im Cache, oder Netzwerkfehler
#   "ungueltig"   RDKit-Sanitize bestanden nicht / kein InChIKey
#
# STEREO-HINWEIS: unsere Molekuele haben keine Stereoinformation,
# der InChIKey ist also der "flache" Schluessel. Gibt es in PubChem
# nur die konkreten Stereoisomere, aber keinen flachen Eintrag, kommt
# faelschlich "unbekannt" heraus. Deshalb steht in jedem Ergebnis
# "stereo_offen" (True = das Molekuel hat unfestgelegte Stereozentren
# oder Doppelbindungs-Stereo). "unbekannt" + stereo_offen=True ist
# ein Kandidat, der noch einmal von Hand oder per
# Verbindungsebenen-Suche gegengeprueft werden sollte.
#
# PubChem-Nutzungsregeln: hoechstens 5 Anfragen/Sekunde und
# 400/Minute. min_intervall=0.25 s ergibt maximal 240/Minute.
# Jede Antwort ("bekannt"/"unbekannt") wird gecacht, Fehler nicht.
# ==========================================================


import json
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request

from datetime import datetime, timezone

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemlabor.chemie import rdkit_bruecke as rb
from chemlabor.pfade import NEUHEIT_CACHE_PFAD


PUBCHEM_BASIS_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey"
)

CACHE_STANDARD = NEUHEIT_CACHE_PFAD


class PubChemFehler(Exception):
    """Netzwerk-/Serverproblem, das auch nach Wiederholungen bleibt."""


# ==========================================================
# Struktur-Schluessel
# ==========================================================


def struktur_schluessel(molekuel):
    """Berechnet aus dem internen Format die Kennungen fuer den
    Abgleich. Rueckgabe None, wenn RDKit das Molekuel ablehnt oder
    kein InChIKey berechenbar ist (kommt bei sehr exotischen Faellen
    wie ungewoehnlichen Radikalen/Ladungen vor)."""

    gueltig, _fehler, m = rb.pruefe(molekuel)

    if not gueltig:
        return None

    try:
        m_ohne_h = Chem.RemoveHs(m)
        smiles = Chem.MolToSmiles(m_ohne_h)
        inchikey = Chem.MolToInchiKey(m_ohne_h)
    except Exception:  # RDKit-Fehlerklassen variieren je Version
        return None

    if not inchikey:
        return None

    return {
        "smiles": smiles,
        "inchikey": inchikey,
        "formel": rdMolDescriptors.CalcMolFormula(m_ohne_h),
        "molmasse": round(rdMolDescriptors.CalcExactMolWt(m_ohne_h), 4),
        "stereo_offen": len(Chem.FindPotentialStereo(m_ohne_h)) > 0,
    }


# ==========================================================
# Pruefer
# ==========================================================


class NeuheitsPruefer:

    def __init__(
            self,
            cache_pfad=CACHE_STANDARD,
            offline=False,
            min_intervall=0.25,
            timeout=15,
            max_wiederholungen=4,
            basis_url=PUBCHEM_BASIS_URL
    ):
        """offline=True: nur der Cache wird befragt, es geht keine
        Anfrage ins Netz. basis_url ist ueberschreibbar (z.B. fuer
        Tests mit einem lokalen Server)."""

        self.offline = offline
        self.min_intervall = min_intervall
        self.timeout = timeout
        self.max_wiederholungen = max_wiederholungen
        self.basis_url = basis_url

        self._letzte_anfrage = 0.0

        self.db = sqlite3.connect(cache_pfad)

        self.db.execute(
            "CREATE TABLE IF NOT EXISTS pubchem ("
            "inchikey TEXT PRIMARY KEY, "
            "status TEXT NOT NULL, "
            "cid INTEGER, "
            "geprueft TEXT NOT NULL)"
        )

        self.db.commit()


    def schliessen(self):
        self.db.close()


    # ------------------------------------------------------
    # Cache
    # ------------------------------------------------------


    def _cache_lesen(self, inchikey):

        zeile = self.db.execute(
            "SELECT status, cid FROM pubchem WHERE inchikey = ?",
            (inchikey,)
        ).fetchone()

        return zeile  # (status, cid) oder None


    def _cache_schreiben(self, inchikey, status, cid):

        self.db.execute(
            "INSERT OR REPLACE INTO pubchem VALUES (?, ?, ?, ?)",
            (
                inchikey,
                status,
                cid,
                datetime.now(timezone.utc).isoformat()
            )
        )

        self.db.commit()


    # ------------------------------------------------------
    # PubChem
    # ------------------------------------------------------


    def _warte(self):
        """Haelt den Mindestabstand zwischen zwei Anfragen ein."""

        rest = self.min_intervall - (
            time.monotonic() - self._letzte_anfrage
        )

        if rest > 0:
            time.sleep(rest)

        self._letzte_anfrage = time.monotonic()


    def _frage_pubchem(self, inchikey):
        """Rueckgabe: ("bekannt", cid) oder ("unbekannt", None).
        Wiederholt bei 429/503/504 und Verbindungsfehlern mit
        wachsender Pause, wirft danach PubChemFehler."""

        url = (
            f"{self.basis_url}/"
            f"{urllib.parse.quote(inchikey)}/cids/JSON"
        )

        anfrage = urllib.request.Request(
            url,
            headers={"User-Agent": "KIGenerator-Neuheitsfilter/0.1"}
        )

        letzter_fehler = None

        for versuch in range(self.max_wiederholungen):

            self._warte()

            try:

                with urllib.request.urlopen(
                        anfrage, timeout=self.timeout
                ) as antwort:
                    daten = json.load(antwort)

                # CID 0 bedeutet bei PubChem "kein Treffer"
                cids = [
                    c for c in
                    daten.get("IdentifierList", {}).get("CID", [])
                    if c
                ]

                if cids:
                    return "bekannt", cids[0]

                return "unbekannt", None

            except urllib.error.HTTPError as fehler:

                # HTTPError ist Unterklasse von URLError - muss
                # deshalb VOR dem URLError-Zweig stehen.

                if fehler.code == 404:
                    return "unbekannt", None

                if fehler.code in (429, 503, 504):
                    letzter_fehler = f"HTTP {fehler.code}"
                    time.sleep(2 ** versuch)
                    continue

                raise PubChemFehler(f"HTTP {fehler.code}") from fehler

            except (urllib.error.URLError, TimeoutError) as fehler:

                letzter_fehler = str(fehler)
                time.sleep(2 ** versuch)

        raise PubChemFehler(
            f"nach {self.max_wiederholungen} Versuchen aufgegeben "
            f"({letzter_fehler})"
        )


    # ------------------------------------------------------
    # Einzelnes Molekuel
    # ------------------------------------------------------


    def pruefe_schluessel(self, schluessel):
        """Ergaenzt das Schluessel-Dict um 'status' und 'cid'."""

        ergebnis = dict(schluessel)
        ergebnis["cid"] = None

        zeile = self._cache_lesen(schluessel["inchikey"])

        if zeile is not None:
            ergebnis["status"], ergebnis["cid"] = zeile
            return ergebnis

        if self.offline:
            ergebnis["status"] = "ungeprueft"
            return ergebnis

        try:
            status, cid = self._frage_pubchem(schluessel["inchikey"])
        except PubChemFehler as fehler:
            ergebnis["status"] = "ungeprueft"
            ergebnis["fehler"] = str(fehler)
            return ergebnis

        self._cache_schreiben(schluessel["inchikey"], status, cid)

        ergebnis["status"] = status
        ergebnis["cid"] = cid

        return ergebnis


    def pruefe_molekuel(self, molekuel):

        schluessel = struktur_schluessel(molekuel)

        if schluessel is None:
            return {"status": "ungueltig"}

        return self.pruefe_schluessel(schluessel)


    # ------------------------------------------------------
    # Ergebnisliste des Labors
    # ------------------------------------------------------


    def pruefe_ergebnisse(
            self,
            ergebnisse,
            min_punkte=90,
            nur_rdkit_gueltig=True,
            limit=None,
            on_fortschritt=print,
            abbruch=None
    ):
        """ergebnisse: Eintraege im Format des Labors, also
        (nummer, phase, molekuel, analyse) - z.B. labor.alle_ergebnisse.

        Filtert nach Kritiker-Punkten (und RDKit-Gueltigkeit), entfernt
        Duplikate ueber den InChIKey (es bleibt der Eintrag mit den
        meisten Punkten) und fragt nur diese eindeutigen Kandidaten ab.
        limit begrenzt die Zahl der Kandidaten (beste zuerst).
        abbruch: optionales threading.Event - wird es gesetzt, endet
        die Pruefung nach der aktuellen Anfrage, und die bis dahin
        geprueften Eintraege werden zurueckgegeben.

        Rueckgabe: Liste von Dicts, 'unbekannt' zuerst, dann nach
        Punkten absteigend."""

        beste = {}

        for nummer, phase, molekuel, analyse in ergebnisse:

            punkte = analyse["bewertung"]["punkte"]

            if punkte < min_punkte:
                continue

            if nur_rdkit_gueltig and not analyse.get("rdkit_gueltig"):
                continue

            schluessel = struktur_schluessel(molekuel)

            if schluessel is None:
                continue

            key = schluessel["inchikey"]

            if key in beste and beste[key]["punkte"] >= punkte:
                continue

            beste[key] = {
                **schluessel,
                "nr": nummer,
                "phase": phase,
                "punkte": punkte,
                "molekuel": molekuel,
            }

        kandidaten = sorted(
            beste.values(),
            key=lambda eintrag: -eintrag["punkte"]
        )

        if limit is not None:
            kandidaten = kandidaten[:limit]

        on_fortschritt(
            f"Neuheitsfilter: {len(kandidaten)} eindeutige Kandidaten "
            f"(>= {min_punkte} Punkte)."
        )

        ergebnis = []

        for i, kandidat in enumerate(kandidaten, start=1):

            if abbruch is not None and abbruch.is_set():
                on_fortschritt(
                    f"Neuheitsfilter abgebrochen nach {i - 1} von "
                    f"{len(kandidaten)} Kandidaten."
                )
                break

            pruefung = self.pruefe_schluessel(kandidat)

            kandidat["status"] = pruefung["status"]
            kandidat["cid"] = pruefung["cid"]

            if "fehler" in pruefung:
                kandidat["fehler"] = pruefung["fehler"]

            ergebnis.append(kandidat)

            if i % 25 == 0:
                on_fortschritt(
                    f"Neuheitsfilter: {i}/{len(kandidaten)} geprueft"
                )

        ergebnis.sort(
            key=lambda eintrag: (
                eintrag["status"] != "unbekannt",
                -eintrag["punkte"]
            )
        )

        zaehler = {}

        for eintrag in ergebnis:
            zaehler[eintrag["status"]] = (
                zaehler.get(eintrag["status"], 0) + 1
            )

        on_fortschritt(f"Neuheitsfilter fertig: {zaehler}")

        return ergebnis


    # ------------------------------------------------------
    # Export
    # ------------------------------------------------------


    @staticmethod
    def exportiere(ergebnisse, pfad, nur_status=("unbekannt",)):
        """Schreibt die Eintraege mit passendem Status als JSON.
        Das interne Molekuel-Dict bleibt dabei erhalten, damit sich
        die Entdeckungen spaeter z.B. ins MoleculeBlueprint laden oder
        erneut pruefen lassen."""

        auswahl = [
            eintrag for eintrag in ergebnisse
            if nur_status is None or eintrag["status"] in nur_status
        ]

        with open(pfad, "w", encoding="utf-8") as datei:
            json.dump(auswahl, datei, indent=2, ensure_ascii=False)

        return len(auswahl)


# ==========================================================
# Schnelltest (live gegen PubChem, wenn online)
# ==========================================================
#
#   python neuheit.py
#
# Erwartung: Wasser und Ethanol "bekannt", ein erfundenes
# Molekuel (hier ein Fluor-Stickstoff-Ring) vermutlich "unbekannt".


if __name__ == "__main__":

    def _mol(elemente, bindungen):
        return {
            "atoms": [
                {"id": i, "element": e}
                for i, e in enumerate(elemente)
            ],
            "bonds": [
                {"atom1": a, "atom2": b, "order": o}
                for a, b, o in bindungen
            ],
        }

    proben = {
        "Wasser": _mol(["O", "H", "H"], [(0, 1, 1), (0, 2, 1)]),
        "Ethanol": _mol(
            ["C", "C", "O", "H", "H", "H", "H", "H", "H"],
            [(0, 1, 1), (1, 2, 1), (0, 3, 1), (0, 4, 1), (0, 5, 1),
             (1, 6, 1), (1, 7, 1), (2, 8, 1)]
        ),
    }

    pruefer = NeuheitsPruefer()

    for name, molekuel in proben.items():
        print(name, "->", pruefer.pruefe_molekuel(molekuel))