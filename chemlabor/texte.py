# ==========================================================
# texte.py
# ==========================================================
#
# Zentrale Uebersetzung aller sichtbaren Texte.
#
#     from chemlabor.texte import t
#     label = ttk.Label(self, text=t("live.experimente"))
#     status.config(text=t("live.trend", text=trend))
#
# Die Sprache wird EINMAL beim Start festgelegt (Umgebungsvariable
# CHEMLABOR_SPRACHE oder setze_sprache() vor dem Bauen des Fensters).
# Standard: Englisch. Ein Umschalten waehrend des Betriebs braeuchte ein
# Neuaufbauen aller Widgets - bewusst nicht enthalten.
#
# Fehlt ein Schluessel in der aktiven Sprache, wird auf Deutsch
# zurueckgefallen, danach auf den Schluessel selbst. tests/test_texte.py
# stellt sicher, dass beide Sprachen vollstaendig sind und dieselben
# Platzhalter benutzen.
#
# Neue Texte: Schluessel in BEIDEN Bloecken unten eintragen.
# Platzhalter in {geschweiften Klammern} muessen in beiden Sprachen
# identisch heissen.
#
# NICHT hier: die deutschen Texte, die der Kritiker erzeugt (Urteil,
# Erklaerung, Probleme) - die werden erst bei der Anzeige uebersetzt,
# siehe chemlabor/gui/kritiker_texte.py.

import os


TEXTE = {

    "de": {

        # ==================================================
        # Tab "Live"
        # ==================================================
        "live.experimente": "Experimente",
        "live.mittel": "Ø Punkte",
        "live.stabil": "Stabil (≥{ab})",
        "live.spanne": "Min – Max",
        "live.phasen": "Beobachtung / Automatik",
        "live.versuche": "({n} Versuche)",
        "live.stabil_wert": "{n} ({p} %)",
        "live.filter": "Verhältnis-Filter: {text}",
        "live.ausschnitt": "(* Durchschnitt über die letzten {n})",
        "live.trend": "Lerntrend: {text}",

        # ==================================================
        # Phasen und Filter (Anzeigetext; interne Werte bleiben)
        # ==================================================
        "phase.alle": "Alle",
        "phase.Beobachtung": "Beobachtung",
        "phase.Automatik": "Automatik",

        "filter.alle": "Alle",
        "filter.top10": "Top 10",
        "filter.top50": "Top 50",
        "filter.stabil": "Nur stabil (≥{ab})",
        "filter.instabil": "Nur instabil (<50)",

        "neufilter.alle": "Alle",
        "neufilter.unbekannt": "Nur unbekannt",
        "neufilter.bekannt": "Nur bekannt",
        "neufilter.ungeprueft": "Nur ungeprüft",

        "status.unbekannt": "unbekannt",
        "status.bekannt": "bekannt",
        "status.ungeprueft": "ungeprüft",

        # ==================================================
        # Tab "Ergebnisse"
        # ==================================================
        "spalte.nr": "Nr",
        "spalte.phase": "Phase",
        "spalte.formel": "Formel",
        "spalte.atome": "Atome",
        "spalte.punkte": "Punkte",
        "spalte.urteil": "Urteil",
        "spalte.rdkit": "RDKit",
        "spalte.pubchem": "PubChem",

        "rdkit.ja_zelle": "ja",
        "rdkit.nein_zelle": "nein",

        "erg.filter": "Filter:",
        "erg.phase": "Phase:",
        "erg.duplikate": "Wiederholungen ausblenden",
        "erg.neuheit_button": "Neuheit prüfen (PubChem)",
        "erg.abbrechen": "■ Abbrechen",
        "erg.offline": "nur Cache (offline)",
        "erg.pubchem_anzeige": "PubChem-Anzeige:",
        "erg.export_button": "Exportieren…",
        "erg.anzahl": "{n} angezeigt · {gesamt} Experimente gesamt",
        "erg.breche_ab": "Breche ab...",
        "erg.keine_zeilen": "Keine Zeilen sichtbar - erst einen Lauf starten oder Filter lockern.",
        "erg.starte": "Starte PubChem-Prüfung...",
        "erg.fortschritt": "PubChem: {i}/{n} Zeilen geprüft...",
        "erg.zusammenfassung": "{u} unbekannt ({s} davon mit offenem Stereo), {b} bekannt, {g} ungeprüft.",
        "erg.abgebrochen": " (abgebrochen)",
        "erg.limit": " – {n} Zeilen über dem Limit von {max} Strukturen nicht geprüft.",
        "erg.log_pruefung": "PubChem-Prüfung: {text}",
        "erg.pruefung_fehler": "PubChem-Prüfung fehlgeschlagen: {fehler}",
        "erg.nichts_export": "Nichts zu exportieren.",
        "erg.export_titel": "Ergebnisse exportieren",
        "erg.exportiert": "{n} Einträge exportiert: {pfad}",
        "erg.export_fehler": "Export fehlgeschlagen: {fehler}",

        # ==================================================
        # Detailansicht
        # ==================================================
        "detail.experiment": "Experiment {nummer}  ({phase})",
        "detail.formel": "Formel: {formel}",
        "detail.punkte": "Punkte: {punkte}  –  {urteil}",
        "detail.rdkit_gueltig": "RDKit-gültig: {antwort}",
        "detail.ja": "Ja",
        "detail.nein": "Nein",
        "detail.rdkit_formel": "RDKit-Summenformel: {formel}",
        "detail.molmasse": "Molmasse: {masse} g/mol",
        "detail.smiles": "Kanonischer SMILES: {smiles}",
        "detail.pubchem": "PubChem: {status}",
        "detail.cid": "CID {cid} (Doppelklick auf die Zeile öffnet den Eintrag)",
        "detail.smiles_ohne_h": "SMILES ohne H: {smiles}",
        "detail.inchikey": "InChIKey: {key}",
        "detail.stereo": "Achtung: Stereo offen – PubChem kennt evtl. nur konkrete Stereoisomere, 'unbekannt' bitte von Hand gegenprüfen.",
        "detail.abfrage_fehler": "Fehler bei der Abfrage: {fehler}",
        "detail.probleme": "Probleme:",
        "detail.keine": " - keine",
        "detail.erklaerung": "Erklärung:",

        # ==================================================
        # Kurzfassung der Einstellungen / Verhaeltnis-Filter
        # ==================================================
        "cfg.beobachtungen": "{n} Beobachtungen",
        "cfg.endlos": "Automatik endlos",
        "cfg.automatik": "{n} Automatik",
        "cfg.max_atome": "max. {n} Atome",
        "cfg.konstruktion": "Konstruktion",
        "cfg.zufall": "Zufall (Best-of-N)",
        "cfg.filter_aus": "Filter aus",
        "cfg.filter": "Filter {n}:1",
        "cfg.autosave": "Autosave",

        "verh.aus": "aus",
        "verh.stand": "{ziel}:1 – {pos} positiv, {neg} negativ, {verw} verworfen",

        # ==================================================
        # Lernkurve
        # ==================================================
        "lern.keine_daten": "Noch keine Daten - die Kurve füllt sich nach den ersten Beobachtungen.",
        "lern.zu_wenig": "Noch zu wenig Daten für eine Kurve ({n} Intervall(e) abgeschlossen).",
        "lern.status": "{n} Intervalle (je {m} Beobachtungen) - {trend}",
        "lern.achse": "Beobachtungen: {a} – {b}",

        # ==================================================
        # Lewis-Ansicht
        # ==================================================
        "lewis.auswaehlen": "Molekül in der Tabelle auswählen.",
        "lewis.keine_bindungen": "Keine Bindungsdaten vorhanden.",
        "lewis.valenz_fehlt": "Valenzdaten fehlen.",
        "lewis.valenz_ungueltig": "Ungültige Valenz bei {element}.",
        "lewis.bindung_unvollstaendig": "Unvollständige Bindungsdaten.",
        "lewis.moeglich": "Darstellung möglich.",
        "lewis.nicht_verfuegbar": "Lewis nicht verfügbar: {grund}",
        "lewis.berechne": "Layout wird berechnet ({n} Atome, {m} Bindungen)...",
        "lewis.fertig": "Lewis-Struktur erzeugt ({n} Atome, {m} Bindungen).",
        "lewis.rohdaten": "Rohdaten:",
        "lewis.rohdaten_atom": "Atom {id}: {element} (Valenz {valenz}, belegt {belegt})",
        "lewis.rohdaten_bindung": "Bindung: {a} - {b} (Ordnung {ordnung})",

        # ==================================================
        # Einstellungen-Dialog
        # ==================================================
        "dlg.titel": "Einstellungen",
        "dlg.beobachtungen": "Beobachtungen",
        "dlg.beobachtungen_hinweis": "Zufallsmoleküle, aus denen die KI zuerst lernt",
        "dlg.automatik": "Automatik-Experimente",
        "dlg.automatik_hinweis": "danach von der KI selbst erzeugt",
        "dlg.kandidaten": "Kandidaten pro Molekül",
        "dlg.kandidaten_hinweis": "Best-of-N im KI-Generator",
        "dlg.max_atome": "Max. Atome",
        "dlg.max_atome_hinweis": "Obergrenze je Molekül",
        "dlg.verhaeltnis": "Verhältnis-Filter",
        "dlg.filter_aktiv": "aktiv, höchstens",
        "dlg.filter_einheit": "schlechte je gutes Molekül",
        "dlg.modus": "Automatik-Modus",
        "dlg.modus_zufall": "Zufall (Best-of-N)",
        "dlg.modus_konstruktion": "Konstruktion (Schritt für Schritt mit RDKit)",
        "dlg.endlos": "Endlos: Automatikphase läuft, bis du auf Stopp drückst",
        "dlg.autosave": "Gedächtnis automatisch sichern",
        "dlg.autosave_pfad": "nach {pfad}",
        "dlg.standardwerte": "Standardwerte",
        "dlg.abbrechen": "Abbrechen",
        "dlg.ok": "OK",
        "dlg.keine_ganzzahl": "{name}: '{text}' ist keine ganze Zahl.",

        # ==================================================
        # Hauptfenster
        # ==================================================
        "fenster.titel": "Chemie-Labor",
        "fenster.start": "▶  Start",
        "fenster.stopp": "■  Stopp",
        "fenster.einstellungen": "⚙  Einstellungen…",
        "fenster.gedaechtnis": "Gedächtnis ▾",
        "fenster.menue_speichern": "Speichern…",
        "fenster.menue_laden": "Laden…",
        "fenster.menue_zusammenfuehren": "Mehrere zusammenführen…",
        "fenster.menue_verwerfen": "Geladenes Gedächtnis verwerfen",
        "fenster.bereit": "Bereit.",
        "fenster.tab_live": "Live",
        "fenster.tab_ergebnisse": "Ergebnisse",
        "fenster.log_offen": "Log ▾",
        "fenster.log_zu": "Log ▸",
        "fenster.log_leeren": "Leeren",
        "fenster.gedaechtnis_geladen": "Gedächtnis: geladen ({n} Beobachtungen)",
        "fenster.gedaechtnis_neu": "Gedächtnis: neu",
        "fenster.labor_fehler_start": "Labor konnte nicht gestartet werden: {fehler}",
        "fenster.gedaechtnis_uebernommen": "Geladenes Gedächtnis übernommen ({n} frühere Beobachtungen).",
        "fenster.autosave_aktiv": "Automatisches Sichern aktiv - alle {n} Experimente nach '{pfad}'.",
        "fenster.laeuft_endlos": "Läuft (Endlos)...",
        "fenster.laeuft": "Läuft...",
        "fenster.fehler_praefix": "Fehler: {fehler}",
        "fenster.stoppe": "Stoppe...",
        "fenster.fertig": "Fertig.",
        "fenster.kein_gedaechtnis": "Kein Gedächtnis vorhanden - erst einen Lauf starten.",
        "fenster.speichern_titel": "Gedächtnis speichern",
        "fenster.gespeichert": "Gedächtnis gespeichert: {pfad}",
        "fenster.speichern_fehler": "Speichern fehlgeschlagen: {fehler}",
        "fenster.laden_erst_stoppen": "Erst den laufenden Lauf stoppen, bevor ein Gedächtnis geladen wird.",
        "fenster.laden_titel": "Gedächtnis laden",
        "fenster.laden_fehler": "Laden fehlgeschlagen: {fehler}",
        "fenster.geladen": "Gedächtnis geladen - wird beim nächsten Start verwendet.",
        "fenster.merge_erst_stoppen": "Erst den laufenden Lauf stoppen, bevor Gedächtnisse zusammengeführt werden.",
        "fenster.merge_titel": "Gedächtnisse zum Zusammenführen auswählen (mind. 2)",
        "fenster.merge_zu_wenig": "Bitte mindestens 2 Dateien auswählen.",
        "fenster.merge_fehler": "Zusammenführen fehlgeschlagen: {fehler}",
        "fenster.merge_speichern_titel": "Zusammengeführtes Gedächtnis speichern (optional)",
        "fenster.merge_speichern_fehler": "Zusammengeführt, aber Speichern fehlgeschlagen: {fehler}",
        "fenster.merge_fertig": "{n} Dateien zusammengeführt ({b} Beobachtungen) - wird beim nächsten Start verwendet.",
        "fenster.verwerfen_erst_stoppen": "Erst den laufenden Lauf stoppen.",
        "fenster.verworfen": "Geladenes Gedächtnis verworfen - der nächste Lauf beginnt neu.",

        "datei.json": "JSON-Dateien",
        "datei.alle": "Alle Dateien",

        # ==================================================
        # KIGedaechtnis (ki/gedaechtnis.py)
        # ==================================================
        "ged.trend_zu_wenig": "Noch zu wenige Daten für einen Lerntrend.",
        "ged.trend_besser": "Verbessert sich (+{d} Punkte im Schnitt seit Beginn).",
        "ged.trend_schlechter": "Verschlechtert sich ({d} Punkte im Schnitt seit Beginn).",
        "ged.trend_stabil": "Stabil, kein deutlicher Trend.",
        "ged.kein_lauf": "Noch kein Lauf verzeichnet.",
        "ged.laeufe": "{n} Lauf/Läufe verzeichnet - {details}",
        "ged.beobachtet": "{n} Molekuele beobachtet.",
        "ged.durchschnitt": "Durchschnittliche Punktzahl aller beobachteten Molekuele: {x}",
        "ged.schwelle": "Adaptive Stabilitaetsschwelle (oberes {p}%-Quartil): {s} Punkte",
        "ged.temperatur": "Softmax-Temperatur: {wert} (Minimum {minimum})",
        "ged.beste_elemente": "Beste Elemente (Bayes-geglaettet): {liste}",
        "ged.bindung_item": "{e1}-{e2} (Ordnung {o}, Ø{x})",
        "ged.beste_bindungen": "Beste Bindungsmuster (Bayes-geglaettet): {liste}",
        "ged.ring_item": "{g}-Ring (Ø{x})",
        "ged.beste_ringe": "Beste Ringgroessen (Bayes-geglaettet): {liste}",
        "ged.lernkurve": "Lernkurve ({n} Punkte, alle {i} Beobachtungen): {trend}",
        "ged.lauf_historie": "Lauf-Historie: {text}",
        "ged.keine_datei": "Mindestens eine Datei wird benötigt.",

        # ==================================================
        # ChemLabor (labor.py) - Meldungen im Log
        # ==================================================
        "lab.konstruktion_gescheitert": "Konstruktion gescheitert (RDKit-Backtracking ohne gültiges Ergebnis)",
        "lab.autosave_ok": "Auto-gespeichert nach {n} Experimenten ({pfad}).",
        "lab.autosave_fehler": "Autosave fehlgeschlagen: {fehler}",
        "lab.fehler": "Fehler ({typ}): {fehler}",
        "lab.weitere_fehler": "Weitere Fehler vom Typ {typ} werden nur noch gezaehlt, nicht mehr einzeln geloggt.",
        "lab.beobachtungsphase": "=== BEOBACHTUNGSPHASE ({n} Durchläufe) ===",
        "lab.beobachtung_fortschritt": "[Beobachtung] {i}/{n} ({p} Punkte)",
        "lab.automatik_endlos": "=== AUTOMATIKPHASE (Endlos-Modus) ===",
        "lab.automatikphase": "=== AUTOMATIKPHASE ({n} Durchläufe) ===",
        "lab.automatik_endlos_fortschritt": "[Automatik/Endlos] {n} Molekuele erzeugt ({p} Punkte)",
        "lab.automatik_fortschritt": "[Automatik] {i}/{n} ({p} Punkte)",
        "lab.titel": "=== KÜNSTLICHES CHEMIE-LABOR v0.7 ===",
        "lab.keine": "keine",
        "lab.beendet": "Labor beendet. Versuche: {v}, Experimente: {e}, Fehler: {f}",
        "lab.verhaeltnis_ende": "Verhaeltnis-Filter (Ziel {ziel}:1 neg:pos): {pos} positiv, {neg} negativ, davon {verw} verworfen.",

        # ==================================================
        # Feldnamen und Fehlermeldungen der Config
        # ==================================================
        "feld.beobachtungen": "Beobachtungen",
        "feld.automatik_experimente": "Automatik-Experimente",
        "feld.kandidaten": "Kandidaten pro Molekuel",
        "feld.max_atome": "Max. Atome",
        "feld.verhaeltnis_sockel": "Verhaeltnis-Sockel",
        "feld.versuche_faktor": "Versuchsfaktor",
        "feld.anzeigeschwelle": "Anzeigeschwelle",
        "feld.fehler_log_limit": "Fehler-Log-Limit",
        "feld.endlos_log": "Endlos-Log-Intervall",
        "feld.autosave_intervall": "Autosave-Intervall",
        "feld.verhaeltnis": "Verhaeltnis Neg:Pos",

        "config.ganzzahl": "{name} muss eine ganze Zahl sein (ist: {wert}).",
        "config.minimum": "{name} muss mindestens {minimum} sein (ist: {wert}).",
        "config.modus": "Automatik-Modus muss einer von {modi} sein (ist: {wert}).",
    },

    "en": {

        # ==================================================
        # "Live" tab
        # ==================================================
        "live.experimente": "Experiments",
        "live.mittel": "Avg. score",
        "live.stabil": "Stable (≥{ab})",
        "live.spanne": "Min – Max",
        "live.phasen": "Observation / Automatic",
        "live.versuche": "({n} attempts)",
        "live.stabil_wert": "{n} ({p}%)",
        "live.filter": "Ratio filter: {text}",
        "live.ausschnitt": "(* average over the last {n})",
        "live.trend": "Learning trend: {text}",

        # ==================================================
        # Phases and filters (display text; internal values stay)
        # ==================================================
        "phase.alle": "All",
        "phase.Beobachtung": "Observation",
        "phase.Automatik": "Automatic",

        "filter.alle": "All",
        "filter.top10": "Top 10",
        "filter.top50": "Top 50",
        "filter.stabil": "Stable only (≥{ab})",
        "filter.instabil": "Unstable only (<50)",

        "neufilter.alle": "All",
        "neufilter.unbekannt": "Unknown only",
        "neufilter.bekannt": "Known only",
        "neufilter.ungeprueft": "Unchecked only",

        "status.unbekannt": "unknown",
        "status.bekannt": "known",
        "status.ungeprueft": "unchecked",

        # ==================================================
        # "Results" tab
        # ==================================================
        "spalte.nr": "No.",
        "spalte.phase": "Phase",
        "spalte.formel": "Formula",
        "spalte.atome": "Atoms",
        "spalte.punkte": "Score",
        "spalte.urteil": "Verdict",
        "spalte.rdkit": "RDKit",
        "spalte.pubchem": "PubChem",

        "rdkit.ja_zelle": "yes",
        "rdkit.nein_zelle": "no",

        "erg.filter": "Filter:",
        "erg.phase": "Phase:",
        "erg.duplikate": "Hide repeats",
        "erg.neuheit_button": "Check novelty (PubChem)",
        "erg.abbrechen": "■ Cancel",
        "erg.offline": "cache only (offline)",
        "erg.pubchem_anzeige": "PubChem view:",
        "erg.export_button": "Export…",
        "erg.anzahl": "{n} shown · {gesamt} experiments in total",
        "erg.breche_ab": "Cancelling...",
        "erg.keine_zeilen": "No rows visible - start a run first or loosen the filters.",
        "erg.starte": "Starting PubChem check...",
        "erg.fortschritt": "PubChem: {i}/{n} rows checked...",
        "erg.zusammenfassung": "{u} unknown ({s} of them with open stereo), {b} known, {g} unchecked.",
        "erg.abgebrochen": " (cancelled)",
        "erg.limit": " – {n} rows beyond the limit of {max} structures were not checked.",
        "erg.log_pruefung": "PubChem check: {text}",
        "erg.pruefung_fehler": "PubChem check failed: {fehler}",
        "erg.nichts_export": "Nothing to export.",
        "erg.export_titel": "Export results",
        "erg.exportiert": "{n} entries exported: {pfad}",
        "erg.export_fehler": "Export failed: {fehler}",

        # ==================================================
        # Detail view
        # ==================================================
        "detail.experiment": "Experiment {nummer}  ({phase})",
        "detail.formel": "Formula: {formel}",
        "detail.punkte": "Score: {punkte}  –  {urteil}",
        "detail.rdkit_gueltig": "RDKit valid: {antwort}",
        "detail.ja": "Yes",
        "detail.nein": "No",
        "detail.rdkit_formel": "RDKit molecular formula: {formel}",
        "detail.molmasse": "Molar mass: {masse} g/mol",
        "detail.smiles": "Canonical SMILES: {smiles}",
        "detail.pubchem": "PubChem: {status}",
        "detail.cid": "CID {cid} (double-clicking the row opens the entry)",
        "detail.smiles_ohne_h": "SMILES without H: {smiles}",
        "detail.inchikey": "InChIKey: {key}",
        "detail.stereo": "Note: stereo open – PubChem may only list specific stereoisomers, please double-check 'unknown' by hand.",
        "detail.abfrage_fehler": "Query error: {fehler}",
        "detail.probleme": "Problems:",
        "detail.keine": " - none",
        "detail.erklaerung": "Explanation:",

        # ==================================================
        # Settings summary / ratio filter
        # ==================================================
        "cfg.beobachtungen": "{n} observations",
        "cfg.endlos": "endless automatic",
        "cfg.automatik": "{n} automatic",
        "cfg.max_atome": "max. {n} atoms",
        "cfg.konstruktion": "construction",
        "cfg.zufall": "random (best-of-N)",
        "cfg.filter_aus": "filter off",
        "cfg.filter": "filter {n}:1",
        "cfg.autosave": "autosave",

        "verh.aus": "off",
        "verh.stand": "{ziel}:1 – {pos} positive, {neg} negative, {verw} discarded",

        # ==================================================
        # Learning curve
        # ==================================================
        "lern.keine_daten": "No data yet - the curve fills up after the first observations.",
        "lern.zu_wenig": "Not enough data for a curve yet ({n} interval(s) completed).",
        "lern.status": "{n} intervals ({m} observations each) - {trend}",
        "lern.achse": "Observations: {a} – {b}",

        # ==================================================
        # Lewis view
        # ==================================================
        "lewis.auswaehlen": "Select a molecule in the table.",
        "lewis.keine_bindungen": "No bond data available.",
        "lewis.valenz_fehlt": "Valence data missing.",
        "lewis.valenz_ungueltig": "Invalid valence on {element}.",
        "lewis.bindung_unvollstaendig": "Incomplete bond data.",
        "lewis.moeglich": "Drawing possible.",
        "lewis.nicht_verfuegbar": "Lewis structure unavailable: {grund}",
        "lewis.berechne": "Calculating layout ({n} atoms, {m} bonds)...",
        "lewis.fertig": "Lewis structure generated ({n} atoms, {m} bonds).",
        "lewis.rohdaten": "Raw data:",
        "lewis.rohdaten_atom": "Atom {id}: {element} (valence {valenz}, used {belegt})",
        "lewis.rohdaten_bindung": "Bond: {a} - {b} (order {ordnung})",

        # ==================================================
        # Settings dialog
        # ==================================================
        "dlg.titel": "Settings",
        "dlg.beobachtungen": "Observations",
        "dlg.beobachtungen_hinweis": "random molecules the AI learns from first",
        "dlg.automatik": "Automatic experiments",
        "dlg.automatik_hinweis": "then generated by the AI itself",
        "dlg.kandidaten": "Candidates per molecule",
        "dlg.kandidaten_hinweis": "best-of-N in the AI generator",
        "dlg.max_atome": "Max. atoms",
        "dlg.max_atome_hinweis": "upper limit per molecule",
        "dlg.verhaeltnis": "Ratio filter",
        "dlg.filter_aktiv": "enabled, at most",
        "dlg.filter_einheit": "bad per good molecule",
        "dlg.modus": "Automatic mode",
        "dlg.modus_zufall": "Random (best-of-N)",
        "dlg.modus_konstruktion": "Construction (step by step with RDKit)",
        "dlg.endlos": "Endless: the automatic phase runs until you press Stop",
        "dlg.autosave": "Save memory automatically",
        "dlg.autosave_pfad": "to {pfad}",
        "dlg.standardwerte": "Defaults",
        "dlg.abbrechen": "Cancel",
        "dlg.ok": "OK",
        "dlg.keine_ganzzahl": "{name}: '{text}' is not a whole number.",

        # ==================================================
        # Main window
        # ==================================================
        "fenster.titel": "Chemistry Lab",
        "fenster.start": "▶  Start",
        "fenster.stopp": "■  Stop",
        "fenster.einstellungen": "⚙  Settings…",
        "fenster.gedaechtnis": "Memory ▾",
        "fenster.menue_speichern": "Save…",
        "fenster.menue_laden": "Load…",
        "fenster.menue_zusammenfuehren": "Merge several…",
        "fenster.menue_verwerfen": "Discard loaded memory",
        "fenster.bereit": "Ready.",
        "fenster.tab_live": "Live",
        "fenster.tab_ergebnisse": "Results",
        "fenster.log_offen": "Log ▾",
        "fenster.log_zu": "Log ▸",
        "fenster.log_leeren": "Clear",
        "fenster.gedaechtnis_geladen": "Memory: loaded ({n} observations)",
        "fenster.gedaechtnis_neu": "Memory: new",
        "fenster.labor_fehler_start": "Could not start the lab: {fehler}",
        "fenster.gedaechtnis_uebernommen": "Loaded memory applied ({n} earlier observations).",
        "fenster.autosave_aktiv": "Automatic saving enabled - every {n} experiments to '{pfad}'.",
        "fenster.laeuft_endlos": "Running (endless)...",
        "fenster.laeuft": "Running...",
        "fenster.fehler_praefix": "Error: {fehler}",
        "fenster.stoppe": "Stopping...",
        "fenster.fertig": "Done.",
        "fenster.kein_gedaechtnis": "No memory available - start a run first.",
        "fenster.speichern_titel": "Save memory",
        "fenster.gespeichert": "Memory saved: {pfad}",
        "fenster.speichern_fehler": "Saving failed: {fehler}",
        "fenster.laden_erst_stoppen": "Stop the running run before loading a memory.",
        "fenster.laden_titel": "Load memory",
        "fenster.laden_fehler": "Loading failed: {fehler}",
        "fenster.geladen": "Memory loaded - will be used on the next start.",
        "fenster.merge_erst_stoppen": "Stop the running run before merging memories.",
        "fenster.merge_titel": "Select memories to merge (at least 2)",
        "fenster.merge_zu_wenig": "Please select at least 2 files.",
        "fenster.merge_fehler": "Merging failed: {fehler}",
        "fenster.merge_speichern_titel": "Save merged memory (optional)",
        "fenster.merge_speichern_fehler": "Merged, but saving failed: {fehler}",
        "fenster.merge_fertig": "{n} files merged ({b} observations) - will be used on the next start.",
        "fenster.verwerfen_erst_stoppen": "Stop the running run first.",
        "fenster.verworfen": "Loaded memory discarded - the next run starts fresh.",

        "datei.json": "JSON files",
        "datei.alle": "All files",

        # ==================================================
        # KIGedaechtnis (ki/gedaechtnis.py)
        # ==================================================
        "ged.trend_zu_wenig": "Not enough data for a learning trend yet.",
        "ged.trend_besser": "Improving (+{d} points on average since the start).",
        "ged.trend_schlechter": "Getting worse ({d} points on average since the start).",
        "ged.trend_stabil": "Stable, no clear trend.",
        "ged.kein_lauf": "No run recorded yet.",
        "ged.laeufe": "{n} run(s) recorded - {details}",
        "ged.beobachtet": "{n} molecules observed.",
        "ged.durchschnitt": "Average score of all observed molecules: {x}",
        "ged.schwelle": "Adaptive stability threshold (upper {p}% quartile): {s} points",
        "ged.temperatur": "Softmax temperature: {wert} (minimum {minimum})",
        "ged.beste_elemente": "Best elements (Bayes-smoothed): {liste}",
        "ged.bindung_item": "{e1}-{e2} (order {o}, Ø{x})",
        "ged.beste_bindungen": "Best bond patterns (Bayes-smoothed): {liste}",
        "ged.ring_item": "{g}-ring (Ø{x})",
        "ged.beste_ringe": "Best ring sizes (Bayes-smoothed): {liste}",
        "ged.lernkurve": "Learning curve ({n} points, every {i} observations): {trend}",
        "ged.lauf_historie": "Run history: {text}",
        "ged.keine_datei": "At least one file is required.",

        # ==================================================
        # ChemLabor (labor.py) - log messages
        # ==================================================
        "lab.konstruktion_gescheitert": "Construction failed (RDKit backtracking without a valid result)",
        "lab.autosave_ok": "Auto-saved after {n} experiments ({pfad}).",
        "lab.autosave_fehler": "Autosave failed: {fehler}",
        "lab.fehler": "Error ({typ}): {fehler}",
        "lab.weitere_fehler": "Further errors of type {typ} are only counted, no longer logged individually.",
        "lab.beobachtungsphase": "=== OBSERVATION PHASE ({n} runs) ===",
        "lab.beobachtung_fortschritt": "[Observation] {i}/{n} ({p} points)",
        "lab.automatik_endlos": "=== AUTOMATIC PHASE (endless mode) ===",
        "lab.automatikphase": "=== AUTOMATIC PHASE ({n} runs) ===",
        "lab.automatik_endlos_fortschritt": "[Automatic/endless] {n} molecules generated ({p} points)",
        "lab.automatik_fortschritt": "[Automatic] {i}/{n} ({p} points)",
        "lab.titel": "=== ARTIFICIAL CHEMISTRY LAB v0.7 ===",
        "lab.keine": "none",
        "lab.beendet": "Lab finished. Attempts: {v}, experiments: {e}, errors: {f}",
        "lab.verhaeltnis_ende": "Ratio filter (target {ziel}:1 neg:pos): {pos} positive, {neg} negative, of which {verw} discarded.",

        # ==================================================
        # Config field names and error messages
        # ==================================================
        "feld.beobachtungen": "Observations",
        "feld.automatik_experimente": "Automatic experiments",
        "feld.kandidaten": "Candidates per molecule",
        "feld.max_atome": "Max. atoms",
        "feld.verhaeltnis_sockel": "Ratio floor",
        "feld.versuche_faktor": "Attempt factor",
        "feld.anzeigeschwelle": "Display threshold",
        "feld.fehler_log_limit": "Error log limit",
        "feld.endlos_log": "Endless-mode log interval",
        "feld.autosave_intervall": "Autosave interval",
        "feld.verhaeltnis": "Ratio neg:pos",

        "config.ganzzahl": "{name} must be an integer (got: {wert}).",
        "config.minimum": "{name} must be at least {minimum} (got: {wert}).",
        "config.modus": "Automatic mode must be one of {modi} (got: {wert}).",
    },
}


def _start_sprache():

    code = os.environ.get("CHEMLABOR_SPRACHE", "en").strip().lower()

    return code if code in TEXTE else "en"


_aktiv = _start_sprache()


def sprache():
    """Aktive Sprache als Code ('de' oder 'en')."""

    return _aktiv


def setze_sprache(code):
    """Aktive Sprache setzen - vor dem Bauen des Fensters aufrufen."""

    global _aktiv

    if code not in TEXTE:
        raise ValueError(
            f"Unknown language {code!r}; available: {sorted(TEXTE)}"
        )

    _aktiv = code


def t_oder(schluessel, standard):
    """Wie t(), liefert aber 'standard' statt des Schluessels, wenn es
    keinen Eintrag gibt - fuer Werte aus der Logik, die nicht
    vollstaendig uebersetzt sein muessen (z.B. ein unbekannter Status).
    Der Text wird NICHT formatiert."""

    text = TEXTE[_aktiv].get(schluessel)

    if text is None:
        text = TEXTE["de"].get(schluessel)

    return standard if text is None else text


def t(schluessel, **werte):
    """Text zum Schluessel in der aktiven Sprache; Platzhalter werden
    mit 'werte' gefuellt."""

    text = t_oder(schluessel, schluessel)

    return text.format(**werte) if werte else text


def zahl(n):
    """Ganzzahl mit Tausendertrennung: 12,345 (en) bzw. 12.345 (de)."""

    text = f"{n:,}"

    return text.replace(",", ".") if _aktiv == "de" else text