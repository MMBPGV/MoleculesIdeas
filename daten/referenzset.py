# ==========================================================
# referenzset.py
# ==========================================================
#
# Kleines Referenzset bekannter Molekuele zum EICHEN des Kritikers:
# Wo vergibt er zu viele, wo zu wenige Punkte? Gedacht als Messlatte fuer
# jede Aenderung am Kritiker (siehe eichung.py im Projektordner) - nicht
# als Trainingsdaten. Es fliesst nirgends in Bewertung oder Lernen ein.
#
# Erwartungsklassen (grobe chemische Einschaetzung, Gedaechtnis-Wissen -
# bitte bei Zweifel korrigieren, das ist die ganze Idee der Liste):
#
#   "stabil"    - isolierbar, alltaeglich (Erwartung: Kritiker >= 75)
#   "grenzfall" - existiert, aber reaktiv/gespannt/metastabil (50-89)
#   "instabil"  - kurzlebig, explosiv oder nur hypothetisch (Kritiker < 50)
#
# Moleküle werden aus SMILES gebaut (RDKit, nur fuer den Bau) und ins
# interne Format mit expliziten H-Atomen uebersetzt.

from rdkit import Chem

# Valenzfeld wie im Generator (P/S mit Maximalvalenz, siehe kritiker.py)
VALENZ = {
    "H": 1, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1, "Si": 4,
    "P": 5, "S": 6, "Cl": 1, "Br": 1, "I": 1,
}

# (Name, SMILES, Erwartung)
REFERENZ = [
    # --- stabil: Grundbausteine ---
    ("Wasser", "O", "stabil"),
    ("Methan", "C", "stabil"),
    ("Ethan", "CC", "stabil"),
    ("Ethen", "C=C", "stabil"),
    ("Ethin", "C#C", "stabil"),
    ("Ammoniak", "N", "stabil"),
    ("Stickstoff N2", "N#N", "stabil"),
    ("Sauerstoff O2", "O=O", "stabil"),
    ("Wasserstoff H2", "[H][H]", "stabil"),
    ("Kohlendioxid", "O=C=O", "stabil"),
    ("Blausaeure", "C#N", "stabil"),
    ("Fluorwasserstoff", "F", "stabil"),
    ("Chlorwasserstoff", "Cl", "stabil"),
    ("Schwefelwasserstoff", "S", "stabil"),
    ("Phosphin", "P", "stabil"),
    ("Silan", "[SiH4]", "stabil"),
    # --- stabil: organisch ---
    ("Formaldehyd", "C=O", "stabil"),
    ("Methanol", "CO", "stabil"),
    ("Ethanol", "CCO", "stabil"),
    ("Essigsaeure", "CC(=O)O", "stabil"),
    ("Ameisensaeure", "C(=O)O", "stabil"),
    ("Aceton", "CC(C)=O", "stabil"),
    ("Acetaldehyd", "CC=O", "stabil"),
    ("Dimethylether", "COC", "stabil"),
    ("Methylamin", "CN", "stabil"),
    ("Harnstoff", "NC(N)=O", "stabil"),
    ("Glycin", "NCC(=O)O", "stabil"),
    ("Acetonitril", "CC#N", "stabil"),
    ("Ethylenglykol", "OCCO", "stabil"),
    ("Propen", "CC=C", "stabil"),
    ("Butadien", "C=CC=C", "stabil"),
    ("Allen", "C=C=C", "stabil"),
    ("Cyclohexan", "C1CCCCC1", "stabil"),
    ("Cyclopentan", "C1CCCC1", "stabil"),
    ("Benzol", "c1ccccc1", "stabil"),
    ("Toluol", "Cc1ccccc1", "stabil"),
    ("Phenol", "Oc1ccccc1", "stabil"),
    ("Anilin", "Nc1ccccc1", "stabil"),
    ("Pyridin", "c1ccncc1", "stabil"),
    ("Furan", "c1ccoc1", "stabil"),
    ("Thiophen", "c1ccsc1", "stabil"),
    ("Naphthalin", "c1ccc2ccccc2c1", "stabil"),
    # --- stabil: Halogene / Hauptgruppen ---
    ("Chlormethan", "CCl", "stabil"),
    ("Chloroform", "ClC(Cl)Cl", "stabil"),
    ("Tetrafluormethan", "FC(F)(F)F", "stabil"),
    ("Stickstofftrifluorid", "FN(F)F", "stabil"),
    ("Bortrifluorid", "FB(F)F", "stabil"),
    ("Schwefeldioxid", "O=S=O", "stabil"),
    ("DMSO", "CS(C)=O", "stabil"),
    ("Schwefelhexafluorid", "FS(F)(F)(F)(F)F", "stabil"),
    ("Phosphorsaeure", "OP(O)(O)=O", "stabil"),
    ("Schwefelsaeure", "OS(O)(=O)=O", "stabil"),
    # --- stabil: Schwefel/Phosphor/N-Heterocyclen (Regressionsschutz) ---
    ("Dimethylsulfon", "CS(C)(=O)=O", "stabil"),
    ("Thionylchlorid", "ClS(Cl)=O", "stabil"),
    ("Trimethylphosphinoxid", "CP(C)(C)=O", "stabil"),
    ("Dimethyldisulfid", "CSSC", "stabil"),
    ("Dimethyltrisulfid", "CSSSC", "stabil"),
    ("Schwefel S8", "S1SSSSSSS1", "stabil"),
    ("Pyrazol", "c1cn[nH]c1", "stabil"),
    ("Imidazol", "c1c[nH]cn1", "stabil"),
    ("1,2,3-Triazol", "c1cn[nH]n1", "stabil"),
    ("Tetrazol", "c1nnn[nH]1", "stabil"),
    ("Azobenzol", "c1ccccc1N=Nc1ccccc1", "stabil"),
    # --- Grenzfaelle ---
    ("Wasserstoffperoxid", "OO", "grenzfall"),
    ("Hydroxylamin", "NO", "grenzfall"),
    ("Hydrazin", "NN", "grenzfall"),
    ("Cyclopropan", "C1CC1", "grenzfall"),
    ("Cyclopropen", "C1=CC1", "grenzfall"),
    ("Cyclobutan", "C1CCC1", "grenzfall"),
    ("Bicyclobutan", "C1C2CC12", "grenzfall"),
    ("Cuban", "C12C3C4C1C5C2C3C45", "grenzfall"),
    ("Diacetylen", "C#CC#C", "grenzfall"),
    ("Dimethylperoxid", "COOC", "grenzfall"),
    ("Fluor F2", "FF", "grenzfall"),
    ("Sauerstoffdifluorid", "FOF", "grenzfall"),
    ("Boran BH3", "[BH3]", "grenzfall"),
    ("Stickstofftrichlorid", "ClN(Cl)Cl", "grenzfall"),
    ("Oxiran", "C1CO1", "grenzfall"),
    ("Dichlormonoxid", "ClOCl", "grenzfall"),
    ("Butatrien", "C=C=C=C", "grenzfall"),
    ("Kohlensuboxid", "O=C=C=C=O", "grenzfall"),
    # --- instabil / kurzlebig / hypothetisch ---
    ("Methylen (Carben)", "[CH2]", "instabil"),
    ("Methylradikal", "[CH3]", "instabil"),
    ("Hydroxylradikal", "[OH]", "instabil"),
    ("Methin CH", "[CH]", "instabil"),
    ("Tetrahedran", "C12C3C1C23", "instabil"),
    ("Cyclopropin", "C1#CC1", "instabil"),
    ("Cyclobutin", "C1#CCC1", "instabil"),
    ("Tetrazet", "N1=NN=N1", "instabil"),
    ("Tetraazatetrahedran N4", "N12N3N1N23", "instabil"),
    ("Triazan", "NNN", "instabil"),
    ("Tetrazan", "NNNN", "instabil"),
    ("Trioxidan", "OOO", "instabil"),
    ("Tetroxid-Kette", "COOOOC", "instabil"),
    ("FOOF", "FOOF", "instabil"),
    ("Pentatetraen", "C=C=C=C=C", "instabil"),
    ("Dichloroxid-Kette", "ClOOCl", "instabil"),
    ("Cumulen im Dreiring", "C1=C=C1", "instabil"),
    ("Ethylidin-Kette", "[CH]=[CH]", "instabil"),
    ("Schwefeltetrahydrid", "[SH4]", "instabil"),
    ("Phosphorpentahydrid", "[PH5]", "instabil"),
    # --- Zusatz: Grenzfaelle und Transiente mit Literaturhinweis ---
    ("Keten", "C=C=O", "grenzfall"),
    ("Diazirin", "C1N=N1", "grenzfall"),
    ("Aziridin", "C1CN1", "grenzfall"),
    ("Peressigsaeure", "CC(=O)OO", "grenzfall"),
    ("Borazin", "B1NBNBN1", "stabil"),
    ("Cyclobutadien", "C1=CC=C1", "instabil"),
    ("Thioformaldehyd", "C=S", "instabil"),
    ("Silanon", "[SiH2]=O", "instabil"),
    ("Disilen", "[SiH2]=[SiH2]", "instabil"),
    ("Silen", "C=[SiH2]", "instabil"),
    ("Phosphaethen", "C=[PH]", "instabil"),
    ("Nitren NH", "[NH]", "instabil"),
    ("Thioaceton", "CC(C)=S", "instabil"),
    ("Thioharnstoff", "NC(N)=S", "stabil"),
    ("Trimethylphosphinsulfid", "CP(C)(C)=S", "stabil"),
    ("Methylentrimethylphosphoran", "C=P(C)(C)C", "grenzfall"),
]

# Kurze Literaturhinweise zur Einteilung - AUS DEM GEDAECHTNIS NIEDERGESCHRIEBEN,
# bitte vor jeder Veroeffentlichung gegen die Originalquelle pruefen
# (z.B. NIST Chemistry WebBook, CRC Handbook, Original-Syntheseliteratur).
HINWEIS = {
    "Wasserstoffperoxid": "isolierbar, Loesungen kaeuflich; zerfaellt exotherm (O-O schwach)",
    "Hydroxylamin": "isolierbar, aber thermisch instabil und explosionsfaehig",
    "Hydrazin": "stabile Fluessigkeit, energiereich (endotherme Bildung), Raketentreibstoff",
    "Cyclopropan": "isolierbar; Ringspannung ca. 27,5 kcal/mol",
    "Cyclopropen": "isolierbar, sehr reaktiv; Doppelbindung im Dreiring",
    "Bicyclobutan": "isolierbar, hohe Spannung (ca. 66 kcal/mol)",
    "Cuban": "1964 synthetisiert (Eaton); kinetisch stabil trotz ca. 166 kcal/mol Spannung",
    "Tetrahedran": "Stammverbindung unbekannt; nur Tetra-tert-butyl-Derivat isoliert (Maier 1978)",
    "Diacetylen": "stabil bei tiefer Temperatur, polymerisiert/explosiv bei Zimmertemperatur",
    "Dimethylperoxid": "isolierbar, aber explosiv",
    "Fluor F2": "stabiles Molekuel, extrem reaktiv",
    "Sauerstoffdifluorid": "isolierbares Gas, starkes Oxidationsmittel",
    "Boran BH3": "Monomer nur transient; dimerisiert zu Diboran (B2H6, 3-Zentren-Bindungen)",
    "Stickstofftrichlorid": "isolierbar, explosiv (endotherm)",
    "Dichlormonoxid": "isolierbar, explosiv zerfallend",
    "Butatrien": "bei tiefer Temperatur handhabbar, polymerisiert",
    "Kohlensuboxid": "isolierbar, polymerisiert bei Zimmertemperatur",
    "Trioxidan": "1990er/2005 bei tiefer Temperatur nachgewiesen, zersetzt sich schnell",
    "FOOF": "Streng 1962; zersetzt sich oberhalb ca. -100 Grad C",
    "Tetrazet": "nicht isoliert, nur theoretisch",
    "Tetraazatetrahedran N4": "nicht isoliert (Tetraeder-N4 nur Rechnungen)",
    "Cyclopropin": "nicht beobachtet",
    "Cyclobutadien": "antiaromatisch, nur in Matrix bei tiefer Temperatur nachgewiesen",
    "Keten": "isolierbar, dimerisiert bei Zimmertemperatur",
    "Diazirin": "existiert, explosionsgefaehrlich; substituierte Diazirine als Photomarker",
    "Aziridin": "isolierbare Fluessigkeit, giftig, Ringspannung ca. 27 kcal/mol",
    "Peressigsaeure": "kaeufliche Loesungen; konzentriert explosiv",
    "Borazin": "isolierbar, 'anorganisches Benzol'",
    "Thioformaldehyd": "nur transient (trimerisiert), im interstellaren Raum nachgewiesen",
    "Silanon": "Stammverbindung nur transient/Matrix",
    "Disilen": "Stammverbindung transient; stabil nur mit sperrigen Resten (West 1981)",
    "Silen": "Stammverbindung transient",
    "Phosphaethen": "Stammverbindung nur bei tiefer Temperatur/transient",
    "Nitren NH": "Radikal, nur transient",
    "Thioaceton": "nur unterhalb ca. -20 Grad C als Monomer; trimerisiert, beruechtigter Geruch",
    "Thioharnstoff": "stabile kristalline Verbindung",
    "Trimethylphosphinsulfid": "stabile kristalline Verbindung",
    "Methylentrimethylphosphoran": "isolierbares Phosphorylid, luft- und feuchtigkeitsempfindlich",
}

# ----------------------------------------------------------
# Fehlerhafte Graphen, wie der Generator sie zufaellig erzeugen kann.
# RDKit lehnt sie ab; sie testen die Pruefmodelle des Kritikers (Valenz,
# Oktett, Halogene, Fragmente ...). Erwartung: alle "instabil".
# ----------------------------------------------------------


def _graph(atome, bindungen):
    """atome: Liste von Elementen; bindungen: Liste (i, j, ordnung)."""

    daten = [
        {"id": i, "element": e, "valence": VALENZ[e], "used_valence": 0}
        for i, e in enumerate(atome)
    ]

    for i, j, ordnung in bindungen:
        daten[i]["used_valence"] += ordnung
        daten[j]["used_valence"] += ordnung

    return {
        "atoms": daten,
        "bonds": [{"atom1": i, "atom2": j, "order": o} for i, j, o in bindungen],
    }


FEHLERHAFT = [
    ("[Fehler] CH5 (5 Bindungen am C)", _graph(["C"] + ["H"] * 5, [(0, i, 1) for i in range(1, 6)])),
    ("[Fehler] NH4 neutral", _graph(["N"] + ["H"] * 4, [(0, i, 1) for i in range(1, 5)])),
    ("[Fehler] NH5", _graph(["N"] + ["H"] * 5, [(0, i, 1) for i in range(1, 6)])),
    ("[Fehler] H3O neutral", _graph(["O"] + ["H"] * 3, [(0, i, 1) for i in range(1, 4)])),
    ("[Fehler] H2F (F mit 2 Bindungen)", _graph(["F", "H", "H"], [(0, 1, 1), (0, 2, 1)])),
    ("[Fehler] H-H-H Kette", _graph(["H", "H", "H"], [(0, 1, 1), (1, 2, 1)])),
    ("[Fehler] H=O (Doppelbindung am H)", _graph(["H", "O", "H"], [(0, 1, 2), (1, 2, 1)])),
    ("[Fehler] CH3-Cl-CH3", _graph(
        ["Cl", "C", "C"] + ["H"] * 6,
        [(0, 1, 1), (0, 2, 1)] + [(1, 3, 1), (1, 4, 1), (1, 5, 1), (2, 6, 1), (2, 7, 1), (2, 8, 1)])),
    ("[Fehler] CH4 + H2O getrennt", _graph(
        ["C", "H", "H", "H", "H", "O", "H", "H"],
        [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (5, 6, 1), (5, 7, 1)])),
    ("[Fehler] einzelnes C-Atom", _graph(["C"], [])),
    ("[Fehler] einzelnes O-Atom", _graph(["O"], [])),
    ("[Fehler] C-C ohne H", _graph(["C", "C"], [(0, 1, 1)])),
    ("[Fehler] O mit 4 Bindungen", _graph(["O"] + ["H"] * 4, [(0, i, 1) for i in range(1, 5)])),
]


def aus_smiles(smiles):
    """SMILES -> Molekuel im internen Format (explizite H-Atome,
    Kekule-Form, used_valence = Summe der Bindungsordnungen)."""

    m = Chem.MolFromSmiles(smiles)

    if m is None:
        raise ValueError(f"SMILES nicht lesbar: {smiles}")

    Chem.Kekulize(m, clearAromaticFlags=True)
    m = Chem.AddHs(m)

    atome = [
        {"id": a.GetIdx(), "element": a.GetSymbol(),
         "valence": VALENZ[a.GetSymbol()], "used_valence": 0}
        for a in m.GetAtoms()
    ]

    bindungen = []

    for b in m.GetBonds():

        ordnung = int(round(b.GetBondTypeAsDouble()))
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()

        bindungen.append({"atom1": i, "atom2": j, "order": ordnung})
        atome[i]["used_valence"] += ordnung
        atome[j]["used_valence"] += ordnung

    return {"atoms": atome, "bonds": bindungen}


def referenz_molekuele():
    """Liste (name, smiles, erwartung, molekuel). Fehlerhafte Graphen haben
    als 'smiles' den Text '(Graph)'."""

    liste = [(n, s, e, aus_smiles(s)) for n, s, e in REFERENZ]
    liste += [(n, "(Graph)", "instabil", m) for n, m in FEHLERHAFT]

    return liste


def hinweis(name):
    """Literaturhinweis zur Einteilung oder ''."""

    return HINWEIS.get(name, "")