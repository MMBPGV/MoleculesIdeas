import random
from collections import Counter, defaultdict

from chemlabor.ki.evolution import Evolution, _zusammenhaengend

VALENZEN = {"H": 1, "C": 4, "N": 3, "O": 2, "F": 1, "Cl": 1, "Si": 4, "S": 2}
BINDUNGSARTEN = {"H": [1], "C": [1, 2, 3], "N": [1, 2, 3], "O": [1, 2], "F": [1],
                 "Cl": [1], "Si": [1, 2], "S": [1, 2]}
ELEMENTE = list(VALENZEN)


def signatur(mol):
    nachbarn = defaultdict(list)
    el = {a["id"]: a["element"] for a in mol["atoms"]}
    for b in mol["bonds"]:
        nachbarn[b["atom1"]].append((el[b["atom2"]], b["order"]))
        nachbarn[b["atom2"]].append((el[b["atom1"]], b["order"]))
    return tuple(sorted((a["element"], tuple(sorted(nachbarn[a["id"]]))) for a in mol["atoms"]))


def zufallsmolekuel(rng, max_atome=10):
    """Einfache Kette aus Zufallselementen (Einfachbindungen)."""
    n = rng.randint(2, max_atome)
    atome = [{"id": i, "element": rng.choice(ELEMENTE), "valence": 0, "used_valence": 0} for i in range(n)]
    for a in atome:
        a["valence"] = VALENZEN[a["element"]]
    bindungen = []
    for i in range(1, n):
        kandidaten = [a for a in atome[:i] if a["valence"] - a["used_valence"] >= 1]
        if not kandidaten or atome[i]["valence"] < 1:
            break
        ziel = rng.choice(kandidaten)
        bindungen.append({"atom1": atome[i]["id"], "atom2": ziel["id"], "order": 1})
        atome[i]["used_valence"] += 1
        ziel["used_valence"] += 1
        if len([b for b in bindungen]) == i and False:
            pass
    verbunden = {0} | {b["atom1"] for b in bindungen} | {b["atom2"] for b in bindungen}
    atome = [a for a in atome if a["id"] in verbunden]
    return {"atoms": atome, "bonds": bindungen}


def pruefe_gueltig(mol, max_atome):
    ids = [a["id"] for a in mol["atoms"]]
    assert len(ids) == len(set(ids)), "doppelte IDs"
    assert 2 <= len(ids) <= max_atome, f"Atomzahl {len(ids)}"
    benutzt = Counter()
    for b in mol["bonds"]:
        assert b["atom1"] in ids and b["atom2"] in ids and b["atom1"] != b["atom2"]
        assert b["order"] in BINDUNGSARTEN[next(a["element"] for a in mol["atoms"] if a["id"] == b["atom1"])]
        assert b["order"] in BINDUNGSARTEN[next(a["element"] for a in mol["atoms"] if a["id"] == b["atom2"])]
        benutzt[b["atom1"]] += b["order"]
        benutzt[b["atom2"]] += b["order"]
    for a in mol["atoms"]:
        assert a["valence"] == VALENZEN[a["element"]]
        assert a["used_valence"] == benutzt[a["id"]], f"used_valence falsch bei {a}"
        assert a["used_valence"] <= a["valence"], f"Valenz ueberschritten bei {a}"
    assert _zusammenhaengend(mol), "nicht zusammenhaengend"
    paare = [frozenset((b["atom1"], b["atom2"])) for b in mol["bonds"]]
    assert len(paare) == len(set(paare)), "Doppelbindung zwischen denselben Atomen"


def bewerte(mol):
    """Testbewertung: viele C, wenige freie Valenzen, Ringe belohnt."""
    c = sum(1 for a in mol["atoms"] if a["element"] == "C")
    frei = sum(a["valence"] - a["used_valence"] for a in mol["atoms"])
    ringe = len(mol["bonds"]) - len(mol["atoms"]) + 1
    return max(0, min(100, 10 * c - 2 * frei + 5 * ringe))


def baue(seed=1, **kwargs):
    rng = random.Random(seed)
    quelle = lambda: zufallsmolekuel(rng)
    return Evolution(ELEMENTE, VALENZEN, BINDUNGSARTEN, [("zufall", quelle, 1.0)], signatur,
                     max_atome=12, rng=random.Random(seed + 1), **kwargs)


def lauf(evo, schritte):
    for _ in range(schritte):
        m = evo.naechstes()
        evo.melde(m, bewerte(m))


def test_mutationen_und_rekombination_erzeugen_gueltige_molekuele():
    evo = baue()
    rng = random.Random(5)
    anzahl = 0
    for _ in range(400):
        m = zufallsmolekuel(rng)
        pruefe_gueltig(m, 12)
        kind = evo.mutiere(m)
        if kind is not None:
            pruefe_gueltig(kind, 12)
            anzahl += 1
        k2 = evo.rekombiniere(m, zufallsmolekuel(rng))
        if k2 is not None:
            pruefe_gueltig(k2, 12)
    assert anzahl > 300


def test_ringoperationen():
    evo = baue()
    rng = random.Random(8)
    m = zufallsmolekuel(rng, 8)
    for _ in range(60):
        n = evo.mutiere(m) or m
        pruefe_gueltig(n, 12)
        m = n


def test_population_verbessert_sich_und_bleibt_gueltig():
    evo = baue(population=20, eliten=3)
    lauf(evo, 20)
    anfang = evo.bestes()[1]
    lauf(evo, 20 * 25)
    ende = evo.bestes()[1]
    assert ende > anfang, (anfang, ende)
    for e in evo.population:
        pruefe_gueltig(e["mol"], 12)
    sigs = [e["sig"] for e in evo.population]
    assert len(sigs) == len(set(sigs)), "doppelte Strukturen in der Population"
    assert len(evo.population) <= 20


def test_sackgasse_loest_einwanderer_aus_und_schutz_wirkt():
    ereignisse = []
    evo = baue(population=16, eliten=2, sackgasse_generationen=3, zufall_grundrate=0.0,
               on_generation=lambda s: ereignisse.append(s))
    # konstante Bewertung -> nie eine Verbesserung -> Sackgasse
    for _ in range(16 * 12):
        m = evo.naechstes()
        evo.melde(m, 50)
    arten = [e["ereignis"] for e in ereignisse]
    assert "sackgasse" in arten, arten
    assert "neustart" in arten, arten
    # in Sackgasse-Generationen kommen Einwanderer
    assert any(s["ereignis"] == "sackgasse" and s["einwanderer"] > 0 for s in ereignisse)


def test_einwanderer_ueberleben_dank_schutz():
    evo = baue(population=12, eliten=2, schutz_generationen=3, zufall_grundrate=0.5)
    evo.besiedle([(zufallsmolekuel(random.Random(i)), 90) for i in range(12)])
    # Einwanderer bekommen schlechte Punkte, Elite hat 90
    for _ in range(12 * 2):
        m = evo.naechstes()
        evo.melde(m, 5)
    schutz_und_schlecht = [e for e in evo.population if e["punkte"] <= 5]
    assert schutz_und_schlecht, "geschuetzte Einwanderer wurden sofort aussortiert"


def test_nicht_gemeldete_kandidaten_zaehlen_nicht():
    evo = baue(population=10, eliten=2)
    for _ in range(40):
        evo.naechstes()           # nie melden (z.B. Bewertung schlug fehl)
    # darf nicht abstuerzen, Population bleibt leer oder klein
    assert len(evo.population) == 0 or all(e["punkte"] is not None for e in evo.population)


def test_besiedle_und_bestes():
    evo = baue(population=5)
    rng = random.Random(3)
    evo.besiedle([(zufallsmolekuel(rng), p) for p in (10, 70, 40)])
    assert evo.bestes()[1] == 70


def test_plateau_keine_dauerhafte_sackgasse_und_umwaelzung():
    """Bei konstant bester Punktzahl (Deckel) muss die Population
    weiterwandern: Ruhmeshalle waechst, mehrere Nischen erreichen das
    Plateau."""

    rng = random.Random(3)

    evo = Evolution(
        ELEMENTE, VALENZEN, BINDUNGSARTEN,
        quellen=[("zufall", lambda: zufallsmolekuel(rng), 1.0)],
        signatur=signatur, population=20, eliten=3, rng=rng,
        plateau=95, max_alter=3,
    )

    for _ in range(600):
        mol = evo.naechstes()
        evo.melde(mol, 100.0)

    assert len(evo.ruhmeshalle) > 5
    assert len(evo._top_nischen) > 10
    assert all(v <= evo.nische_max for v in Counter(
        __import__("chemlabor.ki.evolution", fromlist=["_brutto"])._brutto(e["mol"])
        for e in evo.population).values())