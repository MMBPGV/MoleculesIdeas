"""Eicht den Kritiker an bekannten Molekuelen (chemlabor/chemie/referenzset.py).

Aufruf im Projektordner (dort, wo chemlabor/ liegt):

    python eichung.py                 # Gesamtbericht + eichung.csv
    python eichung.py --kurz          # ohne die lange Molekuel-Tabelle
    python eichung.py --split dev     # nur den Entwicklungsteil bewerten
    python eichung.py --split test    # nur den Testteil (erst am Ende ansehen!)
    python eichung.py --kraftfeld     # Kritiker MIT Kraftfeld-Abzug messen

Ausgabe:
  1. Tabelle: Molekuel, Teil (dev/test), Erwartung, Punkte, Urteil, Treffer
  2. Je Erwartungsklasse: Mittelwert, Spanne, Trefferquote
  3. Trennschaerfe stabil vs. instabil (AUC, 1.0 = perfekt, 0.5 = Zufall),
     getrennt fuer dev / test / alle, mit 95%-Bootstrap-Intervall
  4. Weglass-Analyse: Wie aendern sich AUC und Trefferquote, wenn EIN
     Teilmodell ausgeschaltet wird? (Delta < 0 = das Modell hilft;
     Delta > 0 = es schadet bei diesem Set; 0 = tut hier nichts)
  5. Ausreisser mit den Kritiker-Begruendungen

Entwicklungs- und Testteil:
  Jedes Molekuel gehoert anhand seines SMILES-Texts fest zu 'dev' (ca. 70 %)
  oder 'test' (ca. 30 %) - die Zuordnung aendert sich nicht, wenn das Set
  waechst. Regeln und Abzuege darfst du an 'dev' abstimmen; 'test' ist die
  ehrliche Gegenprobe, ob die Aenderung auch bei Molekuelen traegt, die du
  beim Abstimmen nicht angesehen hast. Wer 'test' beim Abstimmen anschaut,
  hat keine Gegenprobe mehr.

Keine Aenderung an Bewertung oder Gedaechtnis - reine Messung.
"""

import csv
import hashlib
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(".").resolve()))

from chemlabor.chemie.kritiker import Kritiker
from chemlabor.chemie.kritiker_kraftfeld import KritikerKraftfeld
from daten.referenzset import hinweis, referenz_molekuele

# Erwartete Punktbereiche je Klasse
BEREICH = {"stabil": (75, 100), "grenzfall": (50, 89), "instabil": (0, 49)}

TEST_ANTEIL = 0.30
BOOTSTRAP = 1000


def teil(schluessel):
    """'dev' oder 'test' - fest, nur vom SMILES-Text (bei Graphen: vom
    Namen) abhaengig."""

    wert = int(hashlib.md5(schluessel.encode("utf-8")).hexdigest(), 16) % 1000

    return "test" if wert < TEST_ANTEIL * 1000 else "dev"


def begrenze(punkte):

    return max(0, min(100, punkte))


def _kurz(problem):
    """Ein Problem-Eintrag (dict oder Text) als kurze Zeile."""

    if isinstance(problem, dict):
        rest = ", ".join(f"{k}={v}" for k, v in problem.items() if k != "modell")
        return f"{problem.get('modell', '?')}: {rest}"

    return str(problem)


def auc(positiv, negativ):
    """Wahrscheinlichkeit, dass ein zufaelliges 'stabil' mehr Punkte hat
    als ein zufaelliges 'instabil' (Ties zaehlen halb)."""

    if not positiv or not negativ:
        return None

    gewonnen = sum(
        1.0 if p > n else 0.5 if p == n else 0.0
        for p in positiv for n in negativ
    )

    return gewonnen / (len(positiv) * len(negativ))


def auc_intervall(positiv, negativ, rng):
    """95%-Bootstrap-Intervall der AUC (Molekuele werden je Klasse mit
    Zuruecklegen neu gezogen). Bei wenigen Molekuelen breit - und genau
    das soll man sehen."""

    if len(positiv) < 2 or len(negativ) < 2:
        return None

    werte = []

    for _ in range(BOOTSTRAP):
        p = [rng.choice(positiv) for _ in positiv]
        n = [rng.choice(negativ) for _ in negativ]
        werte.append(auc(p, n))

    werte.sort()

    return werte[int(0.025 * BOOTSTRAP)], werte[int(0.975 * BOOTSTRAP)]


def trefferquote(zeilen, punkte_von):
    """Anteil der Molekuele, deren Punkte im Bereich ihrer Klasse liegen."""

    if not zeilen:
        return None

    treffer = sum(
        1 for z in zeilen
        if BEREICH[z["erwartung"]][0] <= punkte_von(z) <= BEREICH[z["erwartung"]][1]
    )

    return treffer / len(zeilen)


def messzahlen(zeilen, punkte_von):
    """(AUC stabil vs. instabil, Trefferquote) fuer eine Zeilenliste."""

    positiv = [punkte_von(z) for z in zeilen if z["erwartung"] == "stabil"]
    negativ = [punkte_von(z) for z in zeilen if z["erwartung"] == "instabil"]

    return auc(positiv, negativ), trefferquote(zeilen, punkte_von)


def formatiere(wert, prozent=False):

    if wert is None:
        return "  –  "

    return f"{wert:5.0%}" if prozent else f"{wert:5.2f}"


def main():

    argumente = sys.argv[1:]
    kurz = "--kurz" in argumente

    auswahl = "alle"

    if "--split" in argumente:
        position = argumente.index("--split") + 1
        auswahl = argumente[position] if position < len(argumente) else "alle"

    if auswahl not in ("dev", "test", "alle"):
        print("--split erwartet dev, test oder alle")
        return 1

    kritiker = KritikerKraftfeld() if "--kraftfeld" in argumente else Kritiker()

    zeilen = []
    abweichungen = 0

    for name, smiles, erwartung, molekuel in referenz_molekuele():

        schluessel = name if smiles == "(Graph)" else smiles

        if auswahl != "alle" and teil(schluessel) != auswahl:
            continue

        daten = kritiker.analysiere(molekuel)

        beitraege = {m.__name__: m(daten)["punkte"] for m in kritiker.modelle}

        endpunkte = daten["bewertung"]["punkte"]
        summe = begrenze(100 + sum(beitraege.values()))

        # Aufschlag/Abzug ausserhalb der Teilmodelle (z.B. Kraftfeld)
        zusatz = endpunkte - summe

        # Pruefung: stimmt die Rekonstruktion aus den Teilmodellen?
        if zusatz != 0 and "--kraftfeld" not in argumente:
            abweichungen += 1

        ohne = {
            modell: begrenze(begrenze(100 + sum(beitraege.values()) - wert) + zusatz)
            for modell, wert in beitraege.items()
        }

        unten, oben = BEREICH[erwartung]

        zeilen.append({
            "name": name, "smiles": smiles, "teil": teil(schluessel),
            "hinweis": hinweis(name),
            "erwartung": erwartung, "punkte": endpunkte,
            "urteil": daten["bewertung"]["urteil"],
            "treffer": unten <= endpunkte <= oben,
            "beitraege": beitraege, "ohne": ohne,
            "probleme": " | ".join(_kurz(p) for p in daten["probleme"][:3]),
        })

    if not zeilen:
        print("Keine Molekuele in dieser Auswahl.")
        return 1

    # 1. Tabelle
    if not kurz:
        print(f"{'Molekuel':26} {'Teil':4} {'Erwartung':10} {'Punkte':>6}  {'Urteil':22} Treffer")
        print("-" * 84)
        for z in zeilen:
            print(f"{z['name'][:26]:26} {z['teil']:4} {z['erwartung']:10} {z['punkte']:6}  "
                  f"{z['urteil']:22} {'ok' if z['treffer'] else 'ABWEICHUNG'}")

    # 2. Klassen
    print("\nJe Erwartungsklasse:")
    for klasse in ("stabil", "grenzfall", "instabil"):
        gruppe = [z for z in zeilen if z["erwartung"] == klasse]
        if not gruppe:
            continue
        p = [z["punkte"] for z in gruppe]
        quote = sum(1 for z in gruppe if z["treffer"]) / len(gruppe)
        print(f"  {klasse:10} n={len(p):3}  Mittel {sum(p)/len(p):5.1f}  "
              f"Spanne {min(p)}-{max(p)}  Trefferquote {quote:4.0%}")

    # 3. Trennschaerfe je Teil
    print("\nTrennschaerfe stabil vs. instabil (AUC, 95%-Intervall durch Bootstrap):")

    rng = random.Random(12345)

    for name_teil in ("dev", "test", "alle"):

        gruppe = zeilen if name_teil == "alle" else [z for z in zeilen if z["teil"] == name_teil]

        positiv = [z["punkte"] for z in gruppe if z["erwartung"] == "stabil"]
        negativ = [z["punkte"] for z in gruppe if z["erwartung"] == "instabil"]

        wert = auc(positiv, negativ)

        if wert is None:
            continue

        spanne = auc_intervall(positiv, negativ, rng)
        text = f"[{spanne[0]:.2f} – {spanne[1]:.2f}]" if spanne else ""

        print(f"  {name_teil:5} n(stabil)={len(positiv):3} n(instabil)={len(negativ):3}  "
              f"AUC {wert:.2f} {text}")

    # 4. Weglass-Analyse
    voll_auc, voll_quote = messzahlen(zeilen, lambda z: z["punkte"])

    print("\nWeglass-Analyse (ein Teilmodell aus, Rest unveraendert):")
    print(f"  {'Teilmodell':34} {'zieht ab':>8} {'AUC':>6} {'dAUC':>6} {'Treffer':>8} {'dTreffer':>8}")
    print(f"  {'(alle Modelle)':34} {'':>8} {formatiere(voll_auc)} {'':>6} {formatiere(voll_quote, True)}")

    for modell in zeilen[0]["beitraege"]:

        wirkt = sum(1 for z in zeilen if z["beitraege"][modell])

        a, q = messzahlen(zeilen, lambda z, m=modell: z["ohne"][m])

        d_auc = None if a is None or voll_auc is None else a - voll_auc
        d_q = None if q is None or voll_quote is None else q - voll_quote

        print(f"  {modell:34} {wirkt:7}x {formatiere(a)} "
              f"{'' if d_auc is None else f'{d_auc:+6.2f}'} {formatiere(q, True)} "
              f"{'' if d_q is None else f'{d_q:+8.0%}'}")

    print("  (dAUC/dTreffer < 0: ohne das Modell wird es schlechter = es hilft)")

    # 5. Ausreisser
    ausreisser = [z for z in zeilen if not z["treffer"]]
    print(f"\nAusreisser: {len(ausreisser)} von {len(zeilen)}")
    for z in ausreisser:
        print(f"  - {z['name']} ({z['smiles']}) [{z['teil']}]: erwartet {z['erwartung']}, "
              f"Kritiker {z['punkte']} ({z['urteil']})")
        if z["probleme"]:
            print(f"      {z['probleme'][:160]}")
        if z["hinweis"]:
            print(f"      Literatur: {z['hinweis']}")

    if abweichungen:
        print(f"\nHinweis: bei {abweichungen} Molekuelen stimmt die Summe der "
              f"Teilmodelle nicht mit der Kritiker-Punktzahl ueberein - die "
              f"Weglass-Analyse ist dort ungenau.")

    felder = ["name", "smiles", "teil", "erwartung", "punkte", "urteil", "treffer", "probleme", "hinweis"]

    with open("eichung.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=felder, extrasaction="ignore")
        w.writeheader()
        w.writerows(zeilen)

    print("\nGespeichert: eichung.csv")

    return 0


if __name__ == "__main__":
    sys.exit(main())