"""Eicht den Kritiker an bekannten Molekuelen (chemlabor/chemie/referenzset.py).

Aufruf im Projektordner (dort, wo chemlabor/ liegt):

    python eichung.py            # Bericht in der Konsole + eichung.csv
    python eichung.py --kurz     # nur Zusammenfassung und Ausreisser
    python eichung.py --kraftfeld  # Kritiker MIT Kraftfeld-Abzug messen

Ausgabe:
  1. Tabelle: Molekuel, Erwartung, Kritiker-Punkte, Urteil, Treffer?
  2. Je Erwartungsklasse: Mittelwert, Spanne, Trefferquote
  3. Trennschaerfe stabil vs. instabil (AUC, 1.0 = perfekt, 0.5 = Zufall)
  4. Welches der 16 Teilmodelle zieht wie oft / wie stark Punkte ab - Modelle,
     die nie etwas abziehen, tragen bei diesen Molekuelen nichts bei.
  5. Ausreisser mit den Kritiker-Begruendungen

Keine Aenderung an Bewertung oder Gedaechtnis - reine Messung.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(".").resolve()))

from chemlabor.chemie.kritiker import Kritiker
from chemlabor.chemie.kritiker_kraftfeld import KritikerKraftfeld
from daten.referenzset import referenz_molekuele

# Erwartete Punktbereiche je Klasse
BEREICH = {"stabil": (75, 100), "grenzfall": (50, 89), "instabil": (0, 49)}


def modellbeitraege(kritiker, molekuel):
    """Punktabzug je Teilmodell (dieselbe Reihenfolge wie im Kritiker)."""

    daten = kritiker.analysiere(molekuel)

    # Die Modelle sind rein lesend auf `daten` - Wiederholung einzeln.
    beitraege = {}

    for modell in kritiker.modelle:
        beitraege[modell.__name__] = modell(daten)["punkte"]

    return daten, beitraege


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


def main():

    kurz = "--kurz" in sys.argv
    # --kraftfeld: Kritiker mit Kraftfeld-Abzug messen
    kritiker = KritikerKraftfeld() if "--kraftfeld" in sys.argv else Kritiker()

    zeilen = []
    je_klasse = defaultdict(list)
    modell_summe = defaultdict(float)
    modell_anzahl = defaultdict(int)

    for name, smiles, erwartung, molekuel in referenz_molekuele():

        daten, beitraege = modellbeitraege(kritiker, molekuel)

        punkte = daten["bewertung"]["punkte"]
        urteil = daten["bewertung"]["urteil"]
        unten, oben = BEREICH[erwartung]
        treffer = unten <= punkte <= oben

        je_klasse[erwartung].append((punkte, treffer))

        for modell, wert in beitraege.items():
            if wert:
                modell_summe[modell] += wert
                modell_anzahl[modell] += 1

        zeilen.append({
            "name": name, "smiles": smiles, "erwartung": erwartung,
            "punkte": punkte, "urteil": urteil, "treffer": treffer,
            "probleme": " | ".join(_kurz(p) for p in daten["probleme"][:3]),
        })

    # 1. Tabelle
    if not kurz:
        print(f"{'Molekuel':26} {'Erwartung':10} {'Punkte':>6}  {'Urteil':22} Treffer")
        print("-" * 78)
        for z in zeilen:
            print(f"{z['name'][:26]:26} {z['erwartung']:10} {z['punkte']:6}  "
                  f"{z['urteil']:22} {'ok' if z['treffer'] else 'ABWEICHUNG'}")

    # 2. Klassen
    print("\nJe Erwartungsklasse:")
    for klasse in ("stabil", "grenzfall", "instabil"):
        werte = je_klasse.get(klasse, [])
        if not werte:
            continue
        p = [w[0] for w in werte]
        quote = sum(1 for w in werte if w[1]) / len(werte)
        print(f"  {klasse:10} n={len(p):3}  Mittel {sum(p)/len(p):5.1f}  "
              f"Spanne {min(p)}-{max(p)}  Trefferquote {quote:4.0%}")

    # 3. Trennschaerfe
    a = auc([w[0] for w in je_klasse["stabil"]], [w[0] for w in je_klasse["instabil"]])
    if a is not None:
        print(f"\nTrennschaerfe stabil vs. instabil (AUC): {a:.2f}")

    # 4. Modelle
    print("\nTeilmodelle (wie oft ziehen sie bei diesem Set Punkte ab):")
    for modell in (m.__name__ for m in kritiker.modelle):
        n = modell_anzahl.get(modell, 0)
        s = modell_summe.get(modell, 0.0)
        print(f"  {modell:36} {n:3}x  Summe {s:7.0f}")

    # 5. Ausreisser
    ausreisser = [z for z in zeilen if not z["treffer"]]
    print(f"\nAusreisser: {len(ausreisser)} von {len(zeilen)}")
    for z in ausreisser:
        print(f"  - {z['name']} ({z['smiles']}): erwartet {z['erwartung']}, "
              f"Kritiker {z['punkte']} ({z['urteil']})")
        if z["probleme"]:
            print(f"      {z['probleme'][:160]}")

    with open("eichung.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(zeilen[0]))
        w.writeheader()
        w.writerows(zeilen)

    print("\nGespeichert: eichung.csv")


if __name__ == "__main__":
    main()