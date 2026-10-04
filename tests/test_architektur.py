# ==========================================================
# test_architektur.py
# ==========================================================
#
# Prueft die Schichtregeln des Projekts (siehe ARCHITEKTUR.md):
#
#   pfade / config  <-  chemie  <-  ki  <-  labor  <-  1
#
# Eine Schicht darf nur Schichten importieren, die in der Reihe
# links von ihr (oder sie selbst) stehen. Ausserdem: tkinter nur in
# 1/, und alle Importe innerhalb des Projekts sind absolut.
#
# Aufruf:   python tests/test_architektur.py      (oder: pytest)

import ast
import pathlib
import sys

WURZEL = pathlib.Path(__file__).resolve().parent.parent
PAKET = WURZEL / "chemlabor"

BASIS = {"pfade", "config"}

ERLAUBT = {
    "pfade": set(),
    "config": {"pfade"},
    "chemie": BASIS | {"chemie"},
    "ki": BASIS | {"chemie", "ki"},
    "labor": BASIS | {"chemie", "ki", "labor"},
    "1": BASIS | {"chemie", "ki", "labor", "1"},
}

# python -m chemlabor (Einstiegspunkt) darf wie die GUI alles importieren
ERLAUBT["__main__"] = ERLAUBT["1"]


def importe(baum):
    """Liefert (zeile, modulname, ebene) fuer jeden Import im Baum."""

    for knoten in ast.walk(baum):

        if isinstance(knoten, ast.Import):
            for alias in knoten.names:
                yield knoten.lineno, alias.name, 0

        elif isinstance(knoten, ast.ImportFrom):

            if knoten.level:
                yield knoten.lineno, knoten.module or "", knoten.level

            elif knoten.module == "chemlabor":
                for alias in knoten.names:
                    yield knoten.lineno, "chemlabor." + alias.name, 0

            else:
                yield knoten.lineno, knoten.module, 0


def verstoesse():

    meldungen = []

    for datei in sorted(PAKET.rglob("*.py")):

        teile = datei.relative_to(PAKET).with_suffix("").parts

        schicht = teile[0]

        if schicht == "__init__":
            continue

        if schicht not in ERLAUBT:
            meldungen.append(
                f"{datei.relative_to(WURZEL)}: unbekannte Schicht "
                f"'{schicht}' - in tests/test_architektur.py (ERLAUBT) eintragen"
            )
            continue

        baum = ast.parse(datei.read_text(encoding="utf-8-sig"))

        for zeile, name, ebene in importe(baum):

            ort = f"{datei.relative_to(WURZEL)}:{zeile}"

            if ebene:
                meldungen.append(f"{ort}: relativer Import - bitte absolut importieren")
                continue

            wurzel_name = name.split(".")[0]

            if wurzel_name == "tkinter" and schicht not in ("1", "__main__"):
                meldungen.append(f"{ort}: tkinter ausserhalb von 1/")

            if wurzel_name == "chemlabor":

                abschnitte = name.split(".")

                ziel = abschnitte[1] if len(abschnitte) > 1 else None

                if ziel is not None and ziel not in ERLAUBT[schicht]:
                    meldungen.append(
                        f"{ort}: Schicht '{schicht}' darf '{ziel}' nicht importieren"
                    )

    return meldungen


def test_schichten():

    assert verstoesse() == []


if __name__ == "__main__":

    probleme = verstoesse()

    for p in probleme:
        print("VERSTOSS:", p)

    print("Architektur ok." if not probleme else f"{len(probleme)} Verstoesse.")

    sys.exit(1 if probleme else 0)
