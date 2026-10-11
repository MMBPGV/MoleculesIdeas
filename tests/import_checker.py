from __future__ import annotations

import ast
import sys
from pathlib import Path
from collections import defaultdict


# ==========================================================
# Einstellungen
# ==========================================================

# Ordner, der untersucht werden soll.
# "." = aktueller Ordner
ROOT = Path("..")


# Ordner, die nicht untersucht werden sollen.
IGNORIEREN = {
    "__pycache__",
    ".git",
    ".venv",
    "venv",
    "env",
    "github_ausschluss",
}


# ==========================================================
# Hilfsfunktionen
# ==========================================================

def ist_python_datei(pfad: Path) -> bool:
    return pfad.suffix.lower() == ".py"


def soll_ignoriert_werden(pfad: Path) -> bool:
    return any(teil in IGNORIEREN for teil in pfad.parts)


def import_name(node: ast.AST) -> str | None:
    """Ermittelt den obersten Paketnamen eines Import-Nodes."""

    if isinstance(node, ast.Import):
        if not node.names:
            return None

        return node.names[0].name.split(".")[0]

    if isinstance(node, ast.ImportFrom):
        if node.level > 0:
            # Relativer Import, z. B.
            # from .foo import bar
            # from ..chemie import xyz
            return "." * node.level + (node.module or "")

        if node.module:
            return node.module.split(".")[0]

    return None


def analysiere_datei(pfad: Path):
    try:
        text = pfad.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        print(f"[WARNUNG] Keine UTF-8-Datei: {pfad}")
        return []

    try:
        baum = ast.parse(text, filename=str(pfad))
    except SyntaxError as e:
        print(f"[FEHLER] Syntaxfehler in {pfad}")
        print(f"         Zeile {e.lineno}: {e.msg}")
        return []

    imports = []

    for node in ast.walk(baum):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            name = import_name(node)

            if name:
                imports.append({
                    "name": name,
                    "zeile": node.lineno,
                    "datei": pfad,
                    "node": node,
                })

    return imports


# ==========================================================
# Hauptanalyse
# ==========================================================

def main():
    print("=" * 70)
    print("Python Import Analyzer")
    print("=" * 70)
    print()

    alle_imports = []
    python_dateien = []

    for pfad in ROOT.rglob("*.py"):
        if soll_ignoriert_werden(pfad):
            continue

        python_dateien.append(pfad)

        imports = analysiere_datei(pfad)
        alle_imports.extend(imports)

    # ------------------------------------------------------
    # Nach Paket sortieren
    # ------------------------------------------------------

    pakete = defaultdict(list)

    for eintrag in alle_imports:
        pakete[eintrag["name"]].append(eintrag)

    # ------------------------------------------------------
    # Ausgabe
    # ------------------------------------------------------

    print(f"Gefundene Python-Dateien: {len(python_dateien)}")
    print(f"Gefundene Import-Anweisungen: {len(alle_imports)}")
    print(f"Verschiedene Import-Namen: {len(pakete)}")
    print()

    print("=" * 70)
    print("IMPORTS")
    print("=" * 70)

    for name in sorted(pakete):
        eintraege = pakete[name]

        print()
        print(f"{name}")

        for eintrag in eintraege:
            print(
                f"    {eintrag['datei']}"
                f": Zeile {eintrag['zeile']}"
            )

    # ------------------------------------------------------
    # Zusammenfassung
    # ------------------------------------------------------

    print()
    print("=" * 70)
    print("ZUSAMMENFASSUNG")
    print("=" * 70)

    for name in sorted(pakete):
        print(f"{name:<30} {len(pakete[name])}x")

    print()
    print("=" * 70)
    print("Hinweis")
    print("=" * 70)
    print()
    print(
        "Diese Liste zeigt alle gefundenen Imports."
    )
    print(
        "Ob ein Import zur Python-Standardbibliothek gehört "
        "oder als externe Abhängigkeit installiert werden muss, "
        "wird bewusst noch nicht automatisch entschieden."
    )


if __name__ == "__main__":
    main()