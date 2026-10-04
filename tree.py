
from pathlib import Path


# Ordner und Dateien, die nicht angezeigt werden sollen
IGNORIEREN = {
    "__pycache__",
    ".git",
    ".venv",
    "venv",
    ".idea",
    ".vscode",
}


def tree(pfad: Path, prefix: str = ""):
    """Gibt den Verzeichnisbaum rekursiv aus."""

    eintraege = sorted(
        [
            p for p in pfad.iterdir()
            if p.name not in IGNORIEREN
        ],
        key=lambda p: (p.is_file(), p.name.lower())
    )

    for index, eintrag in enumerate(eintraege):
        ist_letzter = index == len(eintraege) - 1

        zeichen = "└── " if ist_letzter else "├── "
        print(prefix + zeichen + eintrag.name)

        if eintrag.is_dir():
            neuer_prefix = prefix + ("    " if ist_letzter else "│   ")
            tree(eintrag, neuer_prefix)


if __name__ == "__main__":
    projektordner = Path(__file__).resolve().parent

    print(projektordner.name)
    tree(projektordner)
