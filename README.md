# Chemistry Lab

An experimental Python project for automated molecular generation, chemical rule-based evaluation, and statistical learning.

The project provides an experimental environment for generating molecular structures, evaluating them using simplified chemical rules, and using statistical information from previous experiments to influence future molecule generation.

> **Important:** This project is experimental software. Its chemical evaluations are simplified and must not be interpreted as proof of real-world chemical stability, safety, synthesizability, or existence.

---

## Features

* Automated generation of molecular structures
* Rule-based chemical evaluation
* Statistical learning through a persistent memory system
* RDKit-based molecular validation
* Optional novelty checking against external chemical databases
* Graphical user interface
* Console-based operation
* Learning-curve visualization
* Lewis-structure visualization
* Configurable laboratory parameters
* Automated tests
* Persistent project data

---

## How It Works

The general workflow of the system is:

```text
Molecule Generation
        │
        ▼
Chemical Evaluation
        │
        ▼
Critic / Scoring
        │
        ├───────────────┐
        ▼               ▼
Memory Update       Novelty Check
        │               │
        └───────┬───────┘
                ▼
        Future Generation
```

The system generates molecular structures and evaluates them using a collection of simplified chemical rules.

The resulting evaluation can be stored in the statistical memory. Information from previous experiments can then influence future molecule generation.

The goal is not to simulate chemistry completely, but to experiment with automated molecular generation and adaptive selection.

---

## Chemical Evaluation

The chemical evaluation is performed by the `Kritiker` system.

The relevant modules are located in:

```text
chemlabor/chemie/
├── kritiker.py
├── kritiker_regeln.py
├── neuheit.py
├── rdkit_bruecke.py
└── zufallsgenerator.py
```

### Critic

`kritiker.py` coordinates the evaluation of generated molecules.

### Chemical Rules

`kritiker_regeln.py` contains the individual rules used by the critic.

Depending on the implemented rules, the evaluation can consider aspects such as:

* Valence compatibility
* Bonding constraints
* Element-specific behavior
* Structural consistency
* Electron-related heuristics
* Other project-specific chemical criteria

These rules are intentionally simplified. They provide an experimental scoring system rather than a complete chemical model.

### RDKit

The project can use RDKit to validate molecular structures and obtain additional structural information.

For example, RDKit can be used to determine whether a generated molecular representation is structurally valid and to obtain a molecular formula.

A structure accepted by RDKit is not automatically chemically stable, safe, or experimentally realizable.

---

## Statistical Memory

The statistical system is located in:

```text
chemlabor/ki/
├── gedaechtnis.py
└── konstruktor.py
```

The memory system stores information about previous experiments and can use this information when constructing future molecules.

The system is primarily statistical and rule-based. It is not a large neural-network-based AI model.

This approach is intended to keep the system relatively transparent while still allowing the molecule-generation process to adapt based on previous results.

### Local Memory Data

The directory

```text
daten/gedaechtnis/
```

is reserved for locally generated memory data.

The current public repository does not contain a previously accumulated memory state.

This means that experiments can start with a clean memory and generate their own data locally.

---

## Novelty Checking

Novelty checking is implemented in:

```text
chemlabor/chemie/neuheit.py
```

The system can optionally compare generated molecules or formulas with information from external chemical databases.

A local cache can be created during runtime to reduce repeated external requests.

The cache is intentionally not part of the public repository.

### Important Limitation

A molecule not found in an external database does **not** prove that the molecule is new.

Database searches can be incomplete, compounds can be represented in different ways, and some compounds may not be present in publicly accessible databases.

The novelty system should therefore be interpreted as a database comparison rather than proof of chemical novelty.

---

## Graphical User Interface

The graphical interface is located in:

```text
chemlabor/gui/
```

It provides functionality for interacting with the laboratory, configuring experiments, displaying results, and visualizing molecular information.

Important components include:

```text
app.py
dialog_einstellungen.py
fenster.py
formatierung.py
lernkurve.py
lewis.py
modell.py
tab_ergebnisse.py
tab_live.py
theme.py
```

The GUI is separated from the underlying chemical and statistical logic.

This allows the core laboratory system to be used independently of the graphical interface.

---

## Console Mode

The project also provides a console-based interface.

The main entry points are:

```text
start_gui.py
start_konsole.py
```

Start the graphical interface with:

```bash
python start_gui.py
```

Start the console version with:

```bash
python start_konsole.py
```

The package also contains:

```text
chemlabor/__main__.py
```

which provides a Python module entry point.

---

## Configuration

Configuration is handled by:

```text
chemlabor/config.py
```

The project uses a `LaborConfig` data structure to keep laboratory configuration separate from the implementation.

This allows experiment parameters to be changed without directly modifying the core laboratory logic.

---

## Project Structure

The current public project structure is:

```text
chemistry-lab/
│
├── chemlabor/
│   ├── chemie/
│   │   ├── __init__.py
│   │   ├── kritiker.py
│   │   ├── kritiker_regeln.py
│   │   ├── neuheit.py
│   │   ├── rdkit_bruecke.py
│   │   └── zufallsgenerator.py
│   │
│   ├── gui/
│   │   ├── __init__.py
│   │   ├── app.py
│   │   ├── dialog_einstellungen.py
│   │   ├── fenster.py
│   │   ├── formatierung.py
│   │   ├── lernkurve.py
│   │   ├── lewis.py
│   │   ├── modell.py
│   │   ├── tab_ergebnisse.py
│   │   ├── tab_live.py
│   │   └── theme.py
│   │
│   ├── ki/
│   │   ├── __init__.py
│   │   ├── gedaechtnis.py
│   │   └── konstruktor.py
│   │
│   ├── __init__.py
│   ├── __main__.py
│   ├── config.py
│   ├── labor.py
│   └── pfade.py
│
├── daten/
│   ├── gedaechtnis/
│   └── periodensystem.json
│
├── tests/
│   ├── test_architektur.py
│   ├── test_config.py
│   ├── test_gui_modell.py
│   └── test_gui_smoke.py
│
├── README.md
├── start_gui.py
├── start_konsole.py
└── tree.py
```

---

## Architecture

The project is divided into several major components.

### Chemistry

```text
chemlabor/chemie/
```

Contains:

* Molecular generation
* Chemical evaluation
* Critic rules
* RDKit integration
* Novelty checking

The chemistry package does not depend on the GUI.

### Statistical System

```text
chemlabor/ki/
```

Contains:

* Statistical memory
* Molecule construction
* Learning-related logic

The statistical system interacts with the laboratory and chemistry components without being part of the graphical interface.

### GUI

```text
chemlabor/gui/
```

Contains all graphical user interface components.

The GUI is responsible for presentation and user interaction rather than implementing the underlying chemical rules.

### Laboratory

```text
chemlabor/labor.py
```

Coordinates the different parts of the project and represents the central laboratory logic.

### Configuration

```text
chemlabor/config.py
```

Contains the configuration used by the laboratory.

### Paths and Data

```text
chemlabor/pfade.py
daten/
```

Handle project paths and required project data.

---

## Tests

The project contains tests for important parts of the current implementation.

The available tests are:

```text
tests/
├── test_architektur.py
├── test_config.py
├── test_gui_modell.py
└── test_gui_smoke.py
```

They can be executed individually:

```bash
python tests/test_architektur.py
python tests/test_config.py
python tests/test_gui_modell.py
python tests/test_gui_smoke.py
```

The tests are intended to detect implementation problems and regressions during development.

---

## Project Data

The public project currently contains only the data required by the application.

### Periodic Table

```text
daten/periodensystem.json
```

Contains the element data used by the project.

### Memory Directory

```text
daten/gedaechtnis/
```

The directory is reserved for locally generated memory data.

Generated memory files are intentionally not included in the public repository.

### Runtime Data

Some components can generate additional local runtime data.

For example, the novelty system can create a local database cache.

Runtime-generated data is kept outside the public project.

---

## Development Files

Historical backups, migration files, temporary data, and other development artifacts are kept separately from the public project.

They are not required to run the current version of Chemistry Lab.

The following files and directories are therefore not part of the public repository:

```text
_archiv/
_backup_schritt2/
_backup_vor_migration/
migration.py
aus.py
gedächtnis.json
```

Runtime-generated files such as:

```text
neuheit_cache.sqlite
autosave_gedaechtnis.json
```

are also excluded from the public project.

This separation keeps the repository focused on the current implementation instead of historical development data.

---

## Tree Utility

`tree.py` is a small development utility for displaying the project's directory structure.

It is included in the repository because it can be useful during development and documentation.

It is not required for the operation of Chemistry Lab itself.

---

## Scientific Limitations

This project does not perform full computational chemistry.

The molecular evaluation is based on simplified rules and heuristics. It does not replace:

* Quantum-chemical calculations
* Molecular dynamics
* Experimental chemistry
* Professional chemical databases
* Laboratory validation
* Detailed thermodynamic calculations
* Detailed kinetic calculations

A high score from the critic means only that a molecule performed well according to the implemented rules.

It does not prove that the molecule is:

* Chemically stable
* Safe
* Synthesizable
* Experimentally obtainable
* Thermodynamically favorable
* Kinetically stable
* Previously unknown

The same limitations apply to RDKit validation and database searches.

---

## Future Development

Possible future improvements include:

* More sophisticated chemical rules
* Improved molecular generation
* Improved statistical learning
* More detailed structural analysis
* Improved novelty detection
* More extensive automated tests
* Performance improvements for large experiments
* More advanced molecular visualization
* Improved experiment management
* Additional chemical data sources
The current architecture is intended to make these components independently extensible.

---

## AI-Assisted Development

Parts of this project were developed with assistance from AI-based programming tools.

AI assistance was used for tasks such as:

* Code analysis
* Debugging
* Refactoring
* Architecture planning
* Documentation
* Identifying possible implementation problems

The resulting code is part of an experimental development project and should be reviewed and tested independently.

---

## Requirements

The project requires Python.

Some functionality additionally depends on external Python packages, particularly RDKit.

A typical development environment can be created with:

```bash
python -m venv .venv
```

Activate the environment and install the required dependencies according to the project's current setup.

---

## Running the Project

Start the graphical interface:

```bash
python start_gui.py
```

Start the console version:

```bash
python start_konsole.py
```

For development and testing, execute the individual test files from the project root.

---

## License

This project is released under the MIT License.

See the `LICENSE` file for the complete license text.

---

## Disclaimer

Chemistry Lab is an experimental software project intended for programming, computational experimentation, education, and research into automated molecular generation.

The software must not be used as the sole basis for real-world chemical decisions, laboratory procedures, safety assessments, or claims about the existence or stability of chemical compounds.

Any generated molecular structures should be treated as computational experiments until independently verified using appropriate scientific methods.
