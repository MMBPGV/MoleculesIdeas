# Chemistry Lab

An experimental Python project for automated molecular generation, chemical evaluation, statistical learning, and computational molecular analysis.

Chemistry Lab explores candidate molecular structures, evaluates them with chemical heuristics, and uses information from previous experiments to influence future generations. It also integrates RDKit and neural interatomic potentials for additional structural and energetic analysis.

> **Important:** Chemistry Lab is an experimental research and learning project. A high score, a valid molecular structure, or a calculated energy does not prove that a compound is stable, safe, synthesizable, or experimentally obtainable.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Running Chemistry Lab](#running-chemistry-lab)
- [How It Works](#how-it-works)
- [Chemical Evaluation](#chemical-evaluation)
- [Statistical Learning](#statistical-learning)
- [Computational Molecular Analysis](#computational-molecular-analysis)
- [Novelty Checking](#novelty-checking)
- [Project Structure](#project-structure)
- [Tests](#tests)
- [Scientific Limitations](#scientific-limitations)
- [AI-Assisted Development](#ai-assisted-development)
- [License](#license)

---

## Features

- Automated molecular structure generation
- Rule-based chemical evaluation and scoring
- Statistical memory that influences future experiments
- RDKit-based molecular validation
- Neural-potential-based geometry optimization and energetic estimates
- Automatic selection of supported computational models
- Lewis-structure visualization
- Live experiment monitoring and result analysis
- Learning-curve visualization
- Optional comparison with external chemical databases
- Graphical and console interfaces
- Configurable experiment parameters
- Persistent local experiment data

## Installation

### Requirements

- Python 3
- `pip`
- `venv`
- Tkinter for the graphical interface
- Git, if you want to clone the repository

Python package dependencies are listed in `requirements.txt`.

### Linux

On Debian- or Ubuntu-based systems, install the required system packages:

```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv python3-tk git
```

Clone the repository and enter the project directory:

```bash
git clone https://github.com/MMBPGV/MoleculesIdeas.git
cd MoleculesIdeas
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the Python dependencies:

```bash
python -m pip install -r requirements.txt
```

### Windows

Install Python 3 and make sure it is available through `PATH`.

Clone the repository:

```powershell
git clone https://github.com/MMBPGV/MoleculesIdeas.git
cd MoleculesIdeas
```

Create and activate a virtual environment:

```powershell
python -m venv .venv
.venv\Scripts\activate
```

Install the Python dependencies:

```powershell
python -m pip install -r requirements.txt
```

Tkinter is normally included with the standard Windows Python installation.

### macOS

Install Python 3 and Git, then run:

```bash
git clone https://github.com/MMBPGV/MoleculesIdeas.git
cd MoleculesIdeas

python3 -m venv .venv
source .venv/bin/activate

python -m pip install -r requirements.txt
```

The neural-potential models may download their model weights the first time they are used. This requires an internet connection and additional storage.

## Running Chemistry Lab

Run commands from the repository root with the virtual environment activated.

### Graphical interface

```bash
python start_gui.py
```

On Linux and macOS, you can use:

```bash
python3 start_gui.py
```

### Console interface

```bash
python start_konsole.py
```

On Linux and macOS:

```bash
python3 start_konsole.py
```

The package also provides a module entry point:

```bash
python -m chemlabor
```

## How It Works

The general workflow is:

1. Generate a candidate molecular structure.
2. Evaluate it using the implemented chemical rules.
3. Validate the molecular representation where supported.
4. Store experiment results in the statistical memory.
5. Use accumulated statistical information to influence future generations.
6. Optionally perform further structural, energetic, or database analysis.

The system combines heuristic evaluation with additional computational tools. These components provide different kinds of information and should not be treated as interchangeable forms of evidence.

## Chemical Evaluation

The chemical evaluation system is implemented in `chemlabor/chemie/`.

The critic applies project-specific rules to evaluate generated structures. Depending on the implemented criteria, these can include:

- Valence compatibility
- Bonding constraints
- Element-specific behavior
- Structural consistency
- Electron-related heuristics

The resulting score describes how well a candidate performs according to the implemented evaluation rules. It is not a measured physical property or a direct probability of experimental stability.

### RDKit validation

RDKit is used to validate molecular representations and obtain structural information, such as molecular formulas and molecular masses.

A structure accepted by RDKit is not automatically physically stable or experimentally realizable. The validity of a molecular representation and the real-world stability of a compound are different questions.

## Statistical Learning

The statistical system is located in `chemlabor/ki/`.

It records information from previous experiments and uses that information to influence subsequent molecule construction. This allows the generator to adapt statistically as experiments accumulate.

The memory system is a statistical and rule-based learning mechanism, not a large language model.

Experiment memory and other generated runtime data are stored locally. The public repository does not require a pre-existing memory state.

## Computational Molecular Analysis

The module `chemlabor/chemie/quantenchemie.py` provides additional calculations using neural interatomic potentials.

Depending on the elements present in a candidate, the module can select among supported models, including:

- **ANI-2x** for supported molecules containing H, C, N, O, F, S, and Cl.
- **MACE-OFF23** for supported molecular systems that also include elements such as P, Br, and I.
- **MACE-MP-0** for supported systems including Si and B, with important limitations for molecular accuracy.

The module can perform geometry optimization and estimate formation enthalpies using model energies and calibration. These calculations provide a different source of information from the heuristic critic.

The calculations have defined applicability limits. For example, the models do not support every element, charge state, or electronic configuration. Unsupported cases should be reported as not applicable rather than interpreted as valid numerical predictions.

Model weights may be downloaded on first use.

## Novelty Checking

The module `chemlabor/chemie/neuheit.py` can compare candidates with external chemical database information.

A local cache can reduce repeated requests.

A compound that is not found in a database is **not necessarily new**. Database coverage can be incomplete, and different representations of a compound can complicate comparisons. Novelty checking is a database comparison, not proof that a molecule has never been described.

## Project Structure

The main components are organized as follows:

```text
MoleculesIdeas/
├── chemlabor/
│   ├── chemie/
│   │   ├── kritiker.py
│   │   ├── kritiker_regeln.py
│   │   ├── neuheit.py
│   │   ├── quantenchemie.py
│   │   ├── rdkit_bruecke.py
│   │   └── zufallsgenerator.py
│   ├── gui/
│   │   ├── app.py
│   │   ├── fenster.py
│   │   ├── lernkurve.py
│   │   ├── lewis.py
│   │   ├── tab_ergebnisse.py
│   │   └── tab_live.py
│   ├── ki/
│   │   ├── gedaechtnis.py
│   │   └── konstruktor.py
│   ├── config.py
│   ├── labor.py
│   └── __main__.py
├── daten/
│   └── periodensystem.json
├── tests/
├── start_gui.py
├── start_konsole.py
├── requirements.txt
├── README.md
├── LICENSE
└── .gitignore
```

This is a simplified overview. The repository contains additional modules and files.

Generated memory files, caches, and other runtime data are kept locally rather than committed to the public repository.

## Tests

The project contains automated tests for selected parts of the implementation.

Run individual tests from the repository root, for example:

```bash
python tests/test_architektur.py
python tests/test_config.py
python tests/test_gui_modell.py
python tests/test_gui_smoke.py
```

On Linux and macOS, `python3` can be used instead of `python` when necessary.

The tests help identify implementation errors and regressions. They do not constitute scientific validation of generated compounds or calculated energies.

## Scientific Limitations

Chemistry Lab is an experimental computational project. Its outputs must be interpreted with care.

- Heuristic scores are not direct measurements of stability.
- RDKit validation checks molecular representations, not whether a compound can exist under particular physical conditions.
- Neural potentials have model-specific training domains and applicability limits.
- Estimated energies and formation enthalpies have uncertainties and are not substitutes for validated high-level calculations or experimental measurements.
- Database searches cannot prove that a candidate is novel.
- A favorable calculation does not establish safety, synthesizability, kinetic stability, or experimental existence.

The project is intended for programming, computational experimentation, education, and research into automated molecular generation.

## AI-Assisted Development

AI-based programming tools have been used to support parts of development, including code analysis, debugging, refactoring, architecture planning, and documentation.

The project author defines the goals, evaluates suggestions, makes implementation decisions, and tests the resulting software. The project remains experimental and should be reviewed and validated independently.

## License

Chemistry Lab is released under the MIT License.

See the `LICENSE` file for the complete license text.

## Disclaimer

This software is provided for educational and experimental purposes. Do not use its outputs as the sole basis for real-world chemical decisions, laboratory procedures, or safety assessments.

Generated structures and computational predictions require independent scientific validation before any real-world conclusions can be drawn.

