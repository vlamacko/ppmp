# ppmp - Predicting perturbations of modular proteins

ppmp is a Python module and command-line tool for investigating dynamics of 
modular protein structures. See Documentation section for more information.

### Content

1. Setup
   - prerequisites
1. Usage
   - functions
   - input data structure
1. Documentation
1. Progress

---

## Setup

1. Clone this repo.
      ```bash
     git clone https://github.com/vlamacko/ppmp.git
     cd ppmp
     ```
1. Make sure all prerequisites are installed.
   - Python 3.2+
   - all dependencies:
      - statsmodels (conda, PyPI pip)
      - tqdm (PyPI pip)
      - scipy (conda, PyPI pip)
      - seaborn (conda, PyPI pip)
      - pandas (conda, PyPI pip)
      - numpy (conda, PyPI pip)
      - matplotlib (conda, PyPI pip)
      - biopython (conda, PyPI pip)

    All these packages are available through conda or PyPI repos.

    **pip requirements file is available.** This is useful if we want to use
    the script as a purely command-line tools witout installing ppmp as a package.
    Here is an example using virtualenv (optional):
    ```bash
    virtualenv --python=<path/to/your/python-3.2+> venv    #optional
    source ./venv/bin/activate                             #optional -activate the virtual environment
    pip install -r requirements.txt                        #install required libraries locally
    ```

    If using ppmp command-line only we can export the scripts to the PYTHONPATH 
    for increase ease of use:
    ```bash
    PYTHONPATH="$PWD/ppmp/:$PYTHONPATH"
    export PYTHONPATH
    ```
   **Or use `setup.py` to install the package and its dependencies:**

   Install the package in the environment/system. This should actually install 
   all dependencies automatically.
   ```bash
   python setup.py install
   ```
   If there are some dependencies missing, install them manually from PyPI 
   or conda.

## Usage

- Submodule `ppmp.to_rmsd_csv`
   - Script to calculate RMSD from pdb file and save it in csv files.
   - Can be used after importing the package or as a command-line tool (more
   options in command-line mode planned).
- Submodule `ppmp.model`
   - Script for analysing the RMSD csv files.
   - Can be used after importing the package or as a command-line tool (fully 
   functional).
- Submodule `ppmp.protein`
  - Used by the other submodules to define new class `Protein`.

For more information (docstrings, comments) see the source code. There are also 
few tests available in `/tests/` (both command-line and module use).

### Input data

The default folder structure/naming scheme is the following. It is partially 
adjustable. More options will be added in future. The `data` folder should be 
located in the working directory.
```text
data
├── csv  <------------------>  Directory containing RMSD csv files.
|                                 - form: proteinName.csv
├── json  <----------------->  Directory containg module json files.
|                                 - form: proteinName.csv
├── pdb  <------------------>  Directory containing structural pdb files.
|                                 - unpertubed: proteinName.pdb
|                                 - pertibations x: proteinName_x.pdb
├── test  <----------------->  Directory of test json and pdb files.
|                                 - structure: proteinName.json
|                                 - unpertubed: proteinName.pdb
|                                 - pertibations x: proteinName_x.pdb
└── modules-length  <------->  File defining length of each module.
```
~~See example data in `/tests/data/`.~~ The data will be released later.

### Output data

The output is saved in `out` folder in the current working directory. More 
customisation to be done.

## Documentation

This module has been created as a tool for an undergraduate 4 week project 
about predicting dynamics in modular protein structures. It was conducted by 5 
students of Engineering Mathematics at University of Bristol. The project report
`/docs/protein-dynamics-report.pdf` can be used as an additional documentation.

## Progress

- [ ] finish writing README
- [x] create ./requirements.txt
- [x] create working tests
- [x] relative paths, should work on all OS
- [x] add paper as a doc
- [x] add missing graphs, code
- [x] create setup.py
- [ ] function annotations, type hints
- [ ] test cross-platform compatibility
- [ ] test virtualenv (linux)
- [ ] integrate `./tests/prediction_analysis.py` into the package