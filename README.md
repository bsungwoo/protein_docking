## README: protein_docking

# Protein Docking
An interactive script for testing docking between proteins and small molecules. This pipeline automates the process of downloading receptor and ligand files, converting them to `.pdbqt` format, and running docking simulations using **AutoDock Vina**. The script provides a user-friendly interface for inputting custom parameters and file paths.

---

### Features
1. **Interactive Input**: 
   - Specify input/output directories and docking parameters at runtime.
2. **Ligand & Receptor Handling**:
   - Downloads ligands from **PubChem** and receptors from **AlphaFold**.
   - Converts files to `.pdbqt` format using **Open Babel** and **AutoDockTools**.
3. **Docking Simulation**:
   - Generates configuration files for each ligand-receptor pair.
   - Runs **AutoDock Vina** for docking simulations.
   - Outputs docking results, including scores, to a specified directory.
4. **Parallel Processing**:
   - Utilizes multiple CPU cores for efficient processing.

---

### Setup and Requirements

#### 1. Software Dependencies
- Python 3.9+
- AutoDock Vina (download [here](http://vina.scripps.edu/))
- AlphaFold PDB files (used for receptor download)

#### 2. Python Libraries
Install the required Python libraries:
```bash
pip install pandas requests pubchempy openbabel git+https://github.com/jaimergp/autodocktools-prepare-py3k.git