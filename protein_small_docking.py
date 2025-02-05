import os
import re
import pandas as pd
import requests
import pubchempy as pcp
from openbabel import openbabel
from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

def is_pubchem_cid(s):
    """
    Checks if the input string follows the PubChem CID format:
    - Should be a positive integer.
    - Should not have leading zeros unless it's "0".
    """
    return bool(re.fullmatch(r"[1-9]\d*|0", s))

def download_ligand_sdf(ligand_name, output_sdf):
    """
    Download the specified ligand from PubChem in SDF format using PubChemPy.
    'ligand_name' can be a PubChem-recognized name (e.g. 'aspirin') or CID, etc.
    'output_sdf' is the path where the SDF should be saved.
    """
    try:
        namespace_choice = 'cid' if is_pubchem_cid(str(ligand_name)) else 'name'
        # 'name' namespace means we're looking up by compound name
        # record_type='3d' requests 3D coordinates from PubChem
        pcp.download(
            'SDF',
            output_sdf,
            identifier=ligand_name,
            namespace=namespace_choice,
            record_type='3d',
            overwrite=True  # Overwrite existing file if present
        )
        print(f"[INFO] Successfully downloaded SDF for '{ligand_name}' -> {output_sdf}")
        return True
    except Exception as e:
        print(f"[ERROR] Download failed for '{ligand_name}': {e}")
        return False

def pdb_to_pdbqt_mgltools(receptor_pdb, receptor_pdbqt):
    """
    Convert a PDB receptor file to PDBQT using the same logic as:
        prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt -A hydrogens -U waters

    :param receptor_pdb: Path to the input PDB receptor file
    :param receptor_pdbqt: Path where the PDBQT file should be written
    """
    # Read the PDB
    mols = Read(receptor_pdb)
    if not mols:
        raise ValueError(f"Could not read receptor file '{receptor_pdb}'.")
    mol = mols[0]

    # Prepare receptor as in 'prepare_receptor4.py -A hydrogens -U waters'
    #   - repairs='hydrogens': adds hydrogens if missing
    #   - charges_to_add='Kollman': add Kollman partial charges
    #   - cleanup='waters': remove water residues
    #   - mode='automatic': writes out the file without user interaction
    _ = AD4ReceptorPreparation(
        mol,
        mode='automatic',
        repairs='hydrogens',
        charges_to_add='Kollman',
        cleanup='waters',
        outputfilename=receptor_pdbqt
    )

    print(f"[INFO] Wrote receptor PDBQT: {receptor_pdbqt}")


def convert_to_pdbqt(input_file, output_pdbqt, input_format='sdf'):
    """
    Convert an SDF ligand file to PDBQT using Open Babel.
    """
    assert input_format in ['sdf','pdb'], f"Wrong input for 'input_format': '{input_format}'\nShould be among ['sdf','pdb']"
    try:
        if input_format == 'sdf':
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat(input_format)
            obConversion.SetOutFormat("pdbqt")
        
            mol = openbabel.OBMol()
            # Read the molecule from SDF
            obConversion.ReadFile(mol, input_file)
        
            # Add hydrogens
            mol.AddHydrogens()
        
            # Optional: generate 3D using OBBuilder
            builder = openbabel.OBBuilder()
            builder.Build(mol)
        
            # Write out to PDBQT
            obConversion.WriteFile(mol, output_pdbqt)
        elif input_format == 'pdb':
            pdb_to_pdbqt_mgltools(input_file, output_pdbqt)
        
        print(f"[INFO] Converted {input_format} format file: {input_file} -> {output_pdbqt} using openbabel core classes.")
    except Exception as e:
        print(f"[ERROR] Conversion failed for '{input_file}': {e}")
        return False

def download_alphafold_pdb(uniprot_id, output_pdb):
    """
    Download a PDB file from AlphaFold for a given UniProt ID.
    Uses the standard URL pattern:
      https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb
    Adjust if your data differs.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        r = requests.get(url)
        if r.status_code == 200:
            with open(output_pdb, 'wb') as f:
                f.write(r.content)
            print(f"[INFO] Downloaded AlphaFold PDB: {uniprot_id} -> {output_pdb}")
        else:
            print(f"[ERROR] Failed to download PDB for {uniprot_id}, code: {r.status_code}")
    except Exception as e:
        print(f"[ERROR] Download failed for {uniprot_id}: {e}")


def download_pdbqt_files(csv_file, out_dir):
    """
    Reads a CSV (ligand, receptor),
    Downloads & converts each ligand/receptor,
    Performs docking using the Python 'vina' API,
    Writes out the best scores to a CSV.
    """

    # Create subdirectories for organization
    lig_dir = os.path.join(out_dir, "ligands")
    rec_dir = os.path.join(out_dir, "receptors")
    os.makedirs(lig_dir, exist_ok=True)
    os.makedirs(rec_dir, exist_ok=True)

    # Read the CSV
    df_lr = pd.read_csv(csv_file)
    df_lr['ligand'] = df_lr['ligand'].astype(str)
    df_lr['receptor'] = df_lr['receptor'].astype(str)
    
    for idx in range(len(df_lr)):
        ligand_name = df_lr.iloc[idx]['ligand'].strip()
        receptor_id = df_lr.iloc[idx]['receptor'].strip()

        print(f"\n=== Processing {ligand_name} vs {receptor_id} ===")

        # 1) Download & convert ligand
        ligand_sdf = os.path.join(lig_dir, f"{ligand_name}.sdf")
        ligand_pdbqt = os.path.join(lig_dir, f"{ligand_name}.pdbqt")

        if not os.path.exists(ligand_pdbqt):
            # Download SDF if needed
            if not os.path.exists(ligand_sdf):
                success = download_ligand_sdf(ligand_name, ligand_sdf)
                if not success:
                    continue  # skip if download failed
            convert_to_pdbqt(ligand_sdf, ligand_pdbqt, input_format='sdf')

        # 2) Download & convert receptor
        receptor_pdb = os.path.join(rec_dir, f"{receptor_id}.pdb")
        receptor_pdbqt = os.path.join(rec_dir, f"{receptor_id}.pdbqt")

        if not os.path.exists(receptor_pdbqt):
            if not os.path.exists(receptor_pdb):
                download_alphafold_pdb(receptor_id, receptor_pdb)
            if os.path.exists(receptor_pdb):
                convert_to_pdbqt(receptor_pdb, receptor_pdbqt, input_format='pdb')


if __name__ == "__main__":
    """
    Example usage:
      python vina_python_pipeline.py
    You might call run_vina_docking_with_python() here directly or 
    import it from another script.
    """
    # Example call:
    download_pdbqt_files(
        csv_file="ligand_receptor.csv",
        out_dir="vina_results"
    )