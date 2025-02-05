import os
import re
import pandas as pd
import subprocess
from multiprocessing import Pool, cpu_count

import pandas as pd
import re

def process_vina_results(df):
    """
    Processes an AutoDock Vina output dataframe and extracts docking results, including grid parameters and seed.

    Parameters:
        df (pd.DataFrame): Original dataframe with columns:
            ["Receptor", "Ligand", "Config File", "Output File", "Vina Output"]

    Returns:
        df_mod (pd.DataFrame): Processed dataframe with columns:
            ["Receptor", "Ligand", "Pose", "mode_kcal_mol", "affinity_rmsd_lb", "dist_from_best_mode_rmsd_ub",
             "Grid_X", "Grid_Y", "Grid_Z", "Size_X", "Size_Y", "Size_Z", "Grid_Space", "Exhaustiveness", "Seed"]
    """
    docking_results = []

    for _, row in df.iterrows():
        receptor = row["Receptor"]
        ligand = row["Ligand"]
        vina_output = row["Vina Output"]

        # Pose pattern
        pose_pattern = r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)"
        pose_matches = re.findall(pose_pattern, vina_output, re.MULTILINE)

        # Extract grid center (X, Y, Z)
        grid_center_match = re.search(
            r"Grid center:\s+X\s+(-?\d+\.\d+)\s+Y\s+(-?\d+\.\d+)\s+Z\s+(-?\d+\.\d+)",
            vina_output
        )
        grid_x, grid_y, grid_z = grid_center_match.groups() if grid_center_match else (None, None, None)

        # Extract grid size (X, Y, Z)
        grid_size_match = re.search(
            r"Grid size\s+:\s+X\s+(\d+)\s+Y\s+(\d+)\s+Z\s+(\d+)",
            vina_output
        )
        size_x, size_y, size_z = grid_size_match.groups() if grid_size_match else (None, None, None)

        # Extract grid space
        grid_space_match = re.search(r"Grid space\s+:\s+(\d+\.\d+)", vina_output)
        grid_space = grid_space_match.group(1) if grid_space_match else None

        # Extract exhaustiveness
        exhaustiveness_match = re.search(r"Exhaustiveness:\s+(\d+)", vina_output)
        exhaustiveness = exhaustiveness_match.group(1) if exhaustiveness_match else None

        # Extract seed
        seed_match = re.search(r"random seed:\s+(-?\d+)", vina_output)
        seed = seed_match.group(1) if seed_match else None

        # Populate the results
        for match in pose_matches:
            docking_results.append({
                "Receptor": receptor,
                "Ligand": ligand,
                "Pose": int(match[0]),  # Pose number (1 to n)
                "mode_kcal_mol": float(match[1]),  # Affinity (kcal/mol)
                "affinity_rmsd_lb": float(match[2]),  # RMSD lower bound
                "dist_from_best_mode_rmsd_ub": float(match[3]),  # RMSD upper bound
                "Grid_X": float(grid_x) if grid_x else None,
                "Grid_Y": float(grid_y) if grid_y else None,
                "Grid_Z": float(grid_z) if grid_z else None,
                "Size_X": int(size_x) if size_x else None,
                "Size_Y": int(size_y) if size_y else None,
                "Size_Z": int(size_z) if size_z else None,
                "Grid_Space": float(grid_space) if grid_space else None,
                "Exhaustiveness": int(exhaustiveness) if exhaustiveness else None,
                "Seed": int(seed) if seed else None
            })

    df_mod = pd.DataFrame(docking_results)
    return df_mod


def run_docking(input_csv, vina_exe_path, output_folder, center_x, center_y, center_z, size_x, size_y, size_z, energy_range, exhaustiveness, seed_num):
    """
    Main function to handle docking with user-specified parameters.

    Parameters:
        input_csv (str): Path to the input CSV file with receptor-ligand pairs.
        vina_exe_path (str): Path to the `vina.exe` executable file.
        output_folder (str): Folder to save generated files and results.
        center_x, center_y, center_z (float): Coordinates for the docking center.
        size_x, size_y, size_z (int): Dimensions of the docking box.
        energy_range (int): Energy range parameter for Vina.
        exhaustiveness (int): Exhaustiveness parameter for Vina.
        seed_num (int): Seed number for Vina.
    """
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Load the receptor-ligand pairs from the input CSV
    data = pd.read_csv(input_csv)
    data['ligand'] = data['ligand'].astype(str)
    data['receptor'] = data['receptor'].astype(str)

    # Prepare data for multiprocessing
    rows = [
        (
            row['receptor'],
            row['ligand'],
            vina_exe_path,
            output_folder,
            center_x,
            center_y,
            center_z,
            size_x,
            size_y,
            size_z,
            energy_range,
            exhaustiveness,
            seed_num
        )
        for _, row in data.iterrows()
    ]

    # Use multiprocessing to process docking
    with Pool(cpu_count()) as pool:
        docking_results = pool.map(process_docking, rows)

    # Save the docking results to a CSV file in the output folder
    output_csv_raw = os.path.join(output_folder, "docking_results_raw.csv")
    results_df_raw = pd.DataFrame(docking_results)
    results_df_raw.to_csv(output_csv_raw, index=False)

    output_csv_proc = os.path.join(output_folder, "docking_results_final.csv")
    results_df_proc = process_vina_results(results_df_raw)
    results_df_proc.to_csv(output_csv_proc, index=False)

    print(f"Docking process completed. Results saved in {output_csv_proc}.")


def process_docking(args):
    """
    Processes a single docking task.

    Parameters:
        args (tuple): Contains necessary parameters for docking.

    Returns:
        dict: Result of the docking task.
    """
    (
        receptor,
        ligand,
        vina_exe_path,
        output_folder,
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        energy_range,
        exhaustiveness,
        seed_num
    ) = args

    receptor_ = receptor + ".pdbqt"
    ligand_ = ligand + ".pdbqt"
    receptor_ = os.path.join(output_folder, 'receptors', receptor_)
    ligand_ = os.path.join(output_folder, 'ligands', ligand_)
    
    # Generate the A.txt content
    command_txt_content = f"""
receptor = {receptor_}
ligand = {ligand_}

center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = {size_x}
size_y = {size_y}
size_z = {size_z}

energy_range = {energy_range}

exhaustiveness = {exhaustiveness}
"""
    # Save the command.txt file
    command_txt_path = os.path.join(output_folder, f"command_{ligand}_{receptor}.txt")
    with open(command_txt_path, "w") as f:
        f.write(command_txt_content.strip())

    dock_dir = os.path.join(output_folder, "docking")
    os.makedirs(dock_dir, exist_ok=True)

    # Run the vina command
    vina_output_path = os.path.join(dock_dir, f"dock_{ligand}_{receptor}.pdbqt")
    vina_command = f"{vina_exe_path} --config {command_txt_path} --out {vina_output_path} --seed {seed_num}"
    try:
        # Capture the output of the vina command
        result = subprocess.run(vina_command, shell=True, check=True, capture_output=True, text=True)
        result_status = result.stdout  # Capture detailed command output
    except subprocess.CalledProcessError as e:
        result_status = f"Failed: {e.stderr or str(e)}"

    # Return results for this pair
    return {
        "Receptor": receptor,
        "Ligand": ligand,
        "Config File": command_txt_path,
        "Output File": vina_output_path,
        "Vina Output": result_status  # Include detailed output from Vina
    }
if __name__ == "__main__":
    # Example usage with user-specified parameters
    input_csv = "ligand_receptor.csv"  # Replace with your input CSV file path
    vina_exe_path = "vina.exe"  # Replace with the path to your vina executable
    output_folder = "vina_results"  # Replace with your desired output directory

    # Docking box and Vina parameters
    center_x, center_y, center_z = 10.819, 2.607, -53.797
    size_x, size_y, size_z = 60, 60, 60
    energy_range = 4
    exhaustiveness = 8
    seed_num = 1234

    # Run the docking pipeline
    run_docking(
        input_csv,
        vina_exe_path,
        output_folder,
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        energy_range,
        exhaustiveness,
        seed_num
    )