import os
import pandas as pd
import subprocess
from multiprocessing import Pool, cpu_count
import shutil

def run_docking(input_csv, vina_exe_path, output_folder, center_x, center_y, center_z, size_x, size_y, size_z, energy_range, exhaustiveness):
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
    """
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Load the receptor-ligand pairs from the input CSV
    data = pd.read_csv(input_csv)

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
            exhaustiveness
        )
        for _, row in data.iterrows()
    ]

    # Use multiprocessing to process docking
    with Pool(cpu_count()) as pool:
        docking_results = pool.map(process_docking, rows)

    # Save the docking results to a CSV file in the output folder
    output_csv = os.path.join(output_folder, "docking_results.csv")
    results_df = pd.DataFrame(docking_results)
    results_df.to_csv(output_csv, index=False)

    print(f"Docking process completed. Results saved in {output_csv}.")


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
        exhaustiveness
    ) = args

    receptor += ".pdbqt"
    ligand += ".pdbqt"
    receptor = os.path.join(output_folder, 'receptors', receptor)
    ligand = os.path.join(output_folder, 'ligands', ligand)
    
    # Generate the A.txt content
    command_txt_content = f"""
receptor = {receptor}
ligand = {ligand}

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
    vina_command = f"{vina_exe_path} --config {command_txt_path} --out {vina_output_path}"
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
        exhaustiveness
    )