import os
from protein_small_docking import download_pdbqt_files  # Import from protein_small_docking
from docking_system import run_docking  # Import from docking_system

def interactive_docking_pipeline():
    """
    Interactive docking pipeline for ligand-receptor docking.
    Prompts the user for input and runs the docking workflow.
    """
    # Get user inputs interactively
    input_csv = str(input('Enter the full path to the input CSV file (ligand-receptor pairs): '))
    vina_exe_path = str(input('Enter the full path to vina.exe (default: vina.exe): ') or 'vina.exe')
    outdir = str(input('Enter the directory to save the output files (default: vina_result): ') or 'vina_result')
    center_x = float(input('Enter the docking center X coordinate (default: 10.819): ') or 10.819)
    center_y = float(input('Enter the docking center Y coordinate (default: 2.607): ') or 2.607)
    center_z = float(input('Enter the docking center Z coordinate (default: -53.797): ') or -53.797)
    size_x = float(input('Enter the docking box size X (default: 60): ') or 60)
    size_y = float(input('Enter the docking box size Y (default: 60): ') or 60)
    size_z = float(input('Enter the docking box size Z (default: 60): ') or 60)
    energy_range = int(input('Enter the energy range parameter (default: 4): ') or 4)
    exhaustiveness = int(input('Enter the exhaustiveness parameter (default: 8): ') or 8)

    # Validate input CSV
    if not os.path.exists(input_csv):
        raise ValueError(f"Input CSV file '{input_csv}' does not exist in the specified directory.")

    # Validate Vina executable
    if not os.path.exists(vina_exe_path):
        raise ValueError(f"Vina executable '{vina_exe_path}' does not exist. Please provide a valid path.")

    # Validate or create output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Step 1: Download and Convert Ligands/Receptors to PDBQT
    print("\nStep 1: Downloading and converting ligands and receptors to PDBQT format...")
    download_pdbqt_files(csv_file=input_csv, out_dir=outdir)

    # Step 2: Run Docking
    print("\nStep 2: Running docking procedure...")
    run_docking(
        input_csv=input_csv,
        vina_exe_path=vina_exe_path,
        output_folder=outdir,
        center_x=center_x,
        center_y=center_y,
        center_z=center_z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        energy_range=energy_range,
        exhaustiveness=exhaustiveness
    )

    print("\nDocking pipeline completed successfully!")

if __name__ == "__main__":
    try:
        interactive_docking_pipeline()
    except Exception as e:
        print(f"[ERROR] {e}")
    input()