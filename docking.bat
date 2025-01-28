@echo off
REM Protein Docking Automation Script Using Conda
REM This script creates a Conda environment, installs dependencies, and runs the interactive docking pipeline.

echo ========================================
echo Protein Docking Automation Script
echo ========================================

REM Check if Conda is installed
where conda >nul 2>&1
IF ERRORLEVEL 1 (
    echo [ERROR] Conda is not installed or not in PATH. Please install Miniconda or Anaconda and try again.
    exit /b 1
)

REM Initialize Conda for this script
call conda init >nul 2>&1

REM Check if 'docking' environment exists
echo Checking for Conda environment 'docking'...
conda env list | findstr docking >nul
IF ERRORLEVEL 1 (
    echo [INFO] Conda environment 'docking' does not exist. Creating it now...
    conda create -n docking python=3.9 -y
    call conda activate docking
    echo Installing required Python packages...
    pip install pandas requests pubchempy openbabel git+https://github.com/jaimergp/autodocktools-prepare-py3k.git
) ELSE (
    echo [INFO] Conda environment 'docking' exists. Activating it...
    call conda activate docking
)

REM Run the Python docking pipeline
echo Running interactive docking pipeline...
python interactive_docking_pipeline.py

REM Completion message
echo ========================================
echo Protein Docking Automation Completed!
echo ========================================
pause