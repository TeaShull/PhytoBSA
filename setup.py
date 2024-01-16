#!/usr/bin/env python
import os
import sys
import site
import subprocess

modules_dir = os.path.join(os.getcwd(), 'src', 'modules')
sys.path.append(modules_dir)
from config import (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)


def install_conda_environment(yaml_file):
    try:
        # Check if Mamba is installed
        cmd_result = subprocess.run(
            ['mamba', '--version'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        # If Mamba is installed, create the environment using Mamba
        cmd_result = subprocess.run(
            ['mamba', 'env', 'create', '--file', yaml_file],
            check=True,
            stderr=subprocess.PIPE
        )
        print("Conda environment created using Mamba.")
    except subprocess.CalledProcessError as mamba_error:
        try:
            # Check if Conda is installed
            cmd_result = subprocess.run(
                ['conda', '--version'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            # If Conda is installed, create the environment using Conda
            cmd_result = subprocess.run(
                ['conda', 'env', 'create', '--file', yaml_file],
                check=True,
                stderr=subprocess.PIPE
            )
            print("Conda environment created using Conda.")
            
        except subprocess.CalledProcessError as conda_error:
            error_output = str(mamba_error.stderr, 'utf-8').strip() or str(conda_error.stderr, 'utf-8').strip()
            if "already exists" in error_output:
                print("Conda environment already exists. Proceeding...")
            else:
                print("Error encountered while creating the Conda environment.")
                print(error_output)
            
        except FileNotFoundError:
            print("Conda not found. Please install Conda.")

    except FileNotFoundError:
        print("Mamba or Conda not found. Please install Mamba or Conda.")

def create_directories(directories):
    """
    Check for directories and create them if absent.
    Args:
    directories (list): List of directory paths to be checked/created.
    """
    for directory in directories:
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
                print(f"Directory created: {directory}")
            except OSError as e:
                print(f"Error creating directory {directory}: {e}")

# Install the conda environment using Mamba or conda.
def main():
    conda_yaml = os.path.join(os.getcwd(), "environment.yml")
    install_conda_environment(conda_yaml)
    
    required_directories = (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)
    create_directories(required_directories)

if __name__ == '__main__':
    main()
