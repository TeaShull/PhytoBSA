#!/usr/bin/env python
import os
import subprocess

modules_dir = os.path.join(os.getcwd(), 'src', 'modules')
sys.path.append(modules_dir)

from config import (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)

def install_environment():
    try:
        # Check if mamba is available
        subprocess.check_call(['mamba', '--version'])

        # Use mamba to create/update environment
        print("Detected mamba. Using mamba to install/update environment...")
        subprocess.check_call(['mamba', 'env', 'create', '-f', 'env_thale_BSA.yml'])
        print("Required environment installed/updated successfully.")
    except subprocess.CalledProcessError:
        # mamba not detected, use conda
        try:
            # Check if the required environment is installed with conda
            subprocess.check_call(['conda', 'env', 'create', '-f', 'env_thale_BSA.yml'])

            # Environment already exists
            print("Required environment is already installed.")
        except subprocess.CalledProcessError:
            # Install the environment with conda
            print("Installing required environment using conda...")
            subprocess.check_call(['conda', 'env', 'update', '-f', 'env_thale_BSA.yml'])
            print("Required environment installed/updated successfully.")

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

required_directories = (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)
create_directories(required_directories)
install_dependencies()

setup(
    name="thale_bsa",
    version="1.0.0",
    author="Tim Shull",
    author_email="less.tact_0m@icloud.com",
    description="A python library for bulk segregant analysis of EMS-derived point mutations in Arabidopsis thaliana",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TeaShull/PyAtBSA",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bioinformatics",
    ],
    keywords="Arabidopsis BSA Genomics",
    entry_points={
        "console_scripts": [
            "thale_bsa=src.__main__:main",
        ],
    },
)

