#!/usr/bin/env python
import os
import sys
<<<<<<< HEAD

modules_dir = os.path.join(os.getcwd(), 'src', 'modules')
sys.path.append(modules_dir)
from config import (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)

=======
from setuptools import setup, find_packages
import subprocess
from settings.config import (
    PACKAGE_NAME, BASE_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR
)
>>>>>>> f27469d (major updates to file structure and naming)
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

<<<<<<< HEAD
def main():
    required_directories = (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)
    create_directories(required_directories)

if __name__ == '__main__':
    main()
=======
required_directories = (INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)
create_directories(required_directories)
>>>>>>> f27469d (major updates to file structure and naming)
