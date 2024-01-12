#!/usr/bin/env python
import os
import sys

modules_dir = os.path.join(os.getcwd(), 'src', 'modules')
sys.path.append(modules_dir)
from config import (SRC_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR, REFERENCE_DIR)

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

