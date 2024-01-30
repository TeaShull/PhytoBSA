#!/usr/bin/env python
import os
import sys

from settings.config import (
    BASE_DIR, DATA_DIR, REFERENCE_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR
)

def setup_data_dir():
    """
    Check for directories and create them if absent.
    Args:
    directories (list): List of directory paths to be checked/created.
    """
    required_directories = (
        DATA_DIR, REFERENCE_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR
    )
    for directory in required_directories:
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
                print(f"Directory created: {directory}")
            except OSError as e:
                print(f"Error creating directory {directory}: {e}")