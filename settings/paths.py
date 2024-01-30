#!/usr/bin/env python
import os
import configparser

def check_data_dir(config_ini):
    config = configparser.ConfigParser()
    config.read(config_ini)

    if 'data_dir' in config['Settings']:
        data_dir_prefix = config.get('Settings', 'data_dir')
        
        if data_dir_prefix == 'None':
            raise ValueError("DATA_DIR is not set. Set the data directory using the command: ./phytobsa --data_dir <path-to-dir>")
        elif not os.path.exists(data_dir_prefix):
            raise ValueError(f"Path {data_dir_prefix} does not exist. Set the data directory using the command: ./phytobsa settings --set_data_dir <path-to-dir>")
        else:
            data_dir = os.path.join(data_dir_prefix, 'data')
            return data_dir
    else:
        raise ValueError("DATA_DIR is not set. Set the data directory using the command: ./phytobsa --data_dir <path-to-dir>")
def setup_data_dir():
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


CONFIG_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
CONFIG_INI =  os.path.join(CONFIG_DIR, 'config.ini')
BASE_DIR =  os.path.dirname(CONFIG_DIR)
MODULES_DIR = os.path.join(BASE_DIR, 'modules')
VCF_GEN_SCRIPT = os.path.join(MODULES_DIR, 'subprocess_VCFgen.sh')     


DATA_DIR = check_data_dir(CONFIG_INI)
LOG_DATABASE_NAME = 'phytoBSAlog.db'
LOG_DIR = os.path.join(DATA_DIR, 'logs')
LOG_DATABASE_PATH = os.path.join(LOG_DIR, LOG_DATABASE_NAME)
REFERENCE_DIR = os.path.join(DATA_DIR, 'references')
INPUT_DIR = os.path.join(DATA_DIR, 'input')
OUTPUT_DIR = os.path.join(DATA_DIR, 'output')
setup_data_dir()
