#!/usr/bin/env python
import os
import configparser


def check_data_dir(config):
    '''
    Checks if the data direcory is set. This needs to be done the first time the
    program is run. Will make the conda package a lot easier to make...
    '''

    data_dir_prefix = config.get('SYS', 'data_dir')
        
    if data_dir_prefix == 'None':
        raise ValueError("DATA_DIR is not set. Set the data directory using the command: ./phytobsa settings --set_data_dir <path-to-dir>")
    elif not os.path.exists(data_dir_prefix):
        raise ValueError(f"Path {data_dir_prefix} does not exist. Set the data directory using the command: ./phytobsa settings --set_data_dir <path-to-dir>")
    else:
        data_dir = os.path.join(data_dir_prefix, 'data')
        return data_dir

def setup_data_dir(DATA_DIR, REFERENCE_DIR, INPUT_DIR, OUTPUT_DIR, LOG_DIR): # Doesn't FileUtilities.setup_directory because of circular import...
    '''
    Sets up the data direcotry based on the path set in ./settings/config.ini
    '''
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

def set_threads_limit(config):
    cpus = os.cpu_count()
    try:
        threads_limit = config.get('SYS', 'threads_limit')
        if not threads_limit:
            raise ValueError("Threads limit is not set. You can set it using ./phytobsa settings --set_threads_limit <int>. Proceeding with the(program detected limit - 2 ")
        threads_limit = int(threads_limit)
        if threads_limit > cpus:
            raise ValueError("Invalid threads limit. Please set a threads limit between 1 and the total number of CPUs. Proceeding with the(program detected limit - 2 ")
    except ValueError as e:
        print(e)
        threads_limit = cpus - 2
    
    return threads_limit

# Create the config object and read the configuration file
CONFIG_DIR = os.path.dirname(os.path.realpath(__file__))
CONFIG_INI =  os.path.join(CONFIG_DIR, 'config.ini')

config = configparser.ConfigParser()
config.read(CONFIG_INI)

#Set globals
THREADS_LIMIT = set_threads_limit(config)

BASE_DIR =  os.path.dirname(CONFIG_DIR)
MODULES_DIR = os.path.join(BASE_DIR, 'modules')
VCF_GEN_SCRIPT = os.path.join(MODULES_DIR, 'subprocess_VCFgen.sh')     

DATA_DIR = check_data_dir(config)
REFERENCE_DIR = os.path.join(DATA_DIR, 'references')
INPUT_DIR = os.path.join(DATA_DIR, 'input')
OUTPUT_DIR = os.path.join(DATA_DIR, 'output')
LOG_DIR = os.path.join(DATA_DIR, 'logs')
LOG_DATABASE_NAME = 'phytoBSAlog.db'
LOG_DATABASE_PATH = os.path.join(LOG_DIR, LOG_DATABASE_NAME)



