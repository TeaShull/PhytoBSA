#!/usr/bin/env python
import os

# This file is for configuring paths AND subclasses.  

## PATHS
BASE_DIR = os.getcwd()
SRC_DIR =  os.path.join(BASE_DIR, 'src')
REFERENCE_DIR = os.path.join(BASE_DIR, 'references')
INPUT_DIR = os.path.join(BASE_DIR, 'input')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
LOG_DIR = os.path.join(SRC_DIR, 'logs')
TEMPLATE_DIR = os.path.join(SRC_DIR,'templates')
STATIC_DIR = os.path.join(SRC_DIR,'static')
MODULES_DIR = os.path.join(SRC_DIR,'modules')

