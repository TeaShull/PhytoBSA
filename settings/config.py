#!/usr/bin/env python
import os

# This file is for configuring paths AND subclasses.  
PACKAGE_NAME = "phytobsa"

# PATHS
## I wouldn't mess with these, unless you feel the call to adventure
CONFIG_DIR = os.path.dirname(os.path.realpath(__file__)) #<- YOU ARE HERE
BASE_DIR =  os.path.dirname(CONFIG_DIR)
MODULES_DIR = os.path.join(BASE_DIR, 'modules')
VCF_GEN_SCRIPT = os.path.join(MODULES_DIR, 'subprocess_VCFgen.sh')

"""
 Feel free to edit these to customize your input and output folders. 
just change "DATA_DIR" if you want to change where the data folder is located. 
the rest will autopopulate
"""
s