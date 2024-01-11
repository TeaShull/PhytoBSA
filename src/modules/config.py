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

## SUBCLASSES 
class ExperimentDictionary(dict):
    required_keys = {
        'key':['wt', 'mu', 'allele', 'pairedness', 'vcf_ulid', 'analysis_ulid', 
        'vcf_table_path','output_dir_path', 'analysis_ulid']
    }

    def __setitem__(self, key, value):
        if key not in self.required_keys:
            raise ValueError("Invalid key!")

        if (key == 'wt' or key == 'mu') and not isinstance(value, (list, str)):
            raise ValueError(f"{key} must be a list or a string!")

        if not isinstance(value, str):
            raise ValueError(f"{key} must be a string!")

        super().__setitem__(key, value)