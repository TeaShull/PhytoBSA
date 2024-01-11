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
# #In progress. 
# class ExperimentDictionary:
#     def __init__(self):
#         self.data = 
#             'key': {
#                 'required_keys': ['wt', 'mu', 'allele', 'pairedness', 'vcf_ulid', 'analysis_ulid', 
#                                   'vcf_table_path', 'output_dir_path', 'analysis_ulid'],
#                 'values': {}
#             }

#     def set_value(self, key, item_key, value):

#         if item_key not in self.data[key]['required_keys']:
#             raise ValueError(f"Invalid item key for {key}!")

#         self.data[key]['values'][item_key] = value

#     def get_value(self, key, item_key):
#         if key not in self.data:
#             raise ValueError("Invalid key!")

#         if item_key not in self.data[key]['values']:
#             raise ValueError(f"Item key not found for {key}!")

#         return self.data[key]['values'][item_key]
