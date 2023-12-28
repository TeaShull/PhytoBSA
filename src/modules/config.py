#!/usr/bin/env python
import os
import inspect
from datetime import datetime
import logging

BASE_DIR = os.getcwd()
SRC_DIR =  os.path.join(BASE_DIR, 'src')
INPUT_DIR = os.path.join(BASE_DIR, 'input')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
LOG_DIR = os.path.join(SRC_DIR, 'logs')
TEMPLATE_DIR = os.path.join(SRC_DIR,'templates')
STATIC_DIR = os.path.join(SRC_DIR,'static')
MODULES_DIR = os.path.join(SRC_DIR,'modules')

import logging

def setup_logger(name, log_path):
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)  # Set the default level for the logger

    # Create a FileHandler and set the level to INFO
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.INFO)

    # Create a formatter and add it to the handler
    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)

    # Add the handler to the logger
    logger.addHandler(file_handler)

    return logger

def log_handler(logger, message_type, message):
    log_handler_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
    function_name = inspect.currentframe().f_back.f_code.co_name
    caller_frame = inspect.currentframe().f_back
    script_name = os.path.basename(inspect.getfile(caller_frame))
    
    prefixes = {
        'trigger': '[Flask Trigger]',
        'attempt': '[Attempt]',
        'success': '[Success]',
        'note': '[Note]',
        'fail': '[Fail]',
        'warning': '[WARNING]',
        'delimiter': ''       
    }

    if message_type in prefixes:
        type_prefix = prefixes[message_type]
    else:
        message_type = 'unknown'
        type_prefix = '[Unknown Type]'


    prefix = f"{log_handler_timestamp} {type_prefix}"

    if message_type == 'delimiter':
        message = print_delimiter(message)
    else:        
        message = f"{prefix} (<{script_name}> Function:{function_name}) {message}"

    if message_type == 'note' or 'trigger' or 'attempt' or 'success' or 'unknown':
        logger.info(message)
        print(message)
    elif message_type == 'fail':
        logger.critical(message)
        print(message)
    elif message_type == 'warning':
        logger.warning(message)
        print(message)
    elif message_type =='delimiter':
        logger.info(message)
        print(message)

def print_delimiter(message):
    delimiter_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
    delimiter_message = (f'''
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /||\|||\ /|||\|||\ /||
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||/ \|||\|||/ \|||\|||
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `~   `-~ `-`   `-~ `-`
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
{delimiter_timestamp} {message}
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
        ''')
    return delimiter_message
