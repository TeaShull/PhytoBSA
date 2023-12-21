#!/usr/bin/env python
import os
import inspect
from datetime import datetime

BASE_DIR = os.getcwd()
SRC_DIR =  os.path.join(BASE_DIR, 'src')
INPUT_DIR = os.path.join(BASE_DIR, 'input')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
LOG_DIR = os.path.join(SRC_DIR, 'logs')
TEMPLATE_DIR = os.path.join(SRC_DIR,'templates')
STATIC_DIR = os.path.join(SRC_DIR,'static')
MODULES_DIR = os.path.join(SRC_DIR,'modules')


# General Functions
def error_handler(error_type, message):
    error_handler_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
    function_name = inspect.currentframe().f_back.f_code.co_name
    caller_frame = inspect.currentframe().f_back
    script_name = os.path.basename(inspect.getfile(caller_frame))
    prefixes = {
        'trigger': '[Flask Trigger]',
        'attempt': '[Attempt]',
        'success': '[Success]',
        'fail': '[Fail]'
    }

    if error_type in prefixes:
        type_prefix = prefixes[error_type]
    else:
        type_prefix = '(Unknown Type)'

    prefix = f"{error_handler_timestamp} {type_prefix}"
    if function_name != '<module>':
        error_message = f"{prefix} (Function:{function_name}) {message}"
    else:
        error_message = f"{prefix} (Core logic of:{script_name}) {message}"
    print(error_message)

def print_delimiter(message):
    delimiter_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
    print(f'''
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /||\|||\ /|||\|||\ /||
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||/ \|||\|||/ \|||\|||
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `~   `-~ `-`   `-~ `-`
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
{delimiter_timestamp} {message}
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
        ''')
