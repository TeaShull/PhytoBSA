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

class LogHandler:
    def __init__(self, name, log_filename):
        self.name = name 
        log_dir = LOG_DIR
        self.log_path = os.path.join(log_dir, log_filename)
        self.logger = self.setup_logger()

    def setup_logger(self):
        logger = logging.getLogger(self.name)
        logger.setLevel(logging.INFO)  # Set the default level for the logger

        # Create a FileHandler and set the level to INFO
        file_handler = logging.FileHandler(self.log_path)
        file_handler.setLevel(logging.INFO)

        # Create a formatter and add it to the handler
        formatter = logging.Formatter('%(message)s')
        file_handler.setFormatter(formatter)

        # Add the handler to the logger
        logger.addHandler(file_handler)

        return logger

    def obtain_execution_frames(self):
        current_frame = inspect.currentframe()
        caller_frame = inspect.getouterframes(current_frame)[2]
        
        # Get function name and script name from one level higher
        function_name = caller_frame[0].f_code.co_name
        script_name = os.path.basename(caller_frame[1])
        
        return script_name, function_name

    def construct_message(self, prefix, script_name, function_name, message_in):
        log_handler_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
        if function_name != '<module>':
            message_out = f"{prefix} (<{script_name}> Function:{function_name}) {message_in}"
        else:
            message_out = f"{prefix} (<{script_name}> :module:) {message_in}"
        return message_out 

    def trigger(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[Flask Trigger]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def attempt(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[Attempt]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def success(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[Success]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def note(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[Note]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def fail(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[!-FAIL-!]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def warning(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[!-WARNING-!]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def bash(self, message):
        script_name, function_name = self.obtain_execution_frames()
        log_message = self.construct_message(
            '[sh]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def delimiter(self, message):
        delimiter_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
        log_message = (f'''
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .-. .-.   .
||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /||\|||\ /|||\|||\ /||
|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||/ \|||\|||/ \|||\|||
~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-`   `-~ `-`   `-~ `~   `-~ `-`   `-~ `-`
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
{delimiter_timestamp} {message}
>=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-<
            ''')
        self.logger.info(log_message)
        print(log_message)    

