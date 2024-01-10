import os
import sys
import time
import codecs

import inspect
from datetime import datetime
import logging
import utilities_ulid as ulid
import sqlite3


class ULID:
    # Python class port of https://github.com/alizain/ulid
    # https://github.com/mdipierro/ulid
    # License MIT
    ENCODING = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"
    LENCODING = len(ENCODING)
    
    def __init__(self):
        self.PY3 = sys.version_info[0] == 3
        
    def encode_time_10bytes(self, x):
        s = ''
        while len(s) < 10:
            x, i = divmod(x, self.LENCODING)
            s = self.ENCODING[i] + s
        return s
    
    def encode_random_16bytes(self):
        b = os.urandom(10)
        x = int(codecs.encode(b, 'hex') if self.PY3 else b.encode('hex'), 16)
        s = ''
        while len(s) < 16:
            x, i = divmod(x, self.LENCODING)
            s = self.ENCODING[i] + s
        return s
    
    def convert(self, chars):
        i = 0
        n = len(chars)-1
        for k, c in enumerate(chars):
            i = i + 32**(n-k) * self.ENCODING.index(c)
        return i
    
    def seconds(self, ulid):
        """ return the timestamp from a ulid """
        return 0.001*self.convert(ulid[:10])
    
    def sharding(self, ulid, partitions):
        """ return a sharting partition where to store the ulid"""
        return self.convert(ulid[-16:]) % partitions
    
    def generate_ulid(self):
        timestamp = int(time.time()*1000)
        encoded_time = self.encode_time_10bytes(timestamp)
        encoded_random = self.encode_random_16bytes()
        return encoded_time + encoded_random

class LogHandler:
    def __init__(self, name, db_name="thale_bsa_sqldb.db"):
        # Initialize log file parameters
        self.name = name 
        self.init_timestamp = datetime.now().strftime("%Y.%m.%d-%H.%M")
        
        ulid_generator = ULID()
        self.ulid = ulid_generator.generate_ulid()
        
        # Construct log filename and spin up logger
        self.log_filename = f'{self.init_timestamp}_{self.name}_{self.ulid}.log'
        self.log_path = os.path.join(LOG_DIR, self.log_filename)
        self.logger = self.setup_logger()

        # Initialize log database for storing log parameters between runs
        self.db_name = db_name
        self.conn = sqlite3.connect(self.db_name)
        self._create_tables()

    def setup_logger(self):
        '''Initialize logger, which is designed to be passed to class instances 
        so that logging can be passed around to new functionalities as the program expands.'''
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

    def _obtain_execution_frames(self)->tuple:
        '''
        Obtains the execution frames from the script that the self.log(message) 
        is called. This allows logs to be informative as to what module, and what
        function a message arises from.  

        Args:
        None

        Returns:
        script_name, function_name (tuple)
        '''
        
        current_frame = inspect.currentframe()
        caller_frame = inspect.getouterframes(current_frame)[2]
        
        # Get function name and script name from one level higher
        function_name = caller_frame[0].f_code.co_name
        script_name = os.path.basename(caller_frame[1])
        
        return script_name, function_name

    def _construct_message(
        self, prefix, script_name, function_name, message_in
        )->str:
        '''
        Constructs the log messages. 

        Args:
        prefix(str) - What is the prefix of the message? example [WARNING]
        script_name(str) - what is the script name obtained from obtain_execution_frames? 
        function_name(str) - what is the function name obtained from obtain_execution_frames? 
        message_in(str) - what is the log message? 
            
        Returns: 
        message_out(str)
        '''
        
        log_handler_timestamp = datetime.now().strftime("%Y.%m.%d ~%H:%M")
        if function_name != '<module>':
            message_out = f"{log_handler_timestamp} {prefix} ({script_name}>{function_name}) {message_in}"
        else:
            message_out = f"{log_handler_timestamp} {prefix} ({script_name}> module) {message_in}"
        return message_out 

    def trigger(self, message):
        '''
        log message type. 
        will log and print message [Flask Trigger] "....exc"
        
        all log message types below are similar in function. 
        Code could be condensed. someday 

        Args:
        message(str)

        Returns:
        None. just prints to the log and stdOut
        '''
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[Flask Trigger]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def attempt(self, message):
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[Attempt]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def success(self, message):
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[Success]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def note(self, message):
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[Note]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def fail(self, message):
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[!-FAIL-!]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)
        quit()

    def warning(self, message):
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[!-WARNING-!]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def bash(self, message):
        script_name, function_name = self._obtain_execution_frames()
        log_message = self._construct_message(
            '[sh]', script_name, function_name, message
        )
        self.logger.info(log_message)
        print(log_message)

    def print(self,message):
        self.logger.info(message)
        print(message)

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

    # Log database functions
    def _create_tables(self):
        '''
        Creates the tables for the log database.

        Args:
        None

        Returns:
        None.
        '''
        create_core = '''
            CREATE TABLE IF NOT EXISTS core (
                core_ulid TEXT,
                core_log_path TEXT,
                core_timestamp TEXT,
                PRIMARY KEY(core_ulid)
            )
        '''

        create_vcf = '''
            CREATE TABLE IF NOT EXISTS vcf (
                vcf_ulid TEXT,
                line_name TEXT,
                core_ulid TEXT,
                vcf_log_path TEXT,
                vcf_timestamp TEXT,
                PRIMARY KEY (vcf_ulid, line_name)
            )
        '''

        create_analysis = '''
            CREATE TABLE IF NOT EXISTS analysis (
                analysis_ulid TEXT,
                line_name TEXT,
                core_ulid TEXT,
                vcf_ulid TEXT,
                analysis_log_path TEXT,
                analysis_timestamp TEXT,
                PRIMARY KEY (analysis_ulid, line_name)
            )
        '''

        try:
            self.conn.execute(create_core)
            self.conn.execute(create_vcf)
            self.conn.execute(create_analysis)
            self.conn.commit()
        
        except Exception as e:
            print(f'There was an error during table creation: {e}')


    def add_db_record(
        self, current_line_name=None, core_ulid=None, vcf_ulid=None
    ):
        '''
        Adds records to the log database. 

        Args: 
        for 'core', 
            None  

        for 'vcf', 
            current_line_name(str)
            core_ulid are necissary(str) 

        for 'analysis', 
            current_line_name(str)
            core_ulid(str) 
            vcf_ulid(str) 

        Returns: 
            None
        It just adds logs with the inputted information, along with information
        gathered from the current instance of the LogHandler class
        '''

        try:
            if self.name == 'core':
                add = '''
                    INSERT INTO core (core_ulid, core_log_path, core_timestamp)
                    VALUES(?, ?, ?)
                '''
                values = (self.ulid, self.log_path, self.init_timestamp)
            
            elif self.name == f'vcf_{current_line_name}':
                add = '''
                    INSERT INTO vcf (vcf_ulid, line_name, core_ulid, 
                        vcf_log_path, vcf_timestamp
                    )
                    VALUES (?, ?, ?, ?, ?)
                '''
                values = (self.ulid, current_line_name, core_ulid, 
                    self.log_path, self.init_timestamp
                    )
            
            elif self.name == f'analysis_{current_line_name}':
                add = '''
                    INSERT INTO analysis (
                        analysis_ulid, line_name, core_ulid, 
                        vcf_ulid, analysis_log_path, analysis_timestamp
                    )
                    VALUES (?, ?, ?, ?, ?, ?)
                '''
                values = (
                    self.ulid, current_line_name, core_ulid, 
                    vcf_ulid, self.log_path, self.init_timestamp
                )
            
            else:
                raise ValueError("Invalid table name")

            self.conn.execute(add, values)
            self.conn.commit()
        
        except Exception as e:
            print(f'There was an error adding entries to database:{e}')
