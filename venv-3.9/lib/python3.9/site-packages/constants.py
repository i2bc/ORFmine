"""
constant.Constants - The simple way to deal with environment constants.
"""

import os
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import warnings
import logging
import functools


logger = logging.getLogger(__name__)

def debug(function):
    """
    logging debug decorator
    """
    @functools.wraps(function)
    def wrapper(*args, **kvargs):
        """
        wrap method call
        """
        logger.debug('begin %s %s',
                     args,
                     kvargs,
                     extra={'method': function.__name__})
        result = function(*args, **kvargs)
        logger.debug('end %s', result, extra={'method': function.__name__})
        return result
    return wrapper

VARIABLE = '__CONSTANTS__'
FILENAME = 'constants.ini'

class Constants(object):
    """
    Environement sensitive application constants class
    """

    @debug
    def __init__(self, variable=VARIABLE, filename=FILENAME):
        """
        variable is the name of the environment variable to read the
        environment / config section from default to __CONSTANTS__
        filename is the config filename
        """
        object.__setattr__(self, 'dict', {})
        self.variable = variable
        self.filename = filename
        self.load()

    @debug
    def load(self):
        """
        load the section self.variable from the config file self.filename
        """
        self.get_environment()
        self.read_config()
        self.load_dict()

    @debug
    def get_environment(self):
        """
        returns the value of the environment variable self.variable
        """
        self.environment = os.environ[self.variable]

    @debug
    def read_config(self):
        """
        returns a ConfigParser instance from self.filename
        """
        self.config = configparser.ConfigParser()
        with open(self.filename) as config_file:
            self.config.readfp(config_file)

    @debug
    def load_dict(self):
        """
        load the config items into self.dict
        """
        self.dict = dict (self.config.items(self.environment))
        logger.info('variable: %s, filename: %s, environment: %s, constants: %s',
                    self.variable,
                    self.filename,
                    self.environment,
                    self.dict,
                    extra={'method': 'load'})

    @debug
    def __getitem__(self, item):
        """
        access to environment specific constants in a dictionary manner
        casts to int, float or keep as string
        """
        return self.cast(self.dict[item])

    @debug
    def __getattr__(self, item):
        """
        syntactic sugar, .item rather than ['item']
        """
        return self[item]

    @staticmethod
    @debug
    def cast(constant):
        """
        cast string to int, float, eval or keep as string
        """
        if hasattr(constant, 'startswith') and constant.startswith('0') \
           and '.' not in constant:
            return constant
        try:
            return int(constant)
        except ValueError:
            pass
        try:
            return float(constant)
        except ValueError:
            pass
        try:
            return eval(constant)
        except (NameError, SyntaxError):
            pass
        return constant

    @debug
    def __setitem__(self, item, value):
        """
        dict like assignment - warns when a constant is changed
        """
        if item in self.dict:
            warnings.warn('{0} changed to {1}'.format(item, value))
        self.dict[item] = value

    @debug
    def __setattr__(self, name, value):
        """
        attribute assignment - warns when a constant is changed
        """
        if hasattr(self, 'dict') and name in self.dict:
            warnings.warn('{0} changed to {1}'.format(name, value))
            self.dict[name] = value
        else:
            object.__setattr__(self, name, value)
