# This file was updated from original version in CarveMe by Adam Bates on 7/7/26

import os
from configparser import ConfigParser

__version__ = '2.0.0a'

project_dir = os.path.abspath(os.path.dirname(__file__)) + '/'

config = ConfigParser()
config.read(project_dir + 'config.cfg')
