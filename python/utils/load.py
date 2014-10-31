# -*- coding: UTF-8 -*-
import os, sys

_PIGASUS_PLUGIN_DIR = "PIGASUS_PLUGIN_DIR"

def load(name):
    try:
        path = os.environ[_PIGASUS_PLUGIN_DIR]
    except KeyError:
        print("Error: Please set the environment variable ", _PIGASUS_PLUGIN_DIR)
        sys.exit(1)
    # ... add _PIGASUS_PLUGIN_DIR into PYTHONPATH
    sys.path.insert(0,path)
    # ...
    try:
        mod = __import__("%s" % name)
        return mod
    except:
        print("Error: the module ", name \
                , " is not installed. Please download and save it in " \
                , path \
                , " for any additional information, contact ratnaniahmed@gmail.com")
        sys.exit(1)
