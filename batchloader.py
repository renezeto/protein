#!/usr/bin/python
from __future__ import division
import numpy as np
import sys
import os
import subprocess
import time

#parses sys.argv for key=value arguments
def parseKwargs(stringList):
    stringDict = {}
    for arg in stringList:
        if "=" in arg:
            (key,value) = arg.split("=")
            stringDict[key]=value
    return stringDict




