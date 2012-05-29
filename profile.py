#!/usr/bin/python
import cProfile
from __init__ import *
from numpy import *
from scipy.spatial.distance import squareform
# get "sym_idx" function
from py_symmetric_matrix import *

def main():
  M = random.rand(100,10000)
