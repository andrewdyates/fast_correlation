#!/usr/bin/python
"""Efficiently compute all pairs correlation and save residuals."""
import numpy as np
from py_symmetric_matrix import *
from scipy.spatial.distance import squareform

def load_file(filename, n=None, delimiter="\t"):
  """Load labeled row matrix from file."""
  pass


def correlate_all(M, ranks=False):
  m = np.size(M, 0) # number of rows (variables)
  n = np.size(M, 1) # number of columns (power)

  sums = np.sum(M,1).reshape(m,1)
  stds = np.std(M,1).reshape(m,1) # divided by n

  Dot = squareform(np.dot(M, M.T), checks=False)
  SumProd = squareform(np.dot(sums, sums.T), checks=False)
  StdProd = squareform(np.dot(stds, stds.T), checks=False)

  CorrMatrix = (Dot - (SumProd/n)) / (StdProd*n)
  
  # Correlation Matrix
  return CorrMatrix

