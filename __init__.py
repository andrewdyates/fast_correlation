#!/usr/bin/python
"""Efficiently compute all pairs correlation and save residuals."""
from __future__ import division
import numpy as np
from py_symmetric_matrix import *
from scipy.spatial.distance import squareform


def correlate_all(M):
  """Return all-pairs Pearson's correlation as a squareform matrix. 
  Best on numpy.array(dtype=float)

  Args:
    M: numpy.array row matrix
  Returns:
    squareform top triangle matrix of all-pairs correlation, row order index.

  RUNTIME on random.rand(500,200)
    21.2 ms (improve of 200x over formula)
  RUNTIME on random.rand(15000,250)
  """
  m = np.size(M, 0) # number of rows (variables)
  n = np.size(M, 1) # number of columns (power)

  sums = np.sum(M,1).reshape(m,1)
  stds = np.std(M,1).reshape(m,1) # divided by n

  # TODO: does making this cummlative matter?
  Dot = squareform(np.dot(M, M.T), checks=False)
  SumProd = squareform(np.dot(sums, sums.T), checks=False)
  StdProd = squareform(np.dot(stds, stds.T), checks=False)

  CorrMatrix = (Dot - (SumProd/n)) / (StdProd*n)
  
  # Correlation Matrix
  return CorrMatrix

def cumm_correlate_all(M):
  """Return all-pairs Pearson's correlation as a squareform matrix. 
  Best on numpy.array(dtype=float)
  Is this a performance improvement?

  Args:
    M: numpy.array row matrix
  Returns:
    squareform top triangle matrix of all-pairs correlation, row order index.

  RUNTIME on random.rand(500,200)
    21.2 ms (improve of 200x over formula)
  RUNTIME on random.rand(5000,250)
    3.02 s per loop
  """
  m = np.size(M, 0) # number of rows (variables)
  n = np.size(M, 1) # number of columns (power)

  sums = np.sum(M,1).reshape(m,1)
  stds = np.std(M,1).reshape(m,1) # divided by n

  # TODO: does making this cummlative matter?
  C = squareform(np.dot(M, M.T), checks=False)
  C = C - (squareform(np.dot(sums, sums.T), checks=False)/n)
  C = C / (squareform(np.dot(stds, stds.T), checks=False)*n)
  CorrMatrix = C
  
  # Correlation Matrix
  return CorrMatrix


def get_ranks(M):
  """Convert matrix to ranks. (use in Spearman's)

  Args:
    M: numpy.array row matrix
  Returns:
    numpy.array row matrix of corresponding ranks.

  RUNTIME on random.rand(500,200)
    7.29 ms
  RUNTIME on random.rand(15000,250)
    317 ms 
  """
  return M.argsort(axis=1).argsort(axis=1)


def load_file(filename, n=None, delimiter="\t"):
  """Load labeled row matrix from file."""
  pass

def biologist_correlate(M):
  """Manually compute pearson's correlation for each pair of variables.
  Used for benchmarking purposes only. (this is slow)

  RUNTIME on random.rand(500,200)
    21.7 s per loop
  """ 
  from scipy.stats import pearsonr
  m = np.size(M, 0) # number of rows (variables)

  CorrMatrix = np.zeros(m*(m-1)/2)
  for i in xrange(m):
    for j in xrange(i+1,m):
      idx = sym_idx(i,j,m)
      CorrMatrix[idx] = pearsonr(M[i], M[j])[0]
  return CorrMatrix

def formula_correlate(M):
  """Calculate correlations without optimized dot products.
  Used for benchmarking purposes only.

  Dot product calculation times: 
    16.5ms: reduce(lambda acc,s: acc+s[0]*s[1], zip(M[i],M[j]), 0)
    11.4 ms: np.sum([x*y for x,y in zip(M[i],M[j])])  #numpy sum
    10.6 us: dot(M[i],M[j])

  RUNTIME on random.rand(500,200)
    4.24 s
  """
  m = np.size(M, 0) # number of rows (variables)
  n = np.size(M, 1) # number of columns (power)
  CorrMatrix = np.zeros(m*(m-1)/2)
  sums = np.sum(M,1).reshape(m,1)
  stds = np.std(M,1).reshape(m,1) # divided by n
  
  for i in xrange(m):
    for j in xrange(i+1,m):
      idx = sym_idx(i,j,m)
      dot_product = np.dot(M[i], M[j])
      CorrMatrix[idx] = (dot_product - sums[i]*sums[j]/n) / stds[i]*stds[j]*n
  return CorrMatrix


