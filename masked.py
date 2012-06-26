#!/usr/bin/python
"""Handle masked arrays: import, correlate, other stats.

TODO: 
  make own package
  add symmetric array handling
  break problem into subproblems for parallel processing

SEE:

- import masked array from text
http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html#numpy.genfromtxt

- masked scipy statistics
http://docs.scipy.org/doc/scipy/reference/stats.mstats.html#module-scipy.stats.mstats

- masked numpy functions
http://docs.scipy.org/doc/numpy/reference/maskedarray.baseclass.html#maskedarray-baseclass
http://docs.scipy.org/doc/numpy/reference/routines.ma.html

all pairs spearmanr
all pairs euclidean


USE R!
http://math.furman.edu/~dcs/courses/math47/R/library/Hmisc/html/rcorr.html
"""
from __future__ import division
import os
import numpy as np
import numpy.ma as ma # masked array for missing values
import os
from scipy.stats import mstats


#TEST_FILE=os.path.expanduser("~/Dropbox/biostat/eqtl_data/GSE25935/GSE25935.GPL4133.eQTL.nooutliers.tab")
TEST_FILE=os.path.expanduser("~/Desktop/test.tab")

# this does not (yet) work
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html#numpy.genfromtxt
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html#numpy.loadtxt

# all pairs spearman?
# mstats.spearmanr(M[0], M[1])

# handles comments correctly, but not variable names
def name_iter(fp, varlist):
  """Load a labeled row matrix as a line iterator."""
  for line in fp:
    if line[0] == '#': continue
    name,c,row = line.partition('\t')
    varlist.append(name)
    yield row
  
def correlate_all(M):
  """Return all pairs correlation matrix.

  Note: dot product optimizations cannot be as readily applied
  due to different pairs discarded due to missing values per pair.
  """
  return ma.corrcoef(M)


def all_pairs_spearman(M):
  """This should return a squareform matrix"""
  C = np.zeros((len(M), len(M)))
  for i in xrange(len(M)):
    for j in xrange(i+1, len(M)):
      C[i][j] = mstats.spearmanr(M[i],M[j])[0]
  return C


def all_pairs_pearson(M):
  """This should return a squareform matrix.

  This is about 15% faster than "correlate all", but
  correlate_all calcuates twice the extra values.
  """
  C = np.zeros((len(M), len(M)))
  for i in xrange(len(M)):
    for j in xrange(i+1, len(M)):
      C[i][j] = mstats.pearsonr(M[i],M[j])[0]
  return C

def all_pairs_euclidean(M):
  C = np.zeros((len(M), len(M)))
  for i in xrange(len(M)):
    for j in xrange(i+1, len(M)):
      q=M[i]-M[j]
      C[i][j] = ma.sqrt((q*q.T).sum())
  return C


def run():
  varlist = []
  M = np.genfromtxt(name_iter(open(TEST_FILE), varlist), usemask=True, delimiter='\t')

  return varlist, M




if __name__ == "__main__":
  run()
