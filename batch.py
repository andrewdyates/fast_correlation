#!/usr/bin/python
"""Compute a batch of pairwise relations."""
from __future__ import division
import numpy as np
from scipy.stats import mstats
import numpy.ma as ma # masked array for missing values
import datetime
from py_symmetric_matrix import *
from qsub import *
from dcor import *
import os

LOG_MSG = "#npy_fname=%(npy_fname)s, function=%(function)s, start=%(start)d, end=%(end)d, m=%(m)d, date=%(date)s"
REPORT_N = 1000

def euclidean(x,y):
  q=x-y
  return ma.sqrt((q*q.T).sum())


# this should be in a separate file
FUNCTIONS = {
  'pearson': lambda x, y: mstats.pearsonr(x,y)[0],
  'spearman': lambda x, y: mstats.spearmanr(x,y)[0],
  'euclidean': euclidean,
  'kendalltau': lambda x,y: mstats.kendalltau(x,y)[0],
  'dcor': dcor,
  }

def main(npy_fname=None, function=None, batchname=None, outdir=None, start=None, end=None, m=None):
  """Compute pairs of dependency"""
  assert npy_fname, function
  assert function in FUNCTIONS
  assert os.path.exists(outdir)

  m = int(m)
  assert m > 0

  if end is None: 
    end = m*(m-1) / 2
  else:
    end = int(end)
  if start is None:
    start = 0
  else:
    start = int(start)
  assert start < end, start >= 0

  M = ma.load(npy_fname)

  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npy_fname), function, start, end)

  log_msg = LOG_MSG % {'npy_fname': npy_fname, 'function': function, 'start': start,
  'end': end, 'm': m, 'date': datetime.datetime.now().isoformat(' ')}
  fp_out = open(os.path.join(outdir, batchname+".txt"), 'w')
  fp_out.write(log_msg + "\n")
  print "Started job...", log_msg
  
  f = FUNCTIONS[function]
  for i in xrange(start, end):
    if i % REPORT_N == 0:
      print "Generating pair %d (to %d) in %s..." % \
        (i, end-1, batchname)
    x, y = inv_sym_idx(i, m)
    v = f(M[x], M[y])
    fp_out.write("%.6f\n" % v)
  fp_out.write("#end")

if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
