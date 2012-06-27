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
import errno
import pwd
import shutil

LOG_MSG = "#npy_fname=%(npy_fname)s, function=%(function)s, start=%(start)d, end=%(end)d, m=%(m)d, date=%(date)s"
REPORT_N = 1000
# get username
TMP_DIR = "/tmp/%s" % pwd.getpwuid(os.getuid()).pw_name

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
  assert os.path.isdir(outdir)

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

  # if tmp directory does not exist, create it
  try:
    os.makedirs(TMP_DIR)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

  M = ma.load(npy_fname)

  if batchname is None or batchname in ("None", "NONE", "none"):
    batchname = "%s_%s_%d_%d" % \
      (os.path.basename(npy_fname), function, start, end)

  log_msg = LOG_MSG % {'npy_fname': npy_fname, 'function': function, 'start': start,
  'end': end, 'm': m, 'date': datetime.datetime.now().isoformat(' ')}
  # Create output file in temporary directory
  output_fname = os.path.join(TMP_DIR, batchname+".txt")
  fp_out = open(output_fname, 'w')
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
  fp_out.close() # be sure to close the file to write everything to disk!
  # Copy temporary file to outdir
  shutil.move(output_fname, outdir)

if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
