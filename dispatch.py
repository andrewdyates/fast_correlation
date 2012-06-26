#!/usr/bin/python
"""Dispatch all-pairs relation job batches.

EXAMPLE USE:

python dispatch.py tab_fname=~/Dropbox/biostat/eqtl_data/GSE2034/GSE2034.GPL96.eQTL.normed.tab outdir=~/Desktop/GSE2034_test function=pearson
"""
from __future__ import division
import subprocess
import numpy as np
from scipy.stats import mstats
import numpy.ma as ma # masked array for missing values
import datetime
from py_symmetric_matrix import *
import batch
import os
import sys
import errno
from util import *
from qsub import *
import math


# Assume that batch.py is in the same directory as this script.
DIR = os.path.dirname(os.path.abspath(__file__))
CMD = "time python %s" % os.path.join(DIR, 'batch.py')


def main(tab_fname=None, outdir=None, function=None, k=100000):
  assert all((tab_fname, outdir, function))
  tab_fname = os.path.expanduser(tab_fname)
  outdir = os.path.expanduser(outdir)
  assert os.path.exists(tab_fname)
  if function:
    assert function in batch.FUNCTIONS
  
  k = int(k)
  assert k > 1
  
  # create outdir
  try:
    os.makedirs(outdir)
  except OSError, e:
    if e.errno != errno.EEXIST: raise

  # run for all functions or only a select function
  if function is None:
    ALL_FUNCTIONS = batch.FUNCTIONS
  else:
    ALL_FUNCTIONS = {function: batch.FUNCTIONS[function]}
    
  # for each function, create subdirs
  outdirs = {}
  for function in ALL_FUNCTIONS:
    path = os.path.join(outdir, function)
    try:
      os.makedirs(path)
    except OSError, e:
      if e.errno != errno.EEXIST: raise
    outdirs[function] = os.path.abspath(path)
      
  # import tab
  varlist = []
  # M is masked matrix
  print "Loading %s into masked numpy matrix and varlist..." % tab_fname
  M = np.genfromtxt(name_iter(open(tab_fname), varlist), usemask=True, delimiter='\t')
  m = np.size(M, 0) # number of rows (variables)
  
  # save to file
  print "Saving matrix and varlist..."
  fp = open(tab_fname + ".varlist.txt", "w")
  fp.write('\n'.join(varlist))
  fp.close()
  npy_fname = "%s.npy" % tab_fname
  ma.dump(M, npy_fname)
  
  # dispatch jobs in a loop
  i = 0
  num_pairs = int(m * (m-1) / 2)
  params = {
    'npy_fname': os.path.abspath(npy_fname),
    'batchname': None, # use default batch name
    'm': m,
    }
  print "Changing current working directory to %s." % os.path.abspath(outdir)
  os.chdir(outdir)
  
  while i < num_pairs:
    for function in ALL_FUNCTIONS:
      params.update({
        'function': function, 
        'start': i,
        'end': i+k,
        'outdir': outdirs[function],
        })
      if params['end'] > num_pairs:
        params['end'] = num_pairs
      cmd = CMD + " " + " ".join(["%s=%s" % (key,str(v)) for key,v in params.items()])
      base_job_name = os.path.basename(npy_fname)
      job_name = "%s_%d_to_%d_of_%d" % \
        (base_job_name, params['start'], params['end'], num_pairs)
      # Generate qsub job submission script
      script_txt = fill_template(jobname=job_name, script=cmd, **params)
  
      # Submit job
      print script_txt
      p = subprocess.Popen("qsub", stdin=subprocess.PIPE)
      p.communicate(input=script_txt)
      p.stdin.close()
      
    # increment pairs counter (after submitting for all functions)
    i += k

  print "%d jobs submitted at %s." % (math.ceil(m/i), timestamp())
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
