#!/usr/bin/python
"""Dispatch all-pairs relation job batches.

EXAMPLE USE:

python dispatch.py tab_fname=~/Dropbox/biostat/eqtl_data/GSE2034/GSE2034.GPL96.eQTL.normed.tab outdir=~/Desktop/GSE2034_test function=pearson

Skip first several jobs, dryrun
time python $HOME/fast_correlation/dispatch.py tab_fname=/nfs/01/osu6683/gse2034/GSE2034.GPL96.eQTL.normed.tab outdir=/nfs/01/osu6683/gse2034/gse2034_pairscores k=50000 function=dcor job_offset=1100 dry=True
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


def main(tab_fname=None, outdir=None, function=None, k=500000, dry=False, job_offset=None):
  assert all((tab_fname, outdir, function))
  tab_fname = os.path.expanduser(tab_fname)
  outdir = os.path.expanduser(outdir)
  assert os.path.exists(tab_fname)
  if function:
    assert function in batch.FUNCTIONS
  
  k = int(k)
  assert k > 1
  if job_offset is not None:
    job_offset = int(job_offset)
    assert job_offset >= 0
  
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

  # Only re-create varlist and numpy matrix if they do not yet exist.
  varlist_fname = os.path.join(outdir, os.path.basename(tab_fname) + ".varlist.txt")
  npy_fname = os.path.join(outdir, "%s.npy" % tab_fname)
  if os.path.exists(varlist_fname) and os.path.exists(npy_fname):
    print "Both %s and %s exist, do not recreate varlist and numpy masked matrix files." % \
        (varlist_fname, npy_fname)
    # load numpy matrix to get its size
    m = np.size(ma.load(npy_fname), 0)
  else:
    # import tab
    varlist = []
    # M is masked matrix
    print "Loading %s into masked numpy matrix and varlist..." % tab_fname
    M = np.genfromtxt(name_iter(open(tab_fname), varlist), usemask=True, delimiter='\t')
    m = np.size(M, 0) # number of rows (variables)
  
    # save to file
    print "Saving matrix and varlist..."
    fp = open(varlist_fname, "w")
    fp.write('\n'.join(varlist))
    fp.close()
    ma.dump(M, npy_fname)
  
  # dispatch jobs in a loop
  i = 0
  num_pairs = int(m * (m-1) / 2) # no diagonal: n choose 2
  params = {
    'npy_fname': os.path.abspath(npy_fname),
    'batchname': None, # use default batch name
    'm': m,
    }
  print "Changing current working directory to %s." % os.path.abspath(outdir)
  os.chdir(outdir)

  n_jobs = 0
  n_skipped = 0
  while i < num_pairs:
    
    if job_offset is not None:
      if n_jobs < job_offset:
        print "Skipping job %d... %d offset, %d skipped" % \
            (n_jobs, job_offset, n_skipped)
        n_skipped += 1
        n_jobs += 1
        i += k
        continue

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
      job_name = "%s_%s_%d_to_%d_of_%d" % \
        (base_job_name, function, params['start'], params['end'], num_pairs)
      # Generate qsub job submission script
      script_txt = fill_template(jobname=job_name, script=cmd, **params)
  
      # Submit job
      print script_txt
      if not dry:
        p = subprocess.Popen("qsub", stdin=subprocess.PIPE)
        p.communicate(input=script_txt)
        p.stdin.close()
      
    # increment pairs counter (after submitting for all functions)
    i += k
    n_jobs += 1

  print "%d jobs (%d skipped) submitted at %s." % (n_jobs, n_skipped, timestamp())
  print "expected %d jobs." % math.ceil(num_pairs/k)
  
  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
