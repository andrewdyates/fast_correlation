#!/usr/bin/python
"""Save a copy of a .tab matrix as log scaled.

EXAMPLE USE:
  python log_scale.py mytab.tab
"""
import numpy as np
import numpy.ma as ma
from util import *
import sys

def tostr(x):
  if x == 0 or x:
    return "%f" % x
  else:
    return ""

def main(tab_fname=None):
  varlist = []
  M = np.genfromtxt(name_iter(open(tab_fname), varlist), usemask=True, delimiter='\t')
  Q = ma.MaskedArray(data=M.data, mask=(M.mask|(M.data == 0)))
  Q = ma.log(Q)
  # save matrix back to .tab format
  fp = open(tab_fname + ".logscale.tab", "w")
  for i, row in enumerate(Q):
    fp.write(varlist[i] + '\t')
    fp.write('\t'.join(map(tostr, row)))
    fp.write('\n')
  fp.close()

if __name__ == "__main__":
  main(sys.argv[1])
