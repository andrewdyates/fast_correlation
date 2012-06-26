#!/usr/bin/python
import datetime

def timestamp(date=None):
  if not date:
    date = datetime.datetime.now()
  return date.isoformat(' ')

def name_iter(fp, varlist):
  """Load a labeled row matrix as a line iterator."""
  for line in fp:
    if line[0] == '#': continue
    name,c,row = line.partition('\t')
    varlist.append(name)
    yield row

    
