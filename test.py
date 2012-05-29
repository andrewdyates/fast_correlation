import unittest
from __init__ import *
from numpy import *
from scipy.spatial.distance import squareform
# get "sym_idx" function
from py_symmetric_matrix import *

# Squareform overhead
#%timeit squareform(X, checks=False)
#1 loops, best of 3: 676 ms per loop


# use numpy.ma.MaskedArray for missing values. 1=MASK, 0=TRANSPARENT
# use numpy.ma.MaskedArray for missing values. 1=MASK, 0=TRANSPARENT

# In [12]: Q = arange(10000,dtype=float).reshape(10000,1)
# In [13]: %timeit dot(Q,Q.T)
# 1 loops, best of 3: 1.7 s per loop

# %timeit test(10000)
# 1 loops, best of 3: 3.2 s per loop


class TestMatrix(unittest.TestCase):

  def test_small(self):
    M = arange(8*10, dtype=float).reshape(8,10)
    C = correlate_all(M)
    self.assertTrue(all(C >= 0.999))
    self.assertTrue(all(C < 1.001))

  def test_random(self):
    M = random.rand(100,10000)
    C = correlate_all(M)
    self.assertTrue(all(C**2 >= 0))
    # These may fail... but it's very unlikely for n=10000.
    self.assertTrue(mean(C) <= 0.001)
    self.assertTrue(std(C) <= 0.1)
    self.assertTrue(max(C) <= 1.0)
    self.assertTrue(min(C) >= -1.0)

  def test_big(self):
    M = arange(10000*200, dtype=float).reshape(10000,200)
    C = correlate_all(M)
    self.assertTrue(all(C >= 0.999))
    self.assertTrue(all(C < 1.001))

if __name__ == "__main__":
  unittest.main()

# # sum of x squared, square of sum y (compress these to save memory)
# D_xsq_sqy = np.dot(D_sq_sums, (D_sums.T)**2)
# D_ysq_sqx = np.dot(D_sums**2, D_sq_sums.T)
# D_midterm = (D_xsq_sqy + D_ysq_sqx) / n

# numerator = D_dot - D_sums_prod_n
# denominator = D_sq_sums_prod - D_midterm + D_sums_prod_n**2

# pearsonsr = numerator/denominator**0.5





