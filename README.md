fast_correlation
================

Efficiently compute Pearson's and Spearman's correlations and cache residuals for fast value updates.

Not as fast as Matlab's corr(A), but faster than scipy.spatial.distance.pdist(X, 'correlation') by a few seconds.