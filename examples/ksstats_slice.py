"""Sequential moments matching Appendix B of the manuscript."""
from __future__ import annotations

import numpy as np

from ksstats import KSMoments, KSPairs

rng = np.random.default_rng(0)
n = 8
m = 3
mom = KSMoments(n)
pair = KSPairs(m)
for _ in range(40):
    field = rng.normal(size=n)
    series_at_t = rng.normal(size=m)
    mom.update(field)
    pair.update(series_at_t)
print("std", mom.std())
print("corr", pair.corr())
