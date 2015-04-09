from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose

def test_bad_syntax():
    a = np.array([[1, 1], [1, 1]])
    assert_allclose(a.dot(a), a @ a)
