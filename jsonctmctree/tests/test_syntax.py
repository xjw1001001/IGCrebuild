from __future__ import division, print_function, absolute_import

from numpy.testing import dec

try:
    eval('a @ b')
    has_matmul = True
except:
    has_matmul = False

if has_matmul:
    from jsonctmctree.tests.check_bad_syntax import *

from jsonctmctree.tests.check_good_syntax import *

dec.skipif(True, 'matmul is unavailable')
def test_matmul_stub():
    pass
