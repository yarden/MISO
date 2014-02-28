##
## Run all internals tests
##
import os
import sys
import time

import misopy
import misopy.internal_tests as internal_tests

test_modules = ["test_scores",
                "test_lapack",
                "test_math"]
for test_mod in test_modules:
    print "Testing %s" %(test_mod)
    curr_mod = __import__(test_mod)
    curr_mod.main()
    

