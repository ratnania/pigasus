#!/usr/bin/env python
import sys, os
from time import time

def main():
#    import unittest
    from glob import glob
    testsuitedir = os.path.dirname(__file__)
    sys.path.insert(0, testsuitedir)
    pattern = 'test_*.py'
    wildcard = os.path.join(testsuitedir, pattern)
    testfiles = glob(wildcard)
    testfiles.sort()
    tb = time()
    for testfile in testfiles:
        filename = os.path.basename(testfile)
        testname = os.path.splitext(filename)[0]
        print "======== running ", testname, " ========"
        os.system('python ' + testfile)
    te = time()
    print "-------------------------------------------------------------"
    print "Ran ", len(testfiles) ," tests in ", te-tb, "s"
    print " "
    print "OK"


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    import bootstrap
    main()
