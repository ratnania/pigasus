#!/usr/bin/env python
import sys, os

original_stdout = sys.stdout

def main():
    import unittest
    from glob import glob
    testsuitedir = os.path.dirname(__file__)
    sys.path.insert(0, testsuitedir)
    pattern = 'test_*.py'
    wildcard = os.path.join(testsuitedir, pattern)
    testfiles = glob(wildcard)
    testfiles.sort()
    testsuite = unittest.TestSuite()
    testloader = unittest.TestLoader()
    for testfile in testfiles:
        sys.stdout = original_stdout
        filename = os.path.basename(testfile)
        testname = os.path.splitext(filename)[0]
        if testname in [ \
                        'test_monge_ampere_mesh_picard', \
                        'test_monge_ampere_mesh_picardtwogrids', \
                        'test_monge_ampere_picard', \
                        'test_anisotropicDiffusionMMPDE' \
                       ]:
            continue
        print("======== ", testfile, " =========")
        module = __import__(testname)
        cases = testloader.loadTestsFromModule(module)
        testsuite.addTests(cases)
        cases = []
        for attr in module.__dict__:
            if attr.startswith('test'):
                func = getattr(module, attr)
                case = unittest.FunctionTestCase(func)
                cases.append(case)
        testsuite.addTests(cases)
    runner = unittest.TextTestRunner()
    result = runner.run(testsuite)
    return result.wasSuccessful()

if __name__ == '__main__':
    sys.dont_write_bytecode = True
    import bootstrap
    main()



##!/usr/bin/env python
#import sys, os
#from time import time
#
#def main():
##    import unittest
#    from glob import glob
#    testsuitedir = os.path.dirname(__file__)
#    sys.path.insert(0, testsuitedir)
#    pattern = 'test_*.py'
#    wildcard = os.path.join(testsuitedir, pattern)
#    testfiles = glob(wildcard)
#    testfiles.sort()
#    tb = time()
#    for testfile in testfiles:
#        filename = os.path.basename(testfile)
#        testname = os.path.splitext(filename)[0]
#        if testname in [ 'test_monge_ampere_mesh_MG' \
#                       , 'test_anisotropicDiffusionMMPDE']:
#            continue
#        print "======== running ", testname, " ========"
#        os.system('python ' + testfile)
#    te = time()
#    print "-------------------------------------------------------------"
#    print "Ran ", len(testfiles) ," tests in ", te-tb, "s"
#    print " "
#    print "OK"
#
#
#if __name__ == '__main__':
#    sys.dont_write_bytecode = True
#    import bootstrap
#    main()
