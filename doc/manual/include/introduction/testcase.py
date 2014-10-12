def testcase_Dirichlet():
    # ...
    testcase = {}
    testcase['list_DirFaces'] = [[0,1,2,3]]

    phi = 2 * (1./3) * pi
    eps = 1.0e-6

    c = cos(phi)
    s = sin(phi)

    testcase['A']  = lambda x,y : [  eps * c**2 + s**2 \
                                           , ( eps - 1 ) * c * s \
                                           , ( eps - 1 ) * c * s \
                                           , eps * s**2 + c**2]

    # ...
    # rhs
    # ...
    testcase['f'] = lambda x,y : [1.]
    # ...

    return testcase
