# -*- coding: UTF-8 -*-
import numpy as np

__all__ = ['pigasus']
__author__ = 'ARA'

class pigasus(object):
    """Docstring for class pigasus."""

    #: Doc comment for class attribute pigasus.
    def __init__(self, *args, **kwargs):
        """Creates a new pigasus Object.

        Args:
           geometry (cad_geometry object):  The geometry to bu used for the
           simulation.
           testcase (dictionary):  A dictionary containing data for the PDE.

        Kwargs:
           state (bool): Current state to be in.

        Returns:
           int.  The return code::


        This class is to be used for new models.
        """
        # ...
        try:
            self.fem = kwargs['fem']
        except:
            from . import fem      as fem
            self.fem = fem.fem()
            self.fem.detail = 0
        # ...

        # ...
        try:
            self.geometry = args['geometry']
        except:
            # ...
            try:
                self.geometry = kwargs['geometry']
            except:
                self.geometry = None
        # ...

        # ...
        try:
            self.solverInfo  = kwargs['solverInfo']
        except:
            self.solverInfo = None
        # ...

        # ...
        try:
            self.testcase = kwargs['testcase']
        except:
            self.testcase = None
        # ...

        # ...
        try:
            self.name = kwargs['name']
        except:
            self.name = None
        # ...

    def reset(self):
        for F in self.fields:
            F.reset()

    def initialize(self):
        print("Not yet implemented")

    def assembly(self):
        print("Not yet implemented")

    def solve(self):
        print("Not yet implemented")

    def plot(self):
        print("Not yet implemented")

    def norm(self):
        print("Not yet implemented")

    def __del__(self):
        self.fem.__del__()

if __name__ == '__main__':
    import caid.cad_geometry  as cg
    import pylab                as pl
    import numpy                as np

    # ------------------------------
    def test1():
        from . import fem      as fem
        fe = fem.fem()
        PDE = pigasus(fem=fe)
    # ------------------------------

    # ------------------------------
    def test2():
        PDE = pigasus()
    # ------------------------------

    # ------------------------------
    def test3():
        from . import fem      as fem
        from caid.cad_geometry import line, circle, bilinear
        import caid.cad_geometry as cg

        fe = fem.fem()

        geo1 = cg.cad_geometry(geo=line())
        geo2 = cg.cad_geometry(geo=circle())
        geo3 = cg.cad_geometry(geo=bilinear())

        PDE1 = pigasus(fem=fe, geometry=geo1)
        PDE2 = pigasus(fem=fe, geometry=geo2)
        PDE3 = pigasus(fem=fe, geometry=geo3)
    # ------------------------------

    test1() ; print("test1- OK")
    test2() ; print("test2- OK")
    test3() ; print("test3- OK")
