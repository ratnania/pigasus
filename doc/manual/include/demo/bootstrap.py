#!/usr/bin/env python
import sys, os
import unittest

def getbuilddir():
    from distutils.util import get_platform
    s = os.path.join("build", "lib.%s-%.3s" % (get_platform(), sys.version))
    if (sys.version[:3] >= '2.6' and
        hasattr(sys, 'gettotalrefcount')):
        s += '-pydebug'
    return s

def bootstrap():
    try:
        import pigasus
    except ImportError:
        from os.path import dirname, abspath, join
        top_dir = abspath(join(dirname(__file__), '..'))
        build_dir = join(top_dir, getbuilddir())
        sys.path.insert(0, build_dir)
        import pigasus
        sys.stderr.write(
            "Loaded package 'pigasus' from build/ directory\n")
    return pigasus

bootstrap()
