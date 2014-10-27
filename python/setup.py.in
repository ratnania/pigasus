# -*- coding: UTF-8 -*-
#! /usr/bin/python

import os
from distutils.command.install import install as DistutilsInstall
from distutils.command.build import build as DistutilsBuild
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

NAME    = 'pigasus'
VERSION = '0.9.2'
AUTHOR  = 'Ahmed Ratnani'
EMAIL   = 'ratnaniahmed@gmail.com'
URL     = 'http://www.ratnani.org/pigasus/'
DESCR   = 'Python For IsoGeometric AnalySis and Unified Simulations.'

class MyInstall(DistutilsInstall):
    def run(self):
        DistutilsInstall.run(self)

BUILD_DIR = str(os.environ['PIGASUS_BUILD_DIR'])+"/"
library_dirs = [BUILD_DIR+"lib"]
libraries = [  "tools" \
             , "tracelog" \
             , "geometries" \
             , "blackboxes" \
             , "connectivities" \
             , "grids" \
             , "spm" \
             , "fem" \
             , "assembly"]
include_dirs = []
for lib in libraries:
    include_dirs.append(BUILD_DIR+"fortran/"+lib)

fem_core = Extension('pigasus.fem.core',
                sources = [BUILD_DIR+'python/core/pyfem.f90'],
                f2py_options = ['--quiet'],
                define_macros = [#('F2PY_REPORT_ATEXIT', 0),
#                                 ('F2PY_REPORT_ON_ARRAY_COPY', 0)
                                ],
                include_dirs=include_dirs,
                library_dirs=library_dirs,
                libraries=libraries)

setup(name=NAME,
          version=VERSION,
          author=AUTHOR,
          author_email=EMAIL,
          url=URL,
          description=DESCR,
          license='LICENSE.txt',
          long_description=open('README.md').read(),
    packages=[  'pigasus' \
              , 'pigasus.fem' \
              , 'pigasus.gallery' \
              , 'pigasus.multigrid' \
              , 'pigasus.fit' \
              , 'pigasus.interpolate' \
              , 'pigasus.plugin' \
              , 'pigasus.utils' \
             ],
    package_dir={  'pigasus': 'pigasus'\
                  ,'pigasus.fem': 'fem' \
                  ,'pigasus.gallery': 'gallery' \
                  ,'pigasus.multigrid': 'multigrid'\
                  ,'pigasus.fit': 'fit' \
                  ,'pigasus.interpolate': 'interpolate' \
                  ,'pigasus.plugin': 'plugin' \
                  ,'pigasus.utils':  'utils' \
                  ,},
    ext_modules  = [fem_core],
    cmdclass={'install': MyInstall},
)
