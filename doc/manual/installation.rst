.. role:: envvar(literal)
.. role:: command(literal)
.. role:: file(literal)
.. _installation:


======================
Setup and Installation
======================

Pigasus comes with the following folders:

* **caid** : This is the CAD modeler. More details can be found `here <http://www.ratnani.org/caid>`_.

* **igakit** : contains basic structures and operations for CAD design. The current version was modified to take into account the multi-patch description and its connectivity. More details can be found `here <http://bitbucket.org/dalcinl/igakit>`_.

* **pigasus** : This is the Generic FEM solver. You will find in this website all you need to understand the Pigasus philosophy and write your own scripts.

.. * **spm** : This module is used to encapsulates the treatment of Sparse Matrices. More details can be found `here <http://www.ratnani.org/spm>`_.


Dependecies
***********

You can either install `canopy <https://www.enthought.com/downloads/>`_ or install the following packages **numpy**, **scipy** and **matplotlib**. Details on the installation can be found `here <http://www.scipy.org/install.html>`_.

You will also need a Fortran compiler [#f1]_ and the fortran wrapper *f2py*.

The following packages are optional.

* **petsc4py**: Please download and install `PetSc <http://www.mcs.anl.gov/petsc/documentation/installation.html>`_ then just type::

      sudo pip install petsc4Py

make sure that you defined the **PETSC_DIR** and **PETSC_ARCH** in your *.bashrc* or *.bash_profile*

* **pyAMG**: is an Algebraic MultiGrid solver for Python. You can install it by running::

      sudo pip install pyamg


Installation
************

IGAKIT installation
^^^^^^^^^^^^^^^^^^^

This can be done by the classical setup::

      cd igakit; sudo python setup.py install
      cd ..

PIGASUS installation
^^^^^^^^^^^^^^^^^^^^

The installation has mainly two steps. The first one is to use cmake to compile all shared libraries. The second one, is to use distutils to install the package. Therefor, we need to specify the compiler and its constructor, which is needed for f2py. You must use the same compiler for both cmake and distutils.

The installation needs the following steps::

     cd pigasus
     mkdir build 
     cd build 
     # export the PIGASUS_BUILD_DIR which will be used latter, by the python setup
     export PIGASUS_BUILD_DIR=$PWD
     # choose and export your compiler
     #export FC=ifort ; export FC_CONSTRUCTOR=intelem
     export FC=gfortran  ; export FC_CONSTRUCTOR=gnu95
     # configure with prefix install dir
     cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/installation/pigasus
     make install

.. notes::
   
   This is a local installation. You can also use sudo, but in this case make sure to propagate the environement variables using sudo -E python setup.py install for example. If you are using ifort, please make sure that you have also sourced the compilervars.


.. notes::

   Remember that for f2py you need to give the constructor. For gfortran it's gnu95, for ifort it's intelem

.. notes::

   On ccamu, I had to set up the LD_LIBRARY_PATH, PIGASUS_BUILD_DIR variables before running the python setup installation for pigasus. Otherwise, it can not find and link the libraries compiled by cmake.

Update your *.bashrc* or *.bash_profile*::

     # add PIGASUS variables into .bashrc
     export PIGASUS_BUILD_DIR=path/to/pigasus/build
     export LD_LIBRARY_PATH=$PIGASUS_BUILD_DIR/lib:$LD_LIBRARY_PATH
     # updating the PYTHONPATH environement variable
     # add this line into your .bashrc or .bash_profile
     export PYTHONPATH=/path/to/pigasus/install/directory:$PYTHONPATH
     # define your plugin home
     export PIGASUS_PLUGIN_DIR=path/to/pigasus_plugin

and run on your terminal::

      source $HOME/.bashrc      

In order to check that everything is OK::

     # check that pigasus is installed
     python -c "import pigasus"
     # check that fortran-core is installed
     python -c "import pigasus.fem.core as co; fem = co.pyfem"
     # first test to check if everything is ok
     cd gallery/demo
     python test_poisson_1d.py
     # run all tests (some tests may not run, if you don't specify the right plugin_dir)
     python runtests.py   
         

.. rubric:: Footnotes

.. [#f1] Only  *gfortran* and *ifort* have been tested for the moment. Please inform me if you test any other Compiler.

.. Local Variables:
.. mode: rst
.. End:
