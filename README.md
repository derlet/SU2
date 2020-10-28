# SU2
The underlying theory and examples of the SU(2) analysis framework may be found in the arxiv paper:

P. M. Derlet, Correlated disorder in a well relaxed model binary glass through a local SU(2) bonding topology, arXiv:2007.08878 (2020).

Version 1.0
-----------

Python implementation within OVITO

The python script SU2_v1.0.py is a python modifier for the atomistic visualization software package OVITO (developed and tested on Version 3.2.0), available at www.ovito.org.

Usage within OVITO

The python script may be manually used via the "Python script" modifier. To have it permanently available, place the file (better named as SU2.py) in the Ovito directory share/ovito/scripts/modifiers/

For further information see: http://ovito.org/docs/current/particles.modifiers.python_script.php

Commments:
1) Optimal results are obtained under full periodic boundary conditions
2) Can switch off modified Voronoi procedure via modifiedVoronoi flag
3) Can set maximum iteration number (of modified Voronoi procedure) via maxCount

Fortran90 implementation using voro++

The fortran90 program SU2_v1.0.f90 is a fortran source code that performs a similar SU(2) analysis using the open source Voronoi tessellation package voro++ available at http://math.lbl.gov/voro++

The source code may be compiled using gfortran (https://gcc.gnu.org/fortran/) via:

gfortran -Ofast -fdefault-real-8 -fdefault-integer-8 -ffree-line-length-none -o SU2.x SU2.f90

Comments:
1) The code inputs and outputs data via the *.lammps file format (see https://lammps.sandia.gov)
2) The code uses voro++ via the fortran90 command "execute_command_line", note the correct path must be given within the source code.
3) The code is intended for experienced fortran90 (or any other simple language) users, and can be easily modified for specific file formats.
4) Sample output of the code using the from-the-melt glass sample described in the archive paper is given in "from-the-melt_LJWahnstrom.out". This uses the supplied input lammps file "from-the-melt_LJWahnstrom.lammps". This test run may be executed at the command line via:

./SU2_V1.0.x < from-the-melt_LJWahnstrom.lammps > from-the-melt_LJWahnstrom.out

The outputted atomistic configuration is su2.lammps.
