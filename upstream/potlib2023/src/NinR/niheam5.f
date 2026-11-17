
C   System:          NinR
C   Functional form:
C   Common name:     EAM5 (Ni surface or cluster + H)
C   Number of derivatives: 1
C   Number of bodies: variable
C   Number of electronic surfaces: 1
C   Interface: Section-2
C   Data file: niheam5.tar
C
C   References:      S.E. Wonchoba, W.-P. Hu, and D.G. Truhlar, 
C                    Physical Review B 51, 9985 (1995).
C                    This version of the potential energy surface was provided
C                    to POTLIB by Yao-Yuan Chuang (Nov. 8, 2007).
C
C   Notes:    This potential requires several routines, which are tarred 
C             together as niheam5.tar. This potential is capable of handling
C             a variable number of Ni atoms and up to one H atom.
C


EEEEEEEE         A         M       M    55555555
E               A A        MM     MM    5      
E              A   A       M M   M M    5
EEEEE         AAAAAAA      M  M M  M     555555
E            A       A     M   M   M           5
E           A         A    M       M    5      5
EEEEEEEE   A           A   M       M     555555

****
EAM5
****

This package contains files necessary to run the potential energy function
EAM5 for H interacting with Ni.


REFERENCE:  S.E. Wonchoba, W.-P. Hu, and D.G. Truhlar, Physical Review B 51,
 9985 (1995).

*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************
*******************************************************************************


This package is ready for basic use. 
The GETTING STARTED section familarizes you with the use of the package.
The DESCRIPTION OF FILES section describes in detail the files provided
  in this package.
The ADVANCED FEATURES section provides details for examining your own 
  systems in greater detail.
The ADDING Ni ATOMS section describes how to add Ni atom Cartesians to 
  the input file when their locations are different than those determined
  by the lattice constant.

GETTING STARTED:  EXAMINING THE TESTSUITE
=========================================

To examine the testsuite provided in this package, follow the procedure
outlined below.  EXAMINING THE TESTSUITE IS HIGHLY RECOMMENDED, because
it will familiarize you with the important files, which will be helpful
for running EAM5 on your own systems.  This procedure will also familiarize
you with the use of the package.

1.  Compile the source code, hni5.f, using the fortran compiler on your
    system.

2.  Leave the executable named a.out.
 
3.  Examine the file geom.index, which describes the input data files, 
    geom0 - geom7.  Decide which input geometry file you wish to run
    calculations on.

4.  Run EAM5 for H on Ni by issuing the following command at the prompt:
          driver # 
    where # is the geom# input file you wish to calculate EAM5 energetics
    and derivatives for.  For example, if you wish to calculate EAM5
    energies and derivatives for the input geometry given in the file 
    geom5, type "driver 5" at the prompt.  The EAM5 energy and derivatives
    of the energy with respect to the cartesians will be produced in 
    the output file "fort.9"

5.  Compare your results to the results given in the file test.out.  The
    results should agree.




DESCRIPTION OF THE FILES IN THIS PACKAGE
========================================

README:  This is the file that describes the EAM5 package and how to use
  it.  Evidently you've figured this one out already!

hni5.f:  This is the source code for the EAM5 potential function.  This is
  the fortran file that you need to compile.

driver:  This is the perl script that you run to execute EAM5.  The script
  takes one argument:  a number which is the number describing the geom#
  input file containg the Cartesian coordinates of the gas atoms and any 
  Ni lattice atoms not provided in potfile# (see below).

potfile1:  This is the potential data file which contains the Cartesian
  coordinates of a Ni lattice exposing the (111) crystal face.  The first 
  18 lines of this file are potential parameters that are well-documented
  in the file.  The remaining lines are the Cartesians of the Ni lattice.

potfile2:  Same as potfile1, except the lattice exposes the (100) crystal 
  face.

potfile3:  Same as potfile1, except for use with 3 displaced Ni atoms (see
  ADDING Ni ATOMS section).

geom.index:  This file contains a description of the geom# files (geom0 -
  geom7) included in this package.

geom0 - geom7:  The input files containg the Cartesians (in bohr) of the 
  gas atoms and any Ni lattice atoms not provided in potfile# (see above).

test.out:  These are the results of runs on all geom# input files given in
  this package.  It's a good idea to check your results against these
  results to insure that compilations were correct.

sphere.f:  This is a lattice creator package.  This file is described in 
  the ADVANCED FEATURES section.

SIZES:  This is an include file used with sphere.f.  It is described in the
  ADVANCED FEATURES section.



ADVANCED FEATURES:  RUNNING YOUR OWN INPUT CARTESIANS
===================================================== 

To run your own calculations on your own systems, the input data files 
(geom0 - geom7) and the potential data files (potfile1 - potfile2), and
possibly (but not always) the perl driver script will need to be modified.

1.  CREATE NI LATTICE
        FILES:  sphere.f, SIZES
        The program "sphere.f" is a menu driven code used to create a Ni
lattice, or for that matter any face-centered cubic metal.  The code is well
documented, so it should be easy to modify if necessary.  More than likely
though, only the "SIZES" file will require modification.  In "SIZES",
the user identifies a number of variables --  and these are all explained
within the text of "SIZES".  The most important variables are LATCON,
which is the lattice constant of the metal, and NLAYERS, which is a measure
of the size of the lattice.  MMOVE, MAXNM, and MAXNA are not important
for the Ni lattice.

Adapt "SIZES" file as necessary, compile sphere.f, and run it.  I recommend
the following procedure, which requires running sphere.f multiple times:
        a.  Create lattice (menu selection 1)
        b.  Decide on a point from which a radius will be extended to
              enclose atoms in a sphere which will be included in the
              calculations.  See discussion of Embedded Cluster Method
              in S.E. Wonchoba, W.-P. Hu, and D.G. Truhlar, Phys Rev B
              51, 9985 (1995).
        c.  Use these cartesians in SIZES as XCNTR1, YCNTR1, and ZCNTR1.
              You can also have a second center, whose cartesians are
              XCNTR2, YCNTR2, and ZCNTR2.  If you only have one center,
              make sure both XCNTR2 = XCNTR1, YCNTR2 = YCNTR1.

The program sphere.f has been tested and runs properly on the following
machines:
  IBM RS/6000, AIX versions 3.2.x, 4.0, 4.1, and  4.2
  Cray C90, UNICOS 8.0
It is likely that the program should work properly on most UNIX machines,
but it has been extensively used on the above machines and operating systems.

Potential data files have been provided for Ni with the (111) and (100)
crystal faces exposed, with the Ni atoms spaced with a lattice constant
of 3.5235 Ångstroms.  The sphere.f program also allows for the (110) 
face to be exposed, and for any lattice constant.


2.  CREATE OR MODIFY POTENTIAL DATA FILE (potfile1, potfile2)
        I recommend modifying the existing one.  The first 18 lines
are some potential parameters that are well-documented.
More than likely, the only values that will ever have to be changed are
in the first line (which contains the number of moving metal atoms and
the number of 2ndary zone atoms), the 17th line (the EZERO value), the
18th line (atom types of all moving atoms:  1 = Ni, 2 = H, 3 = C), and
the remaining N lines which are the cartesians of the non-moving lattice
atoms.  Delving deeper into this file will require going through the
potential routine to see what you want to change.
NOTE:  THE CARTESIAN COORDINATES ARE IN ATOMIC UNITS (bohrs).


3.  CREATE OR MODIFY INPUT DATA FILE
        Again, I recommend modifying the existing one.  The first line
of these files is the number of input atoms, and the remaining lines are
the input Cartesians of the input atoms.  The data files are labeled geom#.


4.  CONTINUE
	From here, compile the source code and follow the instructions 
given in the GETTING STARTED section for running the systems.



ADDING Ni ATOMS TO THE INPUT FILE
=================================
	If you want to include Ni atoms in your input file, which you will 
want to do if some of the Ni atoms in your system are offset slightly from 
their pre-defined locations (determined by the lattice constant and
produced from sphere.f), then you need to follow the following procedure:

1.  Decide how many Ni atoms are going to be in the geom# input file.  Call
this number "N".

2.  Edit the potential data file.  
  2a.  Change the second number in the first line of the potential file
         from 0 to "N".
  2b.  Change the third number in the first line of the potential file
         from its current value (call it "M") to ("M" - "N").
  2c.  In line 18, add "N" 1's after the initial 2.  This alerts the 
         program that there are additional Ni atoms to be read in from
         the input data file.
  2d.  Delete from lines 19 on the Cartesian coordinates of the Ni atoms
         which will be altered and moved to the input geom# file.
  2e.  Change the zero of energy in line 17 to the desired value.  Without
          changing this value, the energy will be reported with respect
          to a rigid surface -- the zero of energy will not include the
          "N" new Ni atoms.

3.  Edit the geom# file.
  3a.  Add the Cartesians of the Ni atoms to the geom# file which you 
         are editing.
  3b.  Change the first line of the input file to the total number of 
         input atoms.

Here is an example:

Suppose you wish to determine EAM5 energetics and derivatives for an H 
atom sitting atop a threefold site on Ni(111), and you have exact 
Cartesian coordinates for the three Ni atoms which form the threefold
site above which the H atom is adsorbed (I.e., these Ni atoms are
slightly displaced from their equilibrium lattice positions).

The files which you will alter are geom1 (see geom.index file) and 
potfile1 (which is for Ni(111)).

You have 3 Ni atoms which you want to add to the geom1 input file, so 
the first line of potfile1 should be changed from 
1  0  698    ! # of cov ...
to 
1  3  695    ! # of cov ...

Next, the 18th line of potfile 1 should be changed from
2
to
2  1  1  1 

Next, the 1st, 2nd, and 4th Ni atom Cartesian coordinates should be 
deleted from potfile1 (lines 19, 20, and 22), because these are the
Ni atoms that form the threefold site, and these are the atoms that 
will be placed in geom1.
 
Next, change the 17th line of potfile1 from 178.4662455 to 177.9909896 
(This value is specific for the changes being described here only!)

Finally, add the appropriate Carteisan coordinates describing Ni atoms
1, 2, and 4 (from potfile1) to the file geom1, and change the first line
in geom1 from 1 to 4.

And then compile and run as before.

To make this transition easier to understand, they have been included
in the currecnt package.  Input files geom6 and geom7 and potfile3
are for 3 displaced Ni atoms.
