!********************************************************************** 
!   System:                     O2H3
!   Functional form:            permutationally invariant polynomials-neural network
!   Common name:                O2H3
!   Number of derivatives:      0
!   Number of bodies:           5
!   Number of electronic surfaces: 1
!   Interface: Section-2
!**********************************************************************
!
!   References: 
!   “Quasi-classical trajectory studies on the full-dimensional accurate 
!   potential energy surface for the OH + H2O = H2O + OH reaction”
!   Mengna Bai, Dandan Lu, Jun Li*,  
!   Physical Chemistry Chemical Physics, 2017, 19, 17718-17725
!
!   Notes: see below
!
!**********************************************************************

****************************************************************************
* PIP-NN PES for OH + H2O system at CCSD(T)-F12a/AVTZ level.
* Reference: Physical Chemistry Chemical Physics, 2017, 19, 17718-17725.
****************************************************************************
                            
These files are for the PES of OH + H2O reaction, including
  1. Parameter files: biases*.txt, weights*.txt
  2. Main fortran file containing subroutines: h3o2-pipnn.f

To use this PES, as shown in the main.f, please be sure to do:
  1. Initialization:
       call pes_init() 
     Note: call once and only once in the main code.
  2. Calculate potential energy:
       call h3o2pipNN(c,vpes)
     where,
       c: XYZ cartesian coordinate in HHHOO order and angstrom unit, shown in test-h3o2.xyz
       vpes:  the potential energy(in eV) relative to the reactant asymptote.
