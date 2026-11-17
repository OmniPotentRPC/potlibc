****************************************************************************
* PIP-NN PES for the ground state of HO3 system at MRCI-F12/VDZ-F12 level. *
* Reference: Junxiang Zuo et al., J. Phys. Chem. A, 2020                   *
****************************************************************************
                            
These files are for the PES of H+O3/HO2+O reactions, including
  1. Parameter files: biases*.txt, weights*.txt
  2. Main fortran file containing subroutines: ho3_pes.f

To use this PES, as shown in the test.f, please be sure to do:
  1. Initialization:
       call pes_init 
     Note: call once and only once in the main code.
  2. Calculate potential energy:
       call ho3_pes_interface(q,vpes)
     where,
       q:    XYZ cartesian coordinate in HOOO order and angstrom unit.
       vpes: the potential in eV.
