!********************************************************************** 
!   System:                     FH2O
!   Functional form:            permutationally invariant polynomials-neural network
!   Common name:                FH2O
!   Number of derivatives:      0
!   Number of bodies:           4
!   Number of electronic surfaces: 1
!   Interface: Section-2
!**********************************************************************
!
!   References: J. Li, B. Jiang, and H. Guo “Permutation invariant polynomial
!   neural network approach to fitting potential energy surfaces. II. Four-atom
!   systems”, J. Chem. Phys., 139, 204103 (2013). 
!   J. Li and H. Guo “Quasi-classical trajectory study of F + H2O → HF + OH
!   reaction: influence of barrier height, reactant rotational excitation,
!   and isotopic substitution”, Chin. J. Chem. Phys., 26, 627 (2103). 
!   T. L. Nguyen, J. Li, R. Dawes, J. Stanton, and H. Guo,“Accurate
!   determination of barrier height and kinetics for the F + H2O → HF + OH
!   reaction”,  J. Phys. Chem. A, 117, 8864 (2013).
!   J. Li, B. Jiang, and H. Guo “Spin-orbit corrected full-dimensional 
!   potential energy surfaces for the two lowest-lying electronic states of
!   FH2O and dynamics for the F + H2O → HF + OH reaction”,  J. Chem. Phys.,
!   138, 074309 (2013). 
!   J. Li, R. Dawes, and H. Guo “An ab initio based full-dimensional global
!   potential energy surface for FH2O(X2A') and dynamics for the F + H2O → HF
!   + HO reaction”, J. Chem. Phys., 137, 094304 (2012).
!
!   Notes: see below
!
!   I. Potential energy calculation
!      fh2oNN(ct,vpes,vpesa,vpesb,vpesc) in fh2o-pipnn.f
!
!      H  ct(1,1) ct(2,1) ct(3,1)
!      H  ct(1,2) ct(2,2) ct(3,2)
!      F  ct(1,3) ct(2,3) ct(3,3)
!      O  ct(1,4) ct(2,4) ct(3,4)
!
!      Input: ct(3,4)                      in units of Angstrom
!      Output: vpes,vpesa,vpesb,vpesc      in units of eV 
!              vpes=(vpesa+vpesb+vpesc)/3.0d0
!
!   II. Spin-orbit corrected energy calculation
!      fh2o_soc_nn(xt,vg,socg,vgsoc) in fh2o-pipnn.f
!
!      H  xt(1,1) xt(2,1) xt(3,1)
!      H  xt(1,2) xt(2,2) xt(3,2)
!      F  xt(1,3) xt(2,3) xt(3,3)
!      O  xt(1,4) xt(2,4) xt(3,4)
!
!      Input: xt(3,4)                      in units of Angstrom
!      Output: vg,socg,vgsoc               in units of kcal/mol
!              vg = vpes * 23.0605d0 (see point I.)
!              socg is the spin-orbit correction
!              vgsoc = vg + socg
!
!**********************************************************************


