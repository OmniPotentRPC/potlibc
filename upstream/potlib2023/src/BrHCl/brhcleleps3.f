C
      SUBROUTINE PREPOT
C
C   System:           BrHCl
C   Functional form:  Extended LEPS (London-Erying-Polyani-Sato)
C
C   References for the potential parameters and the potential functional form:
C   Morse Parameters: D. J. Douglas, J. C. Polanyi, and J. J. Sloan
C                     Chem. Phys. 13, 15 (1976).
C   Functional form:  P. J. Kuntz, E. M. Nemth, J. C. Polanyi, 
C                     S. D. Rosner, and C. E. Young 
C                     J. Chem. Phys. 44, 1168 (1966)
C
C   Note: This surface differs from the original BrHCl surface in that 
C         the Sato parameters have been interchanged.  The reference for the
C         original Sato parameters is :
C                     V. K. Babamov, V. Lopez, and R. A. Marcus
C                     J. Chem. Phys. 78, 5621 (1983).
C
C   PREPOT must be called once before any calls to POT.
C   The potential parameters are included in the block data subprogram PPARMS.
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and
C        derivatives of the ground electronic state energy with respect
C        to the coordinates.
C      /PT2CM/ NSURF, NDER, NDUM(21)
C        passes the control flags where
C        NDER  = 0 => no derivatives are computed
C              = 1 => derivatives of the energy for the ground electronic
C                     state with respect to the coordinates are computed
C        NDUM  - not used
C      /PT4CM/ IPTPRT, IDUM(19)
C        IPTPRT passes the FORTRAN unit number used for potential output
C        IDUM   not used
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to
C        Re(AB), the equilibrium bond length for the AB diatomic.
C        Similarly, the terms EASYBC and EASYAC represent the energies
C        in the BC and the AC valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NSURF            0
C      NDER             1
C      IPTPRT           6
C
C*****

         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NFLAG(21)
         COMMON /PT4CM/ IPRT,IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3), 
     *                   DO4Z(3), B(3)
         PARAMETER (CKCAL = 627.5095D0)
         PARAMETER (CANGS =   0.529177106D0)
C
C   Echo the potential parameters
         WRITE (IPRT, 100) DE, RE, BETA, Z
C
100   FORMAT (/, 2X, T5, 'PREPOT has been called for BrHCl ',
     *                   'potential surface no. 3',
     *        /, 2X, T5, 'extended LEPS functional form',
     *        //, 2X, T5, 'Potential energy surface parameters:',
     *        /, 2X, T5, 'Bond', T46, 'Br-H', T58, 'H-Cl', T68, 'Cl-Br',
     *        /, 2X, T5, 'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,/)
C
      DO  10 I = 1,3
C   Convert to atomic units
             DE(I)  = DE(I)/CKCAL
             RE(I)   = RE(I)/CANGS
             BETA(I) = BETA(I)*CANGS
C   Compute useful constants
             ZPO(I)   = 1.0D0 + Z(I)
             OP3Z(I)  = 1.0D0 + 3.0D0*Z(I)
             TOP3Z(I) = 2.0D0*OP3Z(I)
             ZP3(I)   = Z(I) + 3.0D0
             TZP3(I)  = 2.0D0*ZP3(I)
             DO4Z(I)  = DE(I)/4.0D0/ZPO(I)
             B(I)     = BETA(I)*DO4Z(I)*2.0D0
10    CONTINUE
C
C    Set the values of the classical energy in the three asymptotic valleys
             EASYAB = DE(1)
             EASYBC = DE(2)
             EASYAC = DE(3)
C
      RETURN
      END
C*****
C
         SUBROUTINE POT
C
C   System:          ABC
C   Functional form: Extended LEPS (London-Erying-Polyani-Sato)
C   References:      P. J. Kuntz, E. M. Nemth, J. C. Polanyi, S. D. Rosner,
C                    and C. E. Young 
C                    J. Chem. Phys. 44, 1168 (1966)
C
C   The potential parameters must be passed through the common blocks
C   POTCM, ASYCOM, POTCCM, IOCOM, SATOCM, and LEPSCM.  
C   All information passed through the common blocks POTCM, ASYCOM, 
C   SATOCM, and LEPSCM must be in Hartree atomic units.
C
C        For the reaction: A + BC -> AB + C we write:
C                          R(1) = R(A-B)
C                          R(2) = R(B-C)
C                          R(3) = R(C-A)
C
C   NOTE: The potential energy at the reactant asympotote, that is at 
C         A infinitely far from the BC diatomic, BC diatomic at its
C         equilibrium configuration, is set equal to zero.
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3), 
     *                   DO4Z(3), B(3)
         DIMENSION X(3), COUL(3), EXCH(3)
         PARAMETER (R2 = 1.41421356D0)
C
C   Initialize the variable used in the calculation of the energy.
C
         ENERGY = 0.D0
C
C   Check the values of NSURF and NDER for validity.
         IF (NSURF .NE. 0) THEN
             WRITE (IPRT, 900) NSURF
             STOP 'POT 1'
         ENDIF
         IF (NDER .GT. 1) THEN
             WRITE (IPRT, 910) NDER
             STOP 'POT 2'
         ENDIF
C
C   Compute the energy.
C
         DO 10 I = 1,3
               X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
               COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
               EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
               ENERGY  = ENERGY + COUL(I)
10       CONTINUE
C
         RAD = SQRT((EXCH(1)-EXCH(2))**2 + (EXCH(2)-EXCH(3))**2 +
     *              (EXCH(3)-EXCH(1))**2)
C
         ENERGY = ENERGY - RAD/R2 + EASYBC
C
C   Compute the derivatives of the energy with respect to the internal
C   coordinates.
C
         IF (NDER .EQ. 1) THEN
             S = EXCH(1) + EXCH(2) + EXCH(3)
             DO 20 I = 1,3
                   DEDR(I) = B(I)*X(I)*((3.0D0*EXCH(I)-S)/R2*
     *                       (OP3Z(I)*X(I)-ZP3(I))/RAD-
     *                       ZP3(I)*X(I)+OP3Z(I))
20           CONTINUE
         ENDIF
C
900   FORMAT(/,2X,T5,'NSURF has been set equal to ',I5,
     *       /,2X,T5,'This value of NSURF is not allowed for this ',
     *               'potential, ',
     *       /,2X,T5,'only the ground electronic surface, NSURF = 0, ',
     *               'is available')
910   FORMAT(/, 2X,'POT has been called with NDER = ',I5,
     *       /, 2X, 'This value of NDER is not allowed in this ',
     *              'version of the potential.')
C
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPRT /6/
         DATA NSURF /0/
         DATA NDER /1/
         DATA NDUM /21*0/ 
C
C   Initialize the potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA DE/ 90.24D0, 106.4D0, 52.09D0/
         DATA RE/ 1.414D0, 1.275D0, 2.136D0/
         DATA BETA/ 1.851D0, 1.867D0, 1.923D0/
         DATA Z/ 0.02D0, 0.00D0, 0.02D0/
C
         END
C*****
