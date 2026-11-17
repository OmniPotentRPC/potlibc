C
      SUBROUTINE PREPOT
C
C
C   System:           BrH2
C   Functional form:  LEPS (London-Erying-Polyani-Sato)
C   Common name:
C   Interface:        3-XS
C   References:        R. Weston, J. Chem. Phys., Vol. 31, 892 (1959).
C
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Br
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C      Set in the BLOCK DATA subprogram PTPARM
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Br-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Br-H)
C
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
C        In this potential the AB valley represents H infinitely far from
C        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the terms
C        EASYBC and EASYAC represent the energies in the H2 and the BrH 
C        valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      NSURF            0
C      NDER             1 
C      IPTPRT           6
C
C*****
C
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3), 
     *                   DO4Z(3), B(3)
         PARAMETER (CKCAL = 627.5095D0)
         PARAMETER (CANGS =   0.529177106D0)
         save
C
C   Echo the potential parameters
         WRITE (IPTPRT, 100) DE, RE, BETA, Z
C
100   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'BrH2 LEPS potential energy surface',
     *      //,1X,T5,'Parameters:',
     *        /, 2X, T5, 'Bond', T46, 'Br-H', T58, 'H-H', T69, 'H-Br',
     *        /, 2X, T5, 'Dissociation energies (kcal/mol):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Equilibrium bond lengths (Angstroms):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Morse beta parameters (Angstroms**-1):', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,
     *        /, 2X, T5, 'Sato parameters:', 
     *        T44, F10.5, T55, F10.5, T66, F10.5,//,1X,'*****')
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
C   Initialize the asymptotic energy values
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
C   Functional form: LEPS (London-Erying-Polyani-Sato)
C   References:      S. Sato
C                    J. Chem. Phys. 592, 2465 (1955)
C                    
C
C   The potential parameters must be passed through the common blocks
C   PT1CM, PT5CM, PT2CM, PT4CM, SATOCM, and LEPSCM.  
C   All information passed through the common blocks PT1CM, PT5CM, 
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
         COMMON /PT4CM/ IPTPRT
         COMMON /PT1CM/R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /PT5CM/ EASYAB, EASYBC, EASYAC
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
         COMMON /LEPSCM/ ZPO(3), OP3Z(3), ZP3(3), TZP3(3), TOP3Z(3), 
     *                   DO4Z(3), B(3)
         DIMENSION X(3), COUL(3), EXCH(3)
         PARAMETER (R2 = 1.41421356D0)
         save
C
C   Initialize the variable used in the calculation of the energy.
C
         ENERGY = 0.D0
C
C   Check the value of NDER 
         IF (NDER .GT. 1) THEN
             WRITE (IPTPRT, 900) NDER
             STOP 'POT 1'
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
900   FORMAT(/,1X,T5,'Error: POT has been called with NDER = ', I5,
     *       /,1X,T12,'only the first derivatives, NDER = 1, are ',
     *                'coded in this potential')
C
      RETURN
      END
C*****
C
         BLOCK DATA PTPARM
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         COMMON /PT4CM/ IPTPRT, IDUM(19)
         COMMON /PT2CM/ NSURF, NDER, NDUM(21)
         COMMON /SATOCM/ DE(3), RE(3), BETA(3), Z(3) 
C
C   Initialize the flags and the I/O unit numbers for the potential
         DATA IPTPRT, NDER /6, 1/
         DATA NDUM /21*0/
C
C   Initialize the potential parameters; the energy parameters are in
C   kcal/mol, and the lengths are in Angstroms.
         DATA DE/ 90.32790D0, 109.396D0, 90.32790D0/
         DATA RE/ 1.41381D0, 0.74167D0, 1.41381D0/
         DATA BETA/1.81088D0, 1.94198D0, 1.81088D0/
         DATA Z/0.2D0, 0.2D0, 0.2D0/
C
         END
C*****
