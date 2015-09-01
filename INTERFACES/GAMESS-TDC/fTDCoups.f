      PROGRAM fTDCoups
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Fortran part of the TDCoups script for time-derivative couplings for C
C FOMO calculation in GAMESS.                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                17.05.2012                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         DECLARATION
C
      IMPLICIT NONE
C
C         DECLARATION - PARAMETERS
C
      INTEGER MaxOrb
        PARAMETER (MaxOrb = 650)
C
      INTEGER MaxCI
        PARAMETER (MaxCI = 5500)
C
      INTEGER MaxSt
        PARAMETER (MaxSt = 30)
C
C         DECLARATION - FUNCTION
C
      DOUBLE PRECISION FindDet
      EXTERNAL FindDet
C
C         DECLARATION - VARIABLES
C
      INTEGER Indx, VnIndx, Indx2, Indx3, IndIndx, Indx4, Indx5, Indx6, 
     &        Indx7, Indx8, Indx9
C
      INTEGER NOrb, NFunc, NStates, MaxState1, MaxState2
C
      INTEGER NOrbCIMin, NOrbCIMax
C
      INTEGER MODNOrb, OrbNum(MaxOrb), AtmNum(MaxOrb)
C
      INTEGER MODNFunc, PomIndx, FuncNum(5)
C
      INTEGER StAlpha1(MaxSt, MaxCI, MaxOrb), 
     &        StBeta1(MaxSt, MaxCI, MaxOrb),
     &        NExpVect1(MaxSt)
C
      INTEGER StAlpha2(MaxSt, MaxCI, MaxOrb), 
     &        StBeta2(MaxSt, MaxCI, MaxOrb),
     &        NExpVect2(MaxSt)
C
      INTEGER NEleAlpha, NEleBeta
      DATA NEleAlpha, NEleBeta /0, 0/
C
      INTEGER FlipOrbSign
C
      INTEGER NBetter, NSmaller
C
      INTEGER MaxIndx(MaxSt), NAktIndx, SumNOT
C
      INTEGER NStCSF, CSFIde1(MaxSt, MaxCI), CSFIde2(MaxSt, MaxCI), 
     &        NCSF1(MaxSt), NCSF2(MaxSt)
C
      DOUBLE PRECISION OrbCoeff1(MaxOrb, MaxOrb), 
     &                 OrbCoeff2(MaxOrb, MaxOrb)
C
      DOUBLE PRECISION StateEn1(MaxSt), StExpCoeff1(MaxSt, MaxCI), 
     &                 StateSpin1(MaxSt)
C
      DOUBLE PRECISION StateEn2(MaxSt), StExpCoeff2(MaxSt, MaxCI),
     &                 StateSpin2(MaxSt)
C
      DOUBLE PRECISION StCBuff
C
      DOUBLE PRECISION Overlap(2*MaxOrb, 2*MaxOrb)
C
      DOUBLE PRECISION Spq(MaxOrb, MaxOrb)
C
      DOUBLE PRECISION Coupling(MaxSt, MaxSt), NewCoup(MaxSt, MaxSt),
     &                 BestCoup(MaxSt, MaxSt)
C
      DOUBLE PRECISION DesireState
C
      DOUBLE PRECISION TestSignOrb1, TestSignOrb2
C
      DOUBLE PRECISION OrbSum(MaxOrb)
C
      DOUBLE PRECISION OldDiag(MaxSt)
C
      DOUBLE PRECISION AktMax
C
      DOUBLE PRECISION PatVal1(MaxCI, MaxSt), PatVal2(MaxCI, MaxSt)
C
      DOUBLE PRECISION StExpCoeff2B(MaxSt, MaxCI), 
     &                 StExpCoeff2C(MaxSt, MaxCI)
C
      DOUBLE PRECISION SumUp
C
      DOUBLE PRECISION CSFInt1(MaxSt, MaxCI), CSFInt2(MaxSt, MaxCI), 
     &                 SumInt1(MaxSt), SumInt2(MaxSt), CSFDiff(MaxSt)
C
      CHARACTER junk
C
      CHARACTER*64 StateFormat
C
      CHARACTER FRMTPrint*20 
C
      CHARACTER Patern1(MaxSt, MaxCI)*50, Patern2(MaxSt, MaxCI)*50
      CHARACTER AllPatern(MaxCI)*50
C
      LOGICAL IniSigns(MaxSt), IniThisState
C
      LOGICAL ISTHERE, IniBelow, IniNotThere
C
C-----------------------------------------------------------------------
C
C         MAIN 
C     
C     *****  READING INPUT DATA  *****
C       Opening file: orb_tdc1.temp and extracting information from it
      OPEN(11, file='orb_tdc1.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C        Number of orbitals, number of AO functions
      READ (11, *) NOrb, NFunc
C        Preparing to read orbital coefficients
      MODNOrb = MOD(NOrb, 5)
C        Reading orbital coefficients
      DO Indx = 1, NOrb/5
        READ (11, *) 
        READ (11, *) (OrbNum(VnIndx), VnIndx = Indx*5-4, Indx*5)
        READ (11, *) 
        READ (11, *) 
        DO Indx2 = 1, NFunc
          READ (11, *) junk, junk, AtmNum(Indx2), junk, 
     &        (OrbCoeff1(VnIndx, Indx2), VnIndx = Indx*5-4, Indx*5)
        ENDDO
      ENDDO
C        Reading the last non-complete block of orbital coeffcient data
      IF (MODNOrb .NE. 0) THEN
        READ (11, *) 
        READ (11, *) (OrbNum(VnIndx), VnIndx = NOrb-MODNOrb+1, NOrb)
        READ (11, *) 
        READ (11, *) 
        DO Indx2 = 1, NFunc
          READ (11, *) junk, junk, AtmNum(Indx2), junk, 
     &        (OrbCoeff1(VnIndx, Indx2), VnIndx = NOrb-MODNOrb+1,NOrb)
        ENDDO
      ENDIF
C     CCCCCCCCCCCCCCCCCC
CCC   OPEN(21, file='d_orb1.dbg')
CCC   WRITE (21, 9901) (OrbNum(VnIndx), VnIndx = 1, NOrb)
CCC   DO Indx = 1, NFunc
CCC     WRITE (21, 9900) AtmNum(Indx), 
CCC  &               (OrbCoeff1(VnIndx, Indx), VnIndx = 1, NOrb)
CCC   ENDDO
CCC   CLOSE(21)
C        Closing of file orb_tdc1.temp
      CLOSE(11)
C        
C       Opening file: ci_tdc1.temp and extracting information from it
      OPEN(11, file='ci_tdc1.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C        Reading number of states used in FOMO
      READ (11, *) NStates, NOrbCIMin, NOrbCIMax, DesireState
      WRITE (FRMTPrint, *) "(", 2*(NOrbCIMax-NOrbCIMin+1), "I1)"
C        Reading the states information: coeffs, energies and orbital expansion
      DO Indx = 1, NStates
        READ (11, *) junk, StateEn1(Indx), NExpVect1(Indx), 
     &                     StateSpin1(Indx)
        DO Indx2 = 1, NExpVect1(Indx)
          READ (11, *) (StAlpha1(Indx, Indx2, VnIndx), 
     &                           VnIndx = NOrbCIMin, NOrbCIMax), junk,
     &                 (StBeta1(Indx, Indx2, VnIndx), 
     &                           VnIndx = NOrbCIMin, NOrbCIMax), junk,
     &                 StExpCoeff1(Indx, Indx2) 
          DO Indx3 = 1, NOrbCIMin-1
            StAlpha1(Indx, Indx2, Indx3) = 1
            StBeta1(Indx, Indx2, Indx3) = 1
          ENDDO
          WRITE (Patern1(Indx, Indx2), FRMTPrint) 
     &                         (StAlpha1(Indx, Indx2, VnIndx), 
     &                                   VnIndx = NOrbCIMin, NOrbCIMax),
     &                         (StBeta1(Indx, Indx2, VnIndx),
     &                                   VnIndx = NOrbCIMin, NOrbCIMax)
        ENDDO
      ENDDO
      MaxState1 = NStates
C     CCCCCCCCCCCCCCCCCCCCCCC
CCC   OPEN(21, file='d_ci1.dbg')
CCC   DO Indx = 1, MaxState1 
CCC     WRITE (21, 9900) Indx, StateEn1(Indx)
CCC     DO Indx2 = 1, NExpVect1(Indx)
CCC       WRITE (21, '(7I2, "|", 7I2, F14.8)') 
CCC  &                (StAlpha1(Indx, Indx2, VnIndx), 
CCC  &                            VnIndx = NOrbCIMin, NOrbCIMax), 
CCC  &                (StBeta1(Indx, Indx2, VnIndx), 
CCC  &                         VnIndx = NOrbCIMin, NOrbCIMax), 
CCC  &                 StExpCoeff1(Indx, Indx2) 
CCC     ENDDO
CCC   ENDDO
CCC   CLOSE(21)
C        Closing of file ci_tdc1.temp
      CLOSE(11)
C
C       Opening file: csf_tdc1.temp
      OPEN(11, file='csf_tdc1.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C       Number of actual states 
      READ (11, *) NStCSF
C       CSFs intensities and their identification
      DO Indx = 1, NStCSF
        READ (11, *) junk, NCSF1(Indx)
        DO Indx2 = 1, NCSF1(Indx)
          READ (11, *) junk, CSFInt1(Indx, Indx2), CSFIde1(Indx, Indx2)
        ENDDO
      ENDDO
C       Closing of file csf_tdc1.temp
      CLOSE(11)
C
C       Opening file: orb_tdc2.temp and extracting information from it
      OPEN(11, file='orb_tdc2.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C        Number of orbitals, number of AO functions
      READ (11, *) 
C        Reading orbital coefficients
      DO Indx = 1, NOrb/5
        READ (11, *) 
        READ (11, *) 
        READ (11, *) 
        READ (11, *) 
        DO Indx2 = 1, NFunc
          READ (11, *) junk, junk, AtmNum(Indx2), junk, 
     &        (OrbCoeff2(VnIndx, Indx2), VnIndx = Indx*5-4, Indx*5)
        ENDDO
      ENDDO
C        Reading the last non-complete block of orbital coeffcient data
      IF (MODNOrb .NE. 0) THEN
        READ (11, *) 
        READ (11, *) 
        READ (11, *) 
        READ (11, *) 
        DO Indx2 = 1, NFunc
          READ (11, *) junk, junk, AtmNum(Indx2), junk, 
     &        (OrbCoeff2(VnIndx, Indx2), VnIndx = NOrb-MODNOrb+1,NOrb)
        ENDDO
      ENDIF
C     CCCCCCCCCCCCCCCCCCCCCCCCC
CCC   OPEN(21, file='d_orb2.dbg')
CCC   WRITE (21, 9901) (OrbNum(VnIndx), VnIndx = 1, NOrb)
CCC   DO Indx = 1, NFunc
CCC     WRITE (21, 9900) AtmNum(Indx), 
CCC  &               (OrbCoeff2(VnIndx, Indx), VnIndx = 1, NOrb)
CCC   ENDDO
CCC   CLOSE(21)
C        Closing of file orb_tdc2.temp
      CLOSE(11)
C        
C       Opening file: ci_tdc2.temp and extracting information from it
      OPEN(11, file='ci_tdc2.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C        Reading number of states used in FOMO
      READ (11, *) 
C        Reading the states information: coeffs, energies and orbital expansion
      DO Indx = 1, NStates
        READ (11, *) junk, StateEn2(Indx), NExpVect2(Indx), 
     &                     StateSpin2(Indx)
        DO Indx2 = 1, NExpVect2(Indx)
          READ (11, *) (StAlpha2(Indx, Indx2, VnIndx), 
     &                            VnIndx = NOrbCIMin, NOrbCIMax), junk,
     &                 (StBeta2(Indx, Indx2, VnIndx),
     &                           VnIndx = NOrbCIMin, NOrbCIMax), junk,
     &                 StExpCoeff2(Indx, Indx2) 
          DO Indx3 = 1, NOrbCIMin-1
            StAlpha2(Indx, Indx2, Indx3) = 1
            StBeta2(Indx, Indx2, Indx3) = 1
          ENDDO
          WRITE (Patern2(Indx, Indx2), FRMTPrint) 
     &                         (StAlpha2(Indx, Indx2, VnIndx), 
     &                                   VnIndx = NOrbCIMin, NOrbCIMax),
     &                         (StBeta2(Indx, Indx2, VnIndx),
     &                                   VnIndx = NOrbCIMin, NOrbCIMax)
        ENDDO
      ENDDO
      MaxState2 = NStates
C     CCCCCCCCCCCCCCCCCC
CCC   OPEN(21, file='d_ci2.dbg')
CCC   DO Indx = 1, MaxState2
CCC     WRITE (21, 9900) Indx, StateEn2(Indx)
CCC     DO Indx2 = 1, NExpVect2(Indx)
CCC       WRITE (21, '(7I2, "|", 7I2, F14.7)') 
CCC  c                (StAlpha2(Indx, Indx2, VnIndx), 
CCC  &                            VnIndx = NOrbCIMin, NOrbCIMax), 
CCC  &                (StBeta2(Indx, Indx2, VnIndx), 
CCC  &                         VnIndx = NOrbCIMin, NOrbCIMax), 
CCC  &                 StExpCoeff2(Indx, Indx2) 
CCC     ENDDO
CCC   ENDDO
CCC   CLOSE(21)
C        Closing of file ci_tdc2.temp
      CLOSE(11)
C
C       Opening file: csf_tdc2.temp
      OPEN(11, file='csf_tdc2.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C       Number of actual states 
      READ (11, *) NStCSF
C       CSFs intensities and their identification
      DO Indx = 1, NStCSF
        READ (11, *) junk, NCSF2(Indx)
        DO Indx2 = 1, NCSF2(Indx)
          READ (11, *) junk, CSFInt2(Indx, Indx2), CSFIde2(Indx, Indx2)
        ENDDO
      ENDDO
C       Closing of file csf_tdc2.temp
      CLOSE(11)
C
C       Opening file: over_tdc.temp and extracting information from it
      OPEN(11, file='over_tdc.temp', STATUS='OLD', ACCESS='SEQUENTIAL')
C        Preparing to read overlap matrix
      MODNFunc = MOD(NFunc*2, 5)
C        Reading Overlap matrix
      READ (11, *)
      DO Indx = 1, (NFunc*2)/5
        PomIndx = 1
        READ (11, *)
        READ (11, *) (FuncNum(VnIndx), VnIndx = 1, 5)
        READ (11, *) 
        DO Indx2 = Indx*5-4, NFunc*2
          READ (11, *) junk, junk, junk, junk,
     &                 (Overlap(Indx2, (FuncNum(VnIndx))), 
     &                                         VnIndx = 1, PomIndx)
          PomIndx = PomIndx + 1
          IF (PomIndx .GT. 5) THEN
            PomIndx = 5
          ENDIF
        ENDDO
      ENDDO
C        Reading last block of Overlap matrix data
      IF (MODNFunc .NE. 0) THEN
        PomIndx = 1
        READ (11, *)
        READ (11, *) (FuncNum(VnIndx), VnIndx = 1, MODNFunc)
        READ (11, *) 
        DO Indx2 = Indx*5-4, NFunc*2
          READ (11, *) junk, junk, junk, junk,
     &                 (Overlap(Indx2, (FuncNum(VnIndx))), 
     &                                         VnIndx = 1, PomIndx)
          PomIndx = PomIndx + 1
        ENDDO
      ENDIF
C        Filling up the bottom left part of matrix with data
      DO Indx = 1, NFunc*2
        DO Indx2 = Indx, NFunc*2
          Overlap(Indx, Indx2) = Overlap(Indx2, Indx)
        ENDDO
      ENDDO
C
CCC   OPEN(21, file='d_over2.dbg')
CCC   WRITE (21, 9901) (VnIndx, VnIndx = 1, NFunc)
CCC   DO Indx = NFunc + 1, NFunc*2
CCC     WRITE (21, 9900) Indx, (Overlap(Indx, VnIndx), VnIndx=1,NFunc)
CCC   ENDDO
CCC   CLOSE(21)
C
      CLOSE(11)
C
C
      ISTHERE = .FALSE.
      INQUIRE(file='str.state.tmp',EXIST=ISTHERE)
      IF (ISTHERE) THEN
        OPEN(11, file='str.state.tmp')
        READ (11, *) MaxState1
        CLOSE(11)
        MaxState2 = MaxState1
      ENDIF
C
C     *****  END OF READING INPUT DATA  *****
C
C     *****  DATA PROCESSING  *****
C
C       Sort CI vector for time 1 to time 2 CI vector
CCC   DO Indx = 1, MIN(MaxState1, MaxState2)
CCC     DO Indx2 = 1, NExpVect2(Indx)
CCC       IF (Patern2(Indx, Indx2) /= Patern1(Indx, Indx2)) THEN
CCC         DO Indx3 = Indx2, NExpVect1(Indx)
CCC           IF (Patern2(Indx, Indx2) == Patern1(Indx, Indx3)) THEN
CCC             Patern1(Indx, Indx3) = Patern1(Indx, Indx2)
CCC             Patern1(Indx, Indx2) = Patern2(Indx, Indx2)
CCC             StCBuff = StExpCoeff1(Indx, Indx3)
CCC             StExpCoeff1(Indx, Indx3) = StExpCoeff1(Indx, Indx2)
CCC             StExpCoeff1(Indx, Indx2) = StCBuff
CCC             DO Indx4 = NOrbCIMin, NOrbCIMax
CCC               IF (StAlpha1(Indx, Indx2, Indx4) .NE. 
CCC  &                StAlpha1(Indx, Indx3, Indx4)) THEN
CCC                 IF (StAlpha1(Indx, Indx2, Indx4) .EQ. 1) THEN
CCC                   StAlpha1(Indx, Indx2, Indx4) = 0
CCC                   StAlpha1(Indx, Indx3, Indx4) = 1
CCC                 ELSE
CCC                   StAlpha1(Indx, Indx2, Indx4) = 1
CCC                   StAlpha1(Indx, Indx3, Indx4) = 0
CCC                 ENDIF
CCC               ENDIF
CCC               IF (StBeta1(Indx, Indx2, Indx4) .NE. 
CCC  &                StBeta1(Indx, Indx3, Indx4)) THEN
CCC                 IF (StBeta1(Indx, Indx2, Indx4) .EQ. 1) THEN
CCC                   StBeta1(Indx, Indx2, Indx4) = 0
CCC                   StBeta1(Indx, Indx3, Indx4) = 1
CCC                 ELSE
CCC                   StBeta1(Indx, Indx2, Indx4) = 1
CCC                   StBeta1(Indx, Indx3, Indx4) = 0
CCC                 ENDIF
CCC               ENDIF
CCC             ENDDO
CCC           ENDIF
CCC         ENDDO
CCC       ENDIF
CCC     ENDDO
CCC   ENDDO
C
CC    OPEN(21, file='prehoz1.dat')
CC    OPEN(22, file='prehoz2.dat')
CC    DO Indx = 1, MaxState1
CC      WRITE (21, '(I3)') Indx
CC      DO Indx2 = 1, NExpVect1(Indx)
CC        WRITE (21, '(21I1, " | ", 21I1, " | ", F14.7)') 
CC   &                         (StAlpha1(Indx, Indx2, VnIndx), 
CC   &                                   VnIndx = 1, NOrbCIMax),
CC   &                         (StBeta1(Indx, Indx2, VnIndx),
CC   &                                   VnIndx = 1, NOrbCIMax),
CC   &                         StExpCoeff1(Indx, Indx2)
CC      ENDDO
CC    ENDDO
CC    DO Indx = 1, MaxState2
CC      WRITE (22, '(I3)') Indx
CC      DO Indx2 = 1, NExpVect2(Indx)
CC        WRITE (22, '(21I1, " | ", 21I1, " | ", F14.7)') 
CC   &                         (StAlpha2(Indx, Indx2, VnIndx), 
CC   &                                   VnIndx = 1, NOrbCIMax),
CC   &                         (StBeta2(Indx, Indx2, VnIndx),
CC   &                                   VnIndx = 1, NOrbCIMax),
CC   &                         StExpCoeff2(Indx, Indx2)
CC      ENDDO
CC    ENDDO
CC    CLOSE(21)
CC    CLOSE(22)
CCCC  NAktIndx = 1
CCCC  AllPatern(NAktIndx) = Patern2(1, 1)
CCCC  DO Indx = 1, MIN(MaxState1, MaxState2)
CCCC    DO Indx2 = 1, NExpVect2(Indx)
CCCC      SumNOT = 0
CCCC      DO Indx3 = 1, NAktIndx
CCCC        IF (Patern2(Indx, Indx2) .NE. AllPatern(Indx3)) THEN
CCCC          SumNOT = SumNOT+1
CCCC        ENDIF
CCCC        IF (SumNOT .EQ. NAktIndx) THEN
CCCC          NAktIndx = NAktIndx+1
CCCC          AllPatern(NAktIndx) = Patern2(Indx, Indx2)
CCCC        ENDIF
CCCC      ENDDO
CCCC    ENDDO
CCCC    DO Indx2 = 1, NExpVect1(Indx)
CCCC      SumNOT = 0
CCCC      DO Indx3 = 1, NAktIndx
CCCC        IF (Patern1(Indx, Indx2) .NE. AllPatern(Indx3)) THEN
CCCC          SumNOT = SumNOT+1
CCCC        ENDIF
CCCC        IF (SumNOT .EQ. NAktIndx) THEN
CCCC          NAktIndx = NAktIndx+1
CCCC          AllPatern(NAktIndx) = Patern1(Indx, Indx2)
CCCC        ENDIF
CCCC      ENDDO
CCCC    ENDDO
CCCC  ENDDO
C
CCCC  DO Indx = 1, MIN(MaxState1, MaxState2)
CCCC    DO Indx2 = 1, NExpVect2(Indx)
CCCC      DO Indx3 = 1, NAktIndx
CCCC        IF (Patern2(Indx, Indx2) .EQ. AllPatern(Indx3)) THEN
CCCC          PatVal2(Indx3, Indx) = StExpCoeff2(Indx, Indx2)
CCCC        ENDIF
CCCC      ENDDO
CCCC    ENDDO
CCCC    DO Indx2 = 1, NExpVect1(Indx)
CCCC      DO Indx3 = 1, NAktIndx
CCCC        IF (Patern1(Indx, Indx2) .EQ. AllPatern(Indx3)) THEN
CCCC          PatVal1(Indx3, Indx) = StExpCoeff1(Indx, Indx2)
CCCC        ENDIF
CCCC      ENDDO
CCCC    ENDDO
CCCC  ENDDO
C
CCCC  OPEN(21, file='d_seraz1.dbg')
CCCC  DO Indx = 1, NAktIndx
CCCC    WRITE (21, '(I3, 3X, A10, 5000F14.7)') Indx, AllPatern(Indx), 
CCCC &                      (PatVal1(Indx, VnIndx), 
CCCC &                       VnIndx = 1, MIN(MaxState1, MaxState2))
CCCC    WRITE (21, '(I3, 3X, A10, 5000F14.7)') Indx, AllPatern(Indx), 
CCCC &                      (PatVal2(Indx, VnIndx), 
CCCC &                       VnIndx = 1, MIN(MaxState1, MaxState2))
CCCC    WRITE (21, *)
CCCC  ENDDO
CCCC  CLOSE(21)
C
C
C       Erasing history of coupling diagonal 
      DO Indx = 1, MIN(MaxState1, MaxState2)
        OldDiag(Indx) = 1.0D0
      ENDDO
C
C       Calculating couplings between states
C
C
 9500 CONTINUE
CCCCCCCCCxx
C
      CALL CalcSpqMatrix(NOrb, NFunc, OrbCoeff1, OrbCoeff2, Overlap,Spq)
      CALL CalcCoups(MaxState1, MaxState2, NExpVect1, NExpVect2, 
     &               StAlpha1, StAlpha2, StBeta1, StBeta2, NOrbCIMax, 
     &               Spq, StExpCoeff1, StExpCoeff2, Coupling)
C
C       Sometimes GAMESS fails to give same CSFs as calculation before, resulting into non-sense coupling 
C       (part of the CI vector coeffs have opposite signs). If there is that problem, only for the sake of 
C       this step calculation make them same sign (even if they make CIs not orthogonal - that sometimes happens)
C       (NOTE: it doesnt matter if there is state position change, programs finds the state with the most overlap 
C       and tries to change their signs)
C
C       Find the highes coupling
      DO Indx = 1, MIN(MaxState1, MaxState2)
        MaxIndx(Indx) = Indx
        IniSigns(Indx) = .FALSE.
      ENDDO
C
C       Check if there is any below threshold 0.95
      IniBelow = .FALSE.
      DO Indx = 1, MIN(MaxState1, MaxState2)
        IF (Coupling(Indx, MaxIndx(Indx)) .LT. 0.95D0) THEN
          IniBelow = .TRUE.
          IniSigns(Indx) = .TRUE.
        ENDIF
      ENDDO
C
      IF (IniBelow) THEN
C         Check if CSF are OK (sometimes character of state changes, so when that happen we dont do the correction)
        OPEN(21, file='d_csfdiff.tmp')
        DO Indx =1, NStCSF
          SumInt1(Indx) = 0.0D0
          DO Indx2 = 1, NCSF1(Indx)
            IniNotThere = .TRUE.
            DO Indx3 = 1, NCSF2(Indx)
              IF (CSFIde1(Indx, Indx2) .EQ. CSFIde2(Indx, Indx3)) THEN
                IniNotThere = .FALSE.
              ENDIF
            ENDDO
            IF (IniNotThere) THEN
              SumInt1(Indx) = SumInt1(Indx)+CSFInt1(Indx, Indx2)**2
            ENDIF
          ENDDO
          SumInt2(Indx) = 0.0D0
          DO Indx2 = 1, NCSF2(Indx)
            IniNotThere = .TRUE.
            DO Indx3 = 1, NCSF1(Indx)
              IF (CSFIde2(Indx, Indx2) .EQ. CSFIde1(Indx, Indx3)) THEN
                IniNotThere = .FALSE.
              ENDIF
            ENDDO
            IF (IniNotThere) THEN
              SumInt2(Indx) = SumInt2(Indx)+CSFInt2(Indx, Indx2)**2
            ENDIF
          ENDDO
          WRITE (21, '(I3, 2F14.7)') Indx, SumInt1(Indx), SumInt2(Indx) 
        ENDDO
        CLOSE(21)
C         Check the sums of different CSFs intensities, if it is up the limit dont resign CIs
        DO Indx = 1, NStCSF
          IF ((SumInt1(Indx) .GT. 0.15D0) .OR. 
     &        (SumInt2(Indx) .GT. 0.15D0)) THEN
            IniSigns(Indx) = .FALSE.
CCC       ELSE
CCC         IniSigns(Indx) = .TRUE.
          ENDIF
        ENDDO
C         Check for the opposite signs in CSFs if so trigger resign
        DO Indx = 1, NStCSF
          IniThisState = .FALSE.
          SumInt1(Indx) = 0.0D0
          CSFDiff(Indx) = 0.0D0
          DO Indx2 = 1, NCSF2(Indx)
            DO Indx3 = 1, NCSF1(Indx)
              IF (CSFIde2(Indx, Indx2) .EQ. CSFIde1(Indx, Indx3)) THEN
                CSFDiff(Indx) = CSFDiff(Indx)+
     &                          ABS(CSFInt2(Indx, Indx2)**2-
     &                              CSFInt1(Indx, Indx3)**2)
                IF (CSFInt2(Indx, Indx2)*CSFInt1(Indx, Indx3) 
     &              .LT. 0.0D0) THEN
                  IniBelow = .TRUE.
                  IniThisState = .TRUE.
                  SumInt1(Indx) = SumInt1(Indx)+ABS(CSFInt2(Indx,Indx2)*
     &                                              CSFInt1(Indx,Indx3))
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          IF (IniThisState .AND. (SumInt1(Indx) .GT. 0.15D0)) THEN
            IniSigns(Indx) = .TRUE.
          ENDIF
        ENDDO
C
        DO Indx = 1, NStCSF
          IF (CSFDiff(Indx) .GT. 0.2D0) THEN
            IniSigns(Indx) = .FALSE.
          ENDIF
        ENDDO
C
        OPEN(21, file='d_report.dbg')
        WRITE (21, *) 'Below: ', IniBelow
        DO Indx = 1, NStCSF
          WRITE (21, '(I3, L, 2F14.7)') Indx, 
     &               IniSigns(Indx), SumInt1(Indx), CSFDiff(Indx)
        ENDDO
        CLOSE(21)
      ENDIF
C
C       If there is any coupling below the threshold value, try to change CI signs and recalc Coups
      IF (IniBelow) THEN
        DO Indx = 1, MaxState2
          DO Indx2 = 1, NExpVect2(Indx)
            StExpCoeff2B(Indx, Indx2) = StExpCoeff2(Indx, Indx2)
          ENDDO
        ENDDO
        DO Indx = 1, MIN(MaxState1, MaxState2)
          IF (IniSigns(Indx)) THEN
            DO Indx2 = 1, NExpVect2(Indx)
              DO Indx3 = 1, NExpVect1(MaxIndx(Indx))
                IF (Patern2(Indx, Indx2) .EQ. 
     &              Patern1(MaxIndx(Indx), Indx3)) THEN
                  IF (StExpCoeff2B(Indx, Indx2)*
     c                StExpCoeff1(MaxIndx(Indx), Indx3) .LT. 0.0D0) THEN
                    StExpCoeff2B(Indx, Indx2) = -1.0D0*
     &                                         StExpCoeff2B(Indx, Indx2)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CALL CalcCoups(MaxState1, MaxState2, NExpVect1, NExpVect2, 
     &                  StAlpha1, StAlpha2, StBeta1, StBeta2, NOrbCIMax,
     &                  Spq, StExpCoeff1, StExpCoeff2B, NewCoup)
C
CCCCC   WRITE (StateFormat,'("(2I4, A5," I4, "F6.2,  A3," I4, "F6.2)")')
CCCCC&                          MaxState1, MaxState2
CCCCC   WRITE (*, StateFormat) MaxState1, MaxState2, '---',
CCCCC&                      (StateSpin1(VnIndx), VnIndx = 1, MaxState1),
CCCCC&                 '|', (StateSpin2(VnIndx), VnIndx = 1, MaxState2)
CCCCC   WRITE (*, 9906) (VnIndx, VnIndx = 1, MaxState1)
CCCCC   DO Indx = 1, MaxState2
CCCCC     WRITE (*, 9900) Indx, (NewCoup(VnIndx, Indx), 
CCCCC&                                          VnIndx = 1, MaxState1)
CCCCC   ENDDO
CCCCC   WRITE (*, *) 
C
        SumUp = 0.0D0
        DO Indx = 1, MIN(MaxState1, MaxState2)
          SumUp = SumUp+Coupling(Indx, Indx)-NewCoup(Indx, Indx)
        ENDDO
        IF (SumUp .LT. 0.0D0) THEN
          DO Indx = 1, MaxState1
            DO Indx2 = 1, MaxState2
              Coupling(Indx, Indx2) = NewCoup(Indx, Indx2)
            ENDDO
          ENDDO
          IniBelow = .FALSE.
        ELSE
          IniBelow = .TRUE.
        ENDIF
      ENDIF
C
C       Sometimes states change their positions. If the one up there didnt help, lets find the state with highest coupling
C       and switch CI signs by that.
CCx   DO Indx = 1, MIN(MaxState1, MaxState2)
CCx     AktMax = 0.0D0
CCx     DO Indx2 = 1, MIN(MaxState1, MaxState2)
CCx       IF (AktMax .LT. ABS(Coupling(Indx, Indx2))) THEN
CCx         AktMax = ABS(Coupling(Indx, Indx2))
CCx         MaxIndx(Indx) = Indx2
CCx       ENDIF
CCx     ENDDO
CCx   ENDDO
C
C       Check if there is any below threshold 0.98
CCx   IniBelow = .FALSE.
CCx   DO Indx = 1, MIN(MaxState1, MaxState2)
CCx     IF (Coupling(Indx, MaxIndx(Indx)) .LT. 0.98D0) THEN
CCx       IniBelow = .TRUE.
CCx     ENDIF
CCx   ENDDO
CC
CCx   IniBelow = .FALSE.
CC
C
C       If there is any coupling below the threshold value, try to change CI signs and recalc Coups
CCx   IF (IniBelow) THEN
CCx     DO Indx = 1, MaxState2
CCx       DO Indx2 = 1, NExpVect2(Indx)
CCx         StExpCoeff2C(Indx, Indx2) = StExpCoeff2(Indx, Indx2)
CCx       ENDDO
CCx     ENDDO
CCx     DO Indx = 1, MIN(MaxState1, MaxState2)
CCx       DO Indx2 = 1, NExpVect2(Indx)
CCx         DO Indx3 = 1, NExpVect1(MaxIndx(Indx))
CCx           IF (Patern2(Indx, Indx2) .EQ. 
CCx  &            Patern1(MaxIndx(Indx), Indx3)) THEN
CCx             IF (StExpCoeff2C(Indx, Indx2)*
CCx  c              StExpCoeff1(MaxIndx(Indx), Indx3) .LT. 0.0D0) THEN
CCx               StExpCoeff2C(Indx, Indx2) = -1.0D0*
CCx  &                                         StExpCoeff2C(Indx, Indx2)
CCx             ENDIF
CCx           ENDIF
CCx         ENDDO
CCx       ENDDO
CCx     ENDDO
CCx     CALL CalcCoups(MaxState1, MaxState2, NExpVect1, NExpVect2, 
CCx  &                  StAlpha1, StAlpha2, StBeta1, StBeta2, NOrbCIMax,
CCx  &                  Spq, StExpCoeff1, StExpCoeff2C, BestCoup)
C
CCCCC   WRITE (StateFormat,'("(2I4, A5," I4, "F6.2,  A3," I4, "F6.2)")')
CCCCC&                          MaxState1, MaxState2
CCCCC   WRITE (*, StateFormat) MaxState1, MaxState2, '---',
CCCCC&                      (StateSpin1(VnIndx), VnIndx = 1, MaxState1),
CCCCC&                 '|', (StateSpin2(VnIndx), VnIndx = 1, MaxState2)
CCCCC   WRITE (*, 9906) (VnIndx, VnIndx = 1, MaxState1)
CCCCC   DO Indx = 1, MaxState2
CCCCC     WRITE (*, 9900) Indx, (BestCoup(VnIndx, Indx), 
CCCCC&                                          VnIndx = 1, MaxState1)
CCCCC   ENDDO
CCCCC   WRITE (*, *)
C
CCx     SumUp = 0.0D0
CCx     DO Indx = 1, MIN(MaxState1, MaxState2)
CCx       SumUp = SumUp+Coupling(Indx, Indx)-NewCoup(Indx, Indx)
CCx     ENDDO
CCx     IF (SumUp .LT. 0.0D0) THEN
CCx       DO Indx = 1, MaxState1
CCx         DO Indx2 = 1, MaxState2
CCx           Coupling(Indx, Indx2) = BestCoup(Indx, Indx2)
CCx         ENDDO
CCx       ENDDO
CCx     ENDIF
CCx   ENDIF
C
C
C     *****  END OF DATA PROCESSING  *****
C
C       Write down the output TD Couplings
      WRITE (StateFormat,'( "(2I4, A5," I4, "F6.2,  A3," I4, "F6.2)")')
     &                        MaxState1, MaxState2
      WRITE (*, StateFormat) MaxState1, MaxState2, '---',
     &                      (StateSpin1(VnIndx), VnIndx = 1, MaxState1),
     &                 '|', (StateSpin2(VnIndx), VnIndx = 1, MaxState2)
      WRITE (*, 9906) (VnIndx, VnIndx = 1, MaxState1)
      DO Indx = 1, MaxState2
        WRITE (*, 9900) Indx, (Coupling(VnIndx, Indx), 
     &                                        VnIndx = 1, MaxState1)
      ENDDO
      WRITE (*, *) 
C
C      Write down the output: state energies
      OPEN (20, file='fomo_sten.dat', ACCESS='SEQUENTIAL')
      WRITE (20, 9905) MaxState2
      DO Indx = 1, MaxState2
        WRITE (20, 9902) StateEn2(Indx)
      ENDDO
      CLOSE (20)
C
CCCC  DO Indx = 1, MIN(MaxState1, MaxState2)
CCCC    DO Indx2 = 1, NExpVect2(Indx)
CCCC      DO Indx3 = 1, NAktIndx
CCCC        IF (Patern2(Indx, Indx2) .EQ. AllPatern(Indx3)) THEN
CCCC          PatVal2(Indx3, Indx) = StExpCoeff2(Indx, Indx2)
CCCC        ENDIF
CCCC      ENDDO
CCCC    ENDDO
CCCC    DO Indx2 = 1, NExpVect1(Indx)
CCCC      DO Indx3 = 1, NAktIndx
CCCC        IF (Patern1(Indx, Indx2) .EQ. AllPatern(Indx3)) THEN
CCCC          PatVal1(Indx3, Indx) = StExpCoeff1(Indx, Indx2)
CCCC        ENDIF
CCCC      ENDDO
CCCC    ENDDO
CCCC  ENDDO
C
CCCC  OPEN(21, file='d_seraz2.dbg')
CCCC  DO Indx = 1, NAktIndx
CCCC    WRITE (21, '(I3, 3X, A10, 5000F14.7)') Indx, AllPatern(Indx), 
CCCC &                      (PatVal1(Indx, VnIndx), 
CCCC &                       VnIndx = 1, MIN(MaxState1, MaxState2))
CCCC    WRITE (21, '(I3, 3X, A10, 5000F14.7)') Indx, AllPatern(Indx), 
CCCC &                      (PatVal2(Indx, VnIndx), 
CCCC &                       VnIndx = 1, MIN(MaxState1, MaxState2))
CCCC    WRITE (21, *)
CCCC  ENDDO
CCCC  CLOSE(21)
C
C     *****  FORMAT SECTION  *****
C
 9900 FORMAT(I4, 5000F14.8)
C
 9901 FORMAT(5000I14)
C
 9902 FORMAT(5001F14.8)
C
 9903 FORMAT(4I2, "  |  ", 4I2, "  |  ", F14.8)
C
 9904 FORMAT(3I6, 5000I3)
C
 9905 FORMAT(2I4, F14.8)
C
 9906 FORMAT(I13, 5000I14)
C
 9907 FORMAT(2I4, 5000F14.8)
C
99999 CONTINUE
C
C     *****  END OF FORMAT SECTION  *****
C
C     STOP
C
      END
C
C
C     Fortran code for finding determinant of the NxN matrix
C     Stolen from: http://www.dreamincode.net/code/snippet1273.htm
C     Changed a bit...
      DOUBLE PRECISION FUNCTION FindDet(inpMatrix, n, nmax)
C
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(nmax,nmax) :: inpMatrix
      DOUBLE PRECISION, DIMENSION(n,n) :: matrix
      INTEGER, INTENT(IN) :: n, nmax
      DOUBLE PRECISION :: m, temp
      INTEGER :: i, j, k, l
      LOGICAL :: DetExists = .TRUE.
      l = 1
      DO i = 1, n
        DO j = 1, n
          matrix (i, j) = inpMatrix(i, j)
        END DO
      END DO
C       Convert to upper triangular form
      DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
          DetExists = .FALSE.
          DO i = k+1, n
            IF (matrix(i,k) /= 0) THEN
              DO j = 1, n
                temp = matrix(i,j)
                matrix(i,j)= matrix(k,j)
                matrix(k,j) = temp
              END DO
              DetExists = .TRUE.
              l=-l
              EXIT
            ENDIF
          END DO
          IF (DetExists .EQV. .FALSE.) THEN
            FindDet = 0
            return
          END IF
        ENDIF
        DO j = k+1, n
          m = matrix(j,k)/matrix(k,k)
          DO i = k+1, n
            matrix(j,i) = matrix(j,i) - m*matrix(k,i)
          END DO
        END DO
      END DO
C       Calculate determinant by finding product of diagonal elements
      FindDet = l
      DO i = 1, n
        FindDet = FindDet * matrix(i,i)
      END DO
             
      END FUNCTION FindDet 
C
C
C
      SUBROUTINE CalcCoups(MxSt1, MxSt2, NEV1, NEV2, A1, A2, B1, B2, 
     &                     NOMax, S, EC1, EC2, Coups)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Calculates couplings.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  FUNCTIONS:   DOUBLE PRECISION                                       C
C   FindDet(Matrix, Size, MaxSize)                                     C
C      -calculates determinant of square matrix of size "Size"         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Input/Output:                                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                15.03.2013                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         DECLARATION
C
      IMPLICIT NONE
C
C         DECLARATION - PARAMETERS
C
      INTEGER pMxSt
        PARAMETER (pMxSt = 30)
C
      INTEGER pMxO
        PARAMETER (pMxO = 650)
C
      INTEGER pMxCI
        PARAMETER (pMxCI = 5500)
C
C         DECLARATION - FUNCTION
C
      DOUBLE PRECISION FindDet
      EXTERNAL FindDet
C
C         DECLARATION - VARIABLES
C
      INTEGER MxSt1, MxSt2, NOMax
C
      INTEGER NEV1(pMxSt), NEV2(pMxSt)
C
      INTEGER A1(pMxSt, pMxCI, pMxO), A2(pMxSt, pMxCI, pMxO)
C
      INTEGER B1(pMxSt, pMxCI, pMxO), B2(pMxSt, pMxCI, pMxO)
C
      DOUBLE PRECISION S(pMxO, pMxO)
C
      DOUBLE PRECISION EC1(pMxSt, pMxCI), EC2(pMxSt, pMxCI)
C
      DOUBLE PRECISION Coups(pMxSt, pMxSt)
C  ____________________________________________________________________
C
      INTEGER pI, pI2, pI3, pI4, pI5, pI6, pVI
C
      INTEGER II, CountRe, AI
C
      INTEGER OIA1(pMxO), OIA2(pMxO), OIB1(pMxO), OIB2(pMxO)
C
      INTEGER NEA, NEB
C
      DOUBLE PRECISION SOA(pMxO, pMxO), SOB(pMxO, pMxO)
C
      DOUBLE PRECISION DetA, DetB, Det(pMxCI, pMxCI)
C
      DOUBLE PRECISIOn WrDetA(pMxSt, pMxSt), WrDetB(pMxSt, pMxSt)
C
      DOUBLE PRECISION AMax
C
      LOGICAL Re
C
C-----------------------------------------------------------------------
C
C         MAIN  CalcCoups  
C
C       Calculating number of electrons
      NEA = 0
      NEB = 0
      CountRe = 0
      DO pI = 1, NOMax
        NEA = NEA+A1(1, 1, pI)
        NEB = NEB+B1(1, 1, pI)
      ENDDO
C
CCC   OPEN(21, file='zk_ci1.dat')
CCC   DO pI = 1, MxSt1 
CCC     WRITE (21, '(I3, 5000F14.8)') pI
CCC     DO pI2 = 1, NEV1(pI)
CCC       WRITE (21, '(7I2, " |", 7I2, F14.8)') 
CCC  &                (A1(pI, pI2, pVI), 
CCC  &                            pVI = 15, NOMax), 
CCC  &                (B1(pI, pI2, pVI), 
CCC  &                         pVI = 15, NOMax), 
CCC  &                 EC1(pI, pI2) 
CCC     ENDDO
CCC   ENDDO
CCC   CLOSE(21)
CCC   OPEN(21, file='zk_ci2.dat')
CCC   DO pI = 1, MxSt2 
CCC     WRITE (21, '(I3, 5000F14.8)') pI
CCC     DO pI2 = 1, NEV2(pI)
CCC       WRITE (21, '(7I2, " |", 7I2, F14.8)') 
CCC  &                (A2(pI, pI2, pVI), 
CCC  &                            pVI = 15, NOMax), 
CCC  &                (B2(pI, pI2, pVI), 
CCC  &                         pVI = 15, NOMax), 
CCC  &                 EC2(pI, pI2) 
CCC     ENDDO
CCC   ENDDO
CCC   CLOSE(21)
C
CCC   OPEN(22, file='d_detA.dbg')
CCC   OPEN(23, file='d_detB.dbg')
CCC   OPEN(24, file='d_det.dbg')
C
33300 CONTINUE
      Re = .FALSE.
CCCCCCCCCxx
      DO pI = 1, MxSt2
        DO pI2 = 1, MxSt1
CCC       WRITE (22, *) pI, pI2
CCC       WRITE (22, *) " "
CCC       WRITE (23, *) pI, pI2
CCC       WRITE (23, *) " "
CCC       WRITE (24, *) pI, pI2
CCC       WRITE (24, *) " "
C         First we need to calculate the overlap between Slaters of different state
          DO pI3 = 1, NEV2(pI)
C           At first, we need to know which orbitals in expansion vector are occupied
C              time: t, Alpha
            II = 1 
            DO pI5 = 1, NOMax
              IF (A2(pI, pI3, pI5) .EQ. 1) THEN
                OIA2(II) = pI5
                II = II+1
              ENDIF
            ENDDO
C              time: t, Beta 
            II = 1 
            DO pI5 = 1, NOMax
              IF (B2(pI, pI3, pI5) .EQ. 1) THEN
                OIB2(II) = pI5
                II = II+1
              ENDIF
            ENDDO
C
C            Now we look at the t+h part
            DO pI4 = 1, NEV1(pI2)
              IF (ABS(EC1(pI, pI3)*EC2(pI2, pI4)) .LT. 1.0D-5) THEN
                Det(pI3, pI4) = 0.0D0
                CYCLE
              ENDIF
C               Which orbitals are occupied in the t+h ci expansion vector?
C                 time: t+h, Alpha
              II = 1
              DO pI6 = 1, NOMax
                IF (A1(pI2, pI4, pI6) .EQ. 1) THEN
                  OIA1(II) = pI6
                  II = II+1
                ENDIF
              ENDDO
C                time: t+h, Beta 
              II = 1 
              DO pI6 = 1, NOMax
                IF (B1(pI2, pI4, pI6) .EQ. 1) THEN
                  OIB1(II) = pI6
                  II = II+1
                ENDIF
              ENDDO
C
C               Clearing overlap matrix for actual states and state vectors
              DO pI5 = 1, NEA
                DO pI6 = 1, NEA
                  SOA(pI5, pI6) = 0.0D0
                ENDDO
              ENDDO
              DO pI5 = 1, NEB
                DO pI6 = 1, NEB
                  SOB(pI5, pI6) = 0.0D0
                ENDDO
              ENDDO
C
C               Everything is ready to calculate the wf overlap between expansion vectors
C                  Alpha part
              DO pI5 = 1, NEA
                DO pI6 = 1, NEA
                  SOA(pI5, pI6) = S(OIA2(pI5), OIA1(pI6))
                ENDDO
              ENDDO
              DetA = FindDet(SOA, NEA, pMxO)
              WrDetA(pI3, pI4) = DetA
C
C                  Beta part
              DO pI5 = 1, NEB
                DO pI6 = 1, NEB
                  SOB(pI5, pI6) = S(OIB2(pI5), OIB1(pI6))
                ENDDO
              ENDDO
              DetB = FindDet(SOB, NEB, pMxO)
              WrDetB(pI3, pI4) = DetB
C
C                  Making Alpha and Beta part together
              Det(pI3, pI4) = DetA*DetB
C
            ENDDO
CCCCCCCCCxx
CCC         WRITE (22, 1900) pI3, (WrDetA(pI3, pVI), 
CCC  &                               pVI = 1, NEV1(pI2))
CCC         WRITE (23, 1900) pI3, (WrDetB(pI3, pVI), 
CCC  &                               pVI = 1, NEV1(pI2))
CCC         WRITE (24, 1900) pI3, (Det(pI3, pVI), 
CCC  &                               pVI = 1, NEV1(pI2))
CCCCCCCCCxx
          ENDDO
C
C           Now we got overlap matrix between expansion vectors, so we need to combine it with expansion vector weights
          Coups(pI2, pI) = 0.0D0
          DO pI3 = 1, NEV2(pI)
            DO pI4 = 1, NEV1(pI2)
              Coups(pI2, pI) = Coups(pI2, pI)+EC2(pI, pI3)
     &                         *EC1(pI2, pI4)*Det(pI3, pI4)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C       Now check the signs of elemnts which have the highest overlap between all othre states, 
C       if they are negative switch all signs of that CI and recalc that
      DO pI = 1, MIN(MxSt1, MxSt2)
CCCCCC  AMax = 0.0D0
CCCCCC  DO pI2 = 1, Min(MxSt1, MxSt2)
CCCCCC    IF (AMax .LT. ABS(Coups(pI, pI2))) THEN
CCCCCC      AI = pI2
CCCCCC      AMax = ABS(Coups(pI, pI2))
CCCCCC    ENDIF
CCCCCC  ENDDO
CCCCCC  IF (Coups(pI, AI) .LT. 0.0D0) THEN
        IF (Coups(pI, pI) .LT. 0.0D0) THEN
          Re = .TRUE.
CCCCCC    CountRe = CountRe+1
CCCCCC    IF (CountRe .GE. 2) THEN
CCCCCC      Re = .FALSE.
CCCCCC    ELSE
            DO pI2 = 1, NEV2(pI)
              EC2(pI, pI2) = -1.0D0*EC2(pI, pI2)
            ENDDO
CCCCCC    ENDIF
        ENDIF
      ENDDO
      IF (Re) THEN
        GOTO 33300
      ENDIF
C
CCCCCCCCCxx
CCC   CLOSE(22)
CCC   CLOSE(23)
CCC   CLOSE(24)
CCCCCCCCCxx
C
 1900 FORMAT(I4, 5000F14.8)
C
      END
C
C
C
      SUBROUTINE OrbDotProd(N, NF, C1, C2, S)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Calculates dot product of orbitals.                                  C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  FUNCTIONS:   DOUBLE PRECISION                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Input/Output:                                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                15.03.2013                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         DECLARATION
C
      IMPLICIT NONE
C
C         DECLARATION - PARAMETERS
C
      INTEGER pMxO
        PARAMETER (pMxO = 650)
C
C         DECLARATION - VARIABLES
C
      INTEGER N, NF
C
      DOUBLE PRECISION C1(pMxO, pMxO), C2(pMxO, pMxO), S(pMxO)
C  ____________________________________________________________________
C
      INTEGER pI, pI2
C
C-----------------------------------------------------------------------
C
C         MAIN  OrbDotProd 
C
      DO pI = 1, N
        S(pI) = 0.0D0
        DO PI2 = 1, NF
          S(pI) = S(pI)+C1(pI, PI2)*C2(pI, pI2)
        ENDDO
      ENDDO
C
C
      END
C
C
C
      SUBROUTINE CalcSpqMatrix(NO, NF, O1, O2, S, Spq)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Calculates Spq matrix.                                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Input/Output:                                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                15.04.2013                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         DECLARATION
C
      IMPLICIT NONE
C
C         DECLARATION - PARAMETERS
C
      INTEGER pMxO
        PARAMETER (pMxO = 650)
C
C         DECLARATION - VARIABLES
C
      INTEGER NO, NF
C
      DOUBLE PRECISION O1(pMxO, pMxO), O2(pMxO, pMxO), S(2*pMxO, 2*pMxO)
C
      DOUBLE PRECISION Spq(pMxO, pMxO)
C  ____________________________________________________________________
C
      INTEGER pI, pI2, pI3, pI4, pAktI
C
      DOUBLE PRECISION pAktMax
C
      LOGICAL IniAgain
C
C-----------------------------------------------------------------------
C
C         MAIN  CalcSpqMatrix
C
65000 CONTINUE
      IniAgain = .FALSE.
C
C     Calculating the Spq matrix (it should be the same for alpha and beta part of wf)
      DO pI = 1, NO
        DO pI2 = 1, NO
C            Erasing Spq element (to be sure)
          Spq(pI, pI2) = 0.0D0
C             Cycles around AO functions
          DO pI3 = 1, NF
            DO pI4 = 1, NF
              Spq(pI, pI2) = Spq(pI, pI2)+
     &                       O2(pI, pI3)*O1(pI2, pI4)*S(pI3 + NF, pI4)
            ENDDO
          ENDDO
C         WRITE (*, '(2I3, F14.8)') pI, pI2, Spq(pI, pI2)
        ENDDO
      ENDDO
C
CCC   DO pI = 1, NO
CCC     pAktMax = 0.0D0
CCC     DO pI2 = 1, NO
CCC       IF (ABS(Spq(pI, pI2)) .GT. pAktMax) THEN
CCC         pAktI = pI2
CCC         pAktMax = ABS(Spq(pI, pI2))
CCC       ENDIF
CCC     ENDDO
CCC     IF (Spq(pI, pAktI) .LT. 0.0D0) THEN
CCC       IniAgain = .TRUE.
CCC       DO pI2 = 1, NF
CCC         O2(pI, pI2) = -1.0D0*O2(pI, pI2)
CCC       ENDDO
CCC     ENDIF
CCC   ENDDO
CCC   IF (IniAgain) THEN
CCC     GOTO 65000
CCC   ENDIF
C
      OPEN(21,file='d_spqmat.dbg')
      DO pI = 1, NO
        WRITE (21, '(I4, 5000F14.8)') pI, (Spq(pI, pI2), pI2 = 1, NO)
      ENDDO
      CLOSE(21)
C
C
      END
C
C
C
