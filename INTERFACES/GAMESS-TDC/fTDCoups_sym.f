      PROGRAM fTDCoups_sym
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Fortran part of the TDCoups script for time-derivative couplings for C
C FOMO calculation in GAMESS. After calculation of asymetric formula   C
C this program calculates symetric dot products                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                05.02.2013                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         DECLARATION
C
      IMPLICIT NONE
C
C         DECLARATION - PARAMETERS
C
      INTEGER MaxSt
        PARAMETER (MaxSt = 100)
C
C         DECLARATION - FUNCTION
C
C
C         DECLARATION - VARIABLES
C
      INTEGER Indx, VnIndx, Indx2
C
      INTEGER NStates, CountSt2, ChangeStB(MaxSt), NChangesB
C
      INTEGER CountSt1, ChangeStA(MaxSt), NChangesA
C
      INTEGER Sign1(MaxSt, MaxSt), Sign2(MaxSt, MaxSt)
C
      DOUBLE PRECISION TDCoups1(MaxSt, MaxSt), TDCoups2(MaxSt, MaxSt)
C
      DOUBLE PRECISION Sigma1(MaxSt, MaxSt), Sigma2(MaxSt, MaxSt)
C
      DOUBLE PRECISION Sigma(MaxSt, MaxSt)
C
      DOUBLE PRECISION SigDiff1(MaxSt, MaxSt), SigDiff2(MaxSt, MaxSt)
C
      DOUBLE PRECISION OldCoup(MaxSt, MaxSt), Ratio
C
      DOUBLE PRECISION WrSigma(MaxSt, MaxSt)
C
      CHARACTER junk
C
      LOGICAL UNITOK, UNITOP
C
      LOGICAL IniPlus
C
C-----------------------------------------------------------------------
C
C         MAIN 
C     
C     *****  READING INPUT DATA  *****
C       OPEN str.state.tmp provided by r.gamess to read number of states
      OPEN (11, file='str.state.tmp', STATUS='OLD')
C       read number of states required
      READ (11, *) NStates
C
      CLOSE (11)
C
C       OPEN before_asym_tdcoups.dat, tdcoups from previous step     
      OPEN (11, file='before_asym_tdcoups.dat', STATUS='OLD')
C
      READ (11, *) junk, CountSt2
      READ (11, *)
C       read all dot products from previous step (only one that matters)
      DO Indx = 1, NStates
        READ (11, *) junk, (TDCoups1(Indx, VnIndx), VnIndx = 1, NStates)
      ENDDO 
C       if there was sign change in last step 
CCC   DO Indx = 1, CountSt2-NStates
CCC     READ (11, *)
CCC   ENDDO
CCC   READ (11, *) NChangesB
CCC   DO Indx = 1, NChangesB
CCC     READ (11, *) ChangeStB(Indx)
CCC   ENDDO
C
      CLOSE (11)
C
C       OPEN asym_tdcoups.dat, tdcoups from actual step
      OPEN (11, file='asym_tdcoups.dat', STATUS='OLD')
C
      READ (11, *) junk, CountSt1
      READ (11, *)
C       read all dot products from actual step (only one that matters)
      DO Indx = 1, NStates
        READ (11, *) junk, (TDCoups2(Indx, VnIndx), VnIndx = 1, NStates)
      ENDDO 
C
CCC   DO Indx = 1, CountSt1-NStates
CCC     READ (11, *)
CCC   ENDDO
CCC   READ (11, *) NChangesA
CCC   DO Indx = 1, NChangesA
CCC     READ (11, *) ChangeStA(Indx)
CCC   ENDDO
C
      CLOSE (11)
C
C
C     *****  END OF READING INPUT DATA  *****
C
C     *****  DATA PROCESSING  *****
C       Before we calculate sigma we need to flip signs 
CCC   IF (NChangesB .GT. 0) THEN
CCC     DO Indx = 1, NChangesB
CCC       DO Indx2 = 1, NStates
CCC         TDCoups1(Indx2, ChangeStB(Indx)) = 
CCC  &                        -1.0D0*TDCoups1(Indx2, ChangeStB(Indx))
CCC         TDCoups1(ChangeStB(Indx), Indx2) = 
CCC  &                        -1.0D0*TDCoups1(ChangeStB(Indx), Indx2)
CCC       ENDDO
CCC     ENDDO
CCC   ENDIF
CCC   IF (NChangesA .GT. 0) THEN
CCC     DO Indx2 = 1, NStates
CCC       TDCoups2(Indx2, ChangeStA(Indx)) = 
CCC  &                      -1.0D0*TDCoups2(Indx2, ChangeStA(Indx))
CCC       TDCoups2(ChangeStA(Indx), Indx2) = 
CCC  &                      -1.0D0*TDCoups2(ChangeStA(Indx), Indx2)
CCC     ENDDO
CCC   ENDIF
C
CCC   DO Indx = 1, NStates
CCC     WRITE (*, '(I3, 10F14.7)') Indx, (TDCoups1(Indx, VnIndx), 
CCC  &  VnIndx = 1, NStates)
CCC   ENDDO
C       Calculation of sigma_JI(t-dt/2): Sigma2
C                      sigma_JI(t-3dt/2): Sigma1
C          - time factor not incorporated (it is in ABIN)
      INQUIRE(file='sym_data.tmp', EXIST=UNITOK)
      IF (UNITOK) THEN
        OPEN (11, file='sym_data.tmp', STATUS='OLD')
        DO Indx = 1, NStates
          READ (11, *) (OldCoup(Indx, VnIndx), VnIndx = 1, NStates)
        ENDDO
        DO Indx = 1, NStates
          READ (11, *) (Sign1(Indx, VnIndx), VnIndx = Indx, NStates)
        ENDDO
        DO Indx = 1, NStates
          READ (11, *) (Sign2(Indx, VnIndx), VnIndx = Indx, NStates)
        ENDDO
        CLOSE (11)
      ENDIF
C
      DO Indx = 1, NStates
        DO Indx2 = Indx+1, NStates
          SigDiff2(Indx, Indx2) = (TDCoups2(Indx2, Indx)-
     &                             TDCoups2(Indx, Indx2))/2.0D0
          SigDiff1(Indx, Indx2) = (TDCoups1(Indx2, Indx)-
     &                             TDCoups1(Indx, Indx2))/2.0D0
C
          IF (ABS(SigDiff2(Indx, Indx2)-SigDiff1(Indx, Indx2)) .LT. 
     &        ABS(SigDiff2(Indx, Indx2)+SigDiff1(Indx, Indx2))) THEN
            IniPlus = .FALSE.
          ELSE
            IniPlus = .TRUE.
          ENDIF
C
          IF (IniPlus) THEN
            Sigma(Indx, Indx2) = 1.0D0/2.0D0*
     &              (3.0D0*SigDiff2(Indx, Indx2)+SigDiff1(Indx, Indx2))
          ELSE
            Sigma(Indx, Indx2) = 1.0D0/2.0D0*
     &              (3.0D0*SigDiff2(Indx, Indx2)-SigDiff1(Indx, Indx2))
          ENDIF
          WrSigma(Indx, Indx2) = Sigma(Indx, Indx2)
C
          IF (UNITOK) THEN
            Ratio = ABS(Sigma(Indx, Indx2)-OldCoup(Indx, Indx2))/
     &              ABS(Sigma(Indx, Indx2)+OldCoup(Indx, Indx2))
            IF (Ratio .GT. 1.5D0) THEN
              Sign1(Indx, Indx2) = -1.0D0*Sign1(Indx, Indx2)
            ENDIF
            Sigma(Indx, Indx2) = Sigma(Indx, Indx2)*Sign1(Indx, Indx2)
            IF ((.NOT. IniPlus) .AND. (Ratio .GT. 1.5D0)) THEN
              Sign2(Indx, Indx2) = -1.0D0*Sign2(Indx, Indx2)
            ENDIF
            Sigma(Indx, Indx2) = Sigma(Indx, Indx2)*Sign2(Indx, Indx2)
          ELSE
            Sign2(Indx, Indx2) = 1
            IF (SigDiff2(Indx, Indx2) .LE. 0.0D0) THEN
              Sign1(Indx, Indx2) = -1
            ELSE
              Sign1(Indx, Indx2) = 1
            ENDIF
            Sigma(Indx, Indx2) = Sigma(Indx, Indx2)*Sign1(Indx, Indx2)
          ENDIF
        ENDDO
      ENDDO
C
C       Making the Sigma matrix antisymetric
      DO Indx = 1, NStates
        DO Indx2 = Indx+1, NStates
          Sigma(Indx2, Indx) = -1.0D0*Sigma(Indx, Indx2)
        ENDDO
      ENDDO
C
      OPEN (11, file='sym_data.tmp')
      DO Indx = 1, NStates
        WRITE (11, 9910) (WrSigma(Indx, VnIndx), VnIndx = 1, NStates)
      ENDDO
      DO Indx = 1, NStates
        WRITE (11, 9911) (Sign1(Indx, VnIndx), VnIndx = Indx, NStates)
      ENDDO
      DO Indx = 1, NStates
        WRITE (11, 9911) (Sign2(Indx, VnIndx), VnIndx = Indx, NStates)
      ENDDO
      CLOSE (11)
C
CCC   DO Indx = 1, NStates 
CCC     DO Indx2 = Indx+1, NStates
CCC       Sigma1(Indx, Indx2) = (TDCoups1(Indx, Indx2)-
CCC  &                          TDCoups1(Indx2, Indx))/2.0D0
CCC       Sigma2(Indx, Indx2) = (TDCoups2(Indx, Indx2)-
CCC  &                          TDCoups2(Indx2, Indx))/2.0D0
CCC     ENDDO
CCC   ENDDO
C
C       Calculation of final sigma_JI(t)
CCC   DO Indx = 1, NStates
CCC     DO Indx2 = Indx+1, NStates
CCC       Sigma(Indx, Indx2) = (3.0D0*Sigma2(Indx, Indx2)-
CCC  &                                Sigma1(Indx, Indx2))/(-2.0D0)
CCC     ENDDO
CCC   ENDDO
C
C       Sigma matrix need to be multiplied by -1 factor (ABIN business)
CCC   DO Indx = 1, NStates
CCC     DO Indx2 = 1, NStates
CCC       IF (Indx .EQ. Indx2) THEN
CCC         CYCLE
CCC       ELSE
CCC         Sigma(Indx, Indx2) = -1.0D0*Sigma(Indx, Indx2)
CCC       ENDIF
CCC     ENDDO
CCC   ENDDO
C
C       Only for checking save diagonal of TIME2 couplings
      DO Indx = 1, NStates
        Sigma(Indx, Indx) = TDCoups2(Indx, Indx)
      ENDDO
C
C       Write down final results
      WRITE (*, 9902) NStates, NStates
      WRITE (*, 9901) (VnIndx, VnIndx = 1, NStates)
      DO Indx = 1, NStates
        WRITE (*, 9900) Indx, (Sigma(Indx, VnIndx), VnIndx = 1, NStates)
      ENDDO
C
C
C     *****  END OF DATA PROCESSING  *****
C
C     *****  FORMAT SECTION  *****
C
 9900 FORMAT(I4, 5000F14.8)
C
 9901 FORMAT(X, 5000I14)
C
 9902 FORMAT(2I6)
C
 9910 FORMAT(5000F14.7)
C
 9911 FORMAT(5000I3)
C
99999 CONTINUE
C
C     *****  END OF FORMAT SECTION  *****
C
C     STOP
C
      END
