      PROGRAM fNumForce
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C This program calculates numerical force from data presented in the   C
C "force_en.temp" file.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                06.06.2012                 C
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
        PARAMETER (MaxSt = 50)
C
      INTEGER MaxCoord
        PARAMETER (MaxCoord = 500)
C
C         DECLARATION - FUNCTION
C
C
C         DECLARATION - VARIABLES
C
      INTEGER Indx, VnIndx, Indx2
C
      INTEGER NDiff, NState
C
      DOUBLE PRECISION Step
      DATA Step /0.010D0/
C
      DOUBLE PRECISION En1(MaxSt), En2(MaxSt)
C
      DOUBLE PRECISION ForceX(MaxCoord, MaxSt), ForceY(MaxCoord, MaxSt),
     &                 ForceZ(MaxCoord, MaxSt)
C
C         DECLARATION - CONSTANTS
C
      DOUBLE PRECISION Bohr2Ang
      DATA Bohr2Ang /1.889725989D0/
C
C-----------------------------------------------------------------------
C
C         MAIN 
C
      OPEN(11, file='force_en.temp', ACCESS='SEQUENTIAL',STATUS='OLD')
C
      READ (11, *) NDiff 
      READ (11, *) NState
C
      DO Indx = 1, NDiff/3
        READ (11, *) (En1(VnIndx), VnIndx = 1, NState)
        READ (11, *) (En2(VnIndx), VnIndx = 1, NState)
        READ (11, *)
        DO Indx2 = 1, NState
          ForceX(Indx, Indx2) = (En2(Indx2)-En1(Indx2))/Step/2.0D0
        ENDDO
C
        READ (11, *) (En1(VnIndx), VnIndx = 1, NState)
        READ (11, *) (En2(VnIndx), VnIndx = 1, NState)
        READ (11, *)
        DO Indx2 = 1, NState
          ForceY(Indx, Indx2) = (En2(Indx2)-En1(Indx2))/Step/2.0D0
        ENDDO
C
        READ (11, *) (En1(VnIndx), VnIndx = 1, NState)
        READ (11, *) (En2(VnIndx), VnIndx = 1, NState)
        READ (11, *)
        DO Indx2 = 1, NState
          ForceZ(Indx, Indx2) = (En2(Indx2)-En1(Indx2))/Step/2.0D0
        ENDDO
C
      ENDDO
C
      WRITE (*, *) NState
      DO Indx = 1, NState
        WRITE (*, *) "           X             Y             Z"
        DO Indx2 = 1, NDiff/3
CCCCC     WRITE (*, 9900) Indx2, -1.0d0*ForceX(Indx2, Indx)/Bohr2Ang, 
CCCCC&                           -1.0d0*ForceY(Indx2, Indx)/Bohr2Ang, 
CCCCC&                           -1.0d0*ForceZ(Indx2, Indx)/Bohr2Ang
          WRITE (*, 9900) Indx2, ForceX(Indx2, Indx)/Bohr2Ang, 
     &                           ForceY(Indx2, Indx)/Bohr2Ang, 
     &                           ForceZ(Indx2, Indx)/Bohr2Ang
        ENDDO
        WRITE (*, *)
      ENDDO
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
C     *****  END OF FORMAT SECTION  *****
C
C     STOP
C
      END
