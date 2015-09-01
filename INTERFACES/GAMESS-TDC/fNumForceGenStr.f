      PROGRAM fNumForce
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C This program generates structures for numerical force evaluation.    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                05.06.2012                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         DECLARATION
C
      IMPLICIT NONE
C
C         DECLARATION - PARAMETERS
C
      INTEGER MaxAtm
        PARAMETER (MaxAtm = 500)
C
C         DECLARATION - FUNCTION
C
C
C         DECLARATION - VARIABLES
C
      INTEGER Indx, VnIndx, Indx2
C
      INTEGER NAtoms
C
      DOUBLE PRECISION CoordX(MaxAtm), CoordY(MaxAtm), CoordZ(MaxAtm)
C
      DOUBLE PRECISION Step
      DATA Step /0.010D0/
C
      CHARACTER*2 AtmMark(MaxAtm)
C
C-----------------------------------------------------------------------
C
C         MAIN 
C     
      READ (*, *) NAtoms
      READ (*, *)
C
      DO Indx = 1, NAtoms
        READ (*, *) AtmMark(Indx), CoordX(Indx), CoordY(Indx), 
     &              CoordZ(Indx)
      ENDDO
C
      DO Indx = 1, NAtoms
C         X-1
        WRITE (*, *) NAtoms
        WRITE (*, *) Indx, "X - Step"
        DO Indx2 = 1, NAtoms
          IF (Indx .EQ. Indx2) THEN
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2)-Step, 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ELSE
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ENDIF
        ENDDO
C         X+1
        WRITE (*, *) NAtoms
        WRITE (*, *) Indx, "X + Step"
        DO Indx2 = 1, NAtoms
          IF (Indx .EQ. Indx2) THEN
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2)+Step, 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ELSE
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ENDIF
        ENDDO
C         Y-1
        WRITE (*, *) NAtoms
        WRITE (*, *) Indx, "Y - Step"
        DO Indx2 = 1, NAtoms
          IF (Indx .EQ. Indx2) THEN
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2)-Step, CoordZ(Indx2)
          ELSE
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ENDIF
        ENDDO
C         Y+1
        WRITE (*, *) NAtoms
        WRITE (*, *) Indx, "Y + Step"
        DO Indx2 = 1, NAtoms
          IF (Indx .EQ. Indx2) THEN
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2)+Step, CoordZ(Indx2)
          ELSE
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ENDIF
        ENDDO
C         Z-1
        WRITE (*, *) NAtoms
        WRITE (*, *) Indx, "Z - Step"
        DO Indx2 = 1, NAtoms
          IF (Indx .EQ. Indx2) THEN
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)-Step
          ELSE
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ENDIF
        ENDDO
C         Z+1
        WRITE (*, *) NAtoms
        WRITE (*, *) Indx, "Z + Step"
        DO Indx2 = 1, NAtoms
          IF (Indx .EQ. Indx2) THEN
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)+Step
          ELSE
            WRITE (*, 9900) AtmMark(Indx2), CoordX(Indx2), 
     &                      CoordY(Indx2), CoordZ(Indx2)
          ENDIF
        ENDDO
      ENDDO
C
C     *****  FORMAT SECTION  *****
C
 9900 FORMAT(A4, 5000F14.8)
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
