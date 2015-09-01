      PROGRAM fGAMConvStr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Program that converts the xyz file or movie in xyz format to GAMESS  C
C xyz format.                                                          C
C This is a part of the bash script GAMConvStr                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Authors:                                 Lukas Sistik               C
C  Date  PV:                                12.04.2012                 C

C DH TODO: make this work without the bash script
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
        PARAMETER (MaxAtm = 1000)
C
      INTEGER MaxStr
        PARAMETER (MaxStr = 1000)
C
C         DECLARATION - FUNCTION
C
C
C         DECLARATION - VARIABLES
C
      INTEGER Indx, VnIndx, Indx2
C
      INTEGER NAtoms, NStruct
C
      DOUBLE PRECISION CoordX(MaxStr, MaxAtm), CoordY(MaxStr, MaxAtm), 
     &                 CoordZ(MaxStr, MaxAtm), AtmNum(MaxAtm)
C
      CHARACTER*2 AtmMark(MaxAtm)
C
      CHARACTER junk
C
C-----------------------------------------------------------------------
C
C         MAIN 
C
      OPEN(11, file='temp.tmp', STATUS='OLD', ACCESS='SEQUENTIAL')
C         
C         Reading input parameters: 
C            Number of atoms in molecule, number of structures
      READ (11, *) NAtoms, NStruct
C
C         Reading input data:
C            coordinates of atoms
      DO Indx = 1, NStruct
        READ (*, *)
        READ (*, *)
        DO Indx2 = 1, NAtoms
          READ (*, *) AtmMark(Indx2), CoordX(Indx, Indx2), 
     &                CoordY(Indx, Indx2), CoordZ(Indx, Indx2)
        ENDDO
      ENDDO
C
      DO Indx = 1, NAtoms
        SELECT CASE (AtmMark(Indx))
          CASE ("H", "h")
            AtmNum(Indx) = 1.0D0
          CASE ("He", "he", "HE")
            AtmNum(Indx) = 2.0D0
          CASE ("Li", "li", "LI")
            AtmNum(Indx) = 3.0D0
          CASE ("Be", "be", "BE")
            AtmNum(Indx) = 4.0D0
          CASE ("B", "b")
            AtmNum(Indx) = 5.0D0
          CASE ("C", "c")
            AtmNum(Indx) = 6.0D0
          CASE ("N", "n")
            AtmNum(Indx) = 7.0D0
          CASE ("O", "o")
            AtmNum(Indx) = 8.0D0
          CASE ("F", "f")
            AtmNum(Indx) = 9.0D0
          CASE ("Ne", "ne", "NE")
            AtmNum(Indx) = 10.0D0
          CASE ("Na", "na", "NA")
            AtmNum(Indx) = 11.0D0
          CASE ("Mg", "mg", "MG")
            AtmNum(Indx) = 12.0D0
          CASE ("Al", "al", "AL")
            AtmNum(Indx) = 13.0D0
          CASE ("Si", "si", "SI")
            AtmNum(Indx) = 14.0D0
          CASE ("P", "p")
            AtmNum(Indx) = 15.0D0
          CASE ("S", "s")
            AtmNum(Indx) = 16.0D0
          CASE ("Cl", "cl", "CL")
            AtmNum(Indx) = 17.0D0
          CASE ("Ar", "ar", "AR")
            AtmNum(Indx) = 18.0D0
        END SELECT
      ENDDO
C
      DO Indx = 1, NStruct
        WRITE (*, 9901) NAtoms
        WRITE (*, *)
        DO Indx2 = 1, NAtoms
          WRITE (*, 9900) AtmMark(Indx2), AtmNum(Indx2), 
     &                    CoordX(Indx, Indx2), CoordY(Indx, Indx2), 
     &                    CoordZ(Indx, Indx2)
        ENDDO
      ENDDO
C
 9900 FORMAT(A3, F5.1, 3F14.6)
C
 9901 FORMAT(I6)
C
 9902 FORMAT(A2,F15.6, 2F14.6)
C
 9999 FORMAT(2I6, A3, F14.6)
C
99999 CONTINUE
C
C     STOP
C
      END
