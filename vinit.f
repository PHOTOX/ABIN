C------------------------------------------------------------------------
C    
C INITIAL VELOCITY DISTRIBUTION                
C =============================
C
C This subroutine returns a velocity vector VX/VY/VZ for an 
C atom of given MASS chosen randomly from a Gaussian velocity 
C distribution corresponding to a given temperature TEMP. The
C mean of this distribution is zero, the variance is given by
C
C      2     k_B*T
C sigma   = -------
C              m
C
C Original version is code f24 from the book by 
C M. P. Allen & D. J. Tildesley:
C "Computer simulation of Liquids", Appendix G
C
C TEMP      desired temperature (atomic units: E_h / k_B)
C MASS      mass of atom (atomic units: m_e)
C VX,VY,VZ  velocity vector (output)
C NATOM     number of atoms
C NOUT      FORTRAN output channel
C
C Adapted (real*8, no angular velocities, mass)   B. Schmidt, Apr 6, 1995
C
C------------------------------------------------------------------------
        SUBROUTINE vinit(TEMP, MASS, vx,vy,vz,nout,idum)
        use mod_array_size,only:npartmax
        use mod_general, only: natom,pot
        implicit none
        REAL*8,intent(out)    :: VX(npartmax), VY(npartmax),VZ(npartmax)
        real*8,intent(in)     ::  mass(npartmax)
        integer,intent(inout) ::  idum 
        integer,intent(in)    ::  nout
        integer   :: i
        REAL*8    :: TEMP,SIGMA
        REAL*8    :: GAUSSabin
        real*8    :: vcmx,vcmy,vcmz,tm
 
C Initialize the random number generator
        idum = -idum

C Loop over atoms
        do 10 i=1,natom

C Variance of distribution
           sigma = DSQRT ( TEMP/MASS(i) )
 
C Velocity vector of i-th atom
           VX(i) = sigma * GAUSSabin ( IDUM )
           VY(i) = sigma * GAUSSabin ( IDUM )
           VZ(i) = sigma * GAUSSabin ( IDUM )

 10     continue

C Momentum vector of center of mass
        vcmx = 0.d0
        vcmy = 0.d0
        vcmz = 0.d0
        tm   = 0.d0
        do 20 i=1,natom
           vcmx = vcmx + vx(i)*mass(i)
           vcmy = vcmy + vy(i)*mass(i)
           vcmz = vcmz + vz(i)*mass(i)
           tm = tm + mass(i)
 20     continue
        vcmx = vcmx / tm
        vcmy = vcmy / tm
        vcmz = vcmz / tm

C Shift velocities such that momentum of center of mass is zero
        do 30 i=1,natom
           vx(i) = vx(i) - vcmx
           vy(i) = vy(i) - vcmy
           vz(i) = vz(i) - vcmz
 30     continue

C TODO: odstraneni rotace

       if(pot.eq.'2dho')then
        vx(1)=0.0d0
        vy(1)=0.0d0
        vz(1)=0.0d0
       endif

       RETURN
       END
 
C--------------------------------------------------------------------
C                                                                
C  RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         
C                                                                
C  THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.
C                                                                
C  REFERENCE:                                                    
C                                                                
C  KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        
C     ADDISON-WESLEY), 1978                                      
C                                                                
C--------------------------------------------------------------------
        REAL*8 FUNCTION GAUSSabin ( IDUM )
 
 
        REAL*8       A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138d0, A3 = 0.252408784d0 )
        PARAMETER ( A5 = 0.076542912d0, A7 = 0.008355968d0 )
        PARAMETER ( A9 = 0.029899776d0                     )
 
        REAL*8      SUM, R, R2
        REAL*8      RAN1
        INTEGER     I,IDUM
 
C    *******************************************************************
 
        SUM = 0.0d0
 
        DO 10 I = 1, 12
 
           SUM = SUM + RAN1 ( IDUM )
 
10      CONTINUE
 
        R  = ( SUM - 6.0 ) / 4.0
        R2 = R * R
 
        GAUSSabin = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * 
     &         R2 +A1 ) * R
 
        RETURN
        END
 
 
 
