module mod_vinit
   use mod_const, only: DP
   implicit none
   private
   public :: vinit, scalevelocities
   contains
!------------------------------------------------------------------------
!    
! INITIAL VELOCITY DISTRIBUTION                
! =============================
!
! This subroutine returns a velocity vector VX/VY/VZ for an 
! atom of given MASS chosen randomly from a Gaussian velocity 
! distribution corresponding to a given temperature TEMP. The
! mean of this distribution is zero, the variance is given by
!
!      2     k_B*T
! sigma   = -------
!              m
!
! Original version is code f24 from the book by 
! M. P. Allen & D. J. Tildesley:
! "Computer simulation of Liquids", Appendix G
!
! TEMP      desired temperature (atomic units: E_h / k_B)
! MASS      mass of atom (atomic units: m_e)
! VX,VY,VZ  velocity vector (output)
! NATOM     number of atoms
! NOUT      FORTRAN output channel
!
! Adapted (real(DP), no angular velocities, mass)   B. Schmidt, Apr 6, 1995
!
!------------------------------------------------------------------------
SUBROUTINE vinit(TEMP, MASS, vx, vy, vz)
   use mod_general, only: natom, pot, nwalk
   use mod_random,  only: gautrg
   real(DP),intent(out)    :: vx(:,:), vy(:,:), vz(:,:)
   real(DP),intent(in)     :: mass(:)
   real(DP)    :: rans(3*size(mass))
   real(DP)    :: TEMP,SIGMA
   real(DP)    :: vcmx,vcmy,vcmz,tm
   integer     :: iw,iat,pom

   do iw=1,nwalk

      call gautrg(rans,natom*3,0,6)  

      pom=1
      do iat=1,natom
! Variance of distribution
         sigma = SQRT ( TEMP/MASS(iat) )
 
! Velocity vector of i-th atom
         VX(iat,iw) = sigma * rans(pom) 
         VY(iat,iw) = sigma * rans(pom+1) 
         VZ(iat,iw) = sigma * rans(pom+2)
         pom=pom+3

      end do

! Momentum vector of center of mass
      vcmx = 0.d0
      vcmy = 0.d0
      vcmz = 0.d0
      tm   = 0.d0
      do iat=1,natom
         vcmx = vcmx + vx(iat,iw)*mass(iat)
         vcmy = vcmy + vy(iat,iw)*mass(iat)
         vcmz = vcmz + vz(iat,iw)*mass(iat)
         tm = tm + mass(iat)
      end do
      vcmx = vcmx / tm
      vcmy = vcmy / tm
      vcmz = vcmz / tm

! Shift velocities such that momentum of center of mass is zero
      do iat=1,natom
         vx(iat,iw) = vx(iat,iw) - vcmx
         vy(iat,iw) = vy(iat,iw) - vcmy
         vz(iat,iw) = vz(iat,iw) - vcmz
      end do

! TODO: odstraneni rotace

      if(pot.eq.'2dho')then
         vx(1,iw)=0.0d0
         vy(1,iw)=0.0d0
         vz(1,iw)=0.0d0
      endif

   end do !nwalk
   write(*,*)'Removing center of mass velocity'


   RETURN
END subroutine vinit
 
! Scaling of velocities to correct temperature after COM velocity removal           
! or if restarting to different temperature
subroutine ScaleVelocities(vx,vy,vz)
   use mod_const,   only: autok
   use mod_general, only: natom, nwalk
   use mod_system,  only: dime, f, conatom
   use mod_nhc,     only: scaleveloc, temp
   use mod_kinetic, only: ekin_v
   use mod_shake,   only: nshake
   real(DP),intent(out)    :: vx(:,:), vy(:,:), vz(:,:)
   real(DP)  :: ekin_mom, temp_mom, scal

   ekin_mom=ekin_v(vx, vy, vz)

   if(ekin_mom.gt.0.1d-10)then
      temp_mom=2*ekin_mom/(dime*nwalk*natom-nshake-f-dime*conatom*nwalk)
   else
      temp_mom=0.0d0
   end if

   write(*,*)'Initial temperature [K]:',temp_mom*autok

! TODO: pro normal modes nemusi nutne fungovat!
   if(scaleveloc.eq.1.and.temp_mom.gt.0.1e-10)then

      write(*,*)'Scaling velocities to correct temperature.'
      scal = sqrt(temp/temp_mom)
      vx = vx * scal 
      vy = vy * scal
      vz = vz * scal

      ekin_mom=ekin_v(vx, vy, vz)
      temp_mom=2*ekin_mom/(dime*nwalk*natom-nshake-f-dime*conatom*nwalk)
      write(*,*)'Temperature after scaling [K]:',temp_mom*autok
   end if

end subroutine ScaleVelocities


!Following are legacy functions, which are not in use at this point
!--------------------------------------------------------------------
!                                                                
!  RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         
!                                                                
!  THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.
!                                                                
!  REFERENCE:                                                    
!                                                                
!  KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        
!     ADDISON-WESLEY), 1978                                      
!                                                                
!--------------------------------------------------------------------
   FUNCTION GAUSSabin ( IDUM )
      INTEGER,intent(in)  ::   IDUM
      REAL(DP) :: gaussabin
      REAL(DP),parameter :: A1 = 3.949846138d0, A3 = 0.252408784d0 
      real(DP),parameter :: A5 = 0.076542912d0, A7 = 0.008355968d0 
      real(DP),PARAMETER :: A9 = 0.029899776
      REAL(DP) :: SUM, R, R2
      REAL(DP) :: RAN1
      INTEGER     I
!    *******************************************************************
      SUM = 0.0d0
 
      DO I = 1, 12
         SUM = SUM + RAN1 ( IDUM )
      END DO
 
      R  = ( SUM - 6.0d0 ) / 4.0d0
      R2 = R * R
 
      GAUSSabin = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 ) * R
 
      RETURN
   end function gaussabin

!-------------------------------------------------------------------
!     
! (PSEUDO) RANDOM NUMBER GENERATOR     
! ================================
!                            
! RETURNS A UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0. 
! SET IDUM TO ANY NEGATIVE INTEGER VALUE TO INITIALIZE 
! OR REINITIALIZE THE SEQUENCE.
!                                                                      
! FROM "NUMERICAL RECIPES"                                            
!
! SAVE statement for IFF, R, IXn      B. Schmidt, May  31, 1992
! RAN1 returned as REAL(DP)             B. Schmidt, May  31, 1992
!
!-------------------------------------------------------------------
      REAL(DP) FUNCTION RAN1(IDUM)    

      real(DP) r(97)
      integer ix1,ix2,ix3
      integer idum,j,iff

      real(DP), parameter  :: RM1=3.8580247D-6, RM2=7.4373773D-6 
      integer, parameter :: M1=259200,IA1=7141,IC1=54773
      integer, parameter :: M2=134456,IA2=8121,IC2=28411
      integer, parameter :: M3=243000,IA3=4561,IC3=51349                    

      SAVE IFF, R, IX1,IX2,IX3    !?????????
      DATA IFF /0/                                               

      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN                           
         IFF=1                                      
         IX1=MOD(IC1-IDUM,M1)                      
         IX1=MOD(IA1*IX1+IC1,M1)                  
         IX2=MOD(IX1,M2)                         
         IX1=MOD(IA1*IX1+IC1,M1)                
         IX3=MOD(IX1,M3)                       
         DO 11 J=1,97                         
            IX1=MOD(IA1*IX1+IC1,M1)           
            IX2=MOD(IA2*IX2+IC2,M2)          
            R(J)=(dble(IX1)+dble(IX2)*RM2)*RM1
 11      CONTINUE                              
         IDUM=1                               
      ENDIF                                 
      IX1=MOD(IA1*IX1+IC1,M1)              
      IX2=MOD(IA2*IX2+IC2,M2)             
      IX3=MOD(IA3*IX3+IC3,M3)            
      J=1+(97*IX3)/M3                   
      IF(J.GT.97.OR.J.LT.1)STOP      
      RAN1=  R(J) 
      R(J)=(dble(IX1)+dble(IX2)*RM2)*RM1


      RETURN                             
      END                               

end module mod_vinit
