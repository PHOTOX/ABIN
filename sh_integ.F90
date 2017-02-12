
! Routines for handling and propagation of electronic WF during Surface Hopping.
! Currently, the only production method is the Butcher 5th order integrator.
! Euler and 4th-order Runge-Kutta methods are only here for debugging purposes.
! This module is the parent of mod_sh which contains the driver SH routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_sh_integ
   USE mod_const, ONLY: DP
   USE mod_array_size, ONLY: NSTMAX, NTRAJMAX
   implicit none
   private
   public :: sh_integrate_wf, sh_set_initialwf, write_phaserest
   public :: sh_decoherence_correction, check_popsum, sh_TFS_transmat
   public :: sh_write_wf, sh_read_wf, sh_debug_wf
   public :: integ, phase, el_pop
   public :: nstate, gama, dtp, popsumthr

   real(DP) :: cel_re(NSTMAX, NTRAJMAX), cel_im(NSTMAX, NTRAJMAX)
   real(DP) :: el_pop(NSTMAX, NTRAJMAX)
   real(DP) :: gama(NSTMAX, NSTMAX, NTRAJMAX)
   real(DP) :: eshift, dtp, popsumthr=0.001d0
   integer  :: phase = 0, nstate = 1
   character(len=10) :: integ='butcher'

   CONTAINS

   subroutine sh_integrate_wf(en_array_int, en_array_newint, dotproduct_int, dotproduct_newint, itrj)
      real(DP),intent(in) :: en_array_int(NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: en_array_newint(NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_int(NSTMAX,NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_newint(NSTMAX,NSTMAX,NTRAJMAX)
      integer, intent(in) :: itrj

      if(integ.eq.'butcher')then

         call butcherstep(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)

      else if(integ.eq.'rk4')then
        
         call rk4step(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)

      else if(integ.eq.'euler')then

         call eulerstep(en_array_int, en_array_newint, dotproduct_int, itrj)

      end if

      call sh_calc_elpop(itrj)

   end subroutine sh_integrate_wf


   subroutine sh_calc_elpop(itrj)
      integer  :: itrj, ist1
      do ist1 = 1, nstate 
         el_pop(ist1,itrj) = cel_re(ist1,itrj)**2 + cel_im(ist1,itrj)**2
      enddo
   end subroutine sh_calc_elpop


   ! Calculates transitions matrix according to the Tully's fewest switches algorithm
   subroutine sh_TFS_transmat(dotproduct_int, dotproduct_newint, itrj, ist, pop0, t)
      use mod_utils, only: abinerror
      real(DP),intent(in) :: dotproduct_int(NSTMAX,NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_newint(NSTMAX,NSTMAX,NTRAJMAX)
      real(DP),intent(out):: t(NSTMAX,NSTMAX) ! transition matrix
      integer, intent(in) :: itrj, ist ! current state
      real(DP),intent(in) :: pop0 ! population of the current state
      real(DP) :: a_re, a_im
      integer  :: ist1, ist2

      t = 0.00d0

      do ist2=1,nstate 
         a_re = ( cel_re(ist,itrj)*cel_re(ist2,itrj) + cel_im(ist,itrj)*cel_im(ist2,itrj) )
         if (phase.eq.1)then
             a_im=(- cel_im(ist,itrj)*cel_re(ist2,itrj) + cel_re(ist, itrj)*cel_im(ist2, itrj) )
             t(ist,ist2) = a_re*cos(gama(ist, ist2, itrj))-sin(gama(ist, ist2, itrj))*a_im
             t(ist,ist2) = t(ist,ist2)*(dotproduct_int(ist,ist2,itrj) + dotproduct_newint(ist,ist2,itrj) )
         else
!            t(ist,ist2)=2*a_re*dotproduct_int(ist,ist2,itrj)
             t(ist,ist2) = a_re * (dotproduct_int(ist,ist2,itrj) + dotproduct_newint(ist,ist2,itrj) )
         endif
         t(ist,ist2) = t(ist,ist2) * dtp / (pop0 + 1d-20)
      enddo


      do ist2=1,nstate 
         if(t(ist,ist2).gt.1.0_DP)then
            write(*,*)'ERROR: Hopping probability greater than 1.'
            call abinerror('surfacehop')
         end if
         if(t(ist,ist2).lt.0.0d0) t(ist,ist2) = 0.0d0
      enddo


   end subroutine sh_TFS_transmat


   subroutine integstep(k_re, k_im, en, y_re, y_im, dotprod, gam)
   real(DP),intent(out) :: k_re(NSTMAX),k_im(NSTMAX)
   real(DP),intent(in)  :: dotprod(NSTMAX,NSTMAX), gam(NSTMAX,NSTMAX)
   real(DP),intent(in)  :: en(NSTMAX),y_im(NSTMAX),y_re(NSTMAX)
   real(DP) :: g
   integer  :: ist1,ist2


   if (phase.eq.0)then
      do ist1=1,nstate
         k_re(ist1)=en(ist1)*y_im(ist1)
         k_im(ist1)=-en(ist1)*y_re(ist1)
         do ist2=1,nstate
            if(ist1.ne.ist2)then
               k_re(ist1)=k_re(ist1)-y_re(ist2)*dotprod(ist1,ist2)
               k_im(ist1)=k_im(ist1)-y_im(ist2)*dotprod(ist1,ist2)
            end if
         enddo
         k_re(ist1)=dtp*k_re(ist1)
         k_im(ist1)=dtp*k_im(ist1)
      enddo
   end if

   if (phase.eq.1)then
      do ist1=1,nstate
         k_re(ist1) = 0.0_DP
         k_im(ist1) = 0.0_DP
         do ist2 = 1, nstate
            if(ist1.ne.ist2)then
               g = gam(ist1,ist2)
               k_re(ist1)=k_re(ist1)-(y_re(ist2)*cos(g)-y_im(ist2)*sin(g))*dotprod(ist1,ist2)
               k_im(ist1)=k_im(ist1)-(y_im(ist2)*cos(g)+y_re(ist2)*sin(g))*dotprod(ist1,ist2)
            end if
         end do
         k_re(ist1)=dtp*k_re(ist1)
         k_im(ist1)=dtp*k_im(ist1)
      end do
   end if
   end subroutine integstep


   subroutine integ_gama(en_array_int, en_array_newint, itrj)
   real(DP),intent(in) :: en_array_int(:,:)
   real(DP),intent(in) :: en_array_newint(:,:)
   integer, intent(in) :: itrj
   real(DP) :: g
   integer  :: ist1, ist2

   do ist1=1,nstate
      do ist2=1,nstate
         if (ist1.ne.ist2)then
            g = (en_array_int(ist1,itrj)+en_array_newint(ist1,itrj)) / 2
            g = g - (en_array_int(ist2,itrj)+en_array_newint(ist2,itrj)) / 2
            gama(ist1,ist2,itrj) = gama(ist1,ist2,itrj) + g*dtp
         end if
      end do
   end do

   end subroutine integ_gama


   subroutine eulerstep(en_array_int, en_array_newint, dotproduct_int,itrj)
      real(DP),intent(in) :: en_array_int(:,:)
      real(DP),intent(in) :: en_array_newint(:,:)
      real(DP),intent(in) :: dotproduct_int(:,:,:)
      integer,intent(in)  :: itrj
      real(DP) :: dotprod0(NSTMAX,NSTMAX), gam0(NSTMAX, NSTMAX)
      real(DP) :: k1_re(NSTMAX),k1_im(NSTMAX)
      real(DP) :: y_im(NSTMAX),y_re(NSTMAX)
      real(DP) :: en0(NSTMAX)
      integer   :: ist1, ist2

      do ist1=1,nstate
       en0(ist1)=en_array_int(ist1,itrj)+eshift
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
         dotprod0(ist1,ist2)=dotproduct_int(ist1,ist2,itrj)
         gam0(ist1,ist2)=gama(ist1,ist2,itrj)
       enddo
      end do

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0, gam0)

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1)
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1)
      enddo

      call integ_gama(en_array_int, en_array_newint, itrj)
   end subroutine eulerstep


   subroutine rk4step(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)
      real(DP),intent(in) :: en_array_int(NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: en_array_newint(NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_int(NSTMAX,NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_newint(NSTMAX,NSTMAX,NTRAJMAX)
      integer, intent(in) :: itrj
      real(DP) :: dotprod2(NSTMAX,NSTMAX)
      real(DP) :: dotprod0(NSTMAX,NSTMAX),dotprod1(NSTMAX,NSTMAX)
      real(DP) :: k1_re(NSTMAX),k1_im(NSTMAX)
      real(DP) :: k2_re(NSTMAX),k2_im(NSTMAX)
      real(DP) :: k3_re(NSTMAX),k3_im(NSTMAX)
      real(DP) :: k4_re(NSTMAX),k4_im(NSTMAX)
      real(DP) :: y_im(NSTMAX),y_re(NSTMAX)
      real(DP) :: en0(NSTMAX),en1(NSTMAX),en2(NSTMAX)
      real(DP) :: gam0(NSTMAX,NSTMAX),gam1(NSTMAX,NSTMAX),gam2(NSTMAX,NSTMAX)
      integer :: ist1,ist2     

!     initial interpolations
      do ist1=1,nstate
       en0(ist1)=en_array_int(ist1,itrj)+eshift
       en1(ist1)=en_array_newint(ist1,itrj)+eshift
       en2(ist1)=(en0(ist1)+en1(ist1))/2
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
       dotprod0(ist1,ist2)=dotproduct_int(ist1,ist2,itrj)
       dotprod1(ist1,ist2)=dotproduct_newint(ist1,ist2,itrj)
       dotprod2(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod1(ist1,ist2))/2
       if(phase.eq.1) gam0(ist1,ist2)=gama(ist1,ist2,itrj)
       enddo
      enddo

      if (phase.eq.1)then
         call integ_gama(en_array_int, en_array_newint, itrj)
         do ist1=1,nstate
            do ist2=1,nstate
               gam1(ist1,ist2)=gama(ist1,ist2,itrj)
               gam2(ist1,ist2)=(gam0(ist1,ist2)+gam1(ist1,ist2))/2
            end do
         end do 
      end if

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0,gam0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/2
      enddo
      call integstep(k2_re,k2_im,en2,y_re,y_im,dotprod2,gam2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k2_re(ist1)/2
       y_im(ist1)=cel_im(ist1,itrj)+k2_im(ist1)/2
      enddo
      call integstep(k3_re,k3_im,en2,y_re,y_im,dotprod2,gam2)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en1,y_re,y_im,dotprod1,gam1)

      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+k1_re(ist1)/6+k2_re(ist1)/3+k3_re(ist1)/3+k4_re(ist1)/6
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+k1_im(ist1)/6+k2_im(ist1)/3+k3_im(ist1)/3+k4_im(ist1)/6
      enddo

   end subroutine rk4step


   subroutine butcherstep(en_array_int,en_array_newint,dotproduct_int,dotproduct_newint,itrj)
      real(DP),intent(in) :: en_array_int(NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: en_array_newint(NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_int(NSTMAX,NSTMAX,NTRAJMAX)
      real(DP),intent(in) :: dotproduct_newint(NSTMAX,NSTMAX,NTRAJMAX)
      real(DP) :: dotprod2(NSTMAX,NSTMAX),dotprod4(NSTMAX,NSTMAX),dotprod34(NSTMAX,NSTMAX)
      real(DP) :: dotprod0(NSTMAX,NSTMAX),dotprod1(NSTMAX,NSTMAX)
      real(DP) :: k1_re(NSTMAX),k1_im(NSTMAX)
      real(DP) :: k2_re(NSTMAX),k2_im(NSTMAX)
      real(DP) :: k3_re(NSTMAX),k3_im(NSTMAX)
      real(DP) :: k4_re(NSTMAX),k4_im(NSTMAX)
      real(DP) :: k5_re(NSTMAX),k5_im(NSTMAX)
      real(DP) :: k6_re(NSTMAX),k6_im(NSTMAX)
      real(DP) :: y_im(NSTMAX),y_re(NSTMAX)
      real(DP) :: en0(NSTMAX),en1(NSTMAX),en2(NSTMAX),en4(NSTMAX),en34(NSTMAX)
      real(DP) :: gam0(NSTMAX,NSTMAX),gam1(NSTMAX,NSTMAX),gam2(NSTMAX,NSTMAX)
      real(DP) :: gam4(NSTMAX,NSTMAX),gam34(NSTMAX,NSTMAX)
      integer :: ist1,ist2,itrj     !iteration counters

!     initial interpolations
      do ist1=1,nstate
       en0(ist1)=en_array_int(ist1,itrj)+eshift
       en1(ist1)=en_array_newint(ist1,itrj)+eshift
       en2(ist1)=(en0(ist1)+en1(ist1))/2
       en4(ist1)=(en0(ist1)+en2(ist1))/2
       en34(ist1)=(en2(ist1)+en1(ist1))/2
       y_re(ist1)=cel_re(ist1,itrj)
       y_im(ist1)=cel_im(ist1,itrj)
       do ist2=1,nstate
        dotprod0(ist1,ist2)=dotproduct_int(ist1,ist2,itrj)
        dotprod1(ist1,ist2)=dotproduct_newint(ist1,ist2,itrj)
        dotprod2(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod1(ist1,ist2))/2
        dotprod4(ist1,ist2)=(dotprod0(ist1,ist2)+dotprod2(ist1,ist2))/2
        dotprod34(ist1,ist2)=(dotprod2(ist1,ist2)+dotprod1(ist1,ist2))/2
        if (phase.eq.1) gam0(ist1,ist2)=gama(ist1,ist2,itrj)
       enddo
      enddo

      if (phase.eq.1)then
         call integ_gama(en_array_int, en_array_newint, itrj)
         do ist1=1,nstate
            do ist2=1,nstate
               gam1(ist1,ist2)=gama(ist1,ist2,itrj)
               gam2(ist1,ist2)=(gam0(ist1,ist2)+gam1(ist1,ist2))/2
               gam4(ist1,ist2)=(gam0(ist1,ist2)+gam2(ist1,ist2))/2
               gam34(ist1,ist2)=(gam2(ist1,ist2)+gam1(ist1,ist2))/2
            end do
         end do 
      end if

      call integstep(k1_re,k1_im,en0,y_re,y_im,dotprod0,gam0)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/4
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/4
      enddo
      call integstep(k2_re,k2_im,en4,y_re,y_im,dotprod4,gam4)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+k1_re(ist1)/8+k2_re(ist1)/8
       y_im(ist1)=cel_im(ist1,itrj)+k1_im(ist1)/8+k2_im(ist1)/8
      enddo
      call integstep(k3_re,k3_im,en4,y_re,y_im,dotprod4,gam4)
      
      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)-k2_re(ist1)/2+k3_re(ist1)
       y_im(ist1)=cel_im(ist1,itrj)-k2_im(ist1)/2+k3_im(ist1)
      enddo
      call integstep(k4_re,k4_im,en2,y_re,y_im,dotprod2,gam2)
      
      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)+3*k1_re(ist1)/16+9*k4_re(ist1)/16
       y_im(ist1)=cel_im(ist1,itrj)+3*k1_im(ist1)/16+9*k4_im(ist1)/16
      enddo
      call integstep(k5_re,k5_im,en34,y_re,y_im,dotprod34,gam34)

      do ist1=1,nstate
       y_re(ist1)=cel_re(ist1,itrj)-3*k1_re(ist1)/7+2*k2_re(ist1)/7+12*k3_re(ist1)/7 &
                  -12*k4_re(ist1)/7+8*k5_re(ist1)/7
       y_im(ist1)=cel_im(ist1,itrj)-3*k1_im(ist1)/7+2*k2_im(ist1)/7+12*k3_im(ist1)/7 &
                  -12*k4_im(ist1)/7+8*k5_im(ist1)/7
      enddo
      call integstep(k6_re,k6_im,en1,y_re,y_im,dotprod1,gam1)


      do ist1=1,nstate
       cel_re(ist1,itrj)=cel_re(ist1,itrj)+7*k1_re(ist1)/90+32*k3_re(ist1)/90+12*k4_re(ist1)/90 &
                        +32*k5_re(ist1)/90+7*k6_re(ist1)/90
       cel_im(ist1,itrj)=cel_im(ist1,itrj)+7*k1_im(ist1)/90+32*k3_im(ist1)/90+12*k4_im(ist1)/90 &
                        +32*k5_im(ist1)/90+7*k6_im(ist1)/90
      enddo

   end subroutine butcherstep


   subroutine sh_decoherence_correction(en_array_int, alpha, Ekin_mom, istate, itrj)
      use mod_utils, only: abinerror
      real(DP),intent(in) :: en_array_int(NSTMAX,NTRAJMAX)
      real(DP), intent(in)  :: alpha, Ekin_mom
      integer, intent(in)   :: istate, itrj
      integer  :: ist1
      real(DP) :: apom,edif,tau,fact,sum_norm

      do ist1=1,nstate
         if(ist1.ne.istate) then

!           Calculation of the exponential factor               
            edif = abs( en_array_int(ist1,itrj)-en_array_int(istate,itrj) )
            tau = 1 / edif
            apom = alpha / ekin_mom
            tau = tau * (1.0d0+apom)
            fact = dexp(-dtp/tau)
   
            cel_re(ist1,itrj) = cel_re(ist1,itrj) * fact
            cel_im(ist1,itrj) = cel_im(ist1,itrj) * fact
         endif
      enddo

!     RENORMALIZATION OF THE CURRENT STATE     
      sum_norm=1.0d0
      do ist1=1,nstate
         if(ist1.ne.istate) then
            sum_norm=sum_norm-cel_re(ist1,itrj)**2-cel_im(ist1,itrj)**2
         endif
      enddo

      fact = sum_norm / (cel_re(istate,itrj)**2 + cel_im(istate,itrj)**2+1.0d-7)

      !Following should never happen as we check for popsumthr later in this subroutine
      if(fact.lt.0.0d0)then
         write(*,*)'Fatal error in surfacehop during decoherence renormalization.'
         write(*,*)'fact=',fact,'but should be > 0'
         write(*,*)'This usually means inaccurate integration of electronic SE.'
         write(*,*)'Increase number of substeps or use more accurate integrator.'
         call abinerror('surfacehop')
      end if

      fact = sqrt(fact)

      cel_re(istate,itrj) = cel_re(istate,itrj) * fact
      cel_im(istate,itrj) = cel_im(istate,itrj) * fact

      call sh_calc_elpop(itrj)

   end subroutine sh_decoherence_correction


   subroutine write_phaserest(iunit, itrj)
   use mod_general,  only: narchive, it
   use mod_utils,    only: archive_file
   integer, intent(in) :: iunit, itrj
   integer :: ist1, ist2
   ! iunit2 here only temporarily, should write to iunit
   integer :: iunit2
   iunit2 = 601

   open(iunit2,file='phaserest.dat',action="write")

   do ist1=1,nstate
      write(iunit2,*)(gama(ist1,ist2,itrj),ist2=1,nstate)
   end do

   close(iunit2)

   if(modulo(it,narchive).eq.0)then
      call archive_file('phaserest.dat',it)
   end if

   end subroutine write_phaserest


   subroutine sh_set_initialwf(initial_state, initial_poten, itrj)
      use mod_general, only: irest
      real(DP), intent(in) :: initial_poten
      integer, intent(in)  :: initial_state, itrj
      integer  :: ist1

      gama = 0.0d0

      if (irest.eq.0)then
         do ist1=1,nstate
            cel_re(ist1,itrj) = 0.0d0
            cel_im(ist1,itrj) = 0.0d0
         enddo

         cel_re(initial_state,itrj) = 1.0d0
      end if

      ! Eshift hard set as the potential energy of the ground state
      ! Probably should be handled better (and differently for different trajs)
      Eshift = -initial_poten

      call sh_calc_elpop(itrj)

   end subroutine sh_set_initialwf


   subroutine sh_write_wf(outunit,itrj)
      integer, intent(in)  :: outunit, itrj
      integer  :: ist1
      do ist1=1,nstate
         write(outunit,*)cel_re(ist1,itrj),cel_im(ist1,itrj)
      enddo
   end subroutine sh_write_wf


   subroutine sh_read_wf(inunit,itrj)
      integer, intent(in)  :: inunit, itrj
      integer  :: ist1
      do ist1=1,nstate
         read(inunit,*)cel_re(ist1,itrj),cel_im(ist1,itrj)
      enddo
   end subroutine sh_read_wf


   real(DP) function check_popsum(itrj)
      use mod_utils, only: abinerror
      integer,intent(in)  :: itrj
      real(DP)            :: popsum
      integer             :: ist1

      popsum=0.0d0
      do ist1 = 1, nstate 
         popsum = popsum + el_pop(ist1,itrj)
      end do

      if (abs(popsum-1.0d0).gt.popsumthr)then
         write(*,*)'ERROR:Sum of electronic populations = ',popsum
         write(*,*)'which differs from 1.0 by more than popsumthr = ',popsumthr
         write(*,*)'Increase the number of substeps or use more accurate integrator.'
         call abinerror('surfacehop')
      end if

      check_popsum = popsum
   end function check_popsum


   ! This prints detailed debug info about WF during integration
   subroutine sh_debug_wf(ist, itrj, stepfs, t)
      use mod_files, only: UBKL, UWFCOEF, UPHASE 
      integer, intent(in)  :: ist  ! current electronic state
      integer, intent(in)  :: itrj 
      real(DP), intent(in) :: stepfs
      real(DP), intent(in) :: t(NSTMAX,NSTMAX)
      integer  :: ist1, ist2
      character(len=500)  :: formt

      write(formt,'(A7,I3,A7)')'(F15.4,',nstate,'E20.10)'
      write(UBKL,fmt=formt)stepfs,(t(ist,ist1),ist1=1,nstate)

      write(UWFCOEF,fmt=formt,advance="no")stepfs,(cel_re(ist1,itrj),ist1=1,nstate)
      write(formt,'(A1,I3,A7)')'(',nstate,'E20.10)'
      write(UWFCOEF,fmt=formt,advance="no")(cel_im(ist1,itrj),ist1=1,nstate)
      write(UWFCOEF,*)''

      if(phase.eq.1)then
         write(UPHASE,'(F15.2,E20.10)',advance="no")stepfs,gama(2,1,itrj)
         do ist1=3,nstate
            write(formt,'(A1,I3,A7)')'(',ist1-1,'E20.10)'
            write(UPHASE,fmt=formt,advance="no")(gama(ist1,ist2,itrj),ist2=1,ist1-1)
         end do
         write(UPHASE,*)''
      end if

   end subroutine sh_debug_wf

end module mod_sh_integ

