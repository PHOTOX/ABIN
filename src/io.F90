module mod_io
   use mod_const, only: DP
   implicit none
   character(len=50) :: formats(10)
   character(len=50) :: format1 = '(1I8)', format2 = '(1E12.4)', format3 = '(3E12.4)'
   character(len=50) :: format4 = '(1F8.4)', format5 = '(I8,I3)'
contains

   subroutine print_charges(charges, iw_ist)
      use mod_files, only: UCHARGES
      use mod_general, only: natom, it
      use mod_const, only: DP
      integer, intent(in) :: iw_ist
      real(KIND=DP), intent(in) :: charges(:)
      integer :: iat

      write (UCHARGES, format5, advance='no') it, iw_ist
      do iat = 1, natom
         write (UCHARGES, format4, advance='no') charges(iat)
      end do

      write (UCHARGES, *)
   end subroutine print_charges

   subroutine print_dipoles(dipoles, iw, nstates)
      use mod_files, only: UDIP
      use mod_general, only: it
      integer, intent(in) :: iw, nstates
      real(KIND=DP), intent(in) :: dipoles(:)
      real(KIND=DP) :: total_dip
      integer :: ind

      write (UDIP, format1, advance='no') it, iw

      do ind = 1, 3 * nstates, 3
         total_dip = dipoles(ind)**2 + dipoles(ind + 1)**2 + dipoles(ind + 2)**2
         total_dip = dsqrt(total_dip)
         write (UDIP, format2, advance="no") total_dip
      end do

      do ind = 1, 3 * nstates, 3
         write (UDIP, format3, advance="no") dipoles(ind), dipoles(ind + 1), dipoles(ind + 2)
      end do
      write (UDIP, *)

   end subroutine print_dipoles

   subroutine print_transdipoles(tdipoles, istate, ntdip)
      use mod_files, only: UTDIP
      use mod_general, only: it
      integer, intent(in) :: istate, ntdip
      real(KIND=DP), intent(in) :: tdipoles(:)
      real(KIND=DP) :: total_tdip
      integer :: ind

      write (UTDIP, format5, advance='no') it, istate

      do ind = 1, 3 * ntdip, 3
         total_tdip = tdipoles(ind)**2 + tdipoles(ind + 1)**2 + tdipoles(ind + 2)**2
         total_tdip = dsqrt(total_tdip)
         write (UTDIP, format2, advance="no") total_tdip
      end do

      do ind = 1, 3 * ntdip, 3
         write (UTDIP, format3, advance="no") tdipoles(ind), tdipoles(ind + 1), tdipoles(ind + 2)
      end do
      write (UTDIP, *)

   end subroutine print_transdipoles

   function read_forces(fx, fy, fz, num_atom, iw, funit) result(iost)
      real(DP), intent(inout) :: fx(:, :), fy(:, :), fz(:, :)
      integer, intent(in) :: iw, funit, num_atom
      integer :: iat, iost

      ! For SH or EH, when we did not calculate forces...
      ! Needs to be rewritten anyway...
      if (iw < 1) then
         iost = 0
         return
      end if

      ! WARNING: We actually expect gradients in the file !
      do iat = 1, num_atom
         read (funit, *, IOSTAT=iost) fx(iat, iw), fy(iat, iw), fz(iat, iw)
         if (iost /= 0) then
            return
         end if
         ! Convert gradients to forces
         fx(iat, iw) = -fx(iat, iw)
         fy(iat, iw) = -fy(iat, iw)
         fz(iat, iw) = -fz(iat, iw)
      end do

   end function read_forces

end module mod_io
