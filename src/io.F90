module mod_io
   use mod_const, only: DP
   implicit none
   public
   ! TODO: Move these to individial functions
   character(len=50) :: format1 = '(1I8)', format2 = '(1E12.4)', format3 = '(3E12.4)'
   character(len=50) :: format4 = '(1F8.4)', format5 = '(I8,I3)'
contains

   ! TODO: Create a ElectronicStructure type
   ! and add reading and printint routines as type-bound procedures.
   ! The instance of this type should be an array
   ! corresponding to the number of beads.
   subroutine print_charges(charges, iw_ist)
      use mod_files, only: UCHARGES
      use mod_general, only: natom, it
      use mod_const, only: DP
      ! Index iw_ist is indexing either PI beads
      ! or electronic state in Surface Hopping.
      integer, intent(in) :: iw_ist
      real(KIND=DP), intent(in) :: charges(:)
      integer :: iat

      ! TODO: Print time in fs, not time step
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

      ! TODO: Print time in fs, not time step
      write (UDIP, format5, advance='no') it, iw

      ! TODO: Why the heck are we printing total dipmom fist?
      ! Most QM program print Dx, Dy, Dz components first.
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

end module mod_io
