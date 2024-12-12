module mod_water
   implicit none
   private
   public :: watpot, check_water
   integer :: watpot = 1
   save
contains
   subroutine check_water(natom, names)
      use mod_files, only: stdout
      use mod_error, only: fatal_error
      use mod_utils, only: count_atoms_by_name
      integer, intent(in) :: natom
      character(len=2), intent(in) :: names(:)
      integer :: nH, nO
      integer :: error, iat

      write (stdout, *) ''
      select case (watpot)
      case (1)
         write (stdout, *) 'Using water force field q-TIP4P/F : https://dx.doi.org/10.1063/1.3167790'
      case (2)
         write (stdout, *) 'Using water force field TTM2.1-F  : https://dx.doi.org/10.1021/jp056477k'
      case (3)
         write (stdout, *) 'Using water force field TTM3-F : https://dx.doi.org/10.1063/1.2837299'
      case (4)
         write (stdout, *) 'Using water force field TTM4-F : https://dx.doi.org/10.1063/1.2895750'
      case default
         call fatal_error(__FILE__, __LINE__, 'Parameter "watpot" must be 1, 2, 3 or 4')
      end select

      error = 0

      write (stdout, *) 'Checking that input atom types correspond to pure water.'

      ! Note that some of these checks are redundant,
      ! but we're trying to be helpful and provide as much info as we can.

      nO = count_atoms_by_name(names, 'O', natom)
      nH = count_atoms_by_name(names, 'H', natom)
      if (nO * 2 /= nH) then
         write (stdout, *) 'ERROR: oxygen to hydrogen ratio is not 1:2'
         error = 1
      end if

      if (modulo(natom, 3) /= 0) then
         write (stdout, *) 'ERROR: Number of atoms is not divisible by 3.'
         error = 1
      end if

      do iat = 1, natom, 3
         if (names(iat) /= 'O' .or. names(iat + 1) /= 'H' .or. names(iat + 2) /= 'H') then
            error = 1
            write (stdout, *) 'ERROR: Bad element type.'
            write (stdout, *) 'Water atoms must be ordered as "O H H" '
            exit
         end if
      end do

      if (error == 1) then
         call fatal_error(__FILE__, __LINE__, 'This is not pure water.')
      end if

      write (stdout, '(A,I0,A)') 'Detected ', nO, ' water molecules'
   end subroutine check_water

end module mod_water
