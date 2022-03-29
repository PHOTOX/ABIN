module mod_water
   implicit none
   private
   public :: watpot, check_water
   integer :: watpot = 1
   save
contains
   subroutine check_water(natom, names)
      use mod_general, only: my_rank
      use mod_utils, only: abinerror, count_atoms_by_name
      integer, intent(in) :: natom
      character(len=2) :: names(:)
      integer :: nH, nO
      integer :: error, iat

      if (watpot < 1 .or. watpot > 4) then
         write (*, *) 'ERROR: parameter "watpot" must be 1, 2, 3 or 4.'
         write (*, *) '1 - q-TIP4P/F : https://dx.doi.org/10.1063/1.3167790'
         write (*, *) '2 - TTM2.1-F  : https://dx.doi.org/10.1021/jp056477k'
         write (*, *) '3 - TTM3-F    : https://dx.doi.org/10.1063/1.2837299'
         write (*, *) '4 - TTM4-F    : https://dx.doi.org/10.1063/1.2895750'
         call abinerror('check_water')
      end if

      error = 0

      if (my_rank == 0) print*,'Checking that input atom types correspond to pure water.'

      ! Note that some of these checks are redundant,
      ! but we're trying to be helpful and provide as much info as we can.

      nO = count_atoms_by_name(names, 'O', natom)
      nH = count_atoms_by_name(names, 'H', natom)
      if (nO * 2 /= nH) then
         write (*, *) 'ERROR: oxygen to hydrogen ratio is not 1:2'
         error = 1
      end if

      if (modulo(natom, 3) /= 0) then
         write (*, *) 'ERROR: Number of atoms is not divisible by 3.'
         error = 1
      end if

      do iat = 1, natom, 3
         if (names(iat) /= 'O' .or. names(iat + 1) /= 'H' .or. names(iat + 2) /= 'H') then
            error = 1
            write (*, *) 'ERROR: Bad element type.'
            write (*, *) 'Water atoms must be ordered as "O H H" '
         end if
      end do

      if (error == 1) then
         write (*, *) 'This is not pure water.'
         call abinerror('check_water')
      end if

      if (my_rank == 0) print '(A,I0,A)', 'Detected ', nO, ' water molecules'
   end subroutine check_water

end module mod_water
