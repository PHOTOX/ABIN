!-File with core modules                            created by Daniel Hollas,9.2.2012

! We are using modules to initialize some variables and
! for passing global variables to different subroutines.
!------------------------------------------------------------------------------------

! mod_array_size contains various array limits.
! Modify here if you need larger arrays.
! Most of the other arrays are allocated dynamically.
module mod_array_size
   use mod_const, only: DP
   implicit none
   public
   integer, parameter :: MAXCHAIN = 10, MAXTYPES = 10
   integer, parameter :: NBINMAX = 2000, NDISTMAX = 30
   integer, parameter :: NSTMAX = 50, NTRAJMAX = 1
   save
end module mod_array_size

! General simulation parameters
module mod_general
   use mod_const, only: DP
   implicit none
   public
   ! Current time step
   integer :: it = 0
   ! Denotes integrator for equations of motion, see init.F90
   integer :: md = 1
   ! The main switch for the type of dynamics (Clasical MD, PIMD, SH...)
   integer :: ipimd = 0
   ! PIMD parameters, staging transformation, number of beads, NM transform
   integer :: istage = 0, nwalk = 1, inormalmodes = 0
   ! Ab-initio potential
   character(len=100) :: pot = '_none_'
   ! Reference potential for a multiple time step propagator
   character(len=100) :: pot_ref = '_none_'
   ! imini keyword is mostly deprecated
   integer :: imini = 0
   ! number of time steps (length of simulation)
   integer :: nstep = 1
   ! denotes number of internal steps in RESPA algorithm, see mdstep.f90
   integer :: nabin = 50, nstep_ref = 1
   ! output controls (write the propery every nwriteX step):
   ! general output
   integer :: nwrite = 1
   ! output for XYZ coordinates, velocities and forces (Molden format)
   integer :: nwritex = 1, nwritev = 0, nwritef = 0
   integer :: ncalc = 1
   ! How often do we print restart file and how often do we archive it
   integer :: nrest = 1, narchive = 10000
   ! Restart switch, 1 = restart simulation
   integer :: irest = 0
   integer :: icv = 0, anal_ext = 0, idebug = 0
   integer :: ihess
   ! Random number seed
   ! TODO: Default should be set from urandom see:
   ! https://linux.die.net/man/4/urandom
   integer :: irandom = 156873
   ! Number of atoms, taken from XYZ geometry
   integer :: natom = 0
   ! Switch for internal QM/MM, experimental!!
   integer :: iqmmm = 0
   ! Switch for Replica Exchange MD
   integer :: iremd = 0
   ! Internal MPI variables
   integer :: my_rank = 0, mpi_world_size = 1
   ! If you want to set use some exotic settings that we do not normally allow,
   ! set iknow = 1
   integer :: iknow = 0
   ! Linux Process ID, populated automatically for the current ABIN process
   integer :: pid
   ! Future variables for adaptive timestep in SH
   real(DP) :: dt0, sim_time = 0.0D0
   ! Energy restrain MD by Jiri Suchan
   integer :: en_restraint = 0
   save
end module

! Some information about simulated system, especially for distributions and shake
! TODO: Move this to a separate file, and think hard what should be inside this module.
module mod_system
   use mod_const, only: DP
   ! cannot use this
!   use mod_utils, only: abinerror
   implicit none
   public
   real(DP), allocatable :: am(:)
   character(len=2), allocatable :: names(:)
   integer, allocatable :: inames(:)
   integer :: dime = 3 ! dimension of the system
   integer :: f = 3 ! number of constants of motion
   ! (for calculating kinetic temperature)
   ! sensitive to the type of thermostat
   ! currently prob. not implemented correctly
   integer :: conatom = 0 ! number of constrained atoms
   save
contains

   ! Subroutine mass_init() populates the global am() array,
   ! based on the atom names from names() array.
   ! User can also specify non-standard isotopes/elements
   subroutine mass_init(masses, massnames, natom)
      use mod_const, only: AMU
      use mod_general, only: my_rank
      use mod_error, only: fatal_error
      real(DP), intent(in) :: masses(:)
      character(len=2), intent(in) :: massnames(:)
      integer, intent(in) :: natom
      character(len=100) :: error_msg
      integer :: i, j
      ! Accurate values for H1 and H2 taken from:
      ! Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730
      ! Other atomic weights taken from Handbook of Chemistry and Physics, 2013
      ! Original citation: Wieser, M. E., et al., Pure Appl. Chem. 85, 1047, 2013
      allocate (am(natom))
      am = -1.0D0
      do i = 1, natom
         select case (names(i))
         case ('H')
            am(i) = 1.008D0
         case ('H1')
            am(i) = 1.00782503207D0
         case ('H2', 'D')
            am(i) = 2.01410177785D0
         case ('O')
            am(i) = 15.999D0
         case ('S')
            am(i) = 32.06D0
         case ('Se')
            am(i) = 78.971D0
         case ('Te')
            am(i) = 127.60D0
         case ('N')
            am(i) = 14.007D0
         case ('P')
            am(i) = 30.973761998D0
         case ('As')
            am(i) = 74.921595D0
         case ('Sb')
            am(i) = 121.760D0
         case ('Bi')
            am(i) = 208.98040D0
         case ('F')
            am(i) = 18.998403163D0
         case ('Cl')
            am(i) = 35.45D0
         case ('Br')
            am(i) = 79.904D0
         case ('I')
            am(i) = 126.90447D0
         case ('Li')
            am(i) = 6.94D0
         case ('Na')
            am(i) = 22.98976928D0
         case ('K')
            am(i) = 39.0983D0
         case ('Be')
            am(i) = 9.0121831D0
         case ('Mg')
            am(i) = 24.305D0
         case ('Ca')
            am(i) = 40.078D0
         case ('B')
            am(i) = 10.81D0
         case ('Al')
            am(i) = 26.9815385D0
         case ('C')
            am(i) = 12.0110D0
         case ('Si')
            am(i) = 28.085D0
         case ('Ge')
            am(i) = 72.630D0
         case ('Sn')
            am(i) = 118.710D0
         case ('Pb')
            am(i) = 207.2D0
         case ('He')
            am(i) = 4.002602D0
         case ('Ne')
            am(i) = 20.1797D0
         case ('Ar')
            am(i) = 39.948D0
         case ('Kr')
            am(i) = 83.798D0
         case ('Xe')
            am(i) = 131.293D0
         case ('Fe')
            am(i) = 55.845D0
         case ('Ti')
            am(i) = 47.867D0
         case ('V')
            am(i) = 50.9415D0
         case ('Cr')
            am(i) = 51.9961D0
         case ('Mn')
            am(i) = 54.938044D0
         case ('Co')
            am(i) = 58.933194D0
         case ('Ni')
            am(i) = 58.6934D0
         case ('Cu')
            am(i) = 63.546D0
         case ('Zn')
            am(i) = 65.38D0
         case ('Ag')
            am(i) = 107.8682D0
         case ('Au')
            am(i) = 196.966569D0
         case ('Pt')
            am(i) = 195.084D0
         case ('Cd')
            am(i) = 112.414D0
         case ('Hg')
            am(i) = 200.592D0
         case ('U')
            am(i) = 238.02891D0
         case ('Tl')
            am(i) = 204.38D0
         case ('Ba')
            am(i) = 137.327D0
         case ('Ce')
            am(i) = 140.116D0
         case ('Cs')
            am(i) = 132.90545196D0
         case ('Dy')
            am(i) = 162.500D0
         case ('Er')
            am(i) = 167.259D0
         case ('Eu')
            am(i) = 151.964D0
         case ('Gd')
            am(i) = 157.25D0
         case ('Ga')
            am(i) = 69.723D0
         case ('Hf')
            am(i) = 178.49D0
         case ('Ho')
            am(i) = 164.93033D0
         case ('In')
            am(i) = 114.818D0
         case ('Ir')
            am(i) = 192.217D0
         case ('La')
            am(i) = 138.90547D0
         case ('Lu')
            am(i) = 174.9668D0
         case ('Mo')
            am(i) = 95.95D0
         case ('Nd')
            am(i) = 144.242D0
         case ('Nb')
            am(i) = 92.90637D0
         case ('Os')
            am(i) = 190.23D0
         case ('Pd')
            am(i) = 106.42D0
         case ('Pr')
            am(i) = 140.90766D0
         case ('Pa')
            am(i) = 231.03588D0
         case ('Re')
            am(i) = 186.207D0
         case ('Rh')
            am(i) = 102.90550D0
         case ('Rb')
            am(i) = 85.4678D0
         case ('Ru')
            am(i) = 101.07D0
         case ('Sm')
            am(i) = 150.36D0
         case ('Sc')
            am(i) = 44.955908D0
         case ('Sr')
            am(i) = 87.62D0
         case ('Ta')
            am(i) = 180.94788D0
         case ('Tb')
            am(i) = 158.92535D0
         case ('Th')
            am(i) = 232.0377D0
         case ('Tm')
            am(i) = 168.93422D0
         case ('W')
            am(i) = 183.84D0
         case ('Yb')
            am(i) = 173.054D0
         case ('Y')
            am(i) = 88.90584D0
         case ('Zr')
            am(i) = 91.224D0
         case DEFAULT
            if (my_rank == 0) then
               print*,'Atom name ', names(i), ' was not found in the library.'
               print*,'Using user-defined mass from input file'
            end if
         end select
      end do

      ! Check for duplicate user defined atom names
      do i = 1, size(massnames)
         do j = i + 1, size(massnames)
            if (massnames(i) == massnames(j) .and. trim(massnames(i)) /= '') then
               error_msg = 'Duplicate atom names in input array "massnames"'
               call fatal_error(__FILE__, __LINE__, error_msg)
               return
            end if
         end do
      end do

      ! Here we overwrite library values if user provided alternative value
      ! in the input file, or we define new atom names that are not in the library.
      do i = 1, natom
         do j = 1, size(massnames)
            if (names(i) == massnames(j)) then

               if (masses(j) <= 0.0D0) then
                  error_msg = 'Mass cannot be negative. Fix arrays "masses" in the input file.'
                  call fatal_error(__FILE__, __LINE__, error_msg)
                  return
               end if

               if (my_rank == 0) then
                  write (*, *) 'Defining new atom ', names(i), ' with mass=', masses(j)
               end if

               am(i) = masses(j)

            end if
         end do
      end do

      ! Catch any undefined masses
      do i = 1, natom
         if (am(i) <= 0) then
            write (error_msg, '(a,i0,a)') 'Atomic mass for atom '//names(i)//' was not specified'
            call fatal_error(__FILE__, __LINE__, error_msg)
            return
         end if
      end do

      if (natom <= 100 .and. my_rank == 0) then
         print*,''
         print*,'                        ATOMIC MASSES'
         do i = 1, natom
            write (*, '(A2, A1)', advance='no') names(i), ' '
         end do
         write (*, *)
         print*,'The corresponding relative atomic masses are:'
         write (*, *) (am(i), i=1, natom)
         print*,''
      end if

      ! Finally, convert masses to atomic units
      am = am * AMU

   end subroutine mass_init

end module mod_system

module mod_chars
   character(len=*), parameter :: CHKNOW = 'If you know what you are doing, &
    &set iknow=1 (namelist general) to proceed.'
end module mod_chars
