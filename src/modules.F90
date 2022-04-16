! File with global simulation parameters

! We are using modules to initialize some variables and
! for passing global variables to different subroutines.

! mod_array_size contains various array limits.
! Modify here if you need larger arrays.
! Most of the other arrays are allocated dynamically.
module mod_array_size
   use mod_const, only: DP
   implicit none
   public
   integer, parameter :: NSTMAX = 50
   save
end module mod_array_size

! General simulation parameters
module mod_general
   use mod_const, only: DP
   implicit none
   public
   ! Current time step
   integer :: it = 0
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
   ! If you want to set use some exotic settings that we do not normally allow,
   ! set iknow = 1
   integer :: iknow = 0
   ! Initial time step (for future adaptime timestep functionality)
   real(DP) :: dt0
   ! Total simulation time
   real(DP), protected :: sim_time = 0.0D0
   ! Energy restrain MD by Jiri Suchan
   integer :: en_restraint = 0
   save
contains
   subroutine update_simtime(dt)
      use mod_const, only: DP
      real(DP) :: dt
      sim_time = sim_time + dt
   end subroutine
end module

! TODO: Move this to a separate file, and think hard what should be inside this module.
module mod_system
   use mod_const, only: DP
   implicit none
   public
   real(DP), allocatable :: am(:)
   character(len=2), allocatable :: names(:)
   integer :: dime = 3 ! dimension of the system
   integer :: f = 3 ! number of constants of motion
   ! (for calculating kinetic temperature)
   ! sensitive to the type of thermostat
   ! currently prob. not implemented correctly
   integer :: conatom = 0 ! number of constrained atoms
   save
end module mod_system

module mod_chars
   character(len=*), parameter :: CHKNOW = 'If you know what you are doing, &
    &set iknow=1 (namelist general) to proceed.'
end module mod_chars
