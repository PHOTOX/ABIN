! Exact values for ANG,AUTOKCAL and AUTODEBYE,me,AUTOM,AUTOFS,AMU:
! Mohr, Taylor, Newell, Rev. Mod. Phys. 80 (2008) 633-730
! and using the thermochemical calorie (1 cal = 4.184 J):'

! TODO: Update physical constant according to the latest data
! https://physics.nist.gov/cuu/Constants/bibliography.html

! ABIN uses atomic units internally
module mod_const
   implicit none
   ! TODO; Separate DP to its own module, mod_kinds
   ! and rename mod_const to mod_phys_const
   integer, parameter :: DP = kind(1.0D0)
   ! Exact values of Planck and Avogadro constants according
   ! to the SI redefinition of kilogram in 2019.
   ! TODO: Define all other constants based on base SI units.
   real(DP), parameter :: AVOGADRO = 6.02214076D23
   real(DP), parameter :: PLANCK = 6.626070150D-23
   ! Hartree to Joul converstion according to
   ! https://physics.nist.gov/cuu/Constants/Table/allascii.txt
   real(DP), parameter :: AUTOJ = 4.3597447222071D-18 ! hartree -> J
   real(DP), parameter :: AMU = 1822.888484264545D0 ! atomic mass unit
   real(DP), parameter :: ANG = 1.889726132873D0 ! Angstroms -> Bohrs
   real(DP), parameter :: AUTOFS = 0.02418884326505D0 !atomic units to femtosecs
   real(DP), parameter :: PI = 3.14159265358979323846D0
   real(DP), parameter :: AUTOK = 3.1577464D5 ! temperature in a.u. to Kelvins
   real(DP), parameter :: ME = 9.10938215D-31 ! electron mass
   real(DP), parameter :: AUTOM = 5.2917720859D-11 ! atomic length
   real(DP), parameter :: AUTOCM = 2.1947463D5 ! atomic units to cm-1
   real(DP), parameter :: AUTOKCAL = 6.2750946943D2 ! Ha -> kcal/mol
   real(DP), parameter :: CALTOJ = 4.184D0
   real(DP), parameter :: AUTOKK = 3.1577322D5
   real(DP), parameter :: AUTOEV = 27.21138386D0
   real(DP), parameter :: AUTODEBYE = 2.54174623D0
   ! Boltzmann constant in a.u., assumed to be 1.0d0 in the program
   real(DP), parameter :: KBinAU = 0.9999999748284666D0
   save
end module mod_const
