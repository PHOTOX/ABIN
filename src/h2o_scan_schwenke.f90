! To compile:
! gfortran -c h2o_scan_schwenke.f90
! gfortran -fopenmp -L./ -L../water_potentials/ libabin.a ../water_potentials/libwater.a h2o_scan_schwenke.o -labin -lwater -lm -lstdc++ -o schwenke_scan
program schwenke_scan
    use mod_force_h2o
    use mod_const, only: DP, ANG
    use mod_init, only: read_xyz_file
    implicit none

    integer, parameter :: NATOM = 3
    integer, parameter :: NBEADS = 1
    integer, parameter :: WATPOT = 1
    
    ! Displacement for H atom at each timestep
    real(DP) :: delta = 0.0005_DP
    integer, parameter :: NSTEPS = 5000

    real(DP) :: x(3, 1), y(3, 1), z(3, 1)
    real(DP) :: fx(3, 1), fy(3, 1), fz(3, 1)
    real(DP) :: energy, rOH
    integer :: i, funit
    character(len=2) :: atom_names(3)
    character(len=256) :: fname

    atom_names = (/ "O", "H", "H" /)

    fname = "water.xyz"
    open (newunit=funit, file=fname, action="read")
    call read_xyz_file(funit, fname, atom_names, natom, nbeads, x, y, z)
    close (funit)

    ! Convert from Angstromgs to Bohr (atomic units)
    x(:, 1) = x(:, 1) * ANG
    y(:, 1) = y(:, 1) * ANG
    z(:, 1) = z(:, 1) * ANG

    ! print for output
    print*, "Step   ", "rOH   ", "Energy"

    open(99, file="water_geometries.xyz", action="write")
    ! steps for stretching OH bond
    do i = 1, NSTEPS
        ! Calculate the initial OH bond length (rOH) from initial bond positions
        rOH = sqrt((x(2, 1) - x(1, 1))**2 + (y(2, 1) - y(1, 1))**2 + (z(2, 1) - z(1, 1))**2)

        ! Modify the OH bond length by changing the position of a single hydrogen atom
        x(2, 1) = x(2, 1) + delta

        ! Call the force_water function to calculate forces and energy
        energy = 0.0D0
        call force_h2o(x, y, z, fx, fy, fz, energy, natom, nbeads)

        ! Print the current OH bond length and energy
        print '(I4,F16.8,E18.8)', i, rOH, energy

        ! Print XYZ geometry
        write(99, '(I1)') NATOM
        write(99, *)
        write(99, '(A,F12.6,F12.6,F12.6)') "O ", x(1,1) / ANG, y(1,1) / ANG, z(1,1) / ANG
        write(99, '(A,F12.6,F12.6,F12.6)') "H ", x(2,1) / ANG, y(2,1) / ANG, z(2,1) / ANG
        write(99, '(A,F12.6,F12.6,F12.6)') "H ", x(3,1) / ANG, y(3,1) / ANG, z(3,1) / ANG
    end do

    ! close file
    close(99)
end program