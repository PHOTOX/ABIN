-------------------------------
STRUCTURE OF THE SOURCE CODE
-------------------------------

Short descriptions of key source files.

| Path     | Description |
|----------|-------------|
| abin.F90        | Main routine. |
| modules.F90     | Modules containing all basic variables. |
| init.F90        | Big ugly init routine. Reads input parameters, restart files and checks for errors in input. |
| mdstep.f90      | Routines for propagation of equations of motion (velocity Verlet and RESPA). |
| forces.f90      | Main routine for getting forces (force_clas) and quantum forces in PIMD (force_quant). |
| force_abin.f90  | Routine that calls ab initio BASH interfaces. |
| nosehoover.f90  | Nosé-Hoover thermostat. |
| surfacehop.f90  | Surface Hopping dynamics |
| sh_integ.f90    | Propagation of electronic wavefunction in SH dynamics | 
| landau_zener.f90| Landau-Zener dynamics |
| gle.F90         | Generalized Langevin Equation thermostat |
| shake.f90       | Constraints using SHAKE algorithm. |
| analysis.f90    | Driver routine for printing output. |
| arrays.f90      | Module containing all dynamically allocated arrays.| 
| random.f90      | Random number generator. |
| vinit.f90       | Initial velocities. |
| potentials.f90  | Analytical potentials and hessians for harmonic and Morse oscillators and others. |
| density.f90     | Evaluates distance, angle and dihedrals distributions during dynamics. | 
| estimators.f90  | Energy and heat capacity estimators for PIMD. | 
| transform.F90   | Coordinate transformations for PIMD. |
