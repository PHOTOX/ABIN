## ABIN utilities

This file contains a brief description of what utilities and scripts that you can find here.
More complete description is usually given in each script/source code.

1.  create_trajectories.sh, analyzeSH.sh
Scripts capable of launching and analyzing multiple ABIN trajectories.

2. abin-randomint
Fortran code for creating random integers. This one should be in your $PATH.

3. analyze_movie
Fortran code for analysis of xyz movies from abin.
It can analyze bonds, angles and dihedrals and creates histograms.

4. r.abin
This is a simple bash script for launching ABIN on clusters.

5. Exitabin.sh, fetchabin.sh
These should be in your $PATH. 
These scripts expects that ABIN was launched via r.abin
fetchabin.sh copies data from scratch to the current working directory.
Exitabin.sh exits nicely running ABIN job (better than qdel).

6. prum.awk, checkenergy.awk
These are small helper scripts used by other scripts.
They should be in your $PATH.

