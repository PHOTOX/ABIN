! This file must be recompiled at each compilation
! to get the current date and time,
! which are taken from preprocessor constants,
! and the current Git commit, which is passed in from Makefile.
subroutine print_compile_info()
   use iso_fortran_env, only: compiler_version
   use iso_fortran_env, only: compiler_version, compiler_options
   use mod_files, only: stdout
   ! TODO: Decide how we should do versioning
   character(len=*), parameter :: ABIN_VERSION = '1.1'

   write (stdout, *) 'ABIN version '//ABIN_VERSION
   write (stdout, '(a, a, 1x, a)') 'Compiled at ', __TIME__, __DATE__
   write (stdout, *) 'Git commit '//GIT_COMMIT
!$ write (stdout, *) 'Compiled with parallel OpenMP support for PIMD.'
#ifdef USE_FFTW
   write (stdout, *) 'Compiled with FFTW support.'
#endif
#ifdef USE_CP2K
   write (stdout, *) 'Compiled with in-built CP2K interface.'
#endif
#ifdef USE_PLUMED
   write (stdout, *) 'Compiled with PLUMED (static lib).'
#endif
#ifdef USE_MPI
   write (stdout, *) 'Compiled with MPI support.'
   write (stdout, *) '(used for REMD and direct CP2K and TeraChem interfaces.)'
#endif
   write (stdout, *)

   write (stdout, *) 'This program was compiled by '//compiler_version()
   write (stdout, *) 'using the following compiler options:'
   write (stdout, *) compiler_options()
end subroutine print_compile_info
