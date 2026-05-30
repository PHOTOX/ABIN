! This file must be recompiled at each compilation
! to get the current date and time,
! which are taken from preprocessor constants,
! and the current Git commit, which is passed in from Makefile.
! allow(procedure-not-in-module) ! fortitude linter
subroutine print_compile_info()
   use iso_fortran_env, only: compiler_version, compiler_options
   use mod_files, only: stdout
   character(len=*), parameter :: ABIN_VERSION = '1.2.1-dev'

   write (stdout, *) ''
   write (stdout, *) '          COMPILATION INFO'
   write (stdout, *) ''
   write (stdout, '(1x, a, t25, a)') 'ABIN version:', ABIN_VERSION
   write (stdout, '(1x, a, t25, a, 1x, a)') 'Compiled at:', __TIME__, __DATE__
   write (stdout, '(1x, a, t25, a)') 'Git commit:', GIT_COMMIT
   write (stdout, '(1x, a, t25, a)') 'Compiler:', compiler_version()
   write (stdout, *) 'Compiler options:'
   write (stdout, *) compiler_options()
   write (stdout, *)
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
end subroutine print_compile_info
