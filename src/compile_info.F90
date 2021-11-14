! This file must be recompiled at each compilation
! to get the current date and time,
! which are taken from preprocessor constants,
! and the current Git commit, which is passed in from Makefile.
subroutine print_compile_info()
   use iso_fortran_env, only: compiler_version
   use iso_fortran_env, only: compiler_version, compiler_options
   ! TODO: Decide how we should do versioning
   character(len=*), parameter :: ABIN_VERSION = '1.1'

   print '(a)', 'ABIN version '//ABIN_VERSION
   print '(a, a, 1x, a)', 'Compiled at ', __TIME__, __DATE__
   print '(a)', 'Git commit '//GIT_COMMIT
!$ print'(a)', 'Compiled with parallel OpenMP support for PIMD.'
#ifdef USE_FFTW
   print '(a)', 'Compiled with FFTW support.'
#endif
#ifdef USE_CP2K
   print '(a)', 'Compiled with in-built CP2K interface.'
#endif
#ifdef USE_PLUMED
   print '(a)', 'Compiled with PLUMED (static lib).'
#endif
#ifdef USE_MPI
   print '(a)', 'Compiled with MPI support.'
   print '(a)', '(used for REMD and direct CP2K and TeraChem interfaces.)'
#endif
   print '(a)', ' '

   print '(a)', 'This program was compiled by:'
   print '(a)', compiler_version()
   print '(a)', 'using the compiler options:'
   print '(a)', compiler_options()
end subroutine print_compile_info
