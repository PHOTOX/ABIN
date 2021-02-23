! This file must be recompiled at each compilation
! to get the current date and timem
! passed via make parameter (-DCOMPILE_DATE=)
subroutine print_compile_info()
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 6
   use iso_fortran_env, only: compiler_version, compiler_options
#endif

   ! DATE and COMMIT are defined and exported in Makefile
   print*,'Compiled at ', COMPILE_DATE
   print*,'Git commit '//GIT_COMMIT
!$ print*,'Compiled with parallel OpenMP support for PIMD.'
#ifdef USE_FFTW
   write (*, *) 'Compiled with FFTW support.'
#endif
#ifdef USE_CP2K
   write (*, *) 'Compiled with in-built CP2K interface.'
#endif
#ifdef USE_PLUMED
   write (*, *) 'Compiled with PLUMED (static lib).'
#endif
#ifdef USE_MPI
   write (*, *) 'Compiled with MPI support.'
   write (*, *) '(used for REMD and direct CP2K and TeraChem interfaces.)'
#endif
   print*,' '

#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 6
   print*,'This program was compiled by ', &
      compiler_version(), ' using the options: '
   print*,compiler_options()
#endif
end subroutine print_compile_info
