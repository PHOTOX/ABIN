!---This is a template file for user-defined analysis function.
!---This function will be called if anal_ext=1 in input(section general).
!---Should the user need something more then coordinates(velocities,forces),
! he/she must also modify  analysis.f90 and possibly also  abin.f90
module mod_analyze_ext
   use mod_const, only: DP
   use mod_array_size, only: nbinmax
   implicit none
   private
   public :: analyze_ext
   real(DP)  :: spec(nbinmax)
   real(DP)  :: emin=0,emax=0.35
   integer :: nbinen=400
   save
   contains
   subroutine analyze_ext(x,y,z,vx,vy,vz,amt)
   use mod_array_size
   use mod_general
   use mod_system, ONLY: names
   real(DP) x(:,:),y(:,:),z(:,:)
   real(DP) vx(:,:),vy(:,:),vz(:,:)
   real(DP) amt(:,:)


   if(modulo(it,nwrite).eq.0)then
!--print output every nwrite steps              
   endif
   
   end

end module
                                        



