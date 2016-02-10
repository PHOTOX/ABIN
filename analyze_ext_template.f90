!---This is a template file for user-defined analysis function.
!---This function will be called if anal_ext=1 in input (section general).
module mod_analyze_ext
   use mod_const, only: DP
   implicit none
   private
   public    :: analyze_ext
   save
   contains

   subroutine analyze_ext(x,y,z,vx,vy,vz,amt)
   use mod_general, only: it
   real(DP) x(:,:),y(:,:),z(:,:)
   real(DP) vx(:,:),vy(:,:),vz(:,:)
   real(DP) amt(:,:)
   character(len=100)   :: chsystem
   character(len=50)   :: chit


   ! Launch external BASH script.
   ! Can be used to analyze wavefunction on-the-fly
   write(chit,*)it
   chsystem='./analyze_ext.sh '//trim(chit)
   call system(chsystem)

!   if(modulo(it,nwrite).eq.0)then
!--print output every nwrite steps              
!   endif
   
   end subroutine analyze_ext

end module
                                        

