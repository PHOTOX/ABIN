***Abin input
memory, 96, m;
file, 2, input003.wfu, unknown
PUNCH, input003.pun, new
gprint, orbital, civector

Angstrom

geometry=../geom.dat.003

basis=def2-SVP

!-for simple CASSCF, you don't need to modify anything below this

!-we need to get rid of the SAMC records in file 2 (input.wfu,restart file)
!-otherwise, the forces and NACME are wrong for the following reason
!-cpmscf will not rewrite itself it but ather write into following records
!-but the subsequent call to forces would read from old records -> wrong numbers
!-we use file 2 for forces and NACME due to df-casscf
data,truncate,2,3000

if (lastorb.ne.MCSCF)then
   {hf;wf,10,0,0}
   df-casscf;
   occ,7;
   closed,2;
   WF,10,0,0;
   state,4;
   maxiter,40;
   orbital,2101.2
   save,ci=2501.2
endif

data, copy, 2101.2, 3000.2

{df-casscf;
occ,7;
closed,2;
WF,10,0,0;
state,4;
maxiter,2;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1
}

if (status.lt.0) then
   text, MCSCF failed to converge.
   text, Attempting uncoupled iterations.
   text, Enlarging PSPACE.
   {df-casscf;
   occ,7;
   closed,2;
   WF,10,0,0;
   ! Info about pspace: https://www.molpro.net/info/2015.1/doc/manual/node244.html
   ! uncomment in case of convergence difficulties...
   pspace, 2;
   state,4;
   maxiter,40;
   orbital,2101.2
   ciguess,2501.2
   save,ci=2501.2
   diab,3000.2,save=2101.2,method=-1
   {iterations
   do,uncouple,1,to,5}
   }
   if (status.lt.0) then
      text, Uncoupled iterations did not help, exiting.
      STOP
   endif

endif


TEXT, Calculating gradients. 
df-casscf;
occ,7;
closed,2;
WF,10,0,0;
state,4;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1
DM;  ! calculate dipole moments

cpmcscf, grad, 1.1, save=5100.2, accu=1d-6;
cpmcscf, grad, 2.1, save=5101.2, accu=1d-6;
cpmcscf, grad, 3.1, save=5102.2, accu=1d-6;
cpmcscf, grad, 4.1, save=5103.2, accu=1d-6;
forces; samc, 5100.2;
forces; samc, 5101.2;
forces; samc, 5102.2;
forces; samc, 5103.2;
pop; density,2101.2,state=1.1
pop; density,2101.2,state=2.1
pop; density,2101.2,state=3.1
pop; density,2101.2,state=4.1
