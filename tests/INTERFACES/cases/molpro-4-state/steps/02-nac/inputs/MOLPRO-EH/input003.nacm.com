***Abin input
memory,96,m;
file,2,input003.nacm.wfu,old
PUNCH,input003.nacm.pun,new
!gprint,orbital,civector

Angstrom

geometry=../geom.dat.003

basis=def2-SVP

! for simple CASSCF, you don't need to modify anything below this

!data,truncate,2,5100

if (lastorb.ne.MCSCF)then
   TEXT, Whoops, this should not happen! Exiting..
   TEXT, ERROR: Missing CASSCF wavefunction
   stop
endif

df-casscf;
occ,7;
closed,2;
WF,10,0,0;
state,4;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1

cpmcscf, nacm, 1.1, 2.1, save=5001.2, accu=1d-6;
cpmcscf, nacm, 1.1, 3.1, save=5002.2, accu=1d-6;
cpmcscf, nacm, 1.1, 4.1, save=5003.2, accu=1d-6;
cpmcscf, nacm, 2.1, 3.1, save=5004.2, accu=1d-6;
cpmcscf, nacm, 2.1, 4.1, save=5005.2, accu=1d-6;
forces; nacm, 5001.2;
forces; nacm, 5002.2;
forces; nacm, 5003.2;
forces; nacm, 5004.2;
forces; nacm, 5005.2;

df-casscf;
occ,7;
closed,2;
WF,10,0,0;
state,4;
orbital,2101.2
ciguess,2501.2
save,ci=2501.2
diab,3000.2,save=2101.2,method=-1

cpmcscf, nacm, 3.1, 4.1, save=5006.2, accu=1d-6;
forces; nacm, 5006.2;
