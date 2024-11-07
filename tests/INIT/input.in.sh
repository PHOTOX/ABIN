&general
irest=0,
nstep=1,
pot='_mmwater_'
mdtype='sh',
dt=20.,	
/

&nhcopt
inose=0,
/

! Invalid SH parameters
&sh
nstate=100
istate_init=101
integ='invalid'
adjmom=0
couplings='baeck-an'
nac_accu1=7
nac_accu2=8
decoh_alpha=0.0D0
/
