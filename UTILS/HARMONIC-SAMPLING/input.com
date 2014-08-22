if(NPROC_MPPX.gt.1)then
skipped
endif
if(.not.modul_aims)then
skipped
endif

ANGSTROM
SYMMETRY,nosym
ORIENT,noorient
GEOMTYP=xyz
GEOMETRY=./FMSINPOUT/geometry.xyz
BASIS=sto-3g

MPINIT={
hf
}

MPCALC={
 ! Calculate Energy and Forces
 show,itimestep,itrj;
 hf
 ener( 1 )=energy( 1 )
 if(status.lt.0) then
    text,ERROR: could not converge SCF.
    STOP
 endif
 pop;density
 force
}

aims,dir='./FMSINPOUT'
