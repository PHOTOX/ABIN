#!/usr/bin/awk -f

# Tiny ABIN awk script evaluating energy conservation from energies.dat
# In case of multiple trajectories, itrj parameter defines the index of the trajectory.

# Should be called as: 
# checkenergy -v itrj=$i inputfile

function abs(x){return ((x < 0.0) ? -x : x)}
BEGIN{max_de=0.0; en0=0}
{
   if ($1 == "#")
      next
   if ( en0==0 ) { 
      en=$4
      en0=$4
      next
   }
   de=abs(en-$4)
   if (de > max_de ) {
      max_de=de
   }
   en=$4
}

END{
print itrj, max_de*27.2114, 27.2114*(en0-en)
} 
