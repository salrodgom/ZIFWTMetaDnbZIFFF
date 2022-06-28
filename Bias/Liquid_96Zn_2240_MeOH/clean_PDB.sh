#!/bin/bash
timestep=$1
if [ ! -z "$timestep" ] ; then echo "TimeStep=$timestep" ; else timestep=0.0 ; fi
echo "15 " | gmx_mpi traj -f MDNVT.xtc -s p1_fromPDFFile.pdb -oxt traj.pdb -b $timestep -tu ps -nojump yes -n index.ndx
rm -rf \#*
exit 0
