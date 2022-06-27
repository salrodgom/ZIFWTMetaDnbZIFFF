#!/bin/bash
timestep=$1
if [ ! -z "$timestep" ] ; then echo "TimeStep=$timestep" ; else timestep=0.0 ; fi
echo "0 " | gmx traj -f MDNVT.xtc -s p1_fromPDFFile.pdb -oxt traj.pdb -b $timestep -tu ps -nojump yes
rm -rf \#*
exit 0
