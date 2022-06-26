#!/bin/bash
timestep=$2
if [ ! -z "$timestep" ] ; then echo "TimeStep=$timestep" ; else timestep=0.0 ; fi
echo "0 " | gmx traj -f $1 -s p1.pdb -oxt traj.pdb -b $timestep -tu ps -nojump yes -dt 10
sed -i -e '/Ne/d' -e '/H5/d' -e 's/He/Si/g' -e '/Zn/d' traj.pdb
rm -rf \#*
exit 0
