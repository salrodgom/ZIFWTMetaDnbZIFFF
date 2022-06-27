#!/bin/sh -f
plumed sum_hills --bin 512 --min 0 --max 1 --hills $1 --negbias --outfile $2 --kt 0.025679497211149798
