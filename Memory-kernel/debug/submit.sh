#!/bin/bash

for i in 1 2 10 
do 
#	j=$((10**$i))
#        j=$(($i+1))
  	mkdir -p it-$i/
	cd it-$i/
	cp ../dyn.x ../input ../run_multinode.sh ../fort.* .
	sed -i "s/input-itermax/$i/g" input
        sbatch run_multinode.sh
	cd ../
done
