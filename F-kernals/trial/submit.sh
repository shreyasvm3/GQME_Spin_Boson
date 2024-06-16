#!/bin/bash

for i in 5  #{1..4..1} 
do 
#	j=$((10**$i))
        j=$(($i+1))
  	mkdir -p ntraj-10-$j/
	cd ntraj-10-$j/
	cp ../dyn.x ../input ../run_multinode.sh  .
	sed -i "s/input-ntraj/$i/g" input
        sbatch run_multinode.sh
	cd ../
done
