#!/bin/bash

for i in 200 400 1000 2000 4000 
do 
	j=$(($i/200))
#        j=$(($i+1))
  	mkdir -p t-$j/
	cd t-$j/
	cp ../dyn.x ../input ../run_multinode.sh ../fort.* .
	sed -i "s/input-ntime/$i/g" input
        sbatch run_multinode.sh
	cd ../
done
