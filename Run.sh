#!/bin/bash

#SerieNum=1
#QueueType=qbigmem

#rm -rf "Res/Serie$SerieNum"
#mkdir "Res/Serie$SerieNum"

#cp Parameter.py "Res/Serie$SerieNum/Parameter.py"

#./AnnealingLoop.py  $SerieNum
#sbatch MonoAggregateAnnealing.pbs $SerieNum

for SimNum in {0..9}
do
	rm -rf "Res/Serie$SimNum"
	mkdir "Res/Serie$SimNum"
	sed "3s/.*/SimNum =+$SimNum/" Parameter.py > Res/Serie$SimNum/Parameter.py
	sbatch MonoAggregateAnnealing.pbs $SimNum
done
