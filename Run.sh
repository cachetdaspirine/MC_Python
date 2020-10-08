#!/bin/bash

SerieNum=1
QueueType=qbigmem

rm -rf "Res/Serie$SerieNum"
mkdir "Res/Serie$SerieNum"

cp Parameter.py "Res/Serie$SerieNum/Parameter.py"

./AnnealingLoop.py  $SerieNum
#sbatch MonoAggregateAnnealing.pbs $SerieNum
