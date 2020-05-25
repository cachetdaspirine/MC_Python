#!/bin/bash

SerieNum=1

rm -rf "Res/Serie$SerieNum"
mkdir "Res/Serie$SerieNum"

cp Parameter.py "Res/Serie$SerieNum/Parameter.py"

echo "sbatch Run_Annealing.pbs $SerieNum"
sbatch Run_Annealing.pbs $SerieNum
