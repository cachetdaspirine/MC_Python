#!/bin/bash

SerieNum=0

rm -rf "Res/Serie$SerieNum"
mkdir "Res/Serie$SerieNum"

cp Parameter.py "Res/Serie$SerieNum/Parameter.py"

sbatch MonoAggregateAnnealing.pbs $SerieNum
