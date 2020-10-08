#!/home/hugo/anaconda3/bin/python3

import numpy as np
from Anneal import *
import random as rd
import os
import sys

if len(sys.argv)<2:
    print('Number of the serie not specified, please enter a serie name')
    sys.exit()

SerieNum=sys.argv[1]
sys.path.insert(0,'Res/Serie'+str(SerieNum))
from Parameter import *

time_start = time.perf_counter()


for N in range(Npmin,Npmax+1):
    os.system('mkdir Res/Serie'+str(SerieNum)+'/N_'+str(N))

with open('Res/Serie'+str(SerieNum)+'/Energy.out', 'w') as myfile:
    myfile.write('Number_of_particle Energy Final_Acceptance_rate\n')
for S in range(NRepetition):
    #seed=rd.randint(198,1500)
    for N in range(Npmin,Npmax+1):
        Energy,Statfinal=Annealing(Kmain=Kmain,
                                Kcoupling=Kcoupling,
                                Eps=Eps,
                                KVOL=KVOL,
                                J=J,
                                SizeX=2*N+10,
                                SizeY=2*N+10,
                                NumberOfParticle=N,
                                SimNum=S,
                                Path='Res/Serie'+str(SerieNum)+'/N_'+str(N),
                                TimeStepTot=TimeStepTot,
                                Seed=seed)
        with open('Res/Serie'+str(SerieNum)+'/Energy.out', 'a') as myfile:
            myfile.write(str(N)+' '+str(Energy)+' '+str(Statfinal)+'\n')
