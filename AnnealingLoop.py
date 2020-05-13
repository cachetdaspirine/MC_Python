#!/usr/bin/python3
#!/home/hugo/anaconda3/bin/python3

from Anneal import *
import random as rd
import os

time_start = time.perf_counter()

SerieNum=1

os.system('rm -rf Res/Serie'+str(SerieNum))
os.system('mkdir Res/Serie'+str(SerieNum))

Kmain=1.
Kcoupling=1.
Eps=0.1
KVOL=10.
J=0.13
Npmax=140
Npmin=140
NRepetition=5
TimeStepTot=int(2.5*10**4)
seed=None
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
