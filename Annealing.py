#!/home/hugo/anaconda3/bin/python3
#!/usr/bin/python3

from System import* 
from BinarySystem import *
from McMove import *
import time
import os
from pympler import asizeof

time_start = time.perf_counter()
SimNum=2
os.system('rm -rf Res/Sim'+str(SimNum))
os.system('mkdir Res/Sim'+str(SimNum))

Output=True
with open('Res/Sim'+str(SimNum)+'/Energy.out','w') as myfile:
    myfile.write('time ElasticEnergy SurfaceEnergy TotalEnergy \n')
#  ____                                              _                       
# |  _ \    __ _   _ __    __ _   _ __ ___     ___  | |_    ___   _ __   ___ 
# | |_) |  / _` | | '__|  / _` | | '_ ` _ \   / _ \ | __|  / _ \ | '__| / __|
# |  __/  | (_| | | |    | (_| | | | | | | | |  __/ | |_  |  __/ | |    \__ \
# |_|      \__,_| |_|     \__,_| |_| |_| |_|  \___|  \__|  \___| |_|    |___/

TimeStepTot=2000
StatTime=100#TimeStepTot//100
Shaking=100*100
BetaInitial=10
BetaFinal=1.6*10**2
Seed=98987
DEG=0.0125

def CoolDown(time,DE0):
    #return 10**5
    #print(DE0)
    return -6/(7*DE0)*np.log(1-time/TimeStepTot)+1/(7*DE0)
    #return BetaInitial+time/TimeStepTot*(BetaFinal-BetaInitial)
#  ____                  _                      
# / ___|   _   _   ___  | |_    ___   _ __ ___  
# \___ \  | | | | / __| | __|  / _ \ | '_ ` _ \ 
#  ___) | | |_| | \__ \ | |_  |  __/ | | | | | |
# |____/   \__, | |___/  \__|  \___| |_| |_| |_|
#          |___/

Kmain=1.
Kcoupling=0.1
Eps=0.1
KVOL=8.
#----------------
J=0.072#0.135
#----------------
SizeX=100
SizeY=100
NumberOfParticle=24


#   ___    _   _   _____   ____    _   _   _____ 
#  / _ \  | | | | |_   _| |  _ \  | | | | |_   _|
# | | | | | | | |   | |   | |_) | | | | |   | |  
# | |_| | | |_| |   | |   |  __/  | |_| |   | |  
#  \___/   \___/    |_|   |_|      \___/    |_|

with open('Res/Sim'+str(SimNum)+'/Parameter.out','w') as myfile:
    myfile.write('TimeStepTot '+str(TimeStepTot)+'\n')
    myfile.write('StatTime '+str(StatTime)+'\n')
    myfile.write('BetaInitial '+str(BetaInitial)+'\n')
    myfile.write('BetaFinal '+str(BetaFinal)+'\n')
    myfile.write('DE0G '+str(DEG)+'\n')
    myfile.write('Kmain '+str(Kmain)+'\n')
    myfile.write('Kcoupling '+str(Kcoupling)+'\n')
    myfile.write('Eps '+str(Eps)+'\n')
    myfile.write('KVOL '+str(KVOL)+'\n')
    myfile.write('J '+str(J)+'\n')
    myfile.write('SizeX '+str(SizeX)+'\n')
    myfile.write('SizeY '+str(SizeY)+'\n')
    myfile.write('NumberOfParticles '+str(NumberOfParticle)+'\n')
    
#  ___           _   _     _           _   _                  ____                  _                      
# |_ _|  _ __   (_) | |_  (_)   __ _  | | (_)  ____   ___    / ___|   _   _   ___  | |_    ___   _ __ ___  
#  | |  | '_ \  | | | __| | |  / _` | | | | | |_  /  / _ \   \___ \  | | | | / __| | __|  / _ \ | '_ ` _ \ 
#  | |  | | | | | | | |_  | | | (_| | | | | |  / /  |  __/    ___) | | |_| | \__ \ | |_  |  __/ | | | | | |
# |___| |_| |_| |_|  \__| |_|  \__,_| |_| |_| /___|  \___|   |____/   \__, | |___/  \__|  \___| |_| |_| |_|
#                                                                    |___/
rd.seed(Seed)
Beta=BetaInitial
BinSyst=BinarySystem(SizeX,SizeY)
MC=MonteCarlo(NumberOfParticle,SimNum)
for n in range(NumberOfParticle):
    BinSyst.AddMonoAggregateParticle()
    BinSyst.CheckExpansion()
for n in range(Shaking):
    BinSyst.RmRandContiguousParticle()
    BinSyst.AddMonoAggregateParticle()
    Xg,Yg=BinSyst.ComputeCenter()
    BinSyst.TranslateInTheMiddle(Xg,Yg)
system=System(BinSyst.array,eps=Eps,Kmain=Kmain,Kcoupling=Kcoupling,Kvol=KVOL)

print(" __  __           _             _                             ")
print("|  \/  |   __ _  (_)  _ __     | |       ___     ___    _ __  ")
print("| |\/| |  / _` | | | | '_ \    | |      / _ \   / _ \  | '_ \ ")
print("| |  | | | (_| | | | | | | |   | |___  | (_) | | (_) | | |_) |")
print("|_|  |_|  \__,_| |_| |_| |_|   |_____|  \___/   \___/  | .__/ ")
print("                                                       |_|    ")
Beta=0
for t in range(1,TimeStepTot):
    Success=True
    #------Energy before the move---------------------
    Eiel=system.Energy
    Eisurf=J*len(BinSyst.BoundarySite)
    #------Make the move------------------------------    
    MC.McMove(BinSyst)

    #------------------------------------
    #system.SetNodesPosition()
    #system.Evolv(BinSyst.array)
    #system.StoreNodesPosition()
    CopySystem=System(old_system=system)
    CopySystem.Evolv(BinSyst.array)
    #------Store the Energy after the move------------
    Eaftel=CopySystem.Energy
    Eaftsurf=J*len(BinSyst.BoundarySite)

    if((Eaftel+Eaftsurf)-(Eiel+Eisurf)>1):
        print(system.Energy)
        print(CopySystem.Energy)
        #system.PlotPerSite()
        #CopySystem.PlotPerSite()
        del(system)
        system=System(BinSyst.array,eps=Eps,Kmain=Kmain,Kcoupling=Kcoupling,Kvol=KVOL)
        print(system.Energy)
        print("Remake the system")
        #system.PlotPerSite()
        #system.PlotPerSpring()
    
    #------see wether we accept the move or not-------
    if rd.uniform(0,1)>np.exp(-((Eaftel+Eaftsurf)-(Eiel+Eisurf))*Beta) :
        #--Move refused-------------------------------
        MC.Reverse(BinSyst)
        del(CopySystem)
        Success=False
    else:
        del(system)
        system=CopySystem
    #------keep track of the success/fail-------------
    MC.Count(Success,(Eaftel+Eaftsurf)-(Eiel+Eisurf))
    #------Check if the boundary has to be expanded---
    #BinSyst.CheckExpansion()
    
    Xg,Yg=BinSyst.ComputeCenter()
    BinSyst.TranslateInTheMiddle(Xg,Yg)
    
    #------Cool down the system ----------------------
    if t>StatTime:
        Beta=CoolDown(t,MC.avDE/5.)
    #------Make the stats and adapt the McMove--------
    if t%StatTime==0:        
        print("time=",t)
        MC.MakeStat(t,Beta)
        with open('Res/Sim'+str(SimNum)+'/Energy.out','a') as myfile:
            myfile.write(str(t)+" "+str(system.Energy/system.Np)+" "+str(J*len(BinSyst.BoundarySite)/system.Np)+" "+str((system.Energy+J*len(BinSyst.BoundarySite))/system.Np)+"\n")
        if Output:
            system.PrintPerSite('Res/Sim'+str(SimNum)+'/Site_time'+str(t)+'.res')
            system.PrintPerSpring('Res/Sim'+str(SimNum)+'/Spring_time'+str(t)+'.res')
system.PrintPerSite('Res/Sim'+str(SimNum)+'/Site_Final.res')
system.PrintPerSpring('Res/Sim'+str(SimNum)+'/Spring_Final.res')
#system.PlotPerSite()

print(time.perf_counter() - time_start)
