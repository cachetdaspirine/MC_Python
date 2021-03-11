from System import*
from BinarySystem import *
from McMove import *
import time
import os

def CoolDown(time,DE0,Ttot):#,BetaInitial,BetaFinal,DE0):
    return -6/(7*DE0)*np.log(1-time/Ttot)+1/(7*DE0)
    #return BetaInitial+time/TimeStepTot*(BetaFinal-BetaInitial)

def Annealing(
        Kmain=1.,
        Kcoupling=0.1,
        Eps=0.1,
        KVOL=8.,
        J=0.135,
        SizeX=150,
        SizeY=150,
        NumberOfParticle=100,
        SimNum=1,
        Path='',
        BetaInitial=0,
        TimeStepTot=2*10**3,
        Seed=98986,
        ParticleType='Triangle'
        ):
    #os.system('rm -rf '+Path+'Sim'+str(SimNum))
    #os.system('mkdir '+Path+'Sim'+str(SimNum))

    with open(Path+'/Sim'+str(SimNum)+'_Energy.out','w') as myfile:
        myfile.write('time ElasticEnergy SurfaceEnergy TotalEnergy \n')

    StatTime=TimeStepTot//100
    Shaking=NumberOfParticle**2

    with open(Path+'/Sim'+str(SimNum)+'_Parameter.out','w') as myfile:
        myfile.write('ParticleType '+ParticleType+'\n')
        myfile.write('TimeStepTot '+str(TimeStepTot)+'\n')
        myfile.write('StatTime '+str(StatTime)+'\n')
        myfile.write('BetaInitial '+str(BetaInitial)+'\n')
        #myfile.write('BetaFinal '+str(BetaFinal)+'\n')
        #myfile.write('DE0 '+str(DE0)+'\n')
        myfile.write('Kmain '+str(Kmain)+'\n')
        myfile.write('Kcoupling '+str(Kcoupling)+'\n')
        myfile.write('Eps '+str(Eps)+'\n')
        myfile.write('KVOL '+str(KVOL)+'\n')
        myfile.write('J '+str(J)+'\n')
        myfile.write('SizeX '+str(SizeX)+'\n')
        myfile.write('SizeY '+str(SizeY)+'\n')
        myfile.write('NumberOfParticles '+str(NumberOfParticle)+'\n')
    rd.seed(Seed)
    np.random.seed(Seed)
    BinSyst=BinarySystem(SizeX,SizeY,ParticleType=ParticleType)
    MC=MonteCarlo(NumberOfParticle,SimNum,Path=Path)
    for n in range(NumberOfParticle):
        BinSyst.AddMonoAggregateParticle()
        BinSyst.CheckExpansion()
    for n in range(Shaking):
        BinSyst.RmRandContiguousParticle()
        BinSyst.AddMonoAggregateParticle()
        Xg,Yg=BinSyst.ComputeCenter()
        BinSyst.TranslateInTheMiddle(Xg,Yg)
    system=System(BinSyst.array,eps=Eps,Kmain=Kmain,Kcoupling=Kcoupling,Kvol=KVOL,ParticleType=ParticleType)
    print(" __  __           _             _                             ")
    print("|  \/  |   __ _  (_)  _ __     | |       ___     ___    _ __  ")
    print("| |\/| |  / _` | | | | '_ \    | |      / _ \   / _ \  | '_ \ ")
    print("| |  | | | (_| | | | | | | |   | |___  | (_) | | (_) | | |_) |")
    print("|_|  |_|  \__,_| |_| |_| |_|   |_____|  \___/   \___/  | .__/ ")
    print("                                                       |_|    ")
    Beta=BetaInitial
    for t in range(1,TimeStepTot):
        #------Cool down the system ----------------------
        #Beta=CoolDown(t,BetaInitial,BetaFinal,DE0)
        Success=True
        #------Energy before the move---------------------
        Eiel=system.Energy
        Eisurf=J*BinSyst.GetSurface()#len(BinSyst.BoundarySite)
        #------Make the move------------------------------
        MC.McMove(BinSyst)
        #------------------------------------
        CopySystem=System(old_system=system)
        CopySystem.Evolv(BinSyst.array)
        #------Store the Energy after the move------------
        Eaftel=CopySystem.Energy
        Eaftsurf=J*BinSyst.GetSurface()

        # if((Eaftel+Eaftsurf)-(Eiel+Eisurf)>20*J):
            # print(system.Energy)
            # print(CopySystem.Energy)
            # del(system)
            # system=System(BinSyst.array,eps=Eps,Kmain=Kmain,Kcoupling=Kcoupling,Kvol=KVOL,ParticleType=ParticleType)
            # print(system.Energy)

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
        #-----Re-center the system in the middle
        Xg,Yg=BinSyst.ComputeCenter()
        BinSyst.TranslateInTheMiddle(Xg,Yg)
        if t>StatTime:
            if MC.avDE==0:
                    system.PrintPerSite(Path+'/Sim'+str(SimNum)+'_Site_Final.res')
                    system.PrintPerSpring(Path+'/Sim'+str(SimNum)+'_Spring_Final.res')
                    system.PrintSpringPerSite(Path+'/Sim'+str(SimNum)+'_SpringSite_Final.res')
                    return (system.Energy+J*BinSyst.GetSurface())/system.Np, MC.AcceptanceRate
            Beta=CoolDown(t,MC.avDE/8.,TimeStepTot)
        #------Make the stats and adapt the McMove--------
        if t%StatTime==0:
            print("time=",t)
            MC.MakeStat(t,Beta)
            with open(Path+'/Sim'+str(SimNum)+'_Energy.out','a') as myfile:
                myfile.write(str(t)+
                             " "+
                             str(system.Energy/system.Np)+
                             " "+
                             str(J*BinSyst.GetSurface()/system.Np)+
                             " "+
                             str((system.Energy+J*BinSyst.GetSurface())/system.Np)+
                             "\n")
    system.PrintPerSite(Path+'/Sim'+str(SimNum)+'_Site_Final.res')
    system.PrintPerSpring(Path+'/Sim'+str(SimNum)+'_Spring_Final.res')
    system.PrintSpringPerSite(Path+'/Sim'+str(SimNum)+'_SpringSite_Final.res')
    return (system.Energy+J*BinSyst.GetSurface())/system.Np, MC.AcceptanceRate
