import numpy as np
import random as rd
from BinarySystem import *
import os

class MonteCarlo:
    def __init__(self,Np=1,SimNum=0,Path='Res/'):
        self.Path=Path
        self.Success=0
        self.Refuse=0
        self.DEP,self.DEN=0,0
        self.DEPA=0 #number of step with DE positiv and accepted
        self.DE=0
        self.radius=np.inf
        self.Nmove=10
        self.Moved=list()
        self.Np=Np
        self.SimNum=SimNum
        with open(self.Path+'Sim'+str(self.SimNum)+'/Stat.out','w') as myfile:
            myfile.write('time Beta AcceptanceRate RefusalRate Nmove Radius\n')
        with open(self.Path+'Sim'+str(self.SimNum)+'/AdvanceStat.out','w') as myfile:
            myfile.write('time Beta PositiveDERate NegativeDERate AcceptedPositiveDERate average_Positive_DE\n')
    def McMove(self,BinSyst):
        self.Moved.clear()
        for _ in range(1):#self.Nmove):
            I0,J0=BinSyst.RmRandContiguousParticle()            
            I1,J1=BinSyst.AddMonoAggregateParticle(I0,J0,self.radius)
            self.Moved.append((I0,J0,I1,J1))
    def Reverse(self,BinSyst):
        for site in reversed(self.Moved):
            BinSyst.ReverseMove(site[0],site[1],site[2],site[3])
    def Count(self,Success,DE=0):
        #self.DE+=abs(DE)
        if DE>0:
            self.DE+=DE
            self.DEP+=1
        else:
            self.DEN+=1
        if Success:
            if self.Moved[0][0]!=self.Moved[0][2] or self.Moved[0][1]!=self.Moved[0][3]:
                self.Success+=1
            else :
                self.Refuse+=1
        else:
            self.Refuse+=1
        if Success and DE>=0:
            self.DEPA+=1        
    def MakeStat(self,time,Beta):
        Ntot=self.Success+self.Refuse
        #DEPArate=self.DEPA/Ntot
        if self.DEP!=0:
            self.DEPArate=self.DEPA/self.DEP
            self.avDE=self.DE/self.DEP
        else :
            self.DEPArate=0.
            self.avDE=0.
        DEPrate=self.DEP/Ntot
        DENrate=self.DEN/Ntot
        RefusalRate=self.Refuse/Ntot
        self.AcceptanceRate=self.Success/Ntot
        with open(self.Path+'Sim'+str(self.SimNum)+'/Stat.out','a') as myfile:
            myfile.write(str(time)+' '+str(Beta)+' '+str(self.AcceptanceRate)+' '+str(RefusalRate)+' ')
            myfile.write(str(self.Nmove)+' '+str(self.radius)+'\n')
        with open(self.Path+'Sim'+str(self.SimNum)+'/AdvanceStat.out','a') as myfile:
            myfile.write(str(time)+' '+str(Beta)+' '+str(DEPrate)+' '+str(DENrate)+' '+str(self.DEPArate))
            myfile.write(' '+str(self.avDE)+'\n')
        if self.AcceptanceRate > 0.6:
            self.Harder()
        elif self.AcceptanceRate < 0.3 :
            self.Softer()
        self.DE,self.DEN,self.DEP,self.DEPA=0,0,0,0
        self.Success=0
        self.Refuse=0
    def Harder(self):
        # we start by increasing the radius if it's not infinity
        # there are 10 steps of increasment from Np/20  to  Np/2
        # after Np/2 the radius becomes infinity
        if self.radius!=np.inf:
            if self.radius>=self.Np/2:
                self.radius=np.inf
            else :
                self.radius+=self.Np//20
        elif self.Nmove<=self.Np//10:
            #if the radius is already infinity we  multiply  the
            #number of move per step by 2
            self.Nmove+=1
    def Softer(self):
        if self.radius==self.Np//20 and self.Nmove==1:
            return
        if self.radius==np.inf and self.Nmove==1:
            self.radius=self.Np//2
        elif self.radius==np.inf :
            self.Nmove=self.Nmove//2
        elif self.radius>max(self.Np//20,3):
            self.radius-=self.Np//20
        
            
