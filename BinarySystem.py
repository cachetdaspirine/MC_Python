import numpy as np
import random as rd
import copy
TopologieDownHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieUpHex = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
TopologieDownTriangle = [(1,0),(-1,0),(0,1)]
TopologieUpTriangle = [(1,0),(-1,0),(0,-1)]
def distance(i,j,i_,j_):
    return ((i-i_)**2+(j-j_)**2)**0.5

class BinarySystem :
    def __init__(self, sizeX,sizeY,ParticleType='Triangle'):
        if ParticleType == 'Triangle':
            self.TopologieUp = TopologieUpTriangle
            self.TopologieDown = TopologieDownTriangle
        elif ParticleType=='Hexagon':
            self.TopologieUp = TopologieUpHex
            self.TopologieDown = TopologieDownHex
        self.Lx=int(sizeX)
        self.Ly=int(sizeY)
        self.array=np.array([np.zeros(self.Lx,dtype=int) for _ in range(self.Ly)])
        self.Np=0
        self.BoundarySite=set()
        self.OccupiedSite=set()
    def OutputBinary(self,filename='Binary.txt'):
        np.savetxt(filename,self.array)
    def InputBinary(self,filename='Binary.txt'):
        self.array=np.loadtxt(filename,dtype=int)
    def PrintBinary(self):
        for j in reversed(range(self.array.shape[1])):
            for i in range(self.array.shape[0]):
                print(str(self.array[i,j])+" ",end='')
            print('\n',end='')
    def GetSurface(self):
        surf = 0
        #print(self.OccupiedSite.__len__())
        for site in self.OccupiedSite:
            surf+=self.GetBorderNeighbors(site[0],site[1])
        return surf
    def ReverseMove(self,Irm,Jrm,Iadd,Jadd):
        self.RmParticle(Iadd,Jadd)
        self.UpdateAfterRmMono(Iadd,Jadd)
        self.AddParticle(Irm,Jrm)
        self.UpdateAfterAddMono(Irm,Jrm)
    def AddMonoAggregateParticle(self,irm=0,jrm=0,Radius=np.inf):
        #add a particle close to irm,jrm in a given radius
        if self.Np==0:
            self.array[self.Lx//2][self.Ly//2]=1 #if there is no particle, add a newone ine the middle
            self.OccupiedSite.add((self.Lx//2,self.Ly//2))
            i,j=self.Lx//2, self.Ly//2
            for ij in self.GetFreeNeighbors(self.Lx//2, self.Ly//2) :
                self.BoundarySite.add(ij) # we wanna store pair of index in the array boundarySite
            self.Np+=1
        else :
            i,j=self.AddRandomParticle(irm,jrm,Radius)
            self.UpdateAfterAddMono(i,j)
        return i,j
    def AddRandomParticle(self,irm,jrm,Radius):
        if Radius!=np.inf:
            CloseSite=set()
            for BoundSite in self.BoundarySite:
                if distance(irm,jrm,BoundSite[0],BoundSite[1])<Radius:
                    CloseSite.add(BoundSite)
            if len(CloseSite)!=0:
                ij=rd.sample(CloseSite,1)[0]
            else :
                ij=rd.sample(self.BoundarySite,1)[0]
        else:
            ij=rd.sample(self.BoundarySite,1)[0]

        self.array[ij[0],ij[1]]=1
        self.Np+=1
        return ij[0],ij[1]
    def AddParticle(self,i,j):
        if self.array[i,j]==1:
            print('can t add to an occupied site')
        else :
            self.array[i,j]=1
            self.Np+=1
            return i,j
    def RmParticle(self,i,j):
        #Remove a specific particle i,j but don t necessarily keep the aggregate contiguous
        if self.array[i,j]==0:
            print('Can t remove an empty site')
        else :
            self.array[i,j]=0
            self.Np-=1
            return i,j
    def RmRandContiguousParticle(self):
        Fail=True
        while Fail: #0 0 0 0 0  we continue as long as we didn't manage to remove a particle
            i,j=self.RmRandParticle() # remove
            self.UpdateAfterRmMono(i,j) # update
            Fail=self.CheckDiscontiguity(i,j) # check
            if Fail: # if it didn't work, we reverse the move
                self.AddParticle(i,j) # re-add the previously removed particle
                self.UpdateAfterAddMono(i,j) # update
        return i,j
    def CheckDiscontiguity(self,i,j): #return true is it's discontiguous and false if it's contiguous
        Changed=tuple((i,j))
        if len(self.GetOccupiedNeighbors(i,j))==1:
            return False
        for Neigh1 in self.GetOccupiedNeighbors(i,j):
            for Neigh2 in self.GetOccupiedNeighbors(i,j):
                if Neigh1!=Neigh2:
                    if not self.Linked(Changed,Neigh1,Neigh2):
                        return True
        return False
    def Linked(self,Changed,ij1,ij2):
        Clust={ij1}
        VectClust=[ij1]
        k=0
        while k!= len(Clust):
            for Neigh in self.GetOccupiedNeighbors(VectClust[k][0],VectClust[k][1]):
                if Neigh==ij2:
                    return True
                elif Neigh!=Changed and Neigh not in Clust:
                    Clust.add(Neigh)
                    VectClust.append(Neigh)
            k+=1
        return False
    def RmRandParticle(self):
        ij=rd.sample(self.OccupiedSite,1)[0]
        self.array[ij[0],ij[1]]=0
        self.Np-=1
        return ij[0],ij[1]
    def UpdateAfterAddMono(self,i,j):
        self.OccupiedSite.add((i,j))
        self.BoundarySite.remove((i,j))
        for ij in self.GetFreeNeighbors(i,j):
            self.BoundarySite.add(ij)
    def UpdateAfterRmMono(self,i,j):
        #co=copy.copy(self.BoundarySite)
        try:
            self.OccupiedSite.remove((i,j))
            self.BoundarySite.add((i,j))
        except :
            print("try to set unoccupied an already non-occupied site")
            return
        for Neigh in self.GetFreeNeighbors(i,j):
            if len(self.GetOccupiedNeighbors(Neigh[0],Neigh[1]))==0 :
                try:
                    self.BoundarySite.remove((Neigh[0],Neigh[1]))
                except :
                    print((Neigh[0],Neigh[1]))
                    print(self.BoundarySite)
                    print(co)
                    print('ij=('+str(i)+','+str(j)+')')
                    self.PrintBinary()
                    print(self.Lx)
                    print(self.Ly)
                    input()
    def GetOccupiedNeighbors(self,i,j):
        # get the real or window index
        # Choose the topologie to use depending on the up/down
        ij=(i,j)
        if (ij[0]+ij[1])%2==0:
            Res = np.array(self.TopologieDown)+np.array(ij)
        else :
            Res = np.array(self.TopologieUp)+np.array(ij)
        # regularize the result array with only the value that can be inside the state
        Resreg=np.delete(Res,np.argwhere((Res[:,0]>=self.Lx) | (Res[:,0]<0) | (Res[:,1]>=self.Ly) | (Res[:,1]<0)),0)
        #Build a numpy array of tuple
        Resbis=np.empty(Resreg.__len__(),dtype=object)
        Resbis[:] = list(zip(Resreg[:,0],Resreg[:,1]))
        #check the occupancie or not
        Resbis=Resbis[np.array([self.array[r]==1 for r in Resbis ])]
        return set(Resbis)
    def GetFreeNeighbors(self,i,j):
        # get the real or window index
        # Choose the topologie to use depending on the up/down
        ij=(i,j)
        if (ij[0]+ij[1])%2==0:
            Res = np.array(self.TopologieDown)+np.array(ij)
        else :
            Res = np.array(self.TopologieUp)+np.array(ij)
        # regularize the result array with only the value that can be inside the state
        Resreg=np.delete(Res,np.argwhere((Res[:,0]>=self.Lx) | (Res[:,0]<0) | (Res[:,1]>=self.Ly) | (Res[:,1]<0)),0)
        #Build a numpy array of tuple
        Resbis=np.empty(Resreg.__len__(),dtype=object)
        Resbis[:] = list(zip(Resreg[:,0],Resreg[:,1]))
        #check the occupancie or not
        Resbis=Resbis[np.array([self.array[r]==0 for r in Resbis ])]
        return set(Resbis)
    def GetBorderNeighbors(self, i,j):
        # get the real or window index
        # Choose the topologie to use depending on the up/down
        ij=(i,j)
        if (ij[0]+ij[1])%2==0:
            Res = np.array(self.TopologieDown)+np.array(ij)
        else :
            Res = np.array(self.TopologieUp)+np.array(ij)
        # regularize the result array with only the value that can be inside the state
        #Resreg=np.delete(Res,np.argwhere((Res[:,0]>=self.Size) | (Res[:,0]<0) | (Res[:,1]>=self.Size) | (Res[:,1]<0)),0)
        #Build a numpy array of tuple
        Resbis=np.empty(Res.__len__(),dtype=object)
        Resbis[:] = list(zip(Res[:,0],Res[:,1]))
        Resbis = list(Resbis)
        #check the occupancie or not
        for n in reversed(range(Resbis.__len__())):
            if Resbis[n][0]>=0 and Resbis[n][1]>=0 and Resbis[n][0]<self.Lx and Resbis[n][1]<self.Ly:
                if self.array[Resbis[n]]!=0:
                    del Resbis[n]
        #check the occupancie or not
        #Resbis=Resbis[np.array([self.array[r]==0 for r in Resbis ])]
        return set(Resbis).__len__()
    def CheckExpansion(self):
        for ij in self.OccupiedSite:
            if ij[0]>=self.Lx-3 or ij[0]>=self.Ly-3 or ij[0]<=2 or ij[1]<=2:
                self.ExpandSystem(ij[0],ij[1])
                break
    def ExpandSystem(self,i,j): #Expand the system, knowing that (i,j) is touching a boundary
        self.Lx=self.Lx*2
        self.Ly=self.Ly*2
        self.TranslateInTheMiddle(i,j)
    def ComputeCenter(self):
        Xg,Yg=0,0
        for ij in self.OccupiedSite:
            Xg+=ij[0]
            Yg+=ij[1]
        Xg=int(Xg/len(self.OccupiedSite))
        Yg=int(Yg/len(self.OccupiedSite))
        return Xg,Yg
    def TranslateInTheMiddle(self,i,j):
        self.array=np.array([np.zeros(self.Lx,dtype=int) for _ in range(self.Ly)])
        if (i+j)%2==0:
            MiddleX,MiddleY=self.Lx//2, self.Ly//2
        else:
            MiddleX,MiddleY=self.Lx//2, self.Ly//2+1
        NewOccupied=set()
        #NewBoundary=set()
        #co=copy.copy(self.OccupiedSite)
        for ij in self.OccupiedSite:
            NewOccupied.add((ij[0]-i+MiddleX, ij[1]-j+MiddleY))
        self.OccupiedSite=NewOccupied
        #for ij in self.BoundarySite:
        #    NewBoundary.add((ij[0]-i+MiddleX, ij[1]-j+MiddleY))
        #self.BoundarySite=NewBoundary
        for ij in self.OccupiedSite:
            try:
                self.array[ij[0],ij[1]]=1
            except:
                print(self.OccupiedSite)
                self.PrintBinary()
                print(MiddleX)
                print(MiddleY)
                print(i)
                print(j)
                print(self.OccupiedSite)
                input()
        self.BoundarySite.clear()
        for ij in self.OccupiedSite:
            for FreeNeigh in self.GetFreeNeighbors(ij[0],ij[1]):
                self.BoundarySite.add(FreeNeigh)
