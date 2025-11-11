import numpy as np
import math
class BilletBending:    
# CCM: list of sector [curvature [1/m], length [m], diameter of a roll at the end of the sector [m], axial force at the rolls [N]]
# x0 - start [m]
    def __init__(self,CCM,x0,dx=0.05):
        self.vel=1           # casting speed, m/min    
        self.Kapa0=CCM[0][0] # initial curvuture, 1/m
        self.Nf0=0           # Frictional force at the begining
        self.El={}           # list of: curvature; length; x0
        self.BCs=[]          # boundrary conditions: list of (Node, Degree of freedom)  
        self.Nw=4            # number of points for integration
        self.EpsStf=1        # convergence rate for stiffness
        self.Straight_i=0    # number of a node where straight part starts
#---------------------
        self.K0=np.zeros(6)
        self.b0=0
        self.ElDisp=np.zeros(6)
#----------
        self.KK0=np.zeros((3,6))
        self.bb0=np.zeros(3)
        self.BCs.append((0,1)) #mould - angle
        self.BCs.append((0,0)) #mould - disp
        self.Stiff0={}
        ElNum0=0
        for sector in CCM:
            ElN=int(sector[1]/dx)
            for i in range(ElNum0,ElNum0+ElN):
                self.El[i]=[sector[0],sector[1]/ElN,x0]
                x0+=sector[1]/ElN
                self.Stiff0[i]={}
                for j in range(self.Nw):
                    self.Stiff0[i][(0.5+j)/self.Nw]=0
            ElNum0+=ElN
            if sector[0]>0:
                self.StrBC=len(self.BCs)
                self.Straight_i=ElNum0
            self.BCs.append((ElNum0,0))
        self.Forces=np.zeros(len(self.BCs))       # values of force and moments in BCs
        self.AxialForces0=np.zeros(len(self.BCs)) # axial force from ferrostatic
        self.AxialForces=np.zeros(len(self.BCs))  # full axial force
        self.RollDiam=np.zeros(len(self.BCs))     # roll diameter
        for i in range(len(CCM)):
            self.RollDiam[i+2]=CCM[i][2]    #mm
            self.AxialForces0[i+2]=CCM[i][3]
#----------
        self.N=ElNum0+1 # number of nodes (minimum-2)
        #labels of nodes: 0..N-1
        #labels of elements: 0..N-2
        self.Disp=np.zeros(self.N*3)
#----------------------------
#---------FUNCTION-----------
    def TransvToAxialForce(self,Diam,Force):
        return 2*(0.003/Diam+0.35*0.02)*abs(Force)
#----------------------------
#---------FUNCTION-----------
    def Stiff(self,x,kaparate):
        if x<self.Stffx0:
            i=0
            ksi=0
        elif x>=self.StffxN*self.Stffdx+self.Stffx0:
            i=self.StffxN-1
            ksi=1
        else:
            i=int((x-self.Stffx0)/self.Stffdx)
            ksi=(x-self.Stffx0)/self.Stffdx-i
        if abs(kaparate)<self.StffKR0:
            j=0
            nu=0
        elif abs(kaparate)>=self.StffKRN*self.StffdKR+self.StffKR0:
            j=self.StffKRN-1
            nu=1
        else:
            j=int((abs(kaparate)-self.StffKR0)/self.StffdKR)
            nu=(abs(kaparate)-self.StffKR0)/self.StffdKR-j
        return self.StffV[i][j]*(1-ksi)*(1-nu)+self.StffV[i+1][j+1]*ksi*nu+self.StffV[i+1][j]*ksi*(1-nu)+self.StffV[i][j+1]*(1-ksi)*nu # N*m2*sec
#----------------------------
#---------FUNCTION-----------
    def LoadStiffness(self,FileName, separator=';'):
        f=open(FileName,'r')
        Values=f.readline().split(separator)
        self.StffKRN=len(Values)-2
        self.StffKR0=float(Values[1])
        self.StffdKR=(float(Values[-1])-self.StffKR0)/self.StffKRN
        Values=f.readline().split(separator)
        self.Stffx0=float(Values[0])
        self.StffxN=0
        txt=f.readline()
        while txt:
            Values=txt.split(separator)
            self.StffxN+=1
            xk=float(Values[0])
            txt=f.readline()
        self.Stffdx=(xk-self.Stffx0)/self.StffxN
        f.close
        self.StffV=np.zeros((self.StffxN+1,self.StffKRN+1))
        f=open(FileName,'r')
        txt=f.readline()
        txt=f.readline()
        i=0
        while txt:
            Values=txt.split(separator)
            for j in range(0,len(Values)-1):
                self.StffV[i,j]=float(Values[j+1])
            i+=1
            txt=f.readline()
        f.close
#----------------------------
#---------FUNCTION-----------
    def CalcKaparate(self,NEl,ksi):
        self.K0[0]=self.vel/60*(-60+360*ksi-360*ksi**2)/self.El[NEl][1]**3  # 1/(m2*sec)
        self.K0[1]=self.vel/60*(-36+192*ksi-180*ksi**2)/self.El[NEl][1]**2  # 1/(m*sec)
        self.K0[2]=self.vel/60*(-9+36*ksi-30*ksi**2)/self.El[NEl][1]        # 1/sec
        self.K0[3]=self.vel/60*(60-360*ksi+360*ksi**2)/self.El[NEl][1]**3   # 1/(m2*sec)
        self.K0[4]=self.vel/60*(-24+168*ksi-180*ksi**2)/self.El[NEl][1]**2  # 1/(m*sec)
        self.K0[5]=self.vel/60*(3-24*ksi+30*ksi**2)/self.El[NEl][1]         # 1/sec
        self.b0=self.vel/60*(6*(1-2*ksi)*self.El[NEl][0]/self.El[NEl][1])   # 1/(m*sec)
        for Node in range(2):
            for i in range(3):
                if (NEl+Node,i) in self.BCs:
                    self.ElDisp[Node*3+i]=0
                else:
                    self.ElDisp[Node*3+i]=self.Disp[(NEl+Node)*3+i] # {m, 1, 1/m, m, 1, 1/m}
        if NEl==0: self.ElDisp[2]=self.Kapa0
        return (np.dot(self.K0,self.ElDisp)+self.b0) # 1/(m*sec)
#----------------------------
#---------FUNCTION-----------
    def CalcForceLoc_proc(self,NEl,ksi):
        kaparate=self.CalcKaparate(NEl,ksi) # 1/(m*sec)
        EI=self.Stiff(ksi*self.El[NEl][1]+self.El[NEl][2],kaparate) # N*m2*sec
        if self.Stiff0[NEl][ksi]==0:
            self.Stiff0[NEl][ksi]=EI
        else:
            if EI<self.Stiff0[NEl][ksi]:
                EI=EI*self.EpsStf+(1-self.EpsStf)*self.Stiff0[NEl][ksi]
            else:
                EI=EI*self.EpsStf+(1-self.EpsStf)*self.Stiff0[NEl][ksi]
            self.Stiff0[NEl][ksi]=EI
        self.KK0[0][0]=self.vel/60*EI*(-360+720*ksi)/self.El[NEl][1]**4
        self.KK0[0][1]=self.vel/60*EI*(-192+360*ksi)/self.El[NEl][1]**3
        self.KK0[0][2]=self.vel/60*EI*(-36+60*ksi)/self.El[NEl][1]**2
        self.KK0[0][3]=self.vel/60*EI*(360-720*ksi)/self.El[NEl][1]**4
        self.KK0[0][4]=self.vel/60*EI*(-168+360*ksi)/self.El[NEl][1]**3
        self.KK0[0][5]=self.vel/60*EI*(24-60*ksi)/self.El[NEl][1]**2
        for i in range(6):
            self.KK0[1][i]=EI*self.K0[i]
        if self.El[NEl][0]==0:
            self.KK0[2][0]=-self.vel/60*EI*720/self.El[NEl][1]**4
            self.KK0[2][1]=-self.vel/60*EI*360/self.El[NEl][1]**3
            self.KK0[2][2]=-self.vel/60*EI*60/self.El[NEl][1]**2
            self.KK0[2][3]=self.vel/60*EI*720/self.El[NEl][1]**4
            self.KK0[2][4]=-self.vel/60*EI*360/self.El[NEl][1]**3
            self.KK0[2][5]=self.vel/60*EI*60/self.El[NEl][1]**2
        else:
            self.KK0[2][0]=-self.vel/60*EI*720/self.El[NEl][0]/self.El[NEl][1]**5
            self.KK0[2][1]=-self.vel/60*EI*360/self.El[NEl][0]/self.El[NEl][1]**4
            self.KK0[2][2]=-self.vel/60*EI*60/self.El[NEl][0]/self.El[NEl][1]**3
            self.KK0[2][3]=self.vel/60*EI*720/self.El[NEl][0]/self.El[NEl][1]**5
            self.KK0[2][4]=-self.vel/60*EI*360/self.El[NEl][0]/self.El[NEl][1]**4
            self.KK0[2][5]=self.vel/60*EI*60/self.El[NEl][0]/self.El[NEl][1]**3
        #-----
        self.bb0[0]=self.vel/60*EI*12*self.El[NEl][0]/self.El[NEl][1]**2
        self.bb0[1]=EI*self.b0
#----------------------------
#---------FUNCTION-----------
    def CalcForceLoc_func(self,NEl,ksi):
        self.CalcForceLoc_proc(NEl,ksi)
        return (np.dot(self.KK0,self.ElDisp)+self.bb0)
#----------------------------
#---------FUNCTION-----------
    def CalcLocEl(self,NEl,ksi):
        Kpoint=np.zeros((6,6))
        Vpoint=np.zeros(6)
        self.CalcForceLoc_proc(NEl,ksi)
        for i in range(3):
            Vpoint[i]=-self.bb0[i]
            Vpoint[i+3]=-self.bb0[i]
            for j in range(6):
                Kpoint[i][j]=self.KK0[i][j]
                Kpoint[i+3][j]=self.KK0[i][j]
        return Kpoint, Vpoint
#----------------------------
#---------FUNCTION-----------
    def CalcEl(self,NEl):
        dxw=self.El[NEl][1]/self.Nw
        Kint=np.zeros((6,6))
        Vint=np.zeros(6)
        for i in range(self.Nw):        
            Kpoint,Vpoint=self.CalcLocEl(NEl,(0.5+i)/self.Nw)
            Kint+=Kpoint*dxw
            Vint+=Vpoint*dxw
        Mint=np.zeros((6,6))
        if self.El[NEl][0]==0:
            Mint[0][0]=1/self.El[NEl][1]
            Mint[1][0]=1/2
            Mint[1][1]=1/self.El[NEl][1]
            Mint[2][2]=1
            Mint[3][3]=-1/self.El[NEl][1]
            Mint[4][3]=1/2
            Mint[4][4]=-1/self.El[NEl][1]
            Mint[5][5]=-1
        else:
            Mint[0][0]=math.sin(self.El[NEl][1]*self.El[NEl][0])*self.El[NEl][0]/2/(1-math.cos(self.El[NEl][1]*self.El[NEl][0]))
            Mint[0][2]=self.El[NEl][0]/2
            Mint[1][0]=1/2
            Mint[1][1]=1/self.El[NEl][1]
            Mint[1][2]=1/self.El[NEl][0]/self.El[NEl][1]-math.sin(self.El[NEl][1]*self.El[NEl][0])/2/(1-math.cos(self.El[NEl][1]*self.El[NEl][0]))
            Mint[2][0]=-self.El[NEl][0]/2
            Mint[2][2]=math.sin(self.El[NEl][1]*self.El[NEl][0])*self.El[NEl][0]/2/(1-math.cos(self.El[NEl][1]*self.El[NEl][0]))
            Mint[3][3]=-math.sin(self.El[NEl][1]*self.El[NEl][0])*self.El[NEl][0]/2/(1-math.cos(self.El[NEl][1]*self.El[NEl][0]))
            Mint[3][5]=self.El[NEl][0]/2
            Mint[4][3]=1/2
            Mint[4][4]=-1/self.El[NEl][1]
            Mint[4][5]=math.sin(self.El[NEl][1]*self.El[NEl][0])/2/(1-math.cos(self.El[NEl][1]*self.El[NEl][0]))-1/self.El[NEl][0]/self.El[NEl][1]
            Mint[5][3]=-self.El[NEl][0]/2
            Mint[5][5]=-math.sin(self.El[NEl][1]*self.El[NEl][0])*self.El[NEl][0]/2/(1-math.cos(self.El[NEl][1]*self.El[NEl][0]))
        return np.dot(Mint,Kint), np.dot(Mint,Vint)
#----------------------------
#---------FUNCTION-----------
    def Solve(self, Tol=0.01, DispFile='', BCFile=''):
        err=Tol        
        while err>=Tol:
            Kglob=np.zeros((self.N*3,self.N*3))
            Vglob=np.zeros(self.N*3)
            for i in range(self.N-1):
                KEl,VEl=self.CalcEl(i)
                for row in range(6):
                    Vglob[i*3+row]+=VEl[row]
                    for col in range(6):
                        Kglob[i*3+row][i*3+col]+=KEl[row][col]
#-------------Apply BC----------
#    #Degree: 0 - displacement, 1 - rotation, 2 - curvature
            for Bc in self.BCs:
                if Bc[1]==0:
                    for i in range(3*self.N):
                        Kglob[i][Bc[0]*3]=0
                    Kglob[Bc[0]*3][Bc[0]*3]=-1
                if Bc[1]==1:
                    for i in range(3*self.N):
                        Kglob[i][Bc[0]*3+1]=0
                    Kglob[Bc[0]*3+1][Bc[0]*3+1]=-1
#------Apply Axial Forces and Kapa0
            self.AxialForces[1]=self.Nf0+0.6*abs(self.Forces[1])
            Vglob[2]+=self.AxialForces[1]
            for i in range(2,len(self.BCs)):
                self.AxialForces[i]=self.AxialForces0[i]+self.TransvToAxialForce(self.RollDiam[i],self.Forces[i])
                if self.BCs[i][0]<self.Straight_i: Vglob[3*self.BCs[i][0]+2]+=self.AxialForces[i]
            for i in range(3*self.N):
                Vglob[i]-=Kglob[i][2]*self.Kapa0
                Kglob[i][2]=0
            Kglob[3*self.Straight_i+2][2]=1
#------Solve---------------------
            self.Disp=np.linalg.solve(Kglob, Vglob)
            err=0
            for i in range(len(self.BCs)):
                err+=(1-self.Forces[i]/self.Disp[3*self.BCs[i][0]+self.BCs[i][1]])**2
                self.Forces[i]=self.Disp[3*self.BCs[i][0]+self.BCs[i][1]]
            err=err**0.5
            print('Iteration error: '+str(err))
            print(self.Forces)
#-------------Output------------------
        if DispFile!='':
            k=self.StrBC
            f=open(DispFile,'w')
            f.write('Distance [m]; U [mm]; Angle [rad]; Curvature [1/m]; Curvature rate [1/(m*sek)]; Moment [N*m]; Transverse Force [N]; Axial Force [N]\n')
            for i in range(self.N-1):
                f.write(str(self.El[i][2])+';') #x
                if (i,0) in self.BCs: f.write(str('0.0;'))  #U
                else: f.write(str(self.Disp[i*3]*1000)+';') #U
                if (i,1) in self.BCs: f.write(str('0.0;'))  #Angle
                else: f.write(str(self.Disp[i*3+1])+';')    #Angle
                if i==0: f.write(str(self.Kapa0)+';')       #Curvature
                else: f.write(str(self.Disp[i*3+2])+';')    #Curvature
                kaparate=self.CalcKaparate(i,0)
                f.write(str(kaparate)+';')                  #Curvature rate
                f.write(str(kaparate*self.Stiff(self.El[i][2],kaparate))+';') #Moment
                if i==0:
                    Value=self.CalcForceLoc_func(i,0.5/self.Nw)
                else:
                    Value=(self.CalcForceLoc_func(i-1,1-0.5/self.Nw)+self.CalcForceLoc_func(i,0.5/self.Nw))/2
                f.write(str(Value[0])+';')                
                if i<self.Straight_i:
                    Force=Value[2]
                else:
                    if self.BCs[k][0]==i:
                        Force+=self.AxialForces[k]
                        k+=1
                f.write(str(Force)+'\n')    
            f.close()
        if BCFile!='':
            f=open(BCFile,'w')
            Force=self.Disp[2]
            for i in range(self.StrBC,len(self.BCs)):
                Force+=self.AxialForces[i]
            f.write('Straightening force [N]:'+str(Force)+'\n\n')
            f.write('Num; Distance [m]; Transverse Force [N]; Axial Force [N]; Moment [N*m]\n')
            for i in range(1,len(self.BCs)):                
                if self.BCs[i][0]==self.N-1:                 
                    f.write(str(i-1)+';'+str(self.El[self.N-2][1]+self.El[self.N-2][2])+'; '+str(self.Disp[3*(self.N-1)])+'; '+str(self.AxialForces[i]))
                else:              
                    f.write(str(i-1)+';'+str(self.El[self.BCs[i][0]][2])+'; '+str(self.Disp[3*self.BCs[i][0]])+'; '+str(self.AxialForces[i]))
                    if self.BCs[i][0]==0: f.write('; '+str(self.Disp[3*self.BCs[i][0]+1]))
                f.write('\n')
            f.close()
