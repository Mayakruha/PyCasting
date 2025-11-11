mport numpy as np
from math import exp, log, pi
import vtk
from FEMtoolkit import FEMtoolkit, import_abq, FacesNodes
#-----------------------------------------------
#----------------MOULDS-------------------------
#-----------------------------------------------
def Coating(z):# coating thickness, m
    return 0.0005+0.001*z/0.8
class Mould:
    def __init__(self,Height,Walls):
        self.Height=Height
        self.flux_Tmelt=1115     # melting temperature for mould flux, C
        self.flux_Tliq=1145      # melting temperature for mould flux, C
        self.flux_lamda=1.5      # flux conductivity, W/mK
        self.flux_alfa_liq=5900  # HTC for luquid flux, W/m2*K
        self.flux_alfa_sol=1200  # HTC for solid flux, W/m2*K
        self.walls=Walls
        for BcName in self.walls:
            self.walls[BcName].InPort=self.OutPort
            self.walls[BcName].flux_alfa=self.flux_alfa
        self.Zm=0
    def flux_alfa(self,T):
        if T>self.flux_Tliq: return self.flux_alfa_liq
        elif T<self.flux_Tmelt: return self.flux_alfa_sol
        else: return (T-self.flux_Tmelt)/(self.flux_Tliq-self.flux_Tmelt)*(self.flux_alfa_liq-self.flux_alfa_sol)+self.flux_alfa_sol
    def HeatFlow(self,x,y,Ts,BcName):#W/m2
        if BcName=='BROAD': return self.walls[BcName].HeatFlow(x,Ts)
        else: return self.walls[BcName].HeatFlow(y,Ts)
    def set_level(self,z):
        for BcName in self.walls:
            self.walls[BcName].set_level(z)
    def set_Zm(self,z):
        self.Zm=z
        for BcName in self.walls:
            self.walls[BcName].Zm=z
    def InPort(self): #should be replaced by OutPort from model
        return 0, 0, 0
    def OutPort(self):
        MouldLevel, v, Tsol = self.InPort()
        return MouldLevel, v, Tsol, self.Height, self.flux_Tmelt, self.flux_lamda
    def Update(self):
        for BcName in self.walls:
            print('\n** '+BcName+': Heat transfer parameters')
            self.walls[BcName].Update()
class Wall:
    def __init__(self, CoatFunc):                    
        self.CoatFunc=CoatFunc
    def InPort(self): #should be replaced by function outport from mould object
        return 0, 0, 0, 0, 0, 0, 0
    def flux_alfa(self,T): #should be replaced by function of flux htc function
        return 0 
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    # q - heat flux, W
    def HeatUp(self,q):
        self.Twat+=self.Geom_f*q/self.Qwat*60/(4230-3.6562*self.Twat-0.02585*self.Twat**2)
    def Update(self):
        self.Zm, self.v, self.Tsol, self.Height, self.flux_Tmelt, self.flux_lamda = self.InPort()
        self.Twat=self.Twat_in                                       #current water temeprature
        lamda_wat=0.55748+0.0021525*self.Twat-0.0000097*self.Twat**2 #W/m*K
        visc_wat=1.53555258E-06*exp(-0.036*self.Twat)+2.52805091E-07 #m2/sec
        Pr=self.Prandtl(self.Twat)
        v_wat=self.Qwat/60000/self.Wat_sec #m/sec
        Re=v_wat*self.deff/visc_wat
        self.alfa_wat0=0.023*lamda_wat/self.deff*Re**0.8*Pr**0.4*self.geom_coeff(Pr)
        self.set_level(self.Zm)
        self.flux_thick_m=self.flux_lamda*((self.flux_Tmelt-self.Twat)/self.flux_alfa(self.Tsol)/self.v**0.8/(self.Tsol-self.flux_Tmelt)-self.Rm-1/self.alfa_watz)
        if self.flux_thick_m<0.0: self.flux_thick_m=0.0
        print('Prandtl: '+str(Pr))
        print('Water speed, m/sec: '+str(v_wat))
        print('Reynolds: '+str(Re))
        print('Nominal HTC, kW/(m2K): '+str(self.alfa_wat0/1000))
        print('Solid flux thickness at meniscus, mm: '+str(self.flux_thick_m*1000))
    def geom_coeff(self,Pr):
        return 1
    def set_level(self,z):
        self.Coat_thick=self.CoatFunc(z)
        self.Rm=self.Mould_thick/self.Mould_lamda+self.Coat_thick/self.Coat_lamda
        self.z=z
        self.alfa_watz=self.alfa_wat0*(1+self.deff/z/2)
        if self.z>self.Zm:
            self.shrink=self.Thb*0.0085*(z-self.Zm)**0.5
class Plate(Wall):
    def __init__(self,Nch,Sch,Pch,Mould_thick,Width,Thick,CoatFunc):
        #Nch - number of cooling channels
        #Sch - cross section of a cooling channel, m2
        #Pch - perimeter of a cooling channel, m
        self.Twat_in=15              # Temperature of water [Celsius]
        self.Qwat=0                  # Water flow [l/min]
        self.Coat_lamda=80           # coating conductivity, W/mK
        self.Mould_lamda=370         # mould material conductivity, W/mK
        self.Taper=0.001             # taper of a side over height, m
        self.Geom_f=2                # coefficient for calculation of water heating  (1 - actual heat; 2-half of heat from thermal model etc)
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Width=Width             # considered width of the plate, m
        self.Thb=Thick               # size for calculation of shrinkage
        self.deff=4*Sch/Pch          # effective diameter, m
        self.Wat_sec=Nch*Sch
        self.CoatFunc=CoatFunc
    def mould_shape(self,x,z):#m
        return self.Taper*z/self.Height
    def HeatFlow(self,x,Ts):#W/m2
        if self.z==self.Zm:
            if Ts<=self.Tsol:
                self.flux_thickz=self.flux_thick_m
            else:
                self.flux_thickz=self.flux_thick_m*(1+2*(Ts-self.Tsol)/self.Tsol)
        elif self.z>self.Zm:
            self.flux_thickz=self.flux_thick_m+self.shrink-self.mould_shape(x,self.z)+self.mould_shape(x,self.Zm)
            if self.flux_thick_m+self.shrink-self.mould_shape(0,self.z)+self.mould_shape(0,self.Zm)<0.0:self.flux_thickz=self.mould_shape(0,self.z)-self.mould_shape(x,self.z)
        Q=(Ts-self.Twat)/(1/self.alfa_watz+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa(Ts)/self.v**0.8)
        self.TempW=self.Twat+Q/self.alfa_watz
        self.alfa_wat=self.alfa_watz*(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
        return (Ts-self.Twat)/(1/self.alfa_wat+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa(Ts)/self.v**0.8) #mould 
class Plate_faska(Plate):
    def __init__(self,Nch,Sch,Pch,Mould_thick,Width,Thick,CoatFunc):
        #Nch - number of cooling channels
        #Sch - cross section of a cooling channel, m2
        #Pch - perimeter of a cooling channel, m
        self.Twat_in=15              # Temperature of water [Celsius]
        self.Qwat=0                  # Water flow [l/min]
        self.Coat_lamda=80           # coating conductivity, W/mK
        self.Mould_lamda=370         # mould material conductivity, W/mK
        self.Taper=0.001             # taper of a side over height, m
        self.Geom_f=2                # coefficient for calculation of water heating  (1 - actual heat; 2-half of heat from thermal model etc)
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Width=Width             # considered width of the plate, m
        self.Thb=Thick               # size for calculation of shrinkage
        self.deff=4*Sch/Pch          # effective diameter, m
        self.Wat_sec=Nch*Sch
        # geometry of a faska
        self.faska_start=0.15 #m
        self.faska_width=0.0 #m
        self.faska_fi=0.0 #rad
        self.faska_depth=0.004 #m
        self.CoatFunc=CoatFunc    
        self.ms=[1.36,1.62,1.75,1.78,1.75,1.67,1.54,1.35,1.15,0.94,0.72,0.49,0.25,0]
    def mould_shape(self,x,z):#m
        if z<0.25:
            dd=-0.00064+0.008*z
        else:
            i=int((z-0.25)/0.05)
            if i>12:
                dd=0
            else:
                dd=(self.ms[i]+(self.ms[i+1]-self.ms[i])*(z-i*0.05-0.25)/0.05)/1000
        if z>self.faska_start:
            df=(self.Width-x-self.faska_width*(z-self.faska_start)/(self.Height-self.faska_start))*tan(self.faska_fi)
            if df<-self.faska_depth:
                dd=-self.faska_depth
            elif df<dd:
                dd=df
        return self.Taper*z/self.Height+dd
class Circle_tube(Wall):
    #---------- Circle ---------------------------
    # D_billet    billet diameter, m
    # Mould_thick distance between water and mould surface, m
    # Wat_thick   thickness of water layer, m    
    def __init__(self,D_billet,Mould_thick,Wat_thick,Port):
        self.R=D_billet/2            # radius of billet cross section, m
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Wat_thick=Wat_thick     #thickness of water layer, m
        self.deff=2*Wat_thick   # effective diameter, m
        self.Wat_sec=pi*Wat_thick*(2*(D_billet/2+Mould_thick)+Wat_thick)
        self.Port=Port
    def geom_coeff(self, Pr):
        return (1-0.45/(2.4+Pr))*(1+self.Wat_thick/(self.R+self.Mould_thick))**(0.16/Pr**0.15)
    def HeatFlow(self,Ts):#W/m2
        self.alfa_flux=v**0.8*(self.alfa_flux_max*(1-self.Rate_flux_drop/2*(1-tanh(self.Rate_flux_melt*(Ts-self.flux_Tmelt))))+self.Press_coeff*P)      
        self.alfa_wat=self.alfa_wat0*(1+self.Wat_thick/self.z)
        Q=(Ts-self.Twat)/self.R/(1/self.alfa_wat/(self.R+self.Mould_thick+self.Coat_thick)+\
            log(1+self.Mould_thick/(self.R+self.Coat_thick))/self.Mould_lamda+\
            log(1+self.Coat_thick/self.R)/self.Coat_lamda+\
            1/self.alfa_flux/self.R)
        self.TempW=self.Twat+Q*self.R/(self.R+self.Mould_thick)/self.alfa_wat
        self.alfa_wat*=(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
        return (Ts-self.Twat)/self.R/(1/self.alfa_wat/(self.R+self.Mould_thick+self.Coat_thick)+\
                log(1+self.Mould_thick/(self.R+self.Coat_thick))/self.Mould_lamda+\
                log(1+self.Coat_thick/self.R)/self.Coat_lamda+\
                1/self.alfa_flux/self.R) #mould 
class Rect_tube(Wall):
    def __init__(self,W_billet,Th_billet,Mould_thick,Wat_thick,Port):
        self.Wdb=W_billet            # width of billet cross section, m
        self.Thb=Th_billet           # thickness of billet cross section, m
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Wat_thick=Wat_thick     # thickness of water layer, m
        self.deff=2*Wat_thick        # effective diameter, m
        self.Wat_sec=2*Wat_thick*(W_billet+Th_billet+4*Mould_thick+2*Wat_thick)
        # geometry of a faska
        self.faska_start=0.4  #m
        self.faska_width=0.03 #m
        self.faska_depth=0.01 #m
        self.Port=Port
    def mould_shape(self,x,y,z,BcName):#m
        if BcName=='NARROW':
            if z>self.faska_start:
                yf=self.Thb-(z-self.faska_start)/(self.Height-self.faska_start)*self.faska_width
            else:
                yf=self.Thb
            if y>yf:
                return 0.0055*self.Wdb*z/self.Height-(y-yf)/(self.Thb-yf)*(z-self.faska_start)/(self.Height-self.faska_start)*self.faska_depth
            else:
                return 0.0055*self.Wdb*z/self.Height
        elif BcName=='BROAD':
            return 0.0055*self.Thb*z/self.Height
    def HeatFlow(self,x,y,Ts,BcName):#W/m2
        if self.z==self.Zm:
            if Ts<=self.Tsol:
                self.flux_thickz=self.flux_thick_m
            else:
                self.flux_thickz=self.flux_thick_m*(1+2*(Ts-self.Tsol)/self.Tsol)
        elif self.z>self.Zm:
            if BcName=='NARROW':
                self.flux_thickz=self.flux_thick_m+self.shrinkage(self.z)-self.mould_shape(x,y,self.z,BcName)+self.mould_shape(x,y,self.Zm,BcName)
            elif BcName=='BROAD':
                self.flux_thickz=self.flux_thick_m+self.shrinkage(self.z)-self.mould_shape(x,y,self.z,BcName)+self.mould_shape(x,y,self.Zm,BcName)
            if self.flux_thickz<0.0:self.flux_thickz=0.0
        Q=(Ts-self.Twat)/(1/self.alfa_watz+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa/self.v**0.8)
        self.TempW=self.Twat+Q/self.alfa_watz
        self.alfa_wat=self.alfa_watz*(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
        return (Ts-self.Twat)/(1/self.alfa_wat+self.Rm+self.flux_thickz/self.flux_lamda+1/self.flux_alfa/self.v**0.8) #mould 
#-----------------------------------------------
#----------------BILLET-------------------------
#-----------------------------------------------
class Quasi3DTemp:
    def __init__(self,FileName,Mould):
        self.Z_end=0.9  #end of calculation, m        
        #---------- casting parameters---------------------------
        self.v=3.2      #casting speed, m/min
        self.dTemp=10   #temperature of overheating for liquid steel, K
        self.Tsr=15     #ambient temperature, Celsius
        self.MouldLevel=0.1  #mould level,m
        #---------- heat exchange parameters
        self.mould=Mould
        self.mould.InPort=self.OutPort
        self.Zones=[[self.mould.Height,76.74],[1.388,27.54],[3.531,20.82],[5.366,12.54],[8.788,0.3]] #(start of zone, m; water flux, l/m2*min)
        #-----------steel properties----------------------
        self.lamda=29       #steel conductivity, W/m*K    
#        self.alfa_liq=2200  #htc for a border between solid and liquid steel, W/m2K
        self.lamda_liq=190  #conductivity in the area of luquid steel, W/m*K !!!! It doesn't work
        self.L=272E+3       #heat of solidification, J/kg
        self.ro_sol=7300    #solid steel density, kg/m3
        self.ro_liq=7000    #liquid steel density, kg/m3
        self.Cl=795         #heat capacity for liquid steel, J/kg*K
        self.Cr=680         #heat capacity for solid steel, J/kg*K
        self.k=0.2          #power for strain rate in the steel creep equation
        self.Tc=250         #thermal coefficient in the steel creep equation, K
        self.Ac=73000       #coefficient for stress in the steel creep equation, MPa
        #-------------------------------------
        self.wc=(0.125,0.375)       #relative coordinates for flux integretaion
        self.Epsilon=0.001          #tolerance for extraction of temperature from H
        self.KapaN=20               #Number of columns in the output of billet stiffness
        self.PointScreen={'BROAD':12061,'NARROW':201} #Node and BC for data printing
        self.fem=import_abq(FileName)
        self.AREA=np.zeros((self.fem.MaxNodeNum+1))
        self.T=np.zeros(1)
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0,Mn=0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S) #Solidus temperature, Celsius
        self.Tlik=1536-(78*C+7.5*Si+4.9*Mn+1.3*Cr+3.1*Ni+4.7*Cu+3.6*Al+34.4*P+38*S)  #Liquidus temperature, Celsius
        self.Hl=(self.ro_liq*self.Cl+self.ro_sol*self.Cr)*(self.Tlik-self.Tsol)/2+(self.ro_liq+self.ro_sol)*self.L/2 #J/m3
        print('\n***Thermal Properties of steel:')
        print('Solidus temperature, Celsius: '+str(self.Tsol))
        print('Liquidus temperature, Celsius: '+str(self.Tlik))
    def Cef(self,Temp):
        if Temp<self.Tsol: return self.Cr
        elif Temp>self.Tlik: return self.Cl
        else: return (self.L+self.Cr*(self.Tlik-Temp)+self.Cl*(Temp-self.Tsol))/(self.Tlik-self.Tsol)
    def FuncTemp(self,Value):
        if Value<self.Tsol: return (Value-self.Tsol)*self.ro_sol*self.Cr
        elif Value>self.Tlik: return (Value-self.Tlik)*self.ro_liq*self.Cl+self.Hl
        else: return (Value-self.Tsol)/(self.Tlik-self.Tsol)*self.Hl
    def Temperature(self,Value):
        if Value<0: return Value/self.ro_sol/self.Cr+self.Tsol
        elif Value>self.Hl: return (Value-self.Hl)/self.ro_liq/self.Cl+self.Tlik
        else: return self.Tsol+(self.Tlik-self.Tsol)*Value/self.Hl
    def Stiffness(self,KapaRate):
        Value=0
        for Node in range(1,self.fem.MaxNodeNum+1):
            if self.T[Node]<self.Tsol:
                Value+=self.Ac*1E+6*(KapaRate*self.fem.Coord[Node][1])**self.k*np.exp(-(self.T[Node]+273)/self.Tc)*self.AREA[Node]*self.fem.Coord[Node][1] # H*m
        return Value/KapaRate # H*m2*sek
    def Press(self,z):
        return 9.81*self.ro_liq*(z-self.MouldLevel)
    def HeatFlow(self,x,y,z,Ts,Name):#W/m2
        if z<=self.Zones[0][0]: return self.mould.HeatFlow(x,y,Ts,Name) #mould 
        else:
            iz=0
            while iz<len(self.Zones)-1 and self.Zones[iz+1][0]<z:iz+=1
            return 142*(Ts-self.Tsr)*self.Zones[iz][1]**0.55*exp(-0.0012*(Ts+273))+5.670367E-8*((Ts+273)**4-(self.Tsr+273)**4) #spray
    def dxydksi(self,nu,xy):
        return -(1-nu)*xy[0]+(1-nu)*xy[1]+nu*xy[2]-nu*xy[3]
    def dxydnu(self,ksi,xy):
        return -(1-ksi)*xy[0]-ksi*xy[1]+ksi*xy[2]+(1-ksi)*xy[3]
    def OutPort(self):
        return self.MouldLevel, self.v, self.Tsol
#-------------------------------------------------------------------------------
# kj - convergence coefficient (less coefficient, less time step)
# Meshdz - distance between levels for output, [m]
# VTUFileName - name vtu-file for output. There is no output by default
# VTUHFlowFile - name vtu-file for heat flow output. There is no output by default
# MI - name of a set of files for 3D mechanical calculations
# KapaRate - maximum curvature rate for output of section moment inertia. There is no output by default, [1/m*sec]
# Radius - for visualization, if Radius equals to zero the mesh isn't curved, [m]
    def RunCalc(self, kj=0.5, Meshdz=0.01, VTUFileName='', VTUHFlowFile='', MI='', KapaRate=0, Radius=0):
        HeatFlow_level={}
        face=FacesNodes[self.fem.Eltype[1]][0]
        DistMin=np.linalg.norm(self.fem.Coord[self.fem.Elems[1][face[0]]]-self.fem.Coord[self.fem.Elems[1][face[1]]])        
        XEl=np.zeros(4)
        YEl=np.zeros(4)
        for i in range(4):
            XEl[i]=self.fem.Coord[self.fem.Elems[1][i]][0]
            YEl[i]=self.fem.Coord[self.fem.Elems[1][i]][1]        
        if XEl[0]*(YEl[1]-YEl[3])+XEl[1]*(YEl[3]-YEl[0])+XEl[3]*(YEl[0]-YEl[1])>0:
            IndxList=(0,1,2,3)
            IndxList_f=(0,1,2,3)
        else:
            IndxList=(0,3,2,1)
            IndxList_f=(3,2,1,0)
        self.AREA.fill(0)
        for El in range(1,self.fem.MaxElemNum+1):
            for i in range(4):
                XEl[i]=self.fem.Coord[self.fem.Elems[El][i]][0]
                YEl[i]=self.fem.Coord[self.fem.Elems[El][i]][1]
            self.AREA[self.fem.Elems[El][0]]+=abs(2*XEl[1]*YEl[3]-2*XEl[3]*YEl[1]+(3*XEl[0]-XEl[2])*(YEl[1]-YEl[3])+(3*YEl[0]-YEl[2])*(XEl[3]-XEl[1]))/16
            self.AREA[self.fem.Elems[El][1]]+=abs(2*XEl[2]*YEl[0]-2*XEl[0]*YEl[2]+(3*XEl[1]-XEl[3])*(YEl[2]-YEl[0])+(3*YEl[1]-YEl[3])*(XEl[0]-XEl[2]))/16
            self.AREA[self.fem.Elems[El][2]]+=abs(2*XEl[3]*YEl[1]-2*XEl[1]*YEl[3]+(3*XEl[2]-XEl[0])*(YEl[3]-YEl[1])+(3*YEl[2]-YEl[0])*(XEl[1]-XEl[3]))/16
            self.AREA[self.fem.Elems[El][3]]+=abs(2*XEl[0]*YEl[2]-2*XEl[2]*YEl[0]+(3*XEl[3]-XEl[1])*(YEl[0]-YEl[2])+(3*YEl[3]-YEl[1])*(XEl[2]-XEl[0]))/16
            for face in FacesNodes[self.fem.Eltype[El]]:
                dist=np.linalg.norm(self.fem.Coord[self.fem.Elems[El][face[0]]]-self.fem.Coord[self.fem.Elems[El][face[1]]])
                if DistMin>dist: DistMin=dist
        dtau=kj*DistMin**2*min(self.ro_liq*self.Cl,self.ro_sol*self.Cr)/4/self.lamda  #sek
        dZ=self.v*dtau/60   #m
        Count=int(Meshdz/dZ)
        if Count==0:Count=1
        ElRow=int((self.Z_end-self.MouldLevel)/(dZ*Count))
    #------------VTUFileName--------------------
        if VTUFileName!='':
            Points=vtk.vtkPoints()
            mesh=vtk.vtkUnstructuredGrid()
            mesh.Allocate(ElRow*self.fem.MaxElemNum)
            for i in range(ElRow+1):
                for Node in range(1,self.fem.MaxNodeNum+1):
                    if Radius!=0:
                        if i*dZ*Count<Radius*np.pi/2:
                            Points.InsertNextPoint(self.fem.Coord[Node][0],(self.fem.Coord[Node][1]+Radius)*np.cos(i*dZ*Count/Radius),(self.fem.Coord[Node][1]+Radius)*np.sin(i*dZ*Count/Radius))
                        else:
                            Points.InsertNextPoint(self.fem.Coord[Node][0],Radius*np.pi/2-i*dZ*Count,self.fem.Coord[Node][1]+Radius)
                    else:
                        Points.InsertNextPoint(self.fem.Coord[Node][0],self.fem.Coord[Node][1],i*dZ*Count+self.MouldLevel)
            for i in range(ElRow):
                for El in range(1,self.fem.MaxElemNum+1):
                    mesh.InsertNextCell(vtk.VTK_HEXAHEDRON,8,(self.fem.Elems[El][0]+i*self.fem.MaxNodeNum-1,self.fem.Elems[El][1]+i*self.fem.MaxNodeNum-1,
self.fem.Elems[El][2]+i*self.fem.MaxNodeNum-1,self.fem.Elems[El][3]+i*self.fem.MaxNodeNum-1,
self.fem.Elems[El][0]+(i+1)*self.fem.MaxNodeNum-1,self.fem.Elems[El][1]+(i+1)*self.fem.MaxNodeNum-1,
self.fem.Elems[El][2]+(i+1)*self.fem.MaxNodeNum-1,self.fem.Elems[El][3]+(i+1)*self.fem.MaxNodeNum-1))
            mesh.SetPoints(Points)
            Res=vtk.vtkFloatArray()
            Res.SetName('Temp')
            Res.SetNumberOfValues((ElRow+1)*self.fem.MaxNodeNum)
            Res.Fill(self.Tlik+self.dTemp)
    #------------VTUHFlowFile--------------------
        if VTUHFlowFile!='':
            QData_mesh={}
            QData_points={}
            QSurf={}
            SurfNodes={}
            SurfNodeList={}
            SurfNodeNum={}
            for BcName in self.fem.Surfs:
                QData_points[BcName]=vtk.vtkPoints()
                QData_mesh[BcName]=vtk.vtkUnstructuredGrid()
                SurfNodes[BcName]=np.full(self.fem.MaxNodeNum+1,self.fem.MaxNodeNum+1,dtype=np.uint32)
                SurfNodeList[BcName]=[]
                SurfElemNum=0
                SurfNodeNum[BcName]=0
                for SurfEl in self.fem.Surfs[BcName]:
                    SurfElemNum+=len(self.fem.ESets[SurfEl[0]])
                    for El in self.fem.ESets[SurfEl[0]]:
                        for i in range(2):
                            Node=self.fem.Elems[El][FacesNodes[self.fem.Eltype[El]][SurfEl[1]][i]]
                            if SurfNodes[BcName][Node]==self.fem.MaxNodeNum+1:
                                SurfNodes[BcName][Node]=SurfNodeNum[BcName]
                                SurfNodeList[BcName].append(Node)
                                SurfNodeNum[BcName]+=1
                for i in range(ElRow+1):
                    for Node in SurfNodeList[BcName]:
                        QData_points[BcName].InsertNextPoint(self.fem.Coord[Node][0],self.fem.Coord[Node][1],i*dZ*Count+self.MouldLevel)
                QData_mesh[BcName].Allocate(ElRow*2*SurfElemNum)
                for i in range(ElRow):
                    for SurfEl in self.fem.Surfs[BcName]:
                        for El in self.fem.ESets[SurfEl[0]]:
                            Node0=self.fem.Elems[El][FacesNodes[self.fem.Eltype[El]][SurfEl[1]][0]]
                            Node1=self.fem.Elems[El][FacesNodes[self.fem.Eltype[El]][SurfEl[1]][1]]
                            QData_mesh[BcName].InsertNextCell(vtk.VTK_TRIANGLE,3,(i*SurfNodeNum[BcName]+SurfNodes[BcName][Node1],i*SurfNodeNum[BcName]+SurfNodes[BcName][Node0],(i+1)*SurfNodeNum[BcName]+SurfNodes[BcName][Node0]))
                            QData_mesh[BcName].InsertNextCell(vtk.VTK_TRIANGLE,3,(i*SurfNodeNum[BcName]+SurfNodes[BcName][Node1],(i+1)*SurfNodeNum[BcName]+SurfNodes[BcName][Node0],(i+1)*SurfNodeNum[BcName]+SurfNodes[BcName][Node1]))
                QData_mesh[BcName].SetPoints(QData_points[BcName])
                QSurf[BcName]=vtk.vtkFloatArray()
                QSurf[BcName].SetName('Q, [MW/m2]')
                QSurf[BcName].SetNumberOfValues((ElRow+1)*SurfNodeNum[BcName])
    #------------MI--------------------------
        if MI!='':
            MI_model=FEMtoolkit()
            MI_model.MaxNodeNum=self.fem.MaxNodeNum*(ElRow+1)
            MI_model.MaxElemNum=self.fem.MaxElemNum*ElRow
            MI_model.Coord=np.full(MI_model.MaxNodeNum+1,None)
            MI_model.Elems=np.ones(MI_model.MaxElemNum+1,dtype=tuple)
            MI_model.Eltype=np.zeros(MI_model.MaxElemNum+1,dtype=np.int8)
            MI_model.TypeList[6]='C3D8'
            for NsetName in self.fem.NSets:
                MI_model.NSets[NsetName]=[]            
            MI_model.Surfs['InternalSurf']=[]
            MI_model.NodeValue['Temp']={}
            MI_model.FaceLoad['P']={}
            MInode_count=0
            MInodes_0=np.full(self.fem.MaxNodeNum+1,MI_model.MaxNodeNum+1,dtype=np.uint32)
            MInodes_1=np.full(self.fem.MaxNodeNum+1,MI_model.MaxNodeNum+1,dtype=np.uint32)
    #--------------KapaRate------------------
        if KapaRate!=0:
            f=open('BilletStiffness.csv','w')
            for i in range(1,self.KapaN+1):
                f.write(';'+str(i*KapaRate/self.KapaN))
            f.write('\n')        
    #--------------------------------
        ElTemp=np.zeros(4)
        wN=len(self.wc)
        self.T=np.full(self.fem.MaxNodeNum+1,self.Tlik+self.dTemp) #initial temperature, K
        H=np.full(self.fem.MaxNodeNum+1,self.FuncTemp(self.Tlik+self.dTemp))
        Z=self.MouldLevel
        KK=np.zeros((self.fem.MaxElemNum+1,4,4))
        GG=np.zeros((self.fem.MaxElemNum+1,4))
        for El in range(1,self.fem.MaxElemNum+1):
            for i in range(4):
                XEl[i]=self.fem.Coord[self.fem.Elems[El][i]][0]
                YEl[i]=self.fem.Coord[self.fem.Elems[El][i]][1]
            L1=((XEl[0]+XEl[1]-XEl[2]-XEl[3])**2+(YEl[0]+YEl[1]-YEl[2]-YEl[3])**2)**0.5/2
            L2=((XEl[0]+XEl[3]-XEl[1]-XEl[2])**2+(YEl[0]+YEl[3]-YEl[1]-YEl[2])**2)**0.5/2
            GG[El][0]=((XEl[1]-XEl[0])**2+(YEl[1]-YEl[0])**2)**0.5/2 #m
            GG[El][1]=((XEl[2]-XEl[1])**2+(YEl[2]-YEl[1])**2)**0.5/2 #m
            GG[El][2]=((XEl[3]-XEl[2])**2+(YEl[3]-YEl[2])**2)**0.5/2 #m
            GG[El][3]=((XEl[0]-XEl[3])**2+(YEl[0]-YEl[3])**2)**0.5/2 #m
            # n1: ksi=0.5            
            for nu in self.wc:
                KK[El][0]+=dtau*L1/2*np.array((nu-1,1-nu,nu,-nu))/((self.dxydksi(nu,XEl))**2+(self.dxydksi(nu,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][0]]/wN   #sek/m2
                KK[El][1]-=dtau*L1/2*np.array((nu-1,1-nu,nu,-nu))/((self.dxydksi(nu,XEl))**2+(self.dxydksi(nu,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][1]]/wN   #sek/m2
                KK[El][3]+=dtau*L1/2*np.array((nu-0.5,0.5-nu,nu+0.5,-nu-0.5))/((self.dxydksi(nu+0.5,XEl))**2+(self.dxydksi(nu+0.5,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][3]]/wN   #sek/m2
                KK[El][2]-=dtau*L1/2*np.array((nu-0.5,0.5-nu,nu+0.5,-nu-0.5))/((self.dxydksi(nu+0.5,XEl))**2+(self.dxydksi(nu+0.5,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][2]]/wN   #sek/m2
            # n2: nu=0.5            
            for ksi in self.wc:
                KK[El][0]+=dtau*L2/2*np.array((ksi-1,-ksi,ksi,1-ksi))/((self.dxydnu(ksi,XEl))**2+(self.dxydnu(ksi,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][0]]/wN   #sek/m2
                KK[El][3]-=dtau*L2/2*np.array((ksi-1,-ksi,ksi,1-ksi))/((self.dxydnu(ksi,XEl))**2+(self.dxydnu(ksi,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][3]]/wN   #sek/m2
                KK[El][1]+=dtau*L2/2*np.array((ksi-0.5,-ksi-0.5,ksi+0.5,0.5-ksi))/((self.dxydnu(ksi+0.5,XEl))**2+(self.dxydnu(ksi+0.5,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][1]]/wN   #sek/m2
                KK[El][2]-=dtau*L2/2*np.array((ksi-0.5,-ksi-0.5,ksi+0.5,0.5-ksi))/((self.dxydnu(ksi+0.5,XEl))**2+(self.dxydnu(ksi+0.5,YEl))**2)**0.5/self.AREA[self.fem.Elems[El][2]]/wN   #sek/m2
        j=Count-1
        LevelNum=0
        self.mould.Zm=Z
        self.mould.Update()
        print('-------------------------------------')
        print('The calculation uses dZ='+str(dZ)+' m')
        print('Number of steps:'+str(ElRow*Count))
        print('-------------------------------------')
        print(str(ElRow+1)+' levels:')
        print('Level: '+str(self.MouldLevel)+' m')
        ColumnNames=' Axis,m|'
        msgln=' '
        scrln=' {:5.2f}| {:6.1f}| {:6.1f}|{:7.3f}|'
        for BcName in self.mould.walls:
            msgln+=BcName+' | '
            ColumnNames+='Twat,C|Q,kW/m2| Temp,C|Flux,mm|'
        print(msgln)
        print('_________________________________________________________________')
        print(ColumnNames)
    #--------------------------------
        while Z<self.Z_end:
            if Z<=self.Zones[0][0]:
                self.mould.set_level(Z)
                if Z==self.mould.Zm:
                    Flag=False #for check of start of solidification
                    for BcName in self.fem.Surfs:
                        for SurfEl in self.fem.Surfs[BcName]:
                            for El in self.fem.ESets[SurfEl[0]]:
                                for i in range(2):
                                    Node_indx=FacesNodes[self.fem.Eltype[El]][SurfEl[1]][i]
                                    Node=self.fem.Elems[El][Node_indx]
                                    if self.T[Node]>self.Tsol:Flag=True                  
    #--------------------Printing & Output-------------------------------
            j+=1
            if j==Count:
                msgln=' {:6.3f}|'.format(Z)
                for BcName in self.PointScreen:
                    Ts=self.T[self.PointScreen[BcName]]
                    Q=self.HeatFlow(self.fem.Coord[self.PointScreen[BcName]][0],self.fem.Coord[self.PointScreen[BcName]][1],Z,Ts,BcName)
                    msgln+=scrln.format(self.mould.walls[BcName].Twat,Q/1000,Ts,self.mould.walls[BcName].flux_thickz*1000)
                print(msgln)
                if VTUFileName!='':
                    for Node in range(1,self.fem.MaxNodeNum+1): Res.SetValue(self.fem.MaxNodeNum*LevelNum+Node-1,self.T[Node])
                if VTUHFlowFile!='':
                    for BcName in self.fem.Surfs:
                        for Node in SurfNodeList[BcName]:
                            Q=self.HeatFlow(self.fem.Coord[Node][0],self.fem.Coord[Node][1],Z,self.T[Node],BcName)                          
                            QSurf[BcName].SetValue(SurfNodeNum[BcName]*LevelNum+SurfNodes[BcName][Node],Q/1000000)
#                            QSurf[BcName].SetValue(SurfNodeNum[BcName]*LevelNum+SurfNodes[BcName][Node],self.mould.walls[BcName].flux_thickz*1000) #to extract gap
                if MI!='':
                    if LevelNum<ElRow:
                        PFaces={} # the first key - min mumber of nodes; the second key - max mumber of nodes
                                  # list (0-border/1-solid/2-liquid, face index, element number)
                        for El in range(1,self.fem.MaxElemNum+1):
                            MIFlag=True #if element solid
                            GlobEl=LevelNum*self.fem.MaxElemNum+El
                            for i in range(4):
                                if self.T[self.fem.Elems[El][i]]>self.Tsol: MIFlag=False
                            for Face_i in range(len(FacesNodes[self.fem.Eltype[El]])):
                                NodeSet=set()
                                for i in FacesNodes[self.fem.Eltype[El]][Face_i]: NodeSet.add(self.fem.Elems[El][i])
                                minNode=min(NodeSet)
                                maxNode=max(NodeSet)
                                if minNode in PFaces:
                                    if maxNode in PFaces[minNode]:
                                        if MIFlag and PFaces[minNode][maxNode][0]==2:
                                            PFaces[minNode][maxNode]=[0,Face_i,GlobEl]
                                        elif not MIFlag and PFaces[minNode][maxNode][0]==1:
                                            PFaces[minNode][maxNode][0]=0
                                    else:
                                        if MIFlag: PFaces[minNode][maxNode]=[1,Face_i,GlobEl]
                                        else: PFaces[minNode][maxNode]=[2,Face_i,GlobEl]
                                else:
                                    PFaces[minNode]={}
                                    PFaces[minNode][maxNode]=[]
                                    if MIFlag: PFaces[minNode][maxNode]=[1,Face_i,GlobEl]
                                    else: PFaces[minNode][maxNode]=[2,Face_i,GlobEl]
                            if MIFlag:
                                Nodes=[]
                                for i in IndxList:
                                    XEl[i]=self.fem.Coord[self.fem.Elems[El][i]][0]
                                    YEl[i]=self.fem.Coord[self.fem.Elems[El][i]][1]
                                    if MInodes_0[self.fem.Elems[El][i]]==MI_model.MaxNodeNum+1:
                                        MInode_count+=1
                                        MInodes_0[self.fem.Elems[El][i]]=MInode_count
                                        MI_model.Coord[MInode_count]=np.array((XEl[i],YEl[i],Z))
                                    Nodes.append(MInodes_0[self.fem.Elems[El][i]])
                                for i in IndxList:
                                    if MInodes_1[self.fem.Elems[El][i]]==MI_model.MaxNodeNum+1:
                                        MInode_count+=1
                                        MInodes_1[self.fem.Elems[El][i]]=MInode_count
                                        MI_model.Coord[MInode_count]=np.array((XEl[i],YEl[i],Z+dZ*Count))
                                    Nodes.append(MInodes_1[self.fem.Elems[El][i]])
                                MI_model.Elems[GlobEl]=tuple(Nodes)
                                MI_model.Eltype[GlobEl]=6
                        for minNode in PFaces:
                            for maxNode in PFaces[minNode]:
                                if PFaces[minNode][maxNode][0]==0:
                                    if not 'Press_'+str(IndxList_f[PFaces[minNode][maxNode][1]]+3) in MI_model.ESets: MI_model.ESets['Press_'+str(IndxList_f[PFaces[minNode][maxNode][1]]+3)]=[]
                                    MI_model.ESets['Press_'+str(IndxList_f[PFaces[minNode][maxNode][1]]+3)].append(PFaces[minNode][maxNode][2])
                                    if PFaces[minNode][maxNode][2] not in MI_model.FaceLoad['P']: MI_model.FaceLoad['P'][PFaces[minNode][maxNode][2]]=[]
                                    MI_model.FaceLoad['P'][PFaces[minNode][maxNode][2]].append([IndxList_f[PFaces[minNode][maxNode][1]]+2,self.Press(Z+dZ*Count/2)])
                    for Node in range(1,self.fem.MaxNodeNum+1):
                        if self.T[Node]<=self.Tsol: MI_model.NodeValue['Temp'][MInodes_0[Node]]=self.T[Node]
                    for NsetName in self.fem.NSets:
                        for Node in self.fem.NSets[NsetName]:
                            if MInodes_0[Node]<MI_model.MaxNodeNum+1: MI_model.NSets[NsetName].append(MInodes_0[Node])
                    MInodes_0=MInodes_1.copy()
                    MInodes_1.fill(MI_model.MaxNodeNum+1)
                if KapaRate!=0:
                    f.write(str(Z))
                    for i in range(1,self.KapaN+1):
                        f.write(';'+str(self.Stiffness(i*KapaRate/self.KapaN)))
                    f.write('\n')
                LevelNum+=1
                j=0
    #----------H calculation--------------------------------------------
            for BcName in self.fem.Surfs:
                HeatFlow_level[BcName]=0
            for El in range(1,self.fem.MaxElemNum+1):
                for i in range(4): ElTemp[i]=self.T[self.fem.Elems[El][i]]
                for i in range(4): H[self.fem.Elems[El][i]]+=np.dot(KK[El][i],ElTemp)*self.lamda  #J/m3
            for BcName in self.fem.Surfs:
                for SurfEl in self.fem.Surfs[BcName]:
                    for El in self.fem.ESets[SurfEl[0]]:
                        for i in range(2):
                            Node_indx=FacesNodes[self.fem.Eltype[El]][SurfEl[1]][i]
                            Node=self.fem.Elems[El][Node_indx]
                            Q=self.HeatFlow(self.fem.Coord[Node][0],self.fem.Coord[Node][1],Z,self.T[Node],BcName)
                            HeatFlow_level[BcName]+=GG[El][SurfEl[1]]*Q*dZ
                            H[Node]-=dtau*GG[El][SurfEl[1]]*Q/self.AREA[self.fem.Elems[El][Node_indx]] #J/m3
            if Z<=self.mould.Height:
                for BcName in self.fem.Surfs:
                    self.mould.walls[BcName].HeatUp(HeatFlow_level[BcName])
    #----------Temperature calculation---------------            
            for Node in range(1,self.fem.MaxNodeNum+1):self.T[Node]=self.Temperature(H[Node])
    #----------Preparation for a next level---------------
            Z+=dZ
            if Flag: self.mould.set_Zm(Z)
    #----------------------------------------
        if VTUFileName!='':
            mesh.GetPointData().SetScalars(Res)
            output=vtk.vtkXMLUnstructuredGridWriter()
            output.SetInputData(mesh)
            output.SetFileName(VTUFileName)
            output.Write()
        if VTUHFlowFile!='':
            for BcName in self.fem.Surfs:
                QData_mesh[BcName].GetPointData().SetScalars(QSurf[BcName])
                output=vtk.vtkXMLUnstructuredGridWriter()
                output.SetInputData(QData_mesh[BcName])
                output.SetFileName(VTUHFlowFile[:VTUHFlowFile.rfind('.')]+'_'+BcName+'.vtu')                
                output.Write()
        if MI!='':
            for i in range(2,6):
                if 'Press_'+str(i+1) in MI_model.ESets: MI_model.Surfs['InternalSurf'].append(['Press_'+str(i+1),i])
            MI_model.export_abq(MI)
            MI_model.export_ndload(MI[:MI.rfind('.')]+'_Temp.dat','Temp')
            MI_model.export_fcload(MI[:MI.rfind('.')]+'_press.dat','P')
        if KapaRate!=0:
            f.close()
