from math import exp, log, tan, pi
import numpy as np
import sys
#---------------------------------
#-----SOLIDIFICATION CLASS--------
#---------------------------------
class Solidification:
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0001,Mn=0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0,N=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S)                      # Solidus temperatere, Celsius
        self.Tlik=1536-(78*C+7.6*Si+4.9*Mn+1.3*Cr+3.1*Ni+4.7*Cu+3.6*Al+34.4*P+38*S)      # Liquidus temperature, Celsius
        self.beta_liq=0.000031                                                           # thermal expansion of luquid steel, 1/K
        self.beta_sol=(2.7-0.16*C+0.039*Mn-0.1*Si-0.019*Cr-0.016*Ni-0.5*P-0.25*S)/140000 # thermal expansion of solid steel, 1/K
        self.kvis=7.91E-7    # kinematic viscosity of luquid steel, m2/sec
        self.thdif=5.45E-6   # thermal diffusivity of luquid steel, m2/sec
        self.lamda_liq=26.0  # conductivity of liqid steel, W/m*K
        self.lamda=32.0      # conductivity of solid steel, W/m*K
        self.L=272E+3        # heat of solidification, J/kg
        self.ro_sol=7300     # solid steel density, kg/m3
        self.ro_liq=7000     # liquid steel density, kg/m3        
        self.Cl=680          # heat capacity for liquid steel, J/kg*K
        self.Cr=795          # heat capacity for solid steel, J/kg*K
        self.Hl=(self.ro_liq*self.Cl+self.ro_sol*self.Cr)*(self.Tlik-self.Tsol)/2+(self.ro_liq+self.ro_sol)*self.L/2 #J/m3
    def Port_Out(self):
        return self.lamda_liq, self.beta_liq, self.kvis, self.thdif, self.Tsol, self.Tlik, self.Cl, self.ro_liq, self.beta_sol, self.Size
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
    def output_csv(self, vel=0, level=0):
        '''vel - casting speed [m/min], level - mould level [m]/t
        if vel=0, the results are over time'''
        FileName=self.LogFile[:self.LogFile.rfind('.')]+'.csv'
        f=open(FileName,'w')
        if vel==0: f.write(self.ScalarResList[0])
        else:  f.write('Z [m]')
        for i in range(1,len(self.ScalarResList)):
            f.write(';'+ self.ScalarResList[i])
        f.write('\n')
        n=len(self.ScalarResults[0])
        for i in range(n):
            if vel==0: f.write(self.ScalarResults[0])
            else: f.write(str(self.ScalarResults[0]*vel/60+level))
            for j in range(1, len(self.ScalarResList)):
                f.write(';'+str(self.ScalarResults[j][i]))
            f.write('\n')
        f.close()
#----------------------------------------------
#-----HTC CLASSES: ----------------------------
#----------------------------------------------
#----------Parent HTC class--------------------
class HTC: 
    def __init__(self):
        self.LogFile=''
        self.TimePoints=set()
    def Port_In(self):
        '''Import the following:\t
        self.lamda=lamda       - conductivity of liqid steel, W/m*K\t
        self.beta_liq=beta_liq - thermal expansion of luquid steel, 1/K\t
        self.kvis=kvis         - kinematic viscosity of luquid steel, m2/sec\t
        self.thdif=thdif       - thermal diffusivity of luquid steel, m2/sec\t
        self.Tsol=Tsol         - Solidus temperatere, Celsius\t
        self.Tlik=Tlik         - Liqidus temperatere, Celsius\t
        self.Cl=Cl             - heat capacity for liquid steel, J/kg*K\t
        self.ro_liq=ro_liq     - liquid steel density, kg/m3\t
        self.beta_sol=beta_sol - thermal expansion of solid steel, 1/K'''
        pass
    def set_level(self, tm):
        '''tm - time [sec]'''
        pass
    def htc(self, tm, temp, coord, norm):
        '''return: htc [W/m2K], temp [C], Heat flux [W/m2]\t
        tm - time [sec], temp - temperature [C], coord - coordinate [m, m], norm - normal'''
        return 0, 0, 0  
    def HeatUp(self, J, tm, move): 
        '''function to change bulk temperature according to energy\t
        J - average energy density [J/m2], tm - time [sec], move - movement of border [m]'''
#-----------------------------------------------
#----------Simple HTC classes-------------------
#----HTC class to handle a table [time, htc]----
class HTC_tab(HTC):
    def SetParams(self, Tmet, tab):
        '''Tmet - function (tm) returns metal temparature over time\t
        tab - table with two rows: time and htc'''
        self.Tmet=Tmet    # 
        self.htc_tab=tab  # 
    def htc(self, tm, temp, coord, norm):
        '''return: htc [W/m2K], temp [C], Heat flux [W/m2]\t
        tm - time [sec], temp - temperature [C], coord - coordinate [m, m], norm - normal'''
        if self.Tmet(tm)<=temp:
            return 0, self.Tmet(tm), 0
        else:
            alfa=np.interp(tm, self.htc_tab[0], self.htc_tab[1])
            return alfa, self.Tmet(tm), alfa*(self.Tmet(tm)-temp)
#----HTC class to handle natural convection----
class HTC_nat_conv(HTC):    
    def SetParams(self, Tmet, Length):
        '''Tmet - function (tm) returns metal temparature over time\t
        Length - function (tm) returns length over time'''
        self.Tmet=Tmet
        self.Length=Length
        self.lamda_liq, self.beta_liq, self.kvis, self.thdif, self.Tsol, self.Tlik, self.Cl, self.ro_liq, self.beta_sol = self.Port_In()
    def htc(self, tm, temp, coord, norm):
        '''return: htc [W/m2K], temp [C], Heat flux [W/m2]\t
        tm - time [sec], temp - temperature [C], coord - coordinate [m, m], norm - normal'''
        if self.Tmet(tm)<=temp:
            return 0, self.Tmet(tm)
        else:
            Temp=self.Tmet(tm)
            L=self.Length(tm)
            Ra=9.81*self.beta_liq*(Temp-temp)*L**3/self.kvis/self.thdif
            alfa=0.124*Ra**0.309*self.lamda/L
            return alfa, Temp, alfa*(Temp-temp)
#----HTC class to handle  convection in bath----
class HTC_pool(HTC):
    def __init__(self, Radial=False):
        self.kf=1+Radial
        self.LogFile=''
        self.TimePoints=set()
    def SetParams(self, v, Zm, Tbulk, htc_tab=[[0,12.0],[2000.0,2000.0]]):
        '''v    - Casting speed [m/min]\t
        Zm      - Level [m]\t
        Tbulk   - current temeprature\t
        htc_tab - table with two rows [[z,m];[HTC, W/m2K]]'''
        self.v=v
        self.Level=Zm
        self.Tbulk=Tbulk
        self.htc_tab=htc_tab
        self.lamda_liq, self.beta_liq, self.kvis, self.thdif, self.Tsol, self.Tlik, self.Cl, self.ro_liq, self.beta_sol, self.Size = self.Port_In()
        self.X0=self.Size
    def htc(self, tm, temp, coord, norm):
        '''return: htc [W/m2K], temp [C], Heat flux [W/m2]\t
        tm - time [sec], temp - temperature [C], coord - coordinate [m, m], norm - normal'''
        z=self.Level+self.v*tm/60        
        if self.X0<=0:
            return 0, 0
        else:
            alfa=np.interp(z, self.htc_tab[0],self.htc_tab[1])
            Q=alfa*(self.Tbulk - self.Tlik)            
            return alfa, self.Tbulk, Q
    def HeatUp(self, J, tm, move):
        '''function to change bulk temperature according to energy\t
        J - average energy density [J/m2], tm - time [sec], move - accumulated movement of border [m]'''
        self.X0=self.Size+move # current length 
        if self.X0>0: self.Tbulk-=self.kf*J/self.Cl/self.ro_liq/self.X0
#--------------------------------------------------------------------
#- FUNCTIONS for spray htc: HTC[W/m2K]=func(Val, Tspraywat[C], temp[C])---
#--------------------------------------------------------------------
def Nozaki(Val, Tspraywat, temp):
    '''return HTC [W/m2K]\t
    Val       - water surface density [l/min/m2]\t
    Tspraywat - water temperature [C]\t
    temp      - ambient temperature [C]'''
    return 142*(Val**0.55)*(1-7.5E-3*Tspraywat)*exp(-0.0012*(temp+273))
def DirectHTC(Val, Tspraywat, temp): # Val [W/m2K]
    '''return HTC [W/m2K]\t
    Val - HTC [W/m2K]'''
    return Val
#--------------------------------------------------------------------
def SurfCoord(Coord, Norm, Origin_Coord):
    '''Calculate surface coordinate for a node on the surface\t
    Coord        - global coordinates of a node\t
    Norm         - normal to the surface at the node\t
    Origin_Coord - global coordinates of origin point'''
    coord_loc=(Coord[0]-Origin_Coord[0], Coord[1]-Origin_Coord[1], 0)
    s=np.linalg.norm(coord_loc)
    if np.dot((0,0,1),np.cross(coord_loc,(Norm[0],Norm[1],0)))<0: s=-s
    return s
#-----------------------------------------------
#----------CCM HTC classes----------------------
#----General class for CCM----------------------
class CCM(HTC):
    def NozzleDens(self,s, z, z0, W, L, Val_avrg, s0, kw, fi):
        '''return water density, l/m2*min\t
        s, z     - coordinates, m;  s0, z0 - coordinates of nozzel center, m\t
        W, L     - width and length of spray zone, m\t
        Val_avrg - average value\t
        kw, fi   - define water distribution'''
        x=2*abs(s-s0)/W
        if x<1:
            Val0=(1+1/kw)*2.3 / tan(fi/2) / (1-exp(-2.3/tan(fi/2)))*Val_avrg
            return Val0*(1-x**kw)*exp(-4.6*abs(z-z0)/L/tan(fi/2))
        else:
            return 0
    def Gen_Init(self):
        self.TimePoints=set()
        self.Zones={}
        for part in self.SecCool:
            for ZoneName in part[2]:
                if not ZoneName in self.Zones:
                    self.Zones[ZoneName]=[0,0,''] # Number of nozzles; Nozzle flow, l/min or average HTC, W/m2K; Nozzle Name
                self.Zones[ZoneName][0]+=len(part[2][ZoneName])
        for NozzleName in self.Nozzles:
            for ZoneName in self.Nozzles[NozzleName][3]:
                if ZoneName in self.Zones: self.Zones[ZoneName][2]=NozzleName
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    def geom_coeff(self, Prandtl):
        return 1
    def flux_alfa(self,T):
        if T>self.flux_Tliq: return self.flux_alfa_liq
        elif T<self.flux_Tmelt: return self.flux_alfa_sol
        else: return (T-self.flux_Tmelt)/(self.flux_Tliq-self.flux_Tmelt)*(self.flux_alfa_liq-self.flux_alfa_sol)+self.flux_alfa_sol
    def SetParams(self, v, Zm,  Qwat_mould, Taper, Qwat_spray, data_type='flow', Twat_in=20, Tair=20, Tspraywat=20, flux_Tmelt=1115,
                  flux_Tliq=1145, flux_lamda=1.5, flux_alfa_liq=5900, flux_alfa_sol=1200, alfa_roll=500, alfa_nat=15, HTCSpray_func=Nozaki):
        '''v          - Casting speed [m/min]/t
        Zm            - Mould level [m]/t
        Qwat_mould    - Water flow in the mould [l/min]/t
        Taper         - taper of a side over height [m]/t
        Qwat_spray    - Water flow in the spray zone - {zone name: water flux, l/min or htc, W/m2K}/t
        data_type     - 'flow' if Qwat_spray={zone name: water flux, l/min}; 'htc' if Qwat_spray={zone name: htc, W/m2K}/t
        Twat_in       - Temperature of water [C]/t
        Tair          - Temperature of air [C]/t
        Tspraywat     - Temperature of water in the spray system [C]/t
        flux_Tmelt    - melting temperature for mould flux [C]/t
        flux_Tliq     - melting temperature for mould flux [C]/t
        flux_lamda    - flux conductivity [W/mK]/t
        flux_alfa_liq - HTC for luquid flux [W/m2*K]/t
        flux_alfa_sol - HTC for solid flux [W/m2*K]/t
        alfa_roll     - Contact HTC under rolls [W/m2*K]/t
        alfa_nat      - Natural HTC [W/m2*K]/t
        HTCSpray_func - function(Val, Tspraywat[C], temp[C]) [W/m2*K]'''
        self.v=v
        self.Level=Zm
        self.Qwat_mould=Qwat_mould
        self.Taper=Taper
        self.Qwat_spray=Qwat_spray
        self.data_type=data_type
        self.Twat=Twat_in
        self.Tair=Tair
        self.Tspraywat=Tspraywat
        self.flux_Tmelt=flux_Tmelt
        self.flux_Tliq=flux_Tliq
        self.flux_lamda=flux_lamda
        self.flux_alfa_liq=flux_alfa_liq
        self.flux_alfa_sol=flux_alfa_sol
        self.alfa_roll=alfa_roll
        self.alfa_nat=alfa_nat
        if data_type=='flow': self.HTCSpray_func=HTCSpray_func  # Function to provide HTC from spray flow, W/m2*K
        elif data_type=='htc': self.HTCSpray_func=DirectHTC
        self.Zm=Zm
        self.v_wat=self.Qwat_mould/60000/self.Wat_sec # m/sec
        self.lamda_liq, self.beta_liq, self.kvis, self.thdif, self.Tsol, self.Tlik, self.Cl, self.ro_liq, self.beta_sol, self.Size = self.Port_In()
        self.set_level(self.Zm)                       # Calculate: Pr, Re, alfa_wat0, Rm, alfa_watz
        self.flux_thick_m=self.flux_lamda*((self.flux_Tmelt-self.Twat)/self.flux_alfa(self.Tsol)/self.v**0.8/(self.Tsol-self.flux_Tmelt)-self.mould_resist(self.alfa_watz))
        if self.flux_thick_m<0.0: self.flux_thick_m=0.0
        logfile=open(self.LogFile,'a')
        Log_message('\n** Heat transfer parameters in the mould', logfile)
        Log_message('Prandtl: '+str(self.Pr), logfile)
        Log_message('Water speed, m/sec: '+str(self.v_wat), logfile)
        Log_message('Reynolds: '+str(self.Re), logfile)
        Log_message('Nominal HTC, kW/(m2K): '+str(self.alfa_wat0/1000), logfile)
        Log_message('Solid flux thickness at meniscus, mm: '+str(self.flux_thick_m*1000), logfile)
        self.iz=0
        self.RollDist=self.SecCool[self.iz][0]-self.Height-self.SecCool[self.iz][1]/2
        self.Nozzle_z=self.Height+self.RollDist/2
        Log_message('\n** Heat transfer parameters in the second cooling system', logfile)
        for ZoneName in self.Zones:
            Log_message('Zone '+ZoneName+':', logfile)
            Log_message('  Nozzle: '+self.Zones[ZoneName][2], logfile)
            Log_message('  Number of nozzles '+str(self.Zones[ZoneName][0])+':', logfile)
            if self.data_type=='flow':
                self.Zones[ZoneName][1]=self.Qwat_spray[ZoneName]/self.Zones[ZoneName][0]
                Log_message('  Water flow for a nozzle, l/min: '+str(self.Zones[ZoneName][1]), logfile)
            elif self.data_type=='htc':
                self.Zones[ZoneName][1]=self.Qwat_spray[ZoneName]
                Log_message('  Average HTC for a nozzle, W/m2K: '+str(self.Zones[ZoneName][1]), logfile)
            else:
                Log_message('***ERROR: data_type is specified incorrectly', logfile)
                logfile.close()
                sys.exit()
        logfile.close()
        self.TimePoints.clear()
        self.TimePoints.add(0)
        self.TimePoints.add((self.Height-Zm)/v*60)
        for part in self.SecCool:
            self.TimePoints.add((part[0]-Zm)/v*60)
            if part[1]!=0:
                self.TimePoints.add((part[0]-part[1]/2-Zm)/v*60)
                self.TimePoints.add((part[0]+part[1]/2-Zm)/v*60)
    def set_level(self, tm):
        z=self.Level+self.v*tm/60
        if z<=self.Height:
            # -----mould parameters
            lamda_wat=0.55748+0.0021525*self.Twat-0.0000097*self.Twat**2 #W/m*K
            visc_wat=1.53555258E-06*exp(-0.036*self.Twat)+2.52805091E-07 #m2/sec
            self.Pr=self.Prandtl(self.Twat)
            self.Re=self.v_wat*self.deff/visc_wat
            self.alfa_wat0=0.023*lamda_wat/self.deff*self.Re**0.8*self.Pr**0.4*self.geom_coeff(self.Pr)
            self.alfa_watz=self.alfa_wat0*(1+self.deff/z/2)
            self.Coat_thick=np.interp(z, self.Coating[0],self.Coating[1])
        else:
            self.FullDistr={}
            for s in self.PointList:self.FullDistr[s]=0.0
            Flag=False
            while z>self.SecCool[self.iz][0] and self.iz<len(self.SecCool)-1:
                self.iz+=1
                Flag=True
            # -----spray parameters
            if Flag:
                self.RollDist=self.SecCool[self.iz][0]-self.SecCool[self.iz-1][0]-self.SecCool[self.iz][1]/2-self.SecCool[self.iz-1][1]/2
                self.Nozzle_z=self.SecCool[self.iz-1][0]+self.RollDist/2+self.SecCool[self.iz-1][1]/2
            if (self.iz>0 and z-self.SecCool[self.iz-1][0]>self.SecCool[self.iz-1][1]/2) or (self.SecCool[self.iz][0]-z>self.SecCool[self.iz][1]/2):
                for ZoneName in self.SecCool[self.iz][2]:
                    NozzleName=self.Zones[ZoneName][2]
                    if self.data_type=='flow':
                        Val_avrg=self.Zones[ZoneName][1]/self.Nozzles[NozzleName][0]/self.RollDist
                    elif  self.data_type=='htc':
                        Val_avrg=self.Zones[ZoneName][1]
                    for s0 in self.SecCool[self.iz][2][ZoneName]:
                        for i, s in enumerate(self.PointList):
                            self.FullDistr[s]+=self.NozzleDens(s, z, self.Nozzle_z, self.Nozzles[NozzleName][0], self.RollDist, Val_avrg,
                                               s0=s0, kw=self.Nozzles[NozzleName][1], fi=self.Nozzles[NozzleName][2])
    def mould_shape(self,z):#m
        return self.Taper*z/self.Height+np.interp(z,self.Shape[0],self.Shape[1])
    def HeatUp(self, J, tm, move):
        '''function to change bulk temperature according to energy\t
        J - average energy density [J/m2], tm - time [sec], move - movement of border [m]'''
        z=self.Level+self.v*tm/60
        if z<=self.Height:
            self.Twat-=J/60*self.v*self.Perim/self.Qwat_mould*60/(4230-3.6562*self.Twat-0.02585*self.Twat**2)
    def mould_resist(self, alfa_wat):
        return self.Mould_thick/self.Mould_lamda+self.Coat_thick/self.Coat_lamda+1/alfa_wat
    def htc_resist(self, alfa_wat, temp):
        return 1/self.flux_alfa(temp)/self.v**0.8+self.flux_thickz/self.flux_lamda+self.mould_resist(alfa_wat)
    def htc(self, tm, temp, coord, norm):#W/m2
        '''return: htc [W/m2K], temp [C], Heat flux [W/m2]\t
        tm - time [sec], temp - temperature [C], coord - coordinate [m, m], norm - normal'''
        z=self.Level+self.v*tm/60
        s=SurfCoord(coord, norm, self.Origin)
        if z<=self.Height:
            if temp>=self.Tsol:
                self.flux_thickz=self.flux_thick_m*(1+2*(temp-self.Tsol)/self.Tsol)
                self.Zm=z
            if z>self.Zm:
                shrink=450*np.dot(coord,norm)*self.beta_sol*((z-self.Zm)/self.v)**0.5
                self.flux_thickz=self.flux_thick_m+shrink-self.mould_shape(z)+self.mould_shape(self.Zm)
                if self.flux_thickz<0.0:self.flux_thickz=0.0
            Q=(temp-self.Twat)/self.htc_resist(self.alfa_watz,temp)
            self.TempW=self.Twat+Q/self.alfa_watz
            self.alfa_wat=self.alfa_watz*(self.Prandtl(self.Twat)/self.Prandtl(self.TempW))**0.25
            alfa=1/self.htc_resist(self.alfa_wat,temp)
            return alfa, self.Twat, alfa*(self.Twat-temp) #mould
        elif (self.iz>0 and z-self.SecCool[self.iz-1][0]<self.SecCool[self.iz-1][1]/2) or (self.SecCool[self.iz][0]-z<self.SecCool[self.iz][1]/2):
            return self.alfa_roll, self.Tair, self.alfa_roll*(self.Tair-temp)                               # under roll
        else:
            alfa= self.HTCSpray_func(self.FullDistr[s],self.Tspraywat,temp)+self.alfa_nat+ 0.8*5.670367E-8*((temp+273)**4-(self.Tair+273)**4)/(temp-self.Tair)  # spray
            return alfa, self.Tair, alfa*(self.Tair-temp)
#-----------------------------------------------
class CCM_MouldWithChannels(CCM):
    def __init__(self, Height, Perim, Nch, Sch, Pch, Mould_thick, Nozzles, SecCool, Origin=[0,0],
                 Shape=[[0.0,1.2],[0.0,0.0]], Coating=[[0.0,1.2],[0.0,0.0]], Mould_lamda=370, Coat_lamda=80):
        '''Height   - Mould height [m]/t
        Perim       - Width of mould plate [m]/t
        Nch         - number of cooling channels/t
        Sch         - cross section of a cooling channel [m2]/t
        Pch         - perimeter of a cooling channel [m]/t
        Mould_thick - Distance between water and mould surface [m]/t
        Nozzles     - Nozzles parameters - {Nozzle Name:[W, kw, fi, Zones list]}/t
        SecCool     - Second cooling parameters - [[max z-coordinate, contact width, {Zone name:[s0,s1,..],}],...]/t
        Origin      - Origin Point coordinates [m, m]/t
        Shape       - Shape of mould surface - deflection inside, [[axis,m], [deflection,m]]/t
        Coating     - Coating thickness [[axis, m], [thicknes, m]]/t
        Mould_lamda - Mould material conductivity [W/mK]/t
        Coat_lamda  - Coating conductivity [W/mK]'''
        self.Height=Height
        self.Perim=Perim
        self.deff=4*Sch/Pch          # Effective diameter, m
        self.Wat_sec=Nch*Sch         # Summary cross section of cooling channels, m2
        self.Mould_thick=Mould_thick
        self.Nozzles=Nozzles
        self.SecCool=SecCool
        self.Origin=Origin
        self.Shape=Shape
        self.Coating=Coating
        self.Mould_lamda=Mould_lamda
        self.Coat_lamda=Coat_lamda
        self.PointList=[0.0,]        # List of coordinates (s) for htc calculations
        self.LogFile=''
        self.Gen_Init()
#-----------------------------------------------
class Circle_tube(CCM):
    def __init__(self, Height, D_billet, Wat_thick, Mould_thick, Nozzles, SecCool, Origin=[0,0],
                 Shape=[[0.0,1.2],[0.0,0.0]], Coating=[[0.0,1.2],[0.0,0.0]], Mould_lamda=370, Coat_lamda=80):
        '''Height   - Mould height [m]/t
        D_billet    - billet diameter [m]/t
        Wat_thick   - thickness of water layer [m]/t
        Mould_thick - Distance between water and mould surface [m]/t
        Nozzles     - Nozzles parameters - {Nozzle Name:[W, kw, fi, Zones list]}/t
        SecCool     - Second cooling parameters - [[max z-coordinate, contact width, {Zone name:[s0,s1,..],}],...]/t
        Origin      - Origin Point coordinates [m, m]/t
        Shape       - Shape of mould surface - deflection inside, [[axis,m], [deflection,m]]/t
        Coating     - Coating thickness [[axis, m], [thicknes, m]]/t
        Mould_lamda - Mould material conductivity [W/mK]/t
        Coat_lamda  - Coating conductivity [W/mK]'''  
        self.Height=Height
        self.Perim=pi*D_billet                                           # Full perimeter of tube, m
        self.Wat_thick=Wat_thick
        self.Mould_thick=Mould_thick
        self.Nozzles=Nozzles
        self.SecCool=SecCool
        self.Origin=Origin
        self.Shape=Shape
        self.Coating=Coating
        self.Mould_lamda=Mould_lamda
        self.Coat_lamda=Coat_lamda
        self.deff=2*Wat_thick                                            # Effective diameter, m
        self.Wat_sec=pi*Wat_thick*(2*(D_billet/2+Mould_thick)+Wat_thick) # Total cross section of cooling channel, m2
        self.R=D_billet/2                                                # Radius of billet cross section, m
        self.PointList=[0.0,]                                            # List of coordinates (s) for htc calculations
        self.LogFile=''
        self.Gen_Init()
    def geom_coeff(self, Pr):
        return (1-0.45/(2.4+Pr))*(1+self.Wat_thick/(self.R+self.Mould_thick))**(0.16/Pr**0.15)
    def mould_resist(self, alfa_wat):
        return (self.R-self.flux_thickz)*(log(1+self.Coat_thick/self.R)/self.Coat_lamda+\
                log(1+self.Mould_thick/(self.R+self.Coat_thick))/self.Mould_lamda+\
                1/alfa_wat/(self.R+self.Coat_thick+self.Mould_thick))
    def htc_resist(self, alfa_wat, temp):
        return 1/self.flux_alfa(temp)/self.v**0.8+(self.R-self.flux_thickz)*log(1+self.flux_thickz/(self.R-self.flux_thickz))/self.flux_lamda+\
                self.mould_resist(alfa_wat)
#-----------------------------------------------
class Rect_tube(CCM):
    def __init__(self, Height, Perim, Wat_thick, Mould_thick, Nozzles, SecCool, Origin=[0,0],
                 Shape=[[0.0,1.2],[0.0,0.0]], Coating=[[0.0,1.2],[0.0,0.0]], Mould_lamda=370, Coat_lamda=80):
        '''Height   - Mould height [m]/t
        Perim       - Full perimeter of tube [m]/t
        Wat_thick   - thickness of water layer [m]/t
        Mould_thick - Distance between water and mould surface [m]/t
        Nozzles     - Nozzles parameters - {Nozzle Name:[W, kw, fi, Zones list]}/t
        SecCool     - Second cooling parameters - [[max z-coordinate, contact width, {Zone name:[s0,s1,..],}],...]/t
        Origin      - Origin Point coordinates [m, m]/t
        Shape       - Shape of mould surface - deflection inside, [[axis,m], [deflection,m]]/t
        Coating     - Coating thickness [[axis, m], [thicknes, m]]/t
        Mould_lamda - Mould material conductivity [W/mK]/t
        Coat_lamda  - Coating conductivity [W/mK]'''  
        self.Height=Height
        self.Perim=Perim
        self.deff=2*Wat_thick        # Effective diameter, m
        self.Mould_thick=Mould_thick
        self.Nozzles=Nozzles
        self.SecCool=SecCool
        self.Origin=Origin
        self.Shape=Shape
        self.Coating=Coating
        self.Mould_lamda=Mould_lamda
        self.Coat_lamda=Coat_lamda
        self.Wat_sec=Wat_thick*Perim # Total cross section of cooling channel, m2
        self.PointList=[0.0,]        # List of coordinates (s) for htc calculations
        self.LogFile=''
        self.Gen_Init()
#--------------------------------------------------------------------
#-CREEP FUNCTIONS: Stress[MPa]=func(EpsR[1/sec], Temp[C]))-----------
#--------------------------------------------------------------------
def creep_NISK(EpsR,Temp, B=21900E+6, mu=0.2, lmbda=0.004):
    return B*(EpsR**mu)*exp(-(Temp+273)*lmbda)
#---------------------------------
#-------OUTPUT FUNCTIONS----------
#---------------------------------
def Log_message(text, file, Screen=True):
    if Screen: print(text)
    file.write(text+'\n')
