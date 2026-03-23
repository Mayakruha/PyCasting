from math import exp, log, pi
import numpy as np
#---------------------------------
#-----SOLIDIFICATION CLASS--------
#---------------------------------
class Solidification:
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0,Mn =0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0,N=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S)             #Solidus temperatere, Celsius
        self.Tlik=1535.1-(88*C+8*Si+5*Mn+1.5*Cr+4*Ni+5*Cu+3*Al+30*P+25*S+80*N)  #Liquidus temperature, Celsius
        self.beta_liq=0.000031 # thermal expansion of luquid steel, 1/K
        self.beta_sol=(2.7-0.16*C+0.039*Mn-0.1*Si-0.019*Cr-0.016*Ni-0.5*P-0.25*S)/140000 # thermal expansion of solid steel, 1/K
        self.kvis=7.91E-7  # kinematic viscosity of luquid steel, m2/sec
        self.thdif=5.45E-6 # thermal diffusivity of luquid steel, m2/sec
        self.lamda_liq=26.0  # conductivity of liqid steel, W/m*K
        self.lamda=32.0      # conductivity of solid steel, W/m*K
        self.L=272E+3      # heat of solidification, J/kg
        self.ro_sol=7300   # solid steel density, kg/m3
        self.ro_liq=7000   # liquid steel density, kg/m3        
        self.Cl=680        # heat capacity for liquid steel, J/kg*K
        self.Cr=795        # heat capacity for solid steel, J/kg*K
        self.Hl=(self.ro_liq*self.Cl+self.ro_sol*self.Cr)*(self.Tlik-self.Tsol)/2+(self.ro_liq+self.ro_sol)*self.L/2 #J/m3
        self.HTC1.ImportProp(self.lamda_liq, self.beta_liq, self.kvis, self.thdif, self.Tsol, self.Tlik, self.Cl, self.ro_liq, self.beta_sol)
        self.HTC2.ImportProp(self.lamda_liq, self.beta_liq, self.kvis, self.thdif, self.Tsol, self.Tlik, self.Cl, self.ro_liq, self.beta_sol) 
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
#----------------------------------------------
#-----HTC CLASSES: ----------------------------
#----------------------------------------------
#----------Parent HTC class--------------------
class HTC: 
    def __init__(self):
        self.LogFile=''    
    def ImportProp(self, lamda, beta_liq, kvis, thdif, Tsol, Tlik, Cl,ro_liq, beta_sol):
        self.lamda=lamda       # conductivity of liqid steel, W/m*K
        self.beta_liq=beta_liq # thermal expansion of luquid steel, 1/K
        self.kvis=kvis         # kinematic viscosity of luquid steel, m2/sec
        self.thdif=thdif       # thermal diffusivity of luquid steel, m2/sec
        self.Tsol=Tsol         # Solidus temperatere, Celsius
        self.Tlik=Tlik         # Liqidus temperatere, Celsius
        self.Cl=Cl             # heat capacity for liquid steel, J/kg*K
        self.ro_liq=ro_liq     #liquid steel density, kg/m3
        self.beta_sol=beta_sol # thermal expansion of solid steel, 1/K
    def htc(self, tm, temp): # tm - time [sec], temp - temperature [C]
        return 0, 0, 0 # htc [W/m2K], temp [C], Heat flux [W/m2]
    def HeatUp(self, J, tm, move): # J - average energy density [J/m2], tm - time [sec], move - movement of border [m]
        """function to change bulk temperature according to energy"""
#-----------------------------------------------
#----------Simple HTC classes-------------------
#----HTC class to handle a table [time, htc]----
class HTC_tab(HTC):
    def SetParams(self, Tmet, tab):
        self.Tmet=Tmet   # function (tm) returns metal temparature over time
        self.htc_tab=tab    # table with two rows [time, htc]
    def htc(self, tm, temp):
        if self.Tmet(tm)<=temp:
            return 0, self.Tmet(tm), 0
        else:
            alfa=np.interp(tm, self.htc_tab[0], self.htc_tab[1])
            return alfa, self.Tmet(tm), alfa*(self.Tmet(tm)-temp)
#----HTC class to handle natural convection----
class HTC_nat_conv(HTC):    
    def SetParams(self, Tmet, Length):
        self.Tmet=Tmet     # function (tm) returns metal temparature over time
        self.Length=Length # function (tm) returns length over time        
    def htc(self, tm, temp):
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
    def __init__(self, Size, Radial=False):
        self.Size=Size/2  # size of bath to calculate temperature
        self.kf=1+Radial
        self.LogFile=''
    def SetParams(self, v, Zm, Tbulk, htc_tab=[[0,12.0],[2000.0,2000.0]]):
        self.v=v             # Casting speed [m/min]
        self.Level=Zm        # Level [m]
        self.Tbulk=Tbulk     # current temeprature        
        self.htc_tab=htc_tab # table with two rows [[z,m];[HTC, W/m2K]]
        self.X0=self.Size
    def htc(self, tm, temp):
        z=self.Level+self.v*tm/60        
        if self.X0<=0:
            return 0, 0
        else:
            alfa=np.interp(z, self.htc_tab[0],self.htc_tab[1])
            Q=alfa*(self.Tbulk - self.Tlik)            
            return alfa, self.Tbulk, Q
    def HeatUp(self, J, tm, move): # J - average energy density [J/m2], tm - time [sec], move - movement of border [m]
        self.X0=self.Size+move # current length 
        if self.X0>0: self.Tbulk-=self.kf*J/self.Cl/self.ro_liq/self.X0
#-----------------------------------------------
#----------CCM HTC classes----------------------
#----General class for CCM----------------------
class CCM(HTC):
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    def geom_coeff(self,Pr):
        return 1
    def flux_alfa(self,T):
        if T>self.flux_Tliq: return self.flux_alfa_liq
        elif T<self.flux_Tmelt: return self.flux_alfa_sol
        else: return (T-self.flux_Tmelt)/(self.flux_Tliq-self.flux_Tmelt)*(self.flux_alfa_liq-self.flux_alfa_sol)+self.flux_alfa_sol
    def SetParams(self, v, Zm,  Qwat_mould, Taper, Qwat_spray, Twat_in=20, Tair=20, Tspraywat=20,
                  flux_Tmelt=1115, flux_Tliq=1145, flux_lamda=1.5, flux_alfa_liq=5900, flux_alfa_sol=1200):
        self.v=v                          # Casting speed [m/min]
        self.Level=Zm                     # Level [m]
        self.Qwat_mould=Qwat_mould        # Water flow in the mould [l/min]
        self.Taper=Taper                  # taper of a side over height, m
        self.Qwat_spray=Qwat_spray        # Water flow in the spray zone - list of (end of zone, m; water flux, l/m2*min)
        self.Twat=Twat_in              # Temperature of water [Celsius]
        self.Tair=Tair                    # Temperature of air [Celsius]
        self.Tspraywat=Tspraywat          # Temperature of water in the spray system [Celsius]
        self.flux_Tmelt=flux_Tmelt        # melting temperature for mould flux, C
        self.flux_Tliq=flux_Tliq          # melting temperature for mould flux, C
        self.flux_lamda=flux_lamda        # flux conductivity, W/mK
        self.flux_alfa_liq=flux_alfa_liq  # HTC for luquid flux, W/m2*K
        self.flux_alfa_sol=flux_alfa_sol  # HTC for solid flux, W/m2*K
        self.Zm=Zm
        self.v_wat=self.Qwat_mould/60000/self.Wat_sec # m/sec
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
        logfile.close()
    def set_level(self, z):
        lamda_wat=0.55748+0.0021525*self.Twat-0.0000097*self.Twat**2 #W/m*K
        visc_wat=1.53555258E-06*exp(-0.036*self.Twat)+2.52805091E-07 #m2/sec
        self.Pr=self.Prandtl(self.Twat)
        self.Re=self.v_wat*self.deff/visc_wat
        self.alfa_wat0=0.023*lamda_wat/self.deff*self.Re**0.8*self.Pr**0.4*self.geom_coeff(self.Pr)
        self.alfa_watz=self.alfa_wat0*(1+self.deff/z/2)
        self.Coat_thick=np.interp(z, self.Coating[0],self.Coating[1])
    def mould_shape(self,z):#m
        return self.Taper*z/self.Height+np.interp(z,self.Shape[0],self.Shape[1])
    def HeatUp(self, J, tm, move): # J - average energy density [J/m2], tm - time [sec], move - movement of border [m]
        z=self.Level+self.v*tm/60
        if z<=self.Height:
            self.Twat-=J/60*self.v*self.Perim/self.Qwat_mould*60/(4230-3.6562*self.Twat-0.02585*self.Twat**2)
    def mould_resist(self, alfa_wat):
        return self.Mould_thick/self.Mould_lamda+self.Coat_thick/self.Coat_lamda+1/alfa_wat
    def htc_resist(self, alfa_wat, temp):
        return 1/self.flux_alfa(temp)/self.v**0.8+self.flux_thickz/self.flux_lamda+self.mould_resist(alfa_wat)
    def htc(self, tm, temp):#W/m2
        z=self.Level+self.v*tm/60
        if z<=self.Height:
            if z>self.Zm:
                self.shrink=450*self.Thb*self.beta_sol*((z-self.Zm)/self.v)**0.5
            if temp>=self.Tsol:
                self.flux_thickz=self.flux_thick_m*(1+2*(temp-self.Tsol)/self.Tsol)
                self.Zm=z
            self.set_level(z)
            if z>self.Zm:
                self.flux_thickz=self.flux_thick_m+self.shrink-self.mould_shape(z)+self.mould_shape(self.Zm)
                if self.flux_thickz<0.0:self.flux_thickz=0.0
            Q=(temp-self.Twat)/self.htc_resist(self.alfa_watz,temp)
            self.TempW=self.Twat+Q/self.alfa_watz
            self.alfa_wat=self.alfa_watz*(self.Prandtl(self.Twat)/self.Prandtl(self.TempW))**0.25
            alfa=1/self.htc_resist(self.alfa_wat,temp)
            return alfa, self.Twat, alfa*(self.Twat-temp) #mould
        else: 
            iz=0
            while iz<len(self.Qwat_spray):
                if z>self.Qwat_spray[iz][0]:
                    iz+=1
                else:
                    break
            if iz>=len(self.Qwat_spray):
                alfa=5.670367E-8*((temp+273)**4-(self.Tair+273)**4)/(temp-self.Tair)
                return alfa, self.Tair, alfa*(self.Tair-temp)
            else:
                alfa= 142*(self.Qwat_spray[iz][1]**0.55)*(1-7.5E-3*self.Tspraywat)*exp(-0.0012*(temp+273))+5.670367E-8*((temp+273)**4-(self.Tair+273)**4)/(temp-self.Tair) #spray
                return alfa, self.Tair, alfa*(self.Tair-temp)
#-----------------------------------------------
class CCM_MouldWithChannels(CCM):
    def __init__(self, Height, Perim, Nch, Sch, Pch, Mould_thick, Taper_size, Shape=[[0.0,1.2],[0.0,0.0]],
                            Coating=[[0.0,1.2],[0.0,0.0]], Mould_lamda=370, Coat_lamda=80):
        #Nch - number of cooling channels
        #Sch - cross section of a cooling channel, m2
        #Pch - perimeter of a cooling channel, m
        self.Height=Height           # Mould height, m
        self.Perim=Perim             # Width of mould plate, m
        self.deff=4*Sch/Pch          # Effective diameter, m
        self.Wat_sec=Nch*Sch         # Summary cross section of cooling channels, m2
        self.Mould_thick=Mould_thick # Distance between water and mould surface, m
        self.Thb=Taper_size/2        # Size of billet to calculate gap due to shrinkage
        self.Shape=Shape             # Shape of mould surface - deflection inside, [[axis,m], [deflection,m]] 
        self.Coating=Coating         # Coating thickness [[axis, m], [thicknes, m]]
        self.Mould_lamda=Mould_lamda # Mould material conductivity, W/mK
        self.Coat_lamda=Coat_lamda   # Coating conductivity, W/mK
        self.LogFile=''
#-----------------------------------------------
class Circle_tube(CCM):
    def __init__(self, Height, D_billet, Wat_thick, Mould_thick, Shape=[[0.0,1.2],[0.0,0.0]],
                            Coating=[[0.0,1.2],[0.0,0.0]], Mould_lamda=370, Coat_lamda=80):
        # D_billet    billet diameter, m
        # Wat_thick   thickness of water layer, m    
        self.Height=Height                                               # Mould height, m
        self.Perim=pi*D_billet                                           # Full perimeter of tube, m
        self.deff=2*Wat_thick                                            # Effective diameter, m
        self.Wat_sec=pi*Wat_thick*(2*(D_billet/2+Mould_thick)+Wat_thick) # Total cross section of cooling channel, m2
        self.Mould_thick=Mould_thick                                     # Distance between water and mould surface, m
        self.Thb=D_billet/2                                              # Size of billet to calculate gap due to shrinkage
        self.Shape=Shape                                                 # Shape of mould surface - deflection inside, [[axis,m], [deflection,m]] 
        self.Coating=Coating                                             # Coating thickness [[axis, m], [thicknes, m]]
        self.Mould_lamda=Mould_lamda                                     # Mould material conductivity, W/mK
        self.Coat_lamda=Coat_lamda                                       # Coating conductivity, W/mK
        self.LogFile=''
        self.R=D_billet/2                                                # Radius of billet cross section, m
        self.Wat_thick=Wat_thick                                         # Thickness of water layer, m
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
    def __init__(self, Height, W_billet, Wat_thick, Mould_thick, Th_billet, Shape=[[0.0,1.2],[0.0,0.0]],
                            Coating=[[0.0,1.2],[0.0,0.0]], Mould_lamda=370, Coat_lamda=80):
        self.Height=Height                                                      # Mould height, m
        self.Perim=2*(W_billet+Th_billet)                                       # Full perimeter of tube, m
        self.deff=2*Wat_thick                                                   # Effective diameter, m
        self.Wat_sec=2*Wat_thick*(W_billet+Th_billet+4*Mould_thick+2*Wat_thick) # Total cross section of cooling channel, m2
        self.Mould_thick=Mould_thick                                            # Distance between water and mould surface, m
        self.Thb=Th_billet/2                                                    # Size of billet to calculate gap due to shrinkage
        self.Wdb=W_billet                                                       # Width of billet cross section, m
        self.Shape=Shape                                                        # Shape of mould surface - deflection inside, [[axis,m], [deflection,m]] 
        self.Coating=Coating                                                    # Coating thickness [[axis, m], [thicknes, m]]
        self.Mould_lamda=Mould_lamda                                            # Mould material conductivity, W/mK
        self.Coat_lamda=Coat_lamda                                              # Coating conductivity, W/mK
        self.LogFile=''
#--------------------------------------------------------------------
#-CREEP FUNCTIONS: Stress[MPa]=func(EpsR[1/sec], Temp[C]))-----------
#--------------------------------------------------------------------
def creep_NISK(EpsR,Temp, A=21900E+6, k=0.2, lmbda=0.004):
    return A*(EpsR**k)*exp(-(Temp+273)*lmbda)
#---------------------------------
#-------OUTPUT FUNCTIONS----------
#---------------------------------
def Log_message(text, file, Screen=True):
  if Screen: print(text)
  file.write(text+'\n')
#---------------------------------
def output_csv(FileName, model):
    f=open(FileName,'w')
    for Name in model.ScalarResList:
        f.write(Name+';')
    f.write('\n')
    n=len(model.results[0])
    for i in range(n):
        for j in range(len(model.ScalarResList)):
            f.write(str(model.results[j][i])+';')
        f.write('\n')
    f.close()
