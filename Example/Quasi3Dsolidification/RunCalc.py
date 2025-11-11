from Quasi3DTemp import Quasi3DTemp, Mould, Plate
from math import tan, pi
def CoatBroad(z):
    if z<0.218: return 0
    elif z<0.25: return 0.004*(z-0.218)/(0.25-0.218)
    else: return 0.004
def CoatNarrow(z):
    if z<0.19: return 0
    elif z<0.222: return 0.004*(z-0.19)/(0.222-0.19)
    else: return 0.004
BroadPlate=Plate_faska(86,116.137E-6,49.42E-3,0.029,0.675,0.1315,CoatBroad)
BroadPlate.Mould_lamda=370 #W/mK
BroadPlate.Coat_lamda=88   #W/mK
BroadPlate.Taper=0.0005    #m
NarrowPlate=Plate(11,86.137E-6,39.42E-3,0.03,0.1315,0.675,CoatNarrow)
NarrowPlate.Mould_lamda=350    #W/mK
NarrowPlate.Coat_lamda=80      #W/mK
NarrowPlate.Taper=0.675*0.011 #m
mould=Mould(0.9,{'BROAD':BroadPlate,'NARROW':NarrowPlate})
model=Quasi3DTemp('mesh.inp',mould)
model.Zones=[[0.9,90],] #(start of zone [m]; water flux [l/m2*min])
model.Z_end=0.9
#------Casting parameters
model.v=1.3 #m/min
model.MoldLevel=0.08 #m
model.CalcMaterialProp(C=0.14,Mn=0.76,Si=0.011,P=0.0085,S=0.0078,Cr=0.014,Cu=0.011,Ni=0.006)
model.dTemp=25
model.mould.walls['BROAD'].Qwat=4515   # Water flow [l/min]
model.mould.walls['NARROW'].Qwat=400  # Water flow [l/min]
mould.walls['BROAD'].Twat_in=20  # Temperature of water [Celsius]
mould.walls['NARROW'].Twat_in=20 # Temperature of water [Celsius]
#--------Test
model.mould.walls['NARROW'].faska_width=0.027
model.mould.walls['NARROW'].faska_fi=10/180*pi
model.RunCalc(Meshdz=0.02,VTUFileName='Temp.vtu', VTUHFlowFile='Flux.vtu')
