from Usadka1D_slab import Mould, Usadka1D_slab
def CoatBroad(z):
    if z<0.218: return 0
    elif z<0.25: return 0.004*(z-0.218)/(0.25-0.218)
    else: return 0.004
BroadPlate=Mould(0.9,86,116.137E-6,49.42E-3,0.029,0.1315,CoatBroad)
BroadPlate.Mould_lamda=370 # W/mK
BroadPlate.Coat_lamda=88   # W/mK
BroadPlate.Taper=0.0005    # m
BroadPlate.Twat_in=20      # Temperature of water [Celsius]
BroadPlate.Qwat=4515       # Water flow [l/min]
#-----------------------
model=Usadka1D_slab(BroadPlate)
model.Z_list=[0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9]
#------Casting parameters
model.a=0.6         # half of slab width, m
model.v=1.3           # casting speed, m/min
model.MouldLevel=0.08 # m
model.CalcMaterialProp(C=0.14,Mn=0.76,Si=0.011,P=0.0085,S=0.0078,Cr=0.014,Cu=0.011,Ni=0.006)
model.dTemp=25
model.RunCalc()
f=open('Results.csv','w')
f.write('Axis, m; Shrinkage, mm; Surface Temperature, Celcius; Heat Flux, MVt/m2; Solid thickness,mm;'+
'Liquid thickness, mm; Temprature in the middle, Celsius; Water temp, Celsius\n')
for line in model.output:
    for value in line:
        f.write(str(value)+';')
    f.write('\n')
f.close()
