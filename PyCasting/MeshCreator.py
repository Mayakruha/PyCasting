w=1.350/2 #m
a=0.263/2 #m
d=0.03 #m
Nd=30
Nw=170
Na=30
import numpy as np
import sys
sys.path.append('D:\\Scripts\\FEMtoolkit\\')
sys.path.append('C:\\Program Files\\ParaView 5.11.1\\bin\\Lib\\site-packages\\')
from FEMtoolkit import FEMtoolkit
model=FEMtoolkit()
model.MaxNodeNum=(Nw+Nd+1)*(Na+Nd+1)
model.MaxElemNum=(Nw+Nd)*(Na+Nd)
model.Coord=np.full(model.MaxNodeNum+1,None)
model.Elems=np.ones(model.MaxElemNum+1,dtype=tuple)
model.Eltype=np.zeros(model.MaxElemNum+1,dtype=np.int8)
model.TypeList={}
model.TypeList[3]='CPS4'
model.NSets['NARROW']=[]
model.NSets['BROAD']=[]
model.ESets['EAll']=[]
for i in range(Nw+Nd+1):
    model.NSets['BROAD'].append((Nw+Nd+1)*(Na+Nd)+i+1)
    for j in range(Na+Nd+1):
        if i<=Nw: X=(w-d)/Nw*i
        else: X=w-d+d/Nd*(i-Nw)
        if j<=Na: Y=(a-d)/Na*j
        else: Y=a-d+d/Nd*(j-Na)
        model.Coord[(Nw+Nd+1)*j+i+1]=np.array((X,Y,0))
for j in range(Na+Nd+1):
    model.NSets['NARROW'].append((Nw+Nd+1)*(j+1))
for i in range(Nw+Nd):
    for j in range(Na+Nd):
        model.Elems[(Nw+Nd)*j+i+1]=[(Nw+Nd+1)*j+i+1,(Nw+Nd+1)*j+i+2,(Nw+Nd+1)*(j+1)+i+2,(Nw+Nd+1)*(j+1)+i+1]
        model.Eltype[(Nw+Nd)*j+i+1]=3
        model.ESets['EAll'].append((Nw+Nd)*j+i+1)
model.NodeIntoSurf('NARROW')
model.NodeIntoSurf('BROAD')
model.NodeValue['Set']={}
for i in range(1,model.MaxNodeNum+1):
    if (i in model.NSets['NARROW'])and(i in model.NSets['BROAD']): model.NodeValue['Set'][i]=3
    elif i in model.NSets['NARROW']:model.NodeValue['Set'][i]=1
    elif i in model.NSets['BROAD']:model.NodeValue['Set'][i]=2
    else: model.NodeValue['Set'][i]=0
model.export_VTU('D:\\NLMK\\mesh.vtu', 'EAll',Scalars=['Set',])
model.export_abq('D:\\NLMK\\mesh.inp')
