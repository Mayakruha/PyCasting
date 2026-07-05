#Steps for a calculation:
#-Create HTC objects from .Common
#-Create Casting1D object
#-CalcMaterilaProp(C=0.0,....)
#-Apply initial conditions
#-RunCalc()
#-RunShrinkage()
#
import numpy as np
import vtk
from .Common import Log_message, Solidification, creep_NISK
#--------------------------------------------------------
#---Functions applying initial conditions----------------
def meniscus(model, T0):
    for i in range(len(model.T)-1):
        model.T[i]=T0
    model.T[-2]=model.Tlik
    alfa1, T1, Q1 = model.HTC1.htc(0, model.Tlik)
    Tin=model.Tlik-Q1/model.lamda*model.dX
    if Tin<T0:
        model.T[-1]=Tin
    model.k1=model.n-1
    model.k2=model.n
def coldplate(model, Tcld, thick):
    model.k1=int((model.Dim-thick)/2/model.dX)
    model.k2=int((model.Dim+thick)/2/model.dX)
    alfa1, T1 = model.HTC1.htc(0, Tcld)
    alfa2, T2 = model.HTC2.htc(0, Tcld)
    for i in range(len(model.T)):
        if i<model.k1:
            model.T[i]=T1
        elif i>model.k2:
            model.T[i]=T2
        else:
            model.T[i]=Tcld
#------------------------------------------
#-------------Casting class----------------
class Casting1D(Solidification):
    def __init__(self, HTC1, HTC2, Dim, n=100, Radial=False, R0=0, LogFile='Casting1D.log', Rewrite=True):
        self.LogFile=LogFile
        if Rewrite: logfile=open(LogFile,'w')
        else: logfile=open(LogFile,'a')
        Log_message('** Model preparation has started',logfile)
        logfile.close()
        self.HTC1=HTC1       # htc class
        self.HTC1.Port_In=self.Port_Out
        self.HTC1.LogFile=LogFile
        self.HTC2=HTC2       # htc class
        self.HTC2.Port_In=self.Port_Out
        self.HTC2.LogFile=LogFile
        self.Dim=Dim
        self.n=n
        self.R0=R0           # initial radius is used if Radial=True
        self.dX=(Dim-R0)/n
        self.T=np.zeros(n+1)
        self.RadFlg=Radial   # type of analysis. If True - radial
        self.ScalarResList=['Time [sec]', 'minTemp [C]', 'Thickness [mm]', 'BulkTemp1 [C]', 'HTC1 [W/m2K]', 'Flux1 [W/m2]',
               'BulkTemp2 [C]', 'HTC2 [W/m2K]', 'Flux2 [W/m2]']
        self.ArrayResList=['Temp [C]',]
        self.k1=0
        self.k2=0
    def RunCalc(self, FullTime, kj=0.5, Epsilon=0.0001, out_dtau=0.25):
        # kj-convergence coefficient for thermal calculations
        # OUTPUT TIME POINTS
        TPs=list(self.HTC1.TimePoints.union(self.HTC2.TimePoints).union({0,FullTime}))
        TPs.sort()
        for i in range(len(TPs)-1):
            iner_point_num=round((TPs[i+1]-TPs[i])/out_dtau)
            if iner_point_num>1:
                dtime=(TPs[i+1]-TPs[i])/iner_point_num
                for j in range(1, iner_point_num):
                    TPs.append(j*dtime+TPs[i])
        TPs.sort()
        out_iter=0
        # ---
        self.dtau=kj*self.dX*self.dX*min(self.ro_liq*self.Cl,self.ro_sol**self.Cr)/4/self.lamda  #sek
        self.Epsilon=Epsilon
        self.results=[]
        for Name in self.ScalarResList+self.ArrayResList:
            self.results.append([])
        iter_time=0
        minTemp=min(self.T)
        n=len(self.T)                       #!!!!!!!!!
        H=np.zeros(n)
        for i in range(n):
            if i==self.k1 and self.k1!=0:
                alfa1, T1, Q1 = self.HTC1.htc(iter_time,self.T[self.k1])
                H[i]=(self.FuncTemp(self.Tlik-Q1*self.dX/2/self.lamda)+self.Hl)/2
            elif i==self.k2 and self.k2!=n-1:
                alfa2, T2, Q2 = self.HTC2.htc(iter_time,self.T[self.k2])
                H[i]=(self.FuncTemp(self.Tlik-Q2*self.dX/2/self.lamda)+self.Hl)/2
            else:
                H[i]=self.FuncTemp(self.T[i])
        k10=self.k1
        k20=self.k2
        logfile=open(self.LogFile,'a')
        Log_message('\n** Preparation to the calculation',logfile)
        Log_message(' Solidus Temperature, C: {:6.1f}'.format(self.Tsol),logfile)
        Log_message('Liquidus Temperature, C: {:6.1f}'.format(self.Tlik),logfile)
        Log_message('\tINITIAL CONDITIONS: ',logfile)
        Log_message('Maximum Temperature, C: {:6.1f}'.format(max(self.T)),logfile)
        Log_message('Minimum Temperature, C: {:6.1f}'.format(minTemp),logfile)
        Log_message('\tCALCULATION PARAMETERS: ',logfile)
        Log_message('Number of  nodes: '+str(n),logfile)
        Log_message('Distance between nodes, mm: '+str(self.dX*1000),logfile)
        Log_message('Step time, sec: '+str(self.dtau),logfile)
        Log_message('Number of output points: '+str(len(TPs)-1),logfile)
        Log_message('********************************************',logfile)
        Log_message(' Time, sec | Min Temp, C | Thickness, mm | Bulk Temp 1, C | Bulk Temp 2, C | HTC 1, kW/m2K | HTC 2, kW/m2K',logfile)
        logfile.close()
        while iter_time<FullTime and minTemp<=self.Tlik:
            #----------BCs preparation and output-------------------------
            self.HTC1.set_level(iter_time)
            self.HTC2.set_level(iter_time)
            alfa1, T1, Q1 = self.HTC1.htc(iter_time,self.T[self.k1])
            alfa2, T2, Q2 = self.HTC2.htc(iter_time,self.T[self.k2])
            H1=self.FuncTemp(T1)
            H2=self.FuncTemp(T2)
            if iter_time>=TPs[out_iter]:
                self.results[0].append(iter_time)              # Time [sec] 
                self.results[1].append(minTemp)                # minTemp [C]
                self.results[2].append(self.dX*(self.k2-self.k1-1)*1000) # Thickness [mm]
                self.results[3].append(T1)                     # BulkTemp1 [C]
                self.results[4].append(alfa1)                  # HTC1 [W/m2*K]
                self.results[5].append(abs(Q1))                # Flux1 [W/m2]
                self.results[6].append(T2)                     # BulkTemp2 [C]
                self.results[7].append(alfa2)                  # HTC2 [W/m2*K]
                self.results[8].append(abs(Q2))                # Flux2 [W/m2]
                self.results[9].append(self.T.copy())          # Array-Temp [C]
                logfile=open(self.LogFile,'a')
                Log_message('  {:6.1f}   |    {:6.1f}   |    {:6.2f}     |     {:6.1f}     |     {:6.2f}     |    {:6.3f}     |    {:6.3f}'.format(iter_time, minTemp, self.dX*(self.k2-self.k1-1)*1000, T1, T2, alfa1/1000, alfa2/1000),logfile)
                logfile.close()
                out_iter+=1
            #----------Heat calculation-------------------------
            for i in range(n):
                if i==0 and self.k1==0 and self.R0==0:
                    H[i]+=(2+2*self.RadFlg)*self.lamda*self.dtau/self.dX/self.dX*(self.T[i+1]-self.T[i])
                elif i==0 and self.k1==0:
                    H[i]+=2*self.dtau*self.lamda/self.dX/self.dX*(self.T[i+1]-self.T[i])*(1+self.RadFlg/2/(self.R0/self.dX+i))/(1+self.RadFlg/4/(self.R0/self.dX+i))
                elif i<self.k1:
                    H[i]=H1
                elif i==self.k1:
                    H[i]+=2*self.dtau*(self.lamda/self.dX/self.dX*(self.T[i+1]-self.T[i])*(1+self.RadFlg/2/(self.R0/self.dX+i))+Q1/self.dX)/(1+self.RadFlg/4/(self.R0/self.dX+i))
                elif i==self.k2:
                    H[i]+=2*self.dtau*(self.lamda/self.dX/self.dX*(self.T[i-1]-self.T[i])*(1-self.RadFlg/2/(self.R0/self.dX+i))+Q2/self.dX)/(1-self.RadFlg/4/(self.R0/self.dX+i))
                elif i>self.k2:
                    H[i]=H2
                else:
                    H[i]+=self.lamda*self.dtau/self.dX/self.dX*((1-self.RadFlg/2/(self.R0/self.dX+i))*self.T[i-1]+(1+self.RadFlg/2/(self.R0/self.dX+i))*self.T[i+1]-2*self.T[i])
            #----------Solidification boundary-------------------
            if self.k1!=0:
                if (self.Hl-H[self.k1])>2*(H1-self.Hl):
                    kf=self.R0/self.dX+self.k1
                    H[self.k1-1]=self.Hl
                    H[self.k1]=H1*(1-self.RadFlg/2/kf)+(H[self.k1]*(1+self.RadFlg/4/kf)-self.Hl*(1-3*self.RadFlg/4/kf))/2
                    self.k1-=1
                elif (H[self.k1]-self.Hl)>(H1-self.Hl)/2:
                    kf=self.R0/self.dX+self.k1
                    H[self.k1+1]=(H[self.k1]*(1+self.RadFlg/4/kf)+2*H[self.k1+1]*(1+self.RadFlg/kf)-2*H1*(1+self.RadFlg/2/kf))/(1+5*self.RadFlg/4/kf)
                    H[self.k1]=H1
                    self.k1+=1
            if self.k2!=n-1:
                if (self.Hl-H[self.k2])>2*(H2-self.Hl):
                    kf=self.R0/self.dX+self.k2
                    H[self.k2+1]=self.Hl
                    H[self.k2]=H2*(1+self.RadFlg/2/kf)+(H[self.k2]*(1-self.RadFlg/4/kf)-self.Hl*(1+3*self.RadFlg/4/kf))/2
                    self.k2+=1
                elif (H[self.k2]-self.Hl)>(H2-self.Hl)/2:
                    kf=self.R0/self.dX+self.k2
                    H[self.k2-1]=(H[self.k2]*(1-self.RadFlg/4/kf)+2*H[self.k2-1]*(1-self.RadFlg/kf)-2*H2*(1-self.RadFlg/2/kf))/(1-5*self.RadFlg/4/kf)
                    H[self.k2]=H2
                    self.k2-=1
            self.HTC1.HeatUp(Q1*self.dtau, iter_time, self.dX*(self.k1-k10))
            self.HTC2.HeatUp(Q2*self.dtau, iter_time, self.dX*(self.k2-k20))
            #----------Temperature calculation------------------
            for i in range(n):
                if i<self.k1:
                    self.T[i]=T1
                elif i==self.k1 and self.k1!=0:
                    self.T[i]=self.Temperature(H[i])+Q1*self.dX/4/self.lamda
                elif i==self.k2 and self.k2!=n-1:
                    self.T[i]=self.Temperature(H[i])+Q2*self.dX/4/self.lamda
                elif i>self.k2:
                    self.T[i]=T2
                else:
                    self.T[i]=self.Temperature(H[i])
            #----------------------------------------------------
            iter_time+=self.dtau
            minTemp=min(self.T)
        logfile=open(self.LogFile,'a')
        Log_message('\n** The calculation has been finished at '+str(iter_time-self.dtau)+' sec because',logfile)
        if minTemp>self.Tlik:
            Log_message('\tThe minimum temperature ('+str(minTemp)+') exceeds ',logfile)
            Log_message('the liquidus temperature ('+str(self.Tlik)+') ',logfile)
        if iter_time>FullTime:
            Log_message('\tThe time ('+str(iter_time)+') exceeds the specified time ('+str(FullTime)+') ',logfile)
        logfile.close()
    def NormalForce(self, j, n0, ksi2, ksi3, creep_law):
        Nf=[0,0]
        for i in range(n0,self.n+1):
            ksiC2=ksi2-self.SpT[i]
            ksiC3=ksi3-self.SpT[i]
            if self.RadFlg:
                ksiC2+=(ksi2+ksi3/2)*(self.n**2/self.i**2-1)
                for ii in range(i,self.n):
                    ksiC2-=3*(self.SpT[ii]+self.SpT[ii+1])*(ii+0.5)/i/i/2
            self.KsiI[i]=2*((ksiC2**2+ksiC3**2+ksiC2*ksiC3)/3)**0.5
            if self.KsiI[i]==0:
                ValueLoc=0
            else:
                ValueLoc=2*creep_law(self.KsiI[i],self.results[self.Tindx][j][i])/self.KsiI[i]/3*self.dX
            if (i==n0)or(i==self.n):
                Nf[0]+=ValueLoc*(2*ksiC2+ksiC3)/2
                Nf[1]+=ValueLoc*(2*ksiC3+ksiC2)/2
            else:
                Nf[0]+=ValueLoc*(2*ksiC2+ksiC3)
                Nf[1]+=ValueLoc*(2*ksiC3+ksiC2)
        return Nf
    def RunShrinkage(self, Epsilon=0.0001, creep_law=creep_NISK):
        logfile=open(self.LogFile,'a')
        Log_message('\n** Mechanical calculation has started',logfile)
        logfile.close()
        if self.ScalarResList[0]!='Time [sec]':
            logfile=open(self.LogFile,'a')
            Log_message('\n**ERROR: There is no time in the scalar data',logfile)
            logfile.close()
            exit()
        ScalarResList_add=['Shrinkage [-]',]
        ArrayResList_add=['Creep rate [1/sec]','Creep strain [-]']
        for Name in ScalarResList_add:
            if not Name in self.ScalarResList:
                self.ScalarResList.append(Name)
                self.results.insert(len(self.ScalarResList)-1,[])
        for Name in ArrayResList_add:
            if not Name in self.ArrayResList:
                self.ArrayResList.append(Name)
                self.results.append([])
        self.Tindx=len(self.ScalarResList)
        Num=len(self.results[0])
        ksi2=0.0
        ksi3=0.0
        Shrink=0.0
        self.SpT=np.zeros(self.n+1)
        self.KsiI=np.zeros(self.n+1)
        EpsC=np.zeros(self.n+1)
        self.results[self.Tindx-1].append(0.0)
        self.results[-2].append(self.KsiI.copy())
        self.results[-1].append(EpsC.copy())
        for j in range(1,Num):
            dtau=self.results[0][j]-self.results[0][j-1]            
            nk=self.n
            for i in range(self.n+1):
                self.KsiI[i]=0.0
                if nk==self.n and self.results[self.Tindx][j-1][i]<self.Tsol: nk=i
                self.SpT[i]=self.beta_sol*(self.results[self.Tindx][j][i]-self.results[self.Tindx][j-1][i])/dtau
            if self.results[self.Tindx][j-1][self.n]<self.Tsol:
                if ksi2==0: dksi2=-self.beta_sol*100
                else: dksi2=-100*self.Epsilon*ksi2
                if ksi3==0: dksi3=-self.beta_sol*100
                else: dksi3=-100*self.Epsilon*ksi3
                eps=10*self.Epsilon
                while eps>self.Epsilon:
                    F0=self.NormalForce(j-1, nk, ksi2, ksi3, creep_law)
                    F1=self.NormalForce(j-1, nk, ksi2+dksi2, ksi3, creep_law)
                    F2=self.NormalForce(j-1, nk, ksi2, ksi3+dksi3, creep_law)
                    if (F0[0]*F1[0])<0: dksi2=dksi2/2
                    elif ((F1[0]-F0[0])*F0[0])>0: dksi2=-dksi2
                    elif (F0[0]*F1[0])==0: dksi2=0
                    if (F0[1]*F2[1])<0: dksi3=dksi3/2
                    elif ((F2[1]-F0[1])*F0[1])>0: dksi3=-dksi3
                    elif (F0[1]*F2[1])==0: dksi3=0
                    eps=max(abs(2*dksi2/(abs(ksi2)+abs(ksi2+dksi2))),abs(2*dksi3/(abs(ksi3)+abs(ksi3+dksi3))))
                    if (F0[0]*F1[0])<0: ksi2+=dksi2/2
                    elif ((F1[0]-F0[0])*F0[0])>=0: ksi2=ksi2
                    else: ksi2+=dksi2
                    if (F0[1]*F2[1])<0: ksi3+=dksi3/2
                    elif ((F2[1]-F0[1])*F0[1])>=0: ksi3=ksi3
                    else: ksi3+=dksi3
            Shrink+=ksi2*dtau
            self.results[self.Tindx-1].append(-Shrink)
            self.results[-2].append(self.KsiI.copy())
            for i in range(self.n+1):
                if i<nk:
                    EpsC[i]=0.0
                else:
                    EpsC[i]=self.results[-1][-1][i]+dtau*self.KsiI[i]
            self.results[-1].append(EpsC.copy())
        logfile=open(self.LogFile,'a')
        Log_message('\n** Mechanical calculation has finished',logfile)
        logfile.close()
#---------------------------------
#--------FUNCTIONS----------------
#---------------------------------
def output_vtu(FileName, model):
    indx=len(model.ScalarResList)
    Points=vtk.vtkPoints()
    mesh=vtk.vtkUnstructuredGrid()
    Dict_nodes={}
    Dict_nodes[0]={}
    cells=[]
    node=0
    for i in range(len(model.results[indx])-1):
        Dict_nodes[(i+1)]={}
        for j in range(len(model.results[indx][i])-1):
            for i1 in range(2):
                for j1 in range(2):
                    if not j+j1 in Dict_nodes[(i+i1)]:
                        Dict_nodes[(i+i1)][j+j1]=node
                        Points.InsertNextPoint((j+j1)*model.dX,0,model.results[0][i+i1])
                        node+=1
            cells.append((Dict_nodes[i][j],Dict_nodes[(i+1)][j],Dict_nodes[(i+1)][j+1],Dict_nodes[i][j+1]))
    mesh.Allocate(len(cells))
    mesh.SetPoints(Points)
    for cell in cells:
        mesh.InsertNextCell(vtk.VTK_QUAD,4,cell)
    for jj, Name in enumerate(model.ArrayResList):
        Res=vtk.vtkFloatArray()
        Res.SetName(Name)
        Res.SetNumberOfValues(node)
        for i in Dict_nodes:
            for j in Dict_nodes[i]:
                Res.SetValue(Dict_nodes[i][j],model.results[indx+jj][i][j])
        mesh.GetPointData().AddArray(Res)
    output=vtk.vtkXMLUnstructuredGridWriter()
    output.SetInputData(mesh)
    output.SetFileName(FileName)
    output.Write()
