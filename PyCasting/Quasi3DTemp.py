#Steps for a calculation:
#-Create HTC objects from .Common
#-Create Quasi3DTemp object
#-CalcMaterilaProp(C=0.0,....)
#-meniscus()
#-RunThermCalc()
#-Stiffness()
#
import numpy as np
import vtk
from meshio import read, Mesh, CellBlock
from math import log
from FEMtoolkit import FacesNodes, NodeIntoSurf
from .Common import Log_message, Solidification, creep_NISK, SurfCoord
#----------------------------------------------------------------
#---The function creates mesh for a rectangle section------------
def Rect_Sec(FileName, Width, Thickness, MinSize=1, SurfSize=4, ElmSideRatio=5, DepthRatio=0.6):
    points = []
    cells_data = []
    point_sets = {}
    Depth=min(DepthRatio*Width, DepthRatio*Thickness)
    koeff=(Depth-MinSize)/(Depth-SurfSize)
    n_out=int(log(1+Depth/MinSize*(koeff-1), koeff)) 
    MinSizeFact=Depth*(koeff-1)/(koeff**n_out-1)
    nx_cor=int((Width-Depth)/ElmSideRatio/MinSizeFact)
    ny_cor=int((Thickness-Depth)/ElmSideRatio/MinSizeFact)
    for j in range(ny_cor+n_out+1):
        for i in range(nx_cor+n_out+1):
            if i<=nx_cor:
                X=i*(Width-Depth)/nx_cor
            else:
                X=(Width-Depth)+MinSizeFact*(koeff**(i-nx_cor)-1)/(koeff-1)
            if j<=ny_cor:
                Y=j*(Thickness-Depth)/ny_cor
            else:
                Y=(Thickness-Depth)+MinSizeFact*(koeff**(j-ny_cor)-1)/(koeff-1)
            points.append([X,Y,0])
    for j in range(ny_cor+n_out):
        for i in range(nx_cor+n_out):
            cells_data.append([j*(nx_cor+n_out+1)+i,j*(nx_cor+n_out+1)+i+1,
                       (j+1)*(nx_cor+n_out+1)+i+1,(j+1)*(nx_cor+n_out+1)+i])
    point_sets['broad_side']=[]
    point_sets['broad_side_origin']=[]
    for i in range(nx_cor+n_out+1):
        point_sets['broad_side'].append((nx_cor+n_out+1)*(ny_cor+n_out)+i)
    point_sets['broad_side_origin'].append((nx_cor+n_out+1)*(ny_cor+n_out))
    point_sets['narrow_side']=[]
    point_sets['narrow_side_origin']=[]
    for i in range(ny_cor+n_out+1):
        point_sets['narrow_side'].append((nx_cor+n_out+1)*(i+1)-1)
    point_sets['narrow_side_origin'].append(nx_cor+n_out+1-1)
    model=Mesh(points, [CellBlock('quad', np.array(cells_data, dtype=np.int64))], point_sets=point_sets, faces={})
    model.write(FileName)
#--------------------------------------------------------
#---Function applying initial conditions----------------
def meniscus(model, T0):
    for Node in range(len(model.T)):
        model.T[Node]=T0
    for BcName in model.HTC2:
        for Node in model.fem.point_sets[BcName]:
            model.T[Node]=model.Tlik
#-----------------------------------------------
#----------------BILLET-------------------------
#-----------------------------------------------
class Quasi3DTemp(Solidification):
    def __init__(self, FileName, HTC1, HTC2, LogFile='Casting2D.log', Rewrite=True):
        self.LogFile=LogFile
        if Rewrite: logfile=open(LogFile,'w')
        else: logfile=open(LogFile,'a')
        Log_message('** Model preparation has started',logfile)
        logfile.close()
        #-----------Mesh treatment--------------------------
        mesh_orig=read(FileName)
        points=np.zeros((mesh_orig.points.shape[0],2))
        for Node in range(points.shape[0]):
            for j in range(points.shape[1]):
                points[Node][j]=mesh_orig.points[Node][j]/1000 # millimeters -> meters
        cells = []
        self.AREA=np.zeros(points.shape[0])
        self.Size=0
        XEl=np.zeros(4)
        YEl=np.zeros(4)
        self.DistMin=-1
        cell_sum=[0]
        for blck in mesh_orig.cells:
            if blck.type=='quad':
                cell_sum.append(cell_sum[-1]+blck.data.shape[0])
                for El in blck.data:
                    for i, Node in enumerate(El):
                        XEl[i]=points[Node][0]
                        YEl[i]=points[Node][1]
                    self.AREA[El[0]]+=abs(2*XEl[1]*YEl[3]-2*XEl[3]*YEl[1]+(3*XEl[0]-XEl[2])*(YEl[1]-YEl[3])+(3*YEl[0]-YEl[2])*(XEl[3]-XEl[1]))/16
                    self.AREA[El[1]]+=abs(2*XEl[2]*YEl[0]-2*XEl[0]*YEl[2]+(3*XEl[1]-XEl[3])*(YEl[2]-YEl[0])+(3*YEl[1]-YEl[3])*(XEl[0]-XEl[2]))/16
                    self.AREA[El[2]]+=abs(2*XEl[3]*YEl[1]-2*XEl[1]*YEl[3]+(3*XEl[2]-XEl[0])*(YEl[3]-YEl[1])+(3*YEl[2]-YEl[0])*(XEl[1]-XEl[3]))/16
                    self.AREA[El[3]]+=abs(2*XEl[0]*YEl[2]-2*XEl[2]*YEl[0]+(3*XEl[3]-XEl[1])*(YEl[0]-YEl[2])+(3*YEl[3]-YEl[1])*(XEl[2]-XEl[0]))/16
                    for face in FacesNodes[blck.type]:
                        dist=np.linalg.norm(points[El[face[0]]]-points[El[face[1]]])
                        if self.DistMin>dist: self.DistMin=dist
                        elif self.DistMin<0:
                            self.DistMin=dist
                            #-----------Orientation of elemens--------------------------
                            if XEl[0]*(YEl[1]-YEl[3])+XEl[1]*(YEl[3]-YEl[0])+XEl[3]*(YEl[0]-YEl[1])>0:
                                IndxList=(0,1,2,3)
                            else:
                                IndxList=(0,3,2,1)
                    cells.append([El[i] for i in IndxList])
            else:
                cell_sum.append(cell_sum[-1])
        self.fem=Mesh(points, [CellBlock('quad', np.array(cells))], faces={})
        for Node in range(self.fem.points.shape[0]):
            self.Size+=self.AREA[Node]
        cell_sets={}
        for key in mesh_orig.point_sets:
            key_new=key.lower()
            self.fem.point_sets[key_new]=mesh_orig.point_sets[key]
            if not key in mesh_orig.faces and len(mesh_orig.point_sets[key])>1:
                NodeIntoSurf(self.fem, key_new)
        for key in mesh_orig.faces:
            if not key.lower() in self.fem.point_sets:
                key_new=key.lower()
                self.fem.faces[key_new]={}
                SurfNodeSet=set()
                for FaceSet in mesh_orig.faces[key]:
                    self.fem.faces[key_new][FaceSet]=mesh_orig.faces[key][FaceSet]
                    cell_sets[FaceSet]=[]
                    for i in range(len(mesh_orig.cell_sets[FaceSet])):
                        ElType=mesh_orig.cells[i].type
                        if ElType=='quad':
                            for ElNum in mesh_orig.cell_sets[FaceSet][i]:
                                cell_sets[FaceSet].append(cell_sum[i]+ElNum)
                                for Nd_indx in FacesNodes[ElType][mesh_orig.faces[key][FaceSet]]:
                                    SurfNodeSet.add(mesh_orig.cells[i].data[ElNum][Nd_indx])
                self.fem.point_sets[key_new]=list(SurfNodeSet)
                for FaceSet in mesh_orig.faces[key]:
                    self.fem.cell_sets[FaceSet]=[np.array(cell_sets[FaceSet],dtype=np.int64)]
        #-----------Heat transfer--------------------------
        BoundNodes=set()
        self.HTC1=HTC1       # htc class
        self.HTC1.Port_In=self.Port_Out
        self.HTC1.LogFile=LogFile
        self.HTC2={}
        self.Normals={}
        self.HTC2length={}
        self.GG={}
        for key in HTC2:
            self.HTC2[key]=HTC2[key]
            self.HTC2[key].Port_In=self.Port_Out
            self.HTC2[key].LogFile=LogFile
            if key+'_origin' in self.fem.point_sets:
                OriginNode=self.fem.point_sets[key+'_origin'][0]
                self.HTC2[key].Origin=[self.fem.points[OriginNode][0],self.fem.points[OriginNode][1]]
            self.Normals[key]={}
            self.HTC2length[key]=0 # m
            self.GG[key]={}        # m
            for FaceSet in self.fem.faces[key]:
                face_num=self.fem.faces[key][FaceSet]
                for ElNum in self.fem.cell_sets[FaceSet][0]:
                    Node0=self.fem.cells[0].data[ElNum][FacesNodes['quad'][face_num][0]]
                    Node1=self.fem.cells[0].data[ElNum][FacesNodes['quad'][face_num][1]]
                    Norm=[self.fem.points[Node1][1]-self.fem.points[Node0][1],self.fem.points[Node0][0]-self.fem.points[Node1][0]]
                    EdgeLen=np.linalg.norm(Norm)
                    self.HTC2length[key]+=EdgeLen
                    if not Node0 in self.Normals[key]:
                        self.Normals[key][Node0]=Norm
                        self.GG[key][Node0]=EdgeLen/2
                    else:
                        self.Normals[key][Node0] = [ self.Normals[key][Node0][i]+Norm[i] for i in range(2)]
                        self.GG[key][Node0]+=EdgeLen/2
                    if not Node1 in self.Normals[key]:
                        self.Normals[key][Node1]=Norm
                        self.GG[key][Node1]=EdgeLen/2
                    else:
                        self.Normals[key][Node1] = [ self.Normals[key][Node1][i]+Norm[i] for i in range(2)]
                        self.GG[key][Node1]+=EdgeLen/2
            for Node in self.Normals[key]:
                VecLen=np.linalg.norm(self.Normals[key][Node])
                self.Normals[key][Node][:]/=VecLen
            self.HTC2[key].PointList=[]
            for Node in self.fem.point_sets[key]:
                self.HTC2[key].PointList.append(SurfCoord(self.fem.points[Node],self.Normals[key][Node],self.HTC2[key].Origin))
            BoundNodes|=set(self.fem.point_sets[key])
        for Node in list(BoundNodes):
            self.Size-=self.AREA[Node]
        self.T=np.zeros(self.fem.points.shape[0])
    def dxydksi(self,nu,xy):
        return -(1-nu)*xy[0]+(1-nu)*xy[1]+nu*xy[2]-nu*xy[3]
    def dxydnu(self,ksi,xy):
        return -(1-ksi)*xy[0]-ksi*xy[1]+ksi*xy[2]+(1-ksi)*xy[3]
#-------------------------------------------------------------------------------
# FullTime     - Casting time for calculation, [sec]
# kj           - convergence coefficient (less coefficient, less time step)
# out_dtau     - time period between outputs, [sec]
# wc           - relative coordinates for flux integretaion
    def RunThermCalc(self, FullTime, kj=0.5, out_dtau=0.25, wc=(0.125,0.375)):
        HeatFlow_level={}
        XEl=np.zeros(4)
        YEl=np.zeros(4)
        #-----------------------OUTPUT NODES AND TIME POINTS-----------------------
        self.ScalarResList=['Time [sec]', 'Min Temp [C]', 'Max Temp [C]']
        self.BodyResList=['Temp [C]',]
        self.SurfResList=['Flux [MW/m2K]',]
        self.PointScreen={} #Node and BC for data printing
        for key in self.fem.point_sets:
            if len(self.fem.point_sets[key])==1:
                for BcName in self.HTC2:
                    if self.fem.point_sets[key][0] in self.fem.point_sets[BcName]:
                        if not BcName in self.PointScreen: self.PointScreen[BcName]=[]
                        Node=self.fem.point_sets[key][0]
                        self.PointScreen[BcName].append(Node)
                        prefix=BcName+'-Nd'+str(Node)+':'
                        self.ScalarResList+=[prefix+'Bulk Temp [C]', prefix+' Temp [C]']
        self.ScalarResults=[]
        self.BodyResults=[]
        self.SurfResults={}
        for Name in self.ScalarResList:
            self.ScalarResults.append([])
        for Name in self.BodyResList:
            self.BodyResults.append([])
        TPset={0,FullTime}
        for BcName in self.HTC2:
            TPset = TPset | self.HTC2[BcName].TimePoints
            self.SurfResults[BcName]=[]
            for Name in self.SurfResList:
                self.SurfResults[BcName].append([])
        TPs=list(TPset)
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
        iter_time=0
        dtau=kj*self.DistMin**2*min(self.ro_liq*self.Cl,self.ro_sol*self.Cr)/4/self.lamda  #sek
        n=self.fem.points.shape[0]
        minTemp=min(self.T)
        logfile=open(self.LogFile,'a')
        Log_message('\n** Preparation to the calculation',logfile)
        Log_message(' Solidus Temperature, C: {:6.1f}'.format(self.Tsol),logfile)
        Log_message('Liquidus Temperature, C: {:6.1f}'.format(self.Tlik),logfile)
        Log_message('\tINITIAL CONDITIONS: ',logfile)
        Log_message('Maximum Temperature, C: {:6.1f}'.format(max(self.T)),logfile)
        Log_message('Minimum Temperature, C: {:6.1f}'.format(minTemp),logfile)
        Log_message('\tCALCULATION PARAMETERS: ',logfile)
        Log_message('Number of  nodes: '+str(n),logfile)
        Log_message('Min. distance between nodes, mm: '+str(self.DistMin*1000),logfile)
        Log_message('Step time, sec: '+str(dtau),logfile)
        Log_message('Number of output points: '+str(len(TPs)-1),logfile)
        Log_message('********************************************',logfile)
        Log_message(' Time, sec | Min Temp, C | Max Temp, C |',logfile)
        logfile.close()
        #--------------------------------
        ElTemp=np.zeros(4)
        wN=len(wc)
        H=np.zeros(n)
        for i, Temp in enumerate(self.T):
            H[i]=self.FuncTemp(Temp)
        KK=np.zeros((self.fem.cells[0].data.shape[0],4,4))
        for El in range(self.fem.cells[0].data.shape[0]):
            for i in range(4):
                XEl[i]=self.fem.points[self.fem.cells[0].data[El][i]][0]
                YEl[i]=self.fem.points[self.fem.cells[0].data[El][i]][1]
            L1=((XEl[0]+XEl[1]-XEl[2]-XEl[3])**2+(YEl[0]+YEl[1]-YEl[2]-YEl[3])**2)**0.5/2
            L2=((XEl[0]+XEl[3]-XEl[1]-XEl[2])**2+(YEl[0]+YEl[3]-YEl[1]-YEl[2])**2)**0.5/2
            # n1: ksi=0.5            
            for nu in wc:
                KK[El][0]+=dtau*L1/2*np.array((nu-1,1-nu,nu,-nu))/((self.dxydksi(nu,XEl))**2+(self.dxydksi(nu,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][0]]/wN   #sek/m2
                KK[El][1]-=dtau*L1/2*np.array((nu-1,1-nu,nu,-nu))/((self.dxydksi(nu,XEl))**2+(self.dxydksi(nu,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][1]]/wN   #sek/m2
                KK[El][3]+=dtau*L1/2*np.array((nu-0.5,0.5-nu,nu+0.5,-nu-0.5))/((self.dxydksi(nu+0.5,XEl))**2+(self.dxydksi(nu+0.5,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][3]]/wN   #sek/m2
                KK[El][2]-=dtau*L1/2*np.array((nu-0.5,0.5-nu,nu+0.5,-nu-0.5))/((self.dxydksi(nu+0.5,XEl))**2+(self.dxydksi(nu+0.5,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][2]]/wN   #sek/m2
            # n2: nu=0.5            
            for ksi in wc:
                KK[El][0]+=dtau*L2/2*np.array((ksi-1,-ksi,ksi,1-ksi))/((self.dxydnu(ksi,XEl))**2+(self.dxydnu(ksi,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][0]]/wN   #sek/m2
                KK[El][3]-=dtau*L2/2*np.array((ksi-1,-ksi,ksi,1-ksi))/((self.dxydnu(ksi,XEl))**2+(self.dxydnu(ksi,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][3]]/wN   #sek/m2
                KK[El][1]+=dtau*L2/2*np.array((ksi-0.5,-ksi-0.5,ksi+0.5,0.5-ksi))/((self.dxydnu(ksi+0.5,XEl))**2+(self.dxydnu(ksi+0.5,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][1]]/wN   #sek/m2
                KK[El][2]-=dtau*L2/2*np.array((ksi-0.5,-ksi-0.5,ksi+0.5,0.5-ksi))/((self.dxydnu(ksi+0.5,XEl))**2+(self.dxydnu(ksi+0.5,YEl))**2)**0.5/self.AREA[self.fem.cells[0].data[El][2]]/wN   #sek/m2
        #--------------------------------
        while iter_time<FullTime and minTemp<=self.Tlik:
            #----------BCs preparation and output-------------------------
            self.HTC1.set_level(iter_time)
            for BcName in self.HTC2:
                self.HTC2[BcName].set_level(iter_time)             
            #--------------------Printing & Output-------------------------------
            if iter_time>=TPs[out_iter]:
                alfa1, T1, Q1 = self.HTC1.htc(iter_time,self.Tlik, (0.0,0.0), (1.0,0.0)) #just to get T1
                logfile=open(self.LogFile,'a')
                Log_message('  {:6.1f}   |    {:6.1f}   |    {:6.1f}   |'.format(iter_time, minTemp, T1),logfile)
                logfile.close()
                self.ScalarResults[0].append(iter_time) # Time [sec] 
                self.ScalarResults[1].append(minTemp)   # minTemp [C]
                self.ScalarResults[2].append(T1)        # Metal Temp [C]
                Count=0
                for BcName in self.PointScreen:
                    for Node in self.PointScreen[BcName]:
                        alfa2, T2, Q2 = self.HTC2[BcName].htc(iter_time,self.T[Node],self.fem.points[Node],self.Normals[BcName][Node])
                        self.ScalarResults[3+2*Count].append(T2)             # Bulk Temp [C]
                        self.ScalarResults[3+2*Count+1].append(self.T[Node]) # Surface Temp [C]
                        Count+=1
                self.BodyResults[0].append(self.T.copy())                    # Array-Temp [C]                
            #----------Heat calculation-------------------------
            for BcName in self.HTC2:
                HeatFlow_level[BcName]=0
            for El in range(self.fem.cells[0].data.shape[0]):
                for i in range(4): ElTemp[i]=self.T[self.fem.cells[0].data[El][i]]
                for i in range(4): H[self.fem.cells[0].data[El][i]]+=np.dot(KK[El][i],ElTemp)*self.lamda  #J/m3
            for BcName in self.HTC2:
                if iter_time>=TPs[out_iter]:
                    self.SurfResults[BcName][0].append([])
                for i, Node in enumerate(self.fem.point_sets[BcName]):
                    alfa2, T2, Q2 = self.HTC2[BcName].htc(iter_time,self.T[Node],self.fem.points[Node],self.Normals[BcName][Node])
                    Value=self.GG[BcName][Node]*Q2*dtau #J/m 
                    HeatFlow_level[BcName]-=Value
                    H[Node]+=Value/self.AREA[Node]      #J/m3
                    if iter_time>=TPs[out_iter]:
                        self.SurfResults[BcName][0][-1].append(-Q2/1000000)
            if iter_time>=TPs[out_iter]:
                out_iter+=1
#            self.HTC1.HeatUp(Q1*self.dtau, iter_time, self.dX*(self.k1-k10))
            for BcName in self.HTC2:
                self.HTC2[BcName].HeatUp(HeatFlow_level[BcName]/self.HTC2length[BcName], iter_time, 0)
            #----------Temperature calculation---------------            
            for Node in range(n):self.T[Node]=self.Temperature(H[Node])
            #----------Preparation for a next level---------------
            iter_time+=dtau
            minTemp=min(self.T)
#-------------------------------------------------------------------------------
    def output_vtu(self, vel=0, level=0):
        '''vel - casting speed [m/min], level - mould level [m]/t
        if vel=0, the results are over time'''
        VTUFileName=self.LogFile[:self.LogFile.rfind('.')]+'_res.vtu'
        Points=vtk.vtkPoints()
        mesh=vtk.vtkUnstructuredGrid()
        ElRow=len(self.ScalarResults[0])-1
        n=self.fem.points.shape[0]
        Res=vtk.vtkFloatArray()
        Res.SetName('Temp')
        Res.SetNumberOfValues((ElRow+1)*n)
        for i in range(ElRow+1):
            if vel==0: Z=self.ScalarResults[0][i]
            else: Z=self.ScalarResults[0][i]*vel/60+level
            for Node, Coord in enumerate(self.fem.points):
                Points.InsertNextPoint(Coord[0],Coord[1],Z)
                Res.SetValue(n*i+Node,self.BodyResults[0][i][Node])
        mesh.Allocate(ElRow*self.fem.cells[0].data.shape[0])
        for i in range(ElRow):
            for El in self.fem.cells[0].data:
                mesh.InsertNextCell(vtk.VTK_HEXAHEDRON,8,(El[0]+i*n,El[1]+i*n,El[2]+i*n,El[3]+i*n,
                                          El[0]+(i+1)*n,El[1]+(i+1)*n,El[2]+(i+1)*n,El[3]+(i+1)*n))
        mesh.SetPoints(Points)
        mesh.GetPointData().SetScalars(Res)
        output=vtk.vtkXMLUnstructuredGridWriter()
        output.SetInputData(mesh)
        output.SetFileName(VTUFileName)
        output.Write()
#-------------------------------------------------------------------------------
    def output_surf_vtu(self, vel=0, level=0):
        ElRow=len(self.ScalarResults[0])-1
        for BcName in self.HTC2:
            Points=vtk.vtkPoints()
            for i in range(ElRow+1):
                if vel==0: Z=self.ScalarResults[0][i]
                else: Z=self.ScalarResults[0][i]*vel/60+level               
                for Node in self.fem.point_sets[BcName]:
                    Points.InsertNextPoint(self.fem.points[Node][0],self.fem.points[Node][1],Z)
            SurfNodes={}
            for i, Node in enumerate(self.fem.point_sets[BcName]): SurfNodes[Node]=i
            SurfNodeNum=len(self.fem.point_sets[BcName])
            mesh=vtk.vtkUnstructuredGrid()                
            SurfElemNum=0
            for FaceSet in self.fem.faces[BcName]:
                SurfElemNum+=self.fem.cell_sets[FaceSet][0].shape[0]
            mesh.Allocate(ElRow*2*SurfElemNum)
            for i in range(ElRow):
                for FaceSet in self.fem.faces[BcName]:
                    face_num=self.fem.faces[BcName][FaceSet]
                    for El in self.fem.cell_sets[FaceSet][0]:
                        Node0=self.fem.cells[0].data[El][FacesNodes['quad'][face_num][0]]
                        Node1=self.fem.cells[0].data[El][FacesNodes['quad'][face_num][1]]
                        mesh.InsertNextCell(vtk.VTK_TRIANGLE,3,(i*SurfNodeNum+SurfNodes[Node1],i*SurfNodeNum+SurfNodes[Node0],(i+1)*SurfNodeNum+SurfNodes[Node0]))
                        mesh.InsertNextCell(vtk.VTK_TRIANGLE,3,(i*SurfNodeNum+SurfNodes[Node1],(i+1)*SurfNodeNum+SurfNodes[Node0],(i+1)*SurfNodeNum+SurfNodes[Node1]))
            mesh.SetPoints(Points)
            for i, Name in enumerate(self.SurfResList):                
                Res=vtk.vtkFloatArray()
                Res.SetName(Name)
                Res.SetNumberOfValues((ElRow+1)*SurfNodeNum)
                for j in range(ElRow+1):
                    for k in range(SurfNodeNum):
                        Res.SetValue(SurfNodeNum*j+k,self.SurfResults[BcName][i][j][k])
                mesh.GetPointData().SetScalars(Res)
            output=vtk.vtkXMLUnstructuredGridWriter()
            output.SetInputData(mesh)
            output.SetFileName(self.LogFile[:self.LogFile.rfind('.')]+'_'+BcName+'.vtu')
            output.Write()
#-------------------------------------------------------------------------------
    def Stiffness(self, KapaRate, KapaN=20, creep_law=creep_NISK, vel=0, level=0):
        '''
        KapaRate  - maximum curvature rate for output of section moment inertia. There is no output by default, [1/m*sec]/t
        KapaN     - Number of columns in the output of billet stiffness/t
        creep_law - creep function as follows: Stress[MPa]=func(EpsR[1/sec], Temp[C]))/t
        vel - casting speed [m/min], level - mould level [m]/t
        if vel=0, the results are over time        '''
        FileName=self.LogFile[:self.LogFile.rfind('.')]+'_stiffness.csv'
        f=open(FileName,'w')
        if vel==0: f.write(self.ScalarResList[0]+'\\KapaRate [1/m*sec]')
        else: f.write('Z [m]\\KapaRate [1/m*sec]')
        for i in range(1,KapaN+1):
            f.write(';'+str(i*KapaRate/KapaN))
        f.write('\n')
        for j, tm in enumerate(self.ScalarResults[0]):
            if vel==0: f.write(tm)
            else: f.write(tm*vel/60+level)
            for i in range(1,KapaN+1):
                Value=0
                for Node in range(self.fem.points.shape[0]):
                    if self.BodyResults[0][j][Node]<self.Tsol:
                        EpsR=i*KapaRate/KapaN*self.fem.points[Node][1]
                        Value+=1E+6*creep_law(EpsR,self.BodyResults[0][j][Node])*self.AREA[Node]*self.fem.points[Node][1] # H*m
                f.write(';'+str(Value/KapaRate)) # H*m2*sek
            f.write('\n')
        f.close()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
    def Press(self,z):
        return 9.81*self.ro_liq*(z-self.MouldLevel)
    def RunMICalc(self, wc=(0.125,0.375)):
        XEl=np.zeros(4)
        YEl=np.zeros(4)
        ElRow=len(self.results[0])-1
        MI_model=Mesh()
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
        LevelNum=0
        Z=self.results[0][0]
        Count=20
        dZ=0.25
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
                    for i in self.IndxList:
                        XEl[i]=self.fem.Coord[self.fem.Elems[El][i]][0]
                        YEl[i]=self.fem.Coord[self.fem.Elems[El][i]][1]
                        if MInodes_0[self.fem.Elems[El][i]]==MI_model.MaxNodeNum+1:
                            MInode_count+=1
                            MInodes_0[self.fem.Elems[El][i]]=MInode_count
                            MI_model.Coord[MInode_count]=np.array((XEl[i],YEl[i],Z))
                        Nodes.append(MInodes_0[self.fem.Elems[El][i]])
                    for i in self.IndxList:
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
                        if not 'Press_'+str(self.IndxList_f[PFaces[minNode][maxNode][1]]+3) in MI_model.ESets: MI_model.ESets['Press_'+str(self.IndxList_f[PFaces[minNode][maxNode][1]]+3)]=[]
                        MI_model.ESets['Press_'+str(self.IndxList_f[PFaces[minNode][maxNode][1]]+3)].append(PFaces[minNode][maxNode][2])
                        if PFaces[minNode][maxNode][2] not in MI_model.FaceLoad['P']: MI_model.FaceLoad['P'][PFaces[minNode][maxNode][2]]=[]
                        MI_model.FaceLoad['P'][PFaces[minNode][maxNode][2]].append([self.IndxList_f[PFaces[minNode][maxNode][1]]+2,self.Press(Z+dZ*Count/2)])
        for Node in range(1,self.fem.MaxNodeNum+1):
            if self.T[Node]<=self.Tsol: MI_model.NodeValue['Temp'][MInodes_0[Node]]=self.T[Node]
        for NsetName in self.fem.NSets:
            for Node in self.fem.NSets[NsetName]:
                if MInodes_0[Node]<MI_model.MaxNodeNum+1: MI_model.NSets[NsetName].append(MInodes_0[Node])
        MInodes_0=MInodes_1.copy()
        MInodes_1.fill(MI_model.MaxNodeNum+1)
        for i in range(2,6):
            if 'Press_'+str(i+1) in MI_model.ESets: MI_model.Surfs['InternalSurf'].append(['Press_'+str(i+1),i])
        MI='MechFile.inp'
        MI_model.export_abq(MI)
        MI_model.export_ndload(MI[:MI.rfind('.')]+'_Temp.dat','Temp')
        MI_model.export_fcload(MI[:MI.rfind('.')]+'_press.dat','P')
