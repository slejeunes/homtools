#-----------------------------------------------------------------------
#     Plugin for uniform stress condition on the boundary
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#
#-----------------------------------------------------------------------
from abaqus import *
from abaqusConstants import *
from math import ceil, fabs,sqrt
from numpy import array
#from Numeric import array
import __main__
import random

import textRepr
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class HomogeneTraction:
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,RP1=None,RP2=None,RP3=None,set1=None,is_small_strain=False):
        self.RP1=RP1
        self.RP2=RP2
        self.RP3=RP3
        self.set1=set1
        self.is_small_strain=is_small_strain
        pass
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def build_set(self,the_set,a,small_dist):
       nameset='Picked_set'+str(len(a.sets))
       edges=[]
       faces=[]
       nodes=[]
       vert=[]
       refpoints=[]
       for o in the_set:
          if(type(o).__name__=='Edge'): 
             i1=a.instances[o.instanceName]
             edges.append(i1.edges.findAt(o.pointOn))
          elif(type(o).__name__=='Face'): 
             i1=a.instances[o.instanceName]
             faces.append(i1.faces.findAt(o.pointOn))
          elif(type(o).__name__=='ReferencePoint'): 
             refpoints.append(a.referencePoints.findAt(a.getCoordinates(o)))
          elif(type(o).__name__=='MeshNode'): 
             i1=a.instances[o.instanceName]
             nodes.append(i1.nodes.getByBoundingSphere(o.coordinates,small_dist))
          elif(type(o).__name__=='Vertices'): 
             i1=a.instances[o.instanceName]
             vert.append(i1.vertices.findAt(o.pointOn))
          elif(type(o).__name__=='Vertex'): 
             i1=a.instances[o.instanceName]
             vert.append(i1.vertices.findAt(o.pointOn))
          else: print 'unknown type ',type(o).__name__

       a.Set(edges=edges,referencePoints=refpoints,faces=faces,nodes=nodes, name=nameset,vertices=vert)
       return nameset
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def EffectiveTrac(self,is_smallstrain,dim):
       if(not self.CheckInput()): 
         print "Not all the inputs are given !" 
         return 0

       self.dim=dim
       myModel = self.getCurrentModel()
       a = myModel.rootAssembly
       nameS1=self.build_set(self.set1,a,0.00000001)
        
       if(len(a.allSets[nameS1].nodes)==0):
         print "A mesh is needed!" 
         self.suppressPickedSet(a) 
         return 0
       self.is_small_strain=is_smallstrain
       a.Set(referencePoints=(self.RP1,), name="RefMacro1")
       a.Set(referencePoints=(self.RP2,), name="RefMacro2")
       if(dim==3 and (not is_smallstrain)):
           a.Set(referencePoints=(self.RP3,), name="RefMacro3")
       self.MakeSetsandEquations(nameS1,a)
       self.suppressPickedSet(a) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def suppressPickedSet(self,a):
       for key in a.sets.keys():
          if(key.find('Picked_set')!=-1):
            del a.sets[key]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getTotalvolume(self,a): 
       vol=0.
       if(self.dim==2): vol=a.getMassProperties()['area']
       else: vol=a.getMassProperties()['volume']
       return vol
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getGaussScheme(self,face,nbnodes):
       nd=nbnodes
       if(self.dim==2 and  nd==3):
           gpts=(-0.577350269189626,0.577350269189626)
           gw=(1.,1.)
       elif(self.dim==2 and nd==2):
           gpts=(0.,)
           gw=(2.,)
       elif(self.dim==3 and nd==4):
           gpts=((-0.577350269189626,-0.577350269189626),
                  (0.577350269189626,-0.577350269189626),
                  (0.577350269189626,0.577350269189626),
                  (-0.577350269189626,0.577350269189626))
           gw=(1.,1.,1.,1.)
       elif(self.dim==3 and nd==6):
           gpts=((0.445948490915965,0.445948490915965),
                 (0.10810301816807,0.445948490915965),
                 (0.445948490915965,0.10810301816807),
                 (0.091576213509771,0.091576213509771),
                 (0.816847572980458,0.091576213509771),
                 (0.091576213509771,0.816847572980458))
           gw=(0.111690794839005,0.111690794839005,0.111690794839005,
               0.054975871827661,0.054975871827661,0.054975871827661)
       elif(self.dim==3 and nd==3):
           gpts=((1./6.,1./6.),
                 (2./3.,1./6.),
                 (1./6.,2./3.))
           gw=(1./6.,1./6.,1./6.)
       elif(self.dim==3 and nd==8):
           W3_0=5./9.
           W3_1=8./9.

           gpts=((-0.774596669241483,-0.774596669241483),
                 (0.0,-0.774596669241483),
                 (0.774596669241483,-0.774596669241483),
                 (-0.774596669241483,0.0),
                 (0.0,0.0),
                 (0.774596669241483,0.0),
                 (-0.774596669241483,0.774596669241483),
                 (0.0,0.774596669241483),
                 (0.774596669241483,0.774596669241483))
           gw=(W3_0*W3_0,W3_0*W3_1,W3_0*W3_0,
               W3_1*W3_0,W3_1*W3_1,W3_0*W3_1,
               W3_0*W3_0,W3_0*W3_1,W3_0*W3_0)
       return (gw,gpts)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getInterpolation(self,face,chi,nbnodes):
       nd=nbnodes
       if(self.dim==2 and  nd==3):
           Ni=(0.5*(1.-chi)-0.5*(1.-chi*chi),1.-chi*chi,0.5*(1.+chi)-0.5*(1.-chi*chi))
           de=((-0.5+chi,-2.*chi,0.5+chi),)
       elif(self.dim==2 and nd==2):
           Ni=(0.5*(1.-chi),0.5*(1.+chi))
           de=((-0.5,0.5),)
       elif(self.dim==3 and nd==4):
           Ni=((1.0-chi[0]-chi[1]+chi[1]*chi[0])*0.25,
               (1.0+chi[0]-chi[1]-chi[1]*chi[0])*0.25,
               (1.0+chi[0]+chi[1]+chi[1]*chi[0])*0.25,
               (1.0-chi[0]+chi[1]-chi[1]*chi[0])*0.25)
           de=(((-1.0+chi[1])*0.25,-(-1.0+chi[1])*0.25,(1.0+chi[1])*0.25,-(1.0+chi[1])*0.25), 
               ((-1.0+chi[0])*0.25,(-1.0-chi[0])*0.25,(1.0+chi[0])*0.25,(1.0-chi[0])*0.25))
       elif(self.dim==3 and nd==6):
           u=1.0-chi[0]-chi[1]
           Ni=(u*(2.0*u-1.0),
               chi[0]*(2.0*chi[0]-1.0),
               chi[1]*(2.0*chi[1]-1.0),
               4.0*chi[0]*u,               
               4.0*chi[0]*chi[1],
               4.0*chi[1]*u)
           de=((1.0-4.0*u,-1.0+4.0*chi[0],0.0 ,4.0*(u-chi[0]),4.0*chi[1],-4.0*chi[1]), 
               (1.0-4.0*u,0.0,-1.0+4.0*chi[1],-4.0*chi[0],4.0*chi[0],4.0*(u-chi[1])))
       elif(self.dim==3 and nd==3):
           Ni=(1.0-chi[0]-chi[1],
               chi[0],
               chi[1])
           de=((-1.0,1.0,0.0), 
               (-1.0,0.0,1.0))
       elif(self.dim==3 and nd==8):
           ss=chi[0]*chi[0]
           tt=chi[1]*chi[1]
           st=chi[0]*chi[1]
           sst=chi[0]*chi[0]*chi[1]
           stt=chi[0]*chi[1]*chi[1]
           s2=chi[0]*2.0 
           t2=chi[1]*2.0
           st2=chi[0]*chi[1]*2.0
           Ni=((-1.0+st+ss+tt-sst-stt)/4.0,
               (-1.0-st+ss+tt-sst+stt)/4.0,
               (-1.0+st+ss+tt+stt+sst)/4.0, 
               (-1.0-st+ss+tt+sst-stt)/4.0,
               (1.0-chi[1]-ss+sst)/2.0,               
               (1.0+chi[0]-tt-stt)/2.0,               
               (1.0+chi[1]-ss-sst)/2.0,
               (1.0-chi[0]-tt+stt)/2.0)
           de=(((chi[1]+s2-st2-tt)/4.0,(-chi[1]+s2-st2+tt)/4.0,(chi[1]+s2+st2+tt)/4.0,(-chi[1]+s2+st2-tt)/4.0,-chi[0]+st,(1.0-tt)/2.0,-chi[0]-st,(-1.0+tt)/2.0), 
               ((chi[0]+t2-ss-st2)/4.0,(-chi[0]+t2-ss+st2)/4.0,(chi[0]+t2+ss+st2)/4.0,(-chi[0]+t2+ss-st2)/4.0,(-1.0+ss)/2.0,-chi[1]-st,(1.0-ss)/2.0,-chi[1]+st))
       return (Ni,de)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getNormal(self,face,nodes,de):
       elems=face.getElements()
       if(len(elems)>0):
         center=[0.]*3
         count=0
         #calculate center of element belong to the face   
         for nod in elems[0].getNodes():
           count=count+1 
           center[0]=nod.coordinates[0]+center[0]            
           center[1]=nod.coordinates[1]+center[1]
           center[2]=nod.coordinates[2]+center[2]
         if(count!=0):
           for i in range(0,3):
             center[i]= center[i]/(count*1.0)
         #evaluate the normal
         i=0
         dx_dchi=[0.]*3
         dx_deta=[0.]*3
         for nod in nodes:
            dx_dchi[0]=dx_dchi[0]+de[0][i]*nod.coordinates[0]
            dx_dchi[1]=dx_dchi[1]+de[0][i]*nod.coordinates[1]
            dx_dchi[2]=dx_dchi[2]+de[0][i]*nod.coordinates[2]
            i=i+1
         if(self.dim==3):
           i=0
           for nod in nodes:
             dx_deta[0]=dx_deta[0]+de[1][i]*nod.coordinates[0]
             dx_deta[1]=dx_deta[1]+de[1][i]*nod.coordinates[1]
             dx_deta[2]=dx_deta[2]+de[1][i]*nod.coordinates[2]
             i=i+1
         else: dx_deta[2]=1.0
         ndS=[0.]*3
         ndS[0]=dx_dchi[1]*dx_deta[2]-dx_dchi[2]*dx_deta[1]
         ndS[1]=-dx_dchi[0]*dx_deta[2]+dx_dchi[2]*dx_deta[0]
         ndS[2]=dx_dchi[0]*dx_deta[1]-dx_dchi[1]*dx_deta[0]
         # check the normal orientation
         prodsca=(nodes[0].coordinates[0]-center[0])*ndS[0]+(nodes[0].coordinates[1]-center[1])*ndS[1]+(nodes[0].coordinates[2]-center[2])*ndS[2]
         if(prodsca<0): 
           ndS[0]=-1.*ndS[0]
           ndS[1]=-1.*ndS[1]  
           ndS[2]=-1.*ndS[2]  
         return ndS            
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getBoundaryElements(self,setnodes,a,s1):
       elementfaces=[]
       if(self.dim==3):
          setface=a.allSets[s1].faces
          print setface
          for face in setface:
              elementfaces.extend(face.getElementFaces())
       if(self.dim==2):       
        for n1 in setnodes:
          edgesorfaces=n1.getElemEdges()
          for edf in edgesorfaces:
             isin=True
             for nd in edf.getNodes():
               if not nd in setnodes: 
                    isin=False
                    break
             if(isin and edf not in elementfaces): elementfaces.append(edf)         
       return elementfaces  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def MakeSetsandEquations(self,s1,a):
       setnodes=a.allSets[s1].nodes
       nodesRef=a.allSets["RefMacro1"].referencePoints
       refcoord=a.getCoordinates(nodesRef[0])
       mod = self.getCurrentModel()
       vol=self.getTotalvolume(a)
       #build an array of elements face or edges on the boundary
       elementfaces=self.getBoundaryElements(setnodes,a,s1)

       print 'fin elements faces'
       # build the assembly of nodes contributions
       nodevalue={}
       for face in elementfaces:
          nodes=face.getNodes()
          (gw,gpt)=self.getGaussScheme(face,len(nodes))
          dof=[0.0]*len(nodes)
          ddof=array([dof]*self.dim) 
          for np in range(0,len(gw)):
             (Ni,de)=self.getInterpolation(face,gpt[np],len(nodes)) 
             normal=self.getNormal(face,nodes,de)
             for nd in range(0,len(Ni)): 
                for j in range(0,self.dim):
                   ddof[j][nd]=ddof[j][nd]+gw[np]*Ni[nd]*normal[j]
          i=0
          for n in nodes:
            if(not nodevalue.has_key(n.label)): nodevalue[n.label]=[0.0]*self.dim
            for j in range(0,self.dim):
               nodevalue[n.label][j]=nodevalue[n.label][j]+ddof[j][i]/vol
            i=i+1
       print 'fin node assembly'
       # build the original set of equations
       ntens=4
       dofnum=[1,1,2,2]
       nvalue=[1,2,1,2]
       macrodof=[(1,"RefMacro1",-1),(1,"RefMacro2",-0.5),(1,"RefMacro2",-0.5),(2,"RefMacro1",-1)]
       if(self.dim==2 and not self.is_small_strain):
          macrodof=[(1,"RefMacro1",-1),(1,"RefMacro2",-1),(2,"RefMacro2",-1),(2,"RefMacro1",-1)]
       elif(self.dim==3):
          ntens=9
          dofnum=[1,1,1,2,2,2,3,3,3]
          nvalue=[1,2,3,1,2,3,1,2,3]
          macrodof=[(1,"RefMacro1",-1),(2,"RefMacro1",-1),(3,"RefMacro1",-1),(1,"RefMacro2",-1),(2,"RefMacro2",-1),(3,"RefMacro2",-1),(1,"RefMacro3",-1),(2,"RefMacro3",-1),(3,"RefMacro3",-1)]
          if(self.is_small_strain):
              macrodof=[(1,"RefMacro1",-1),(1,"RefMacro2",-0.5),(2,"RefMacro2",-0.5),(1,"RefMacro2",-0.5),(2,"RefMacro1",-1),(3,"RefMacro2",-0.5),(2,"RefMacro2",-0.5),(3,"RefMacro2",-0.5),(3,"RefMacro1",-1)]

       nodesterms=[[] for i in range(ntens)]
       i=0
       for n in setnodes:
         name1="Num"+str(n.label)
         a.Set(nodes=setnodes[i:i+1],name=name1)         
         for j in range(0,ntens):
            nodesterms[j].append((nodevalue[n.label][nvalue[j]-1],name1,dofnum[j]))
         i=i+1
       #add macro nodes to equations
       for j in range(0,ntens):
         nodesterms[j].append((macrodof[j][2],macrodof[j][1],macrodof[j][0])) 

       #re-organize the equations to remove duplicated dofs (simple Gauss elimination)
       if(self.dim==2):
         for i in range(0,ntens,self.dim):
           if(fabs(nodesterms[i][0][0])>1.e-12):
             pivot=(nodesterms[i+1][0][0]/nodesterms[i][0][0])
             for j in range(1,len(nodesterms[i+1])-1):
                newval=nodesterms[i+1][j][0]-pivot*nodesterms[i][j][0]
                nodesterms[i+1][j]=(newval,nodesterms[i+1][j][1],nodesterms[i+1][j][2])
             j=len(nodesterms[i])-1
             newval=-pivot*nodesterms[i][j][0]
             nodesterms[i+1].append((newval,nodesterms[i][j][1],nodesterms[i][j][2]))
             del nodesterms[i+1][0]
           else: 
             del nodesterms[i][0]
       else:
         # first pass
         for i in range(0,ntens,self.dim):
           if(fabs(nodesterms[i][0][0])>1.e-12):
             pivot1=(nodesterms[i+1][0][0]/nodesterms[i][0][0])
             pivot2=(nodesterms[i+2][0][0]/nodesterms[i][0][0])
             for j in range(1,len(nodesterms[i+1])-1):
                newval1=nodesterms[i+1][j][0]-pivot1*nodesterms[i][j][0]
                newval2=nodesterms[i+2][j][0]-pivot2*nodesterms[i][j][0]
                nodesterms[i+1][j]=(newval1,nodesterms[i+1][j][1],nodesterms[i+1][j][2])
                nodesterms[i+2][j]=(newval2,nodesterms[i+2][j][1],nodesterms[i+2][j][2])
             j=len(nodesterms[i])-1
             newval1=-pivot1*nodesterms[i][j][0]
             newval2=-pivot2*nodesterms[i][j][0]
             nodesterms[i+1].append((newval1,nodesterms[i][j][1],nodesterms[i][j][2]))
             nodesterms[i+2].append((newval2,nodesterms[i][j][1],nodesterms[i][j][2]))
             del nodesterms[i+1][0]
             del nodesterms[i+2][0]
           else: 
             del nodesterms[i][0]
         #second pass    
         for i in range(0,ntens,self.dim):
           if(fabs(nodesterms[i+1][0][0])>1.e-12):
             pivot=(nodesterms[i+2][0][0]/nodesterms[i+1][0][0])
             for j in range(1,len(nodesterms[i+2])-2):
                newval=nodesterms[i+2][j][0]-pivot*nodesterms[i+1][j][0]
                nodesterms[i+2][j]=(newval,nodesterms[i+2][j][1],nodesterms[i+2][j][2])
             j=len(nodesterms[i+1])-1
             newval=nodesterms[i+2][j][0]-pivot*nodesterms[i+1][j][0]
             nodesterms[i+2][j]=(newval,nodesterms[i+2][j][1],nodesterms[i+2][j][2])
             j=len(nodesterms[i+1])-2
             newval=-pivot*nodesterms[i+1][j][0]
             nodesterms[i+2].append((newval,nodesterms[i+1][j][1],nodesterms[i+1][j][2]))
             del nodesterms[i+2][0]
           else: 
             del nodesterms[i+1][0]

       #remove 0 ...
       newterms=[[] for i in range(ntens)]
       for i in range(0,ntens):
         newterms[i]=[term for term in nodesterms[i] if fabs(term[0])>1.e-8]

       #write equations to odb model
       for i in range(0,ntens):
          namec="Constraint-trac"+str(i+1)
          mod.Equation(name=namec, terms=(tuple(newterms[i])))

 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def CheckInput(self):
       if(self.RP1==None): return False
       elif(self.set1==None): return False
       return True  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setRefpoint1(self,RP1): 
       self.RP1=RP1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setRefpoint2(self,RP2): 
       self.RP2=RP2
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setRefpoint3(self,RP3): 
       self.RP3=RP3
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setGroup1(self,set1):
       self.set1=set1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getCurrentModel(self):
      vpName = session.currentViewportName
      modelName = session.sessionState[vpName]['modelName']
      return mdb.models[modelName]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getCurrentViewport(self):
      vpName = session.currentViewportName
      return session.viewports[vpName]
