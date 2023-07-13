#-----------------------------------------------------------------------
#     Plugin for periodic condition on the boundary for plates and beams
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 14/04/11
#
#
#-----------------------------------------------------------------------
import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../Materials'))
from periodicBoundary import PeriodicBoundary,Tree
from abaqus import *
from abaqusConstants import *
from math import ceil, fabs,sqrt
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
from copy import deepcopy
from numpy import *


class PeriodicBoundaryStruct(PeriodicBoundary):
    "generation of boundary periodic condition for Structures" 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,RP1=None,RP2=None,set1=None,set2=None,n1=None,n2=None,refplane=None):
        PeriodicBoundary.__init__(self,RP1=RP1,RP2=RP2,set1=set1,set2=set2)    
        self.n1=n1
        self.n2=n2
        self.refplane=refplane
        self.is_small_strain=True
        self.dim=3
        self.refplane=None
        self.CSYS=None
        self.set1=[]
        self.set2=[]
        pass

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def PeriodicStruct(self,Vx,Vy,Vz,is_plate):
       self.Vx=Vx
       self.Vy=Vy
       self.Vz=Vz
       self.is_plate=is_plate
       # check the local csys
       ret=self.verifLOCALCSYS()
       if(ret!=""):
         print ret
         return 0  
       small=[Vx,Vy,Vz]
       small.sort()
       if(not self.CheckInput()): 
         print "ERROR: Not all the inputs are given !" 
         return 0

       myModel = self.getCurrentModel()
       a = myModel.rootAssembly
       # build all the sets needed
       nameS1=[]
       nameS2=[]
       for i in range(0,len(self.set1)):
         nameS1.append(self.build_set(self.set1[i],a,small[0]*0.00000001))
         nameS2.append(self.build_set(self.set2[i],a,small[0]*0.00000001))
       a.Set(referencePoints=(self.RP1,), name="RefMacro1")
       a.Set(referencePoints=(self.RP2,), name="RefMacro2")
       self.is_small_strain=True
       # check if we have a mesh
       if(len(a.allSets[nameS1[0]].nodes)==0):
         print "ERROR: A mesh is needed!" 
         self.SuppressPickedSet(a) 
         return 0
       # check the mesh periodicity
       for i in range(0,len(self.set1)): 
         if(not self.VerifNodeCoord(nameS1[i],nameS2[i],a)):
           self.SuppressPickedSet(a) 
           return 0
         self.MakeNodeSetsandEquations(nameS1[i],nameS2[i],a)

       self.SuppressPickedSet(a) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def MakeNodeSetsandEquations(self,s1,s2,a):
       dim=3
       nodes1=a.allSets[s1].nodes
       nodes2=a.allSets[s2].nodes

       nodesref1=a.allSets[s1].referencePoints
       nodesref2=a.allSets[s2].referencePoints

       nodes2elem=range(len(nodes2))
       
       mod = self.getCurrentModel()
       i=0

       A=matrix([[self.n1[0],self.n1[1],self.n1[2]],
                 [self.n2[0],self.n2[1],self.n2[2]],
                 [self.n3[0],self.n3[1],self.n3[2]]])
       Ainv=A.T
       ori=self.CSYS.origin.pointOn
       lori=inner(Ainv.T,ori)
       mytree=Tree()
       mytree.getCurrentGraph('Const_')
       for n1 in nodes1:
         name1="Num"+str(n1.label)
         if(a.sets.has_key(name1)==0): a.Set(nodes=nodes1[i:i+1],name=name1)
         for j in nodes2elem:
             c2=(nodes2[j].coordinates[0]-self.Vx,nodes2[j].coordinates[1]-self.Vy,nodes2[j].coordinates[2]-self.Vz)
             dist=sqrt((n1.coordinates[0]-c2[0])**2+(n1.coordinates[1]-c2[1])**2+(n1.coordinates[2]-c2[2])**2)
             if(dist<0.1*self.carac):  
               name2="Num"+str(nodes2[j].label)
               nname2=nodes2[j].label
               a.Set(nodes=nodes2[j:j+1],name=name2)
               x1=squeeze(asarray(inner(A,n1.coordinates)-lori))
               x2=squeeze(asarray(inner(A,nodes2[j].coordinates)-lori))
               nodes2elem.remove(j)
               break
         # Tree construction and checking of cycling DOF's
         branch=[int(n1.label),int(nname2)]
         if(mytree.addBranch(branch) and not mytree.iscycle()):
           cname='Const_'+str(n1.label)+'_'+str(nname2)
           self.MakeSingleEquation(mod,cname,name1,name2,x1,x2,Ainv)
         else:
           mytree.restoreTree()
         i=i+1
       # checking of DOF ordering to avoid undesired DOF suppression...
       self.checkMasterAndSlaves(mytree)

       #Same treatement for RefPoints 
     #  myreftree=Tree()
     #  myreftree.getCurrentGraph('ConstRef_')
     #  nodes2elem=range(len(nodesref2))

     #  for n1 in nodesref1:
     #    name1="Ref_"+str(len(a.sets))
     #    nname1=find(n1,a.referencePoints.items())
     #    n1coord=a.getCoordinates(n1) 
     #    if(a.sets.has_key(name1)==0): a.Set(referencePoints=(n1,),name=name1)
     #    for j in nodes2elem:
      #     n2coord=a.getCoordinates(nodesref2[j])
      #     c2=(n2coord[0]-self.Vx,n2coord[1]-self.Vy,n2coord[2]-self.Vz)
      #     dist=sqrt((n1coord[0]-c2[0])**2+(n1coord[1]-c2[1])**2+(n1coord[2]-c2[2])**2)
      #     if(dist<1e-5):  
      #       name2="Ref_"+str(len(a.sets))
      #       nname2=find(nodesref2[j],a.referencePoints.items())
      #       a.Set(referencePoints=(nodesref2[j],),name=name2) 
      #       x1_x2[0]=n1coord[0]-n2coord[0]
      #       x1_x2[1]=n1coord[1]-n2coord[1]
      #       x1_x2[2]=n1coord[2]-n2coord[2]
      #       nodes2elem.remove(j)
      #       break
           # Tree construction and checking of cycling DOF's
      #   branch=[int(nname1),int(nname2)]
      #   if(myreftree.addBranch(branch) and not myreftree.iscycle()):
      #       cname='ConstRef_'+str(nname1)+'_'+str(nname2)
      #       self.MakeSingleEquation(mod,cname,name1,name2,x1_x2,dim,True)
      #   else:
      #       myreftree.restoreTree()
      # self.checkMasterAndSlaves(myreftree)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def MakeSingleEquation(self,mod,cname,name1,name2,x1,x2,A):

       if(self.is_small_strain and self.is_plate):
            H=matrix([[x1[0]-x2[0],0,(x1[1]-x2[1])/2.,x1[2]*x1[0]-x2[2]*x2[0],0,(x1[2]*x1[1]-x2[2]*x2[1])/2.],
                 [0,x1[1]-x2[1],(x1[0]-x2[0])/2.,0,x1[2]*x1[1]-x2[2]*x2[1],(x1[2]*x1[0]-x2[2]*x2[0])/2.],
                 [0,0,0,-0.5*(x1[0]*x1[0]-x2[0]*x2[0]),-0.5*(x1[1]*x1[1]-x2[1]*x2[1]),-0.5*(x1[0]*x1[1]-x2[0]*x2[1])]])
            H=(A)*H
            mod.Equation(name=cname+'_1', terms=((1.0, name1, 1),(-1.0, name2, 1),(-H[0,0],"RefMacro1",1),
                                                 (-H[0,1],"RefMacro1",2),(-H[0,2],"RefMacro1",3),(-H[0,3],"RefMacro2",1),(-H[0,4],"RefMacro2",2),(-H[0,5],"RefMacro2",3)))
            mod.Equation(name=cname+'_2', terms=((1.0, name1, 2),(-1.0, name2, 2),(-H[1,0],"RefMacro1",1),
                                                 (-H[1,1],"RefMacro1",2),(-H[1,2],"RefMacro1",3),(-H[1,3],"RefMacro2",1),(-H[1,4],"RefMacro2",2),(-H[1,5],"RefMacro2",3)))
            mod.Equation(name=cname+'_3', terms=((1.0, name1, 3),(-1.0, name2, 3),(-H[2,0],"RefMacro1",1),
                                                 (-H[2,1],"RefMacro1",2),(-H[2,2],"RefMacro1",3),(-H[2,3],"RefMacro2",1),(-H[2,4],"RefMacro2",2),(-H[2,5],"RefMacro2",3)))
       elif(self.is_small_strain and not self.is_plate):
            H=matrix([[0,0,0,0,(x1[2]*x1[2]-x2[2]*x2[2])/2.,-(x1[2]*x1[1]-x2[2]*x2[1])],
                 [0,0,0,-(x1[2]*x1[2]-x2[2]*x2[2])/2.,0.,(x1[2]*x1[0]-x2[2]*x2[0])],
                 [x1[2]-x2[2],0,0,(x1[2]*x1[1]-x2[2]*x2[1]),-(x1[2]*x1[0]-x2[2]*x2[0]),0.]])
            H=(A)*H
            mod.Equation(name=cname+'_1', terms=((1.0, name1, 1),(-1.0, name2, 1),(-H[0,0],"RefMacro1",1),
                                                 (-H[0,1],"RefMacro1",2),(-H[0,2],"RefMacro1",3),(-H[0,3],"RefMacro2",1),(-H[0,4],"RefMacro2",2),(-H[0,5],"RefMacro2",3)))
            mod.Equation(name=cname+'_2', terms=((1.0, name1, 2),(-1.0, name2, 2),(-H[1,0],"RefMacro1",1),
                                                 (-H[1,1],"RefMacro1",2),(-H[1,2],"RefMacro1",3),(-H[1,3],"RefMacro2",1),(-H[1,4],"RefMacro2",2),(-H[1,5],"RefMacro2",3)))
            mod.Equation(name=cname+'_3', terms=((1.0, name1, 3),(-1.0, name2, 3),(-H[2,0],"RefMacro1",1),
                                                 (-H[2,1],"RefMacro1",2),(-H[2,2],"RefMacro1",3),(-H[2,3],"RefMacro2",1),(-H[2,4],"RefMacro2",2),(-H[2,5],"RefMacro2",3)))
       else:

            print "Not yet implemented"
            return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setShellRefplane(self,PL1):
       self.refplane=PL1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setCSYS(self,CSYS):
       self.CSYS=CSYS
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def verifLOCALCSYS(self):
       self.n1=self.CSYS.axis1.direction
       self.n2=self.CSYS.axis2.direction
       self.n3=self.CSYS.axis3.direction

#       plN=self.refplane.normal
#       pla1=self.CSYS.axis1.direction
#       pla2=self.CSYS.axis2.direction
#       pla3=self.CSYS.axis3.direction
#       scal1=pla1[0]*plN[0]+pla1[1]*plN[1]+pla1[2]*plN[2]
 #      scal2=pla2[0]*plN[0]+pla2[1]*plN[1]+pla2[2]*plN[2]
 #      scal3=pla3[0]*plN[0]+pla3[1]*plN[1]+pla3[2]*plN[2]
 #      is_col=False
 #      if(scal1==0 and scal2==0 and scal3!=0):
  #          self.n1=self.CSYS.axis1.direction
  #          self.n2=self.CSYS.axis2.direction
  #          self.n3=self.CSYS.axis3.direction
  #          is_col=True
  #     elif(scal1!=0 and scal2==0 and scal3==0):
  #          self.n1=self.CSYS.axis2.direction
  #          self.n2=self.CSYS.axis3.direction
  #          self.n3=self.CSYS.axis1.direction            
  #          is_col=True
  #     elif(scal1==0 and scal2!=0 and scal3==0):
  #          self.n1=self.CSYS.axis3.direction
  #          self.n2=self.CSYS.axis1.direction
  #          self.n3=self.CSYS.axis2.direction
  #          is_col=True
  #     if(not is_col): return "ERROR: The normal of the datum plane is not include in the local CSYS"
  #     apoint=self.refplane.pointOn
  #     d=-self.n3[0]*apoint[0]-self.n3[1]*apoint[1]-self.n3[2]*apoint[2]                     
  #     origin=self.CSYS.origin.pointOn
   #    if(fabs(self.n3[0]*origin[0]+self.n3[1]*origin[1]+self.n3[2]*origin[2]+d)>1e-5): return "ERROR: the local CSYS does not have his origin on the plane"  
       return ""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setGroup2(self,set2,ind): 
       if(len(self.set2)==ind): 
          self.set2.append(set2)
       else:
          self.set2[ind]=set2
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setGroup1(self,set1,ind):
       if(len(self.set1)==ind): 
          self.set1.append(set1)
       else:
          self.set1[ind]=set1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def showMidPlane(self):
       if(self.CSYS):
	 highlight(self.CSYS)
         myModel = self.getCurrentModel()
         a = myModel.rootAssembly
         self.plane=a.DatumPlaneByPointNormal(point=self.CSYS.origin.pointOn,normal=self.CSYS.axis3)  
         highlight(self.plane)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def removeMidPlane(self):
       if(self.CSYS):
	 unhighlight(self.CSYS)
         myModel = self.getCurrentModel()
         a = myModel.rootAssembly
         a.deleteFeatures(featureNames=(self.plane.name,))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def showRefAxis(self):
       if(self.CSYS):
	 highlight(self.CSYS)
         myModel = self.getCurrentModel()
         a = myModel.rootAssembly
         self.axis=a.DatumAxisByThruEdge(edge=self.CSYS.axis3)  
         highlight(self.axis)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def removeRefAxis(self):
       if(self.CSYS):
	 unhighlight(self.CSYS)
         myModel = self.getCurrentModel()
         a = myModel.rootAssembly
         a.deleteFeatures(featureNames=(self.axis.name,))

