#-----------------------------------------------------------------------
#     Plugin for linear deformation on the boundary
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#     
#-----------------------------------------------------------------------
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class HomogeneBoundary:
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
    def Effective(self,is_smallstrain,dim):
       if(not self.CheckInput()): 
         print "Not all the inputs are given !" 
         return 0

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
       self.MakeSetsandEquations(nameS1,a,dim)
       self.suppressPickedSet(a) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def suppressPickedSet(self,a):
       for key in a.sets.keys():
          if(key.find('Picked_set')!=-1):
            del a.sets[key]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def MakeSetsandEquations(self,s1,a,dim):
       nodes1=a.allSets[s1].nodes
       
       nodesRef=a.allSets["RefMacro1"].referencePoints
       refcoord=a.getCoordinates(nodesRef[0])
       x1_x2=[]
       mod = self.getCurrentModel()
       i=0
       x1_x2=[0.]*3
       for n1 in nodes1:
         name1="Num"+str(n1.label)
         a.Set(nodes=nodes1[i:i+1],name=name1)
         for k in range(0,3):
           x1_x2[k]=n1.coordinates[k]-refcoord[k]  
         if(self.is_small_strain==True):
           if(dim==2):
             mod.Equation(name='Constraint-1-'+name1, terms=((1.0, name1, 1),(-x1_x2[0],"RefMacro1",1),(-x1_x2[1]/2,"RefMacro2",1)))
             mod.Equation(name='Constraint-2-'+name1, terms=((1.0, name1, 2),(-x1_x2[0]/2,"RefMacro2",1),(-x1_x2[1],"RefMacro1",2)))
           else: 
             mod.Equation(name='Constraint-1-'+name1, terms=((1.0, name1, 1),(-x1_x2[0],"RefMacro1",1),(-x1_x2[1]/2,"RefMacro2",1),(-x1_x2[2]/2,"RefMacro2",2)))
             mod.Equation(name='Constraint-2-'+name1, terms=((1.0, name1, 2),(-x1_x2[0]/2,"RefMacro2",1),(-x1_x2[1],"RefMacro1",2),(-x1_x2[2]/2,"RefMacro2",3)))
             mod.Equation(name='Constraint-3-'+name1, terms=((1.0, name1, 3),(-x1_x2[0]/2,"RefMacro2",2),(-x1_x2[1]/2,"RefMacro2",3),(-x1_x2[2],"RefMacro1",3)))
         else:
           if(dim==2):
             mod.Equation(name='Constraint-1-'+name1, terms=((1.0, name1, 1),(-x1_x2[0],"RefMacro1",1),(-x1_x2[1],"RefMacro1",2)))
             mod.Equation(name='Constraint-2-'+name1, terms=((1.0, name1, 2),(-x1_x2[0],"RefMacro2",1),(-x1_x2[1],"RefMacro2",2))) 
           else:
             mod.Equation(name='Constraint-1-'+name1, terms=((1.0, name1, 1),(-x1_x2[0],"RefMacro1",1),(-x1_x2[1],"RefMacro1",2),(-x1_x2[2],"RefMacro1",3)))
             mod.Equation(name='Constraint-2-'+name1, terms=((1.0, name1, 2),(-x1_x2[0],"RefMacro2",1),(-x1_x2[1],"RefMacro2",2),(-x1_x2[2],"RefMacro2",3)))
             mod.Equation(name='Constraint-3-'+name1, terms=((1.0, name1, 3),(-x1_x2[0],"RefMacro3",1),(-x1_x2[1],"RefMacro3",2),(-x1_x2[2],"RefMacro3",3)))
         milestone("Process equation"," ",i,len(nodes1))
         i=i+1

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




