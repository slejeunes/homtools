#-----------------------------------------------------------------------
#     Plugin for periodic condition on the boundary (2D & 3D)
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
from copy import deepcopy



def getCurrentModel():
   vpName = session.currentViewportName
   modelName = session.sessionState[vpName]['modelName']
   return mdb.models[modelName]
def switch(a,b):
   tmp=b
   b=a
   a=tmp
def find(f,seq):
   for item in seq:
      (key,value)=item
      if (value==f): return key
   return -1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      A class to describe a graph of type tree: without any cycles
#         designed to avoid duplicated DOF in constraint equations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Tree:
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self):
       self.thetree=[]
       self.oldtree=[]
       self.keyword=''
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getCurrentGraph(self,keyword): 
       self.keyword=keyword
       mod = getCurrentModel()
       conts=mod.constraints
       for i in range(0,len(conts.keys())):
          name=conts.keys()[i]
          if(name.find(keyword)>-1):
             (title,n1,n2,ddl)=name.split('_')      
             if(int(ddl)==1):
               pos=-1
               for j in range(0,len(self.thetree)):
                  if(n1 in self.thetree[j]):
                    twig=[n1,n2]
                    self.insertTwigInBranch(twig,self.thetree[j])
                    self.reduceTree(twig,self.thetree[j])
                    pos=j                     
                    break
                  elif(n2 in self.thetree[j]): 
                    twig=[n2,n1]
                    self.insertTwigInBranch(twig,self.thetree[j])
                    self.reduceTree(twig,self.thetree[j])
                    pos=j
                    break
               if(pos==-1): self.thetree.append([n1,n2])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def insertTwigInBranch(self,twig,branch):
       if(twig[0] not in branch):
           print 'bug in insert twig! twig=',twig,' branch=',branch
           return False
       ind=branch.index(twig[0])
       if(ind==len(branch)-1): branch.append(twig[1]) 
       elif(ind==0): branch.insert(ind,twig[1])
       else: return False
       return True
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def reduceTree(self,twig,curbranch):
       k=-1
       for branch in self.thetree:
          k=k+1
          if((branch is not curbranch) and (twig[1] in branch) and (twig[1] in curbranch)):
              ind0=curbranch.index(twig[1])
              ind1=branch.index(twig[1])
              if(ind0==0 and ind1==len(branch)-1): 
                  for j in range(ind1-1,-1,-1):
                     curbranch.insert(0,branch[j])
                  del self.thetree[k]
                  k=k-1
#                  return True
              elif(ind0==0 and ind1==0):
                  for j in range(1,len(branch)):
                     curbranch.insert(0,branch[j])                   
                  del self.thetree[k]
                  k=k-1
#                  return True
              elif(ind0==len(curbranch)-1 and ind1==len(branch)-1):
                  for j in range(ind1-1,-1,-1):
                     curbranch.append(branch[j])
                  del self.thetree[k]
                  k=k-1
#                  return True
              elif(ind0==len(curbranch)-1 and ind1==0):
                  for j in range(1,len(branch)):
                     curbranch.append(branch[j])
                  del self.thetree[k]
                  k=k-1
#                  return True
              else: 
                  print 'Tree reduction error! (ind0=',ind0,' ind1=',ind1,')','branch',branch,'curbranch',curbranch,'twig',twig
                  return False
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def addBranch(self,branch): 
       n1=branch[0]  
       n2=branch[1]
       abranch=None
       self.oldtree=deepcopy(self.thetree)
       for curbranch in self.thetree:
         if(n1 in curbranch and n2 in curbranch): 
            return False
         elif(n1 in curbranch and abranch==None):
            twig=[n1,n2]
            if(not self.insertTwigInBranch(twig,curbranch)): return False
            self.reduceTree(twig,curbranch)
            abranch=curbranch
            return True
         elif(n2 in curbranch and abranch==None):
            twig=[n2,n1]
            if(not self.insertTwigInBranch(twig,curbranch)): return False
            self.reduceTree(twig,curbranch)
            abranch=curbranch
            return True
       if(abranch==None): self.thetree.append(branch)
       return True
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def restoreTree(self):
       del self.thetree
       self.thetree=deepcopy(self.oldtree)       
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def iscycle(self):
       for curbranch in self.thetree:
         for node in curbranch:
            if(curbranch.count(node)>1): return True
       return False 
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class PeriodicBoundary:
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,RP1=None,RP2=None,RP3=None,set1=None,set2=None,set3=None):
        self.RP1=RP1
        self.RP2=RP2
        self.RP3=RP3
        self.set1=set1
        self.set2=set2
        self.set3=set3
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
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def Periodic(self,Vx,Vy,Vz,is_smallstrain,dim):
       self.Vx=Vx
       self.Vy=Vy
       self.Vz=Vz
       self.carac=1.e-4
       small=[Vx,Vy,Vz]
       small.sort()
       self.is_small_strain=is_smallstrain
       if(not self.CheckInput()): 
         print "Not all the inputs are given !" 
         return 0

       myModel = self.getCurrentModel()
       a = myModel.rootAssembly
       nameS3=''
       # build all the sets needed
       nameS1=self.build_set(self.set1,a,small[0]*0.00000001)
       nameS2=self.build_set(self.set2,a,small[0]*0.00000001)
       if(self.set3!=None):
         nameS3=self.build_set(self.set3,a,small[0]*0.00000001)
       a.Set(referencePoints=(self.RP1,), name="RefMacro1")
       a.Set(referencePoints=(self.RP2,), name="RefMacro2")
       if(dim==3 and (not is_smallstrain)):
           a.Set(referencePoints=(self.RP3,), name="RefMacro3")
       
       # check if we have a mesh
       if(len(a.allSets[nameS1].nodes)==0 or len(a.allSets[nameS2].nodes)==0):
         print "A mesh is needed!" 
         self.SuppressPickedSet(a) 
         return 0
       # check the mesh periodicity
       if(not self.VerifNodeCoord(nameS1,nameS2,a)):
           # self.SuppressPickedSet(a) 
           return 0
       self.MakeNodeSetsandEquations(nameS1,nameS2,nameS3,a,dim)

       self.SuppressPickedSet(a)    

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def SuppressPickedSet(self,a):
       for key in a.sets.keys():
          if(key.find('Picked_set')!=-1):
            del a.sets[key]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def NodeDist(self,n1,n2): 
        dist=sqrt((n1.coordinates[0]-n2.coordinates[0])**2
             +(n1.coordinates[1]-n2.coordinates[1])**2+(n1.coordinates[2]-n2.coordinates[2])**2)
        return dist
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def MakeSingleEquation(self,mod,cname,name1,name2,x1_x2,dim,is_ref):
       if(not self.is_small_strain):
          if(dim==2):
            mod.Equation(name=cname+'_1', terms=((1.0, name1, 1),(-1.0, name2, 1),(-x1_x2[0],"RefMacro1",1),
                                                 (-x1_x2[1],"RefMacro1",2)))
            mod.Equation(name=cname+'_2', terms=((1.0, name1, 2),(-1.0, name2, 2),(-x1_x2[0],"RefMacro2",1),
                                                 (-x1_x2[1],"RefMacro2",2)))
          else:
            mod.Equation(name=cname+'_1', terms=((1.0, name1, 1),(-1.0, name2, 1),(-x1_x2[0],"RefMacro1",1),
                                                 (-x1_x2[1],"RefMacro1",2),(-x1_x2[2],"RefMacro1",3)))
            mod.Equation(name=cname+'_2', terms=((1.0, name1, 2),(-1.0, name2, 2),(-x1_x2[0],"RefMacro2",1),
                                                 (-x1_x2[1],"RefMacro2",2),(-x1_x2[2],"RefMacro2",3)))
            mod.Equation(name=cname+'_3', terms=((1.0, name1, 3),(-1.0, name2, 3),(-x1_x2[0],"RefMacro3",1),
                                                 (-x1_x2[1],"RefMacro3",2),(-x1_x2[2],"RefMacro3",3)))
       else:
          if(dim==2):
            mod.Equation(name=cname+'_1', terms=((1.0, name1, 1),(-1.0, name2, 1),(-x1_x2[1]*0.5,"RefMacro2",1),
                                                 (-x1_x2[0],"RefMacro1",1)))
            mod.Equation(name=cname+'_2', terms=((1.0, name1, 2),(-1.0, name2, 2),(-x1_x2[0]*0.5,"RefMacro2",1),
                                                 (-x1_x2[1],"RefMacro1",2)))
          else:
            mod.Equation(name=cname+'_1', terms=((1.0, name1, 1),(-1.0, name2, 1),(-x1_x2[0],"RefMacro1",1), 
                                                 (-x1_x2[1]*0.5,"RefMacro2",1),(-x1_x2[2]*0.5,"RefMacro2",2)))
            mod.Equation(name=cname+'_2', terms=((1.0, name1, 2),(-1.0, name2, 2),(-x1_x2[0]*0.5,"RefMacro2",1),
                                                 (-x1_x2[1],"RefMacro1",2),(-x1_x2[2]*0.5,"RefMacro2",3)))
            mod.Equation(name=cname+'_3', terms=((1.0, name1, 3),(-1.0, name2, 3),(-x1_x2[0]*0.5,"RefMacro2",2),
                                                 (-x1_x2[1]*0.5,"RefMacro2",3),(-x1_x2[2],"RefMacro1",3)))
       if(is_ref): mod.Equation(name=cname+'_6', terms=((1.0, name1, 6),(-1.0, name2, 6)))
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def MakeNodeSetsandEquations(self,s1,s2,s3,a,dim):
       nodes1=a.allSets[s1].nodes
       nodes2=a.allSets[s2].nodes

       nodesref1=a.allSets[s1].referencePoints
       nodesref2=a.allSets[s2].referencePoints

       rem=[0.]*len(nodes1)
       nodes2elem=range(len(nodes2))
       if(s3!=""):
         for n3 in a.allSets[s3].nodes:
           j=0
           for n1 in nodes1:
             if(self.NodeDist(n1,n3)<1.e-8): 
                 rem[j]=1
             j=j+1  

       x1_x2=[0.]*3
       mod = self.getCurrentModel()
       i=0

       mytree=Tree()
       mytree.getCurrentGraph('Const_')
       for n1 in nodes1:
         name1="Num"+str(n1.label)+n1.instanceName
         if(rem[i]==0):
           if(a.sets.has_key(name1)==0): a.Set(nodes=nodes1[i:i+1],name=name1)
           for j in nodes2elem:
             c2=(nodes2[j].coordinates[0]-self.Vx,nodes2[j].coordinates[1]-self.Vy,nodes2[j].coordinates[2]-self.Vz)
             dist=sqrt((n1.coordinates[0]-c2[0])**2+(n1.coordinates[1]-c2[1])**2+(n1.coordinates[2]-c2[2])**2)
             if(dist<0.1*self.carac):  
               name2="Num"+str(nodes2[j].label)+nodes2[j].instanceName
               nname2=str(nodes2[j].label)+nodes2[j].instanceName
               a.Set(nodes=nodes2[j:j+1],name=name2)
               x1_x2[0]=n1.coordinates[0]-nodes2[j].coordinates[0]
               x1_x2[1]=n1.coordinates[1]-nodes2[j].coordinates[1]
               x1_x2[2]=n1.coordinates[2]-nodes2[j].coordinates[2]
               nodes2elem.remove(j)
               break
           # Tree construction and checking of cycling DOF's
           branch=[str(n1.label)+n1.instanceName,nname2]
           if(mytree.addBranch(branch) and not mytree.iscycle()):
             cname='Const_'+str(n1.label)+n1.instanceName+'_'+str(nname2)
             self.MakeSingleEquation(mod,cname,name1,name2,x1_x2,dim,False)
           else:
             mytree.restoreTree()
         i=i+1
       # checking of DOF ordering to avoid undesired DOF suppression...
       self.checkMasterAndSlaves(mytree)

       #Same treatement for RefPoints 
       myreftree=Tree()
       myreftree.getCurrentGraph('ConstRef_')
       nodes2elem=range(len(nodesref2))

       for n1 in nodesref1:
         name1="Ref_"+str(len(a.sets))
         nname1=find(n1,a.referencePoints.items())
         n1coord=a.getCoordinates(n1) 
         if(a.sets.has_key(name1)==0): a.Set(referencePoints=(n1,),name=name1)
         for j in nodes2elem:
           n2coord=a.getCoordinates(nodesref2[j])
           c2=(n2coord[0]-self.Vx,n2coord[1]-self.Vy,n2coord[2]-self.Vz)
           dist=sqrt((n1coord[0]-c2[0])**2+(n1coord[1]-c2[1])**2+(n1coord[2]-c2[2])**2)
           if(dist<0.1*self.carac):  
             name2="Ref_"+str(len(a.sets))
             nname2=find(nodesref2[j],a.referencePoints.items())
             a.Set(referencePoints=(nodesref2[j],),name=name2) 
             x1_x2[0]=n1coord[0]-n2coord[0]
             x1_x2[1]=n1coord[1]-n2coord[1]
             x1_x2[2]=n1coord[2]-n2coord[2]
             nodes2elem.remove(j)
             break
           # Tree construction and checking of cycling DOF's
         branch=[int(nname1),int(nname2)]
         if(myreftree.addBranch(branch) and not myreftree.iscycle()):
             cname='ConstRef_'+str(nname1)+'_'+str(nname2)
             self.MakeSingleEquation(mod,cname,name1,name2,x1_x2,dim,True)
         else:
             myreftree.restoreTree()
       self.checkMasterAndSlaves(myreftree)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def checkMasterAndSlaves(self,the_tree):
       mod = getCurrentModel()
       conts=mod.constraints
       constname=conts.keys()
       for branch in the_tree.thetree:
          for i in range(1,len(branch)):
             cname1=the_tree.keyword+str(branch[i])+'_'+str(branch[i-1])+'_1'
             cname2=the_tree.keyword+str(branch[i])+'_'+str(branch[i-1])+'_2'
             cname3=the_tree.keyword+str(branch[i])+'_'+str(branch[i-1])+'_3'
             cname6=the_tree.keyword+str(branch[i])+'_'+str(branch[i-1])+'_6'
             
             if(cname1 in constname):
                  theterms=conts[cname1].terms
                  if(len(theterms)==4): newterms=(theterms[1],theterms[0],theterms[2],theterms[3])
                  elif(len(theterms)>4): newterms=(theterms[1],theterms[0])+theterms[2:len(theterms)] 
                  else: newterms=(theterms[1],theterms[0],theterms[2],theterms[3],theterms[4])
                  conts[cname1].setValues(newterms)   
                  conts.changeKey(fromName=cname1,toName=the_tree.keyword+str(branch[i-1])+'_'+str(branch[i])+'_1')
             if(cname2 in constname):
                  theterms=conts[cname2].terms
                  if(len(theterms)==4): newterms=(theterms[1],theterms[0],theterms[2],theterms[3])
                  elif(len(theterms)>4): newterms=(theterms[1],theterms[0])+theterms[2:len(theterms)] 
                  else: newterms=(theterms[1],theterms[0],theterms[2],theterms[3],theterms[4])
                  conts[cname2].setValues(newterms)   
                  conts.changeKey(fromName=cname2,toName=the_tree.keyword+str(branch[i-1])+'_'+str(branch[i])+'_2')
             if(cname3 in constname):
                  theterms=conts[cname3].terms
                  if(len(theterms)==4): newterms=(theterms[1],theterms[0],theterms[2],theterms[3])
                  elif(len(theterms)>4): newterms=(theterms[1],theterms[0])+theterms[2:len(theterms)]  
                  else: newterms=(theterms[1],theterms[0],theterms[2],theterms[3],theterms[4])
                  conts[cname3].setValues(newterms)   
                  conts.changeKey(fromName=cname3,toName=the_tree.keyword+str(branch[i-1])+'_'+str(branch[i])+'_3')
             if(cname6 in constname):
                  theterms=conts[cname6].terms
                  newterms=(theterms[1],theterms[0])
                  conts[cname6].setValues(newterms)   
                  conts.changeKey(fromName=cname6,toName=the_tree.keyword+str(branch[i-1])+'_'+str(branch[i])+'_6')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def GetCarcLength(self,s1,a):
         small=1.e8
         el1=a.allSets[s1].elements        
         for el in el1:
            tabed=el.getElemEdges()
            for ed in tabed:
               nds=ed.getNodes()
               x1=nds[0].coordinates
               x2=nds[1].coordinates
               dist=sqrt((x1[0]-x2[0])**2.+(x1[1]-x2[1])**2.+(x1[2]-x2[2])**2.)
               if(dist>0 and small>dist): small=dist
         return small
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def VerifNodeCoord(self,s1,s2,a):
         nodes1=a.allSets[s1].nodes
         nodes2=a.allSets[s2].nodes
         self.carac=self.GetCarcLength(s1,a)
         nodes2elem=range(len(nodes2))
         for n1 in nodes1:
           found=False
           for j in nodes2elem:
             c2=(nodes2[j].coordinates[0]-self.Vx,nodes2[j].coordinates[1]-self.Vy,nodes2[j].coordinates[2]-self.Vz)
             dist=sqrt((n1.coordinates[0]-c2[0])**2.+(n1.coordinates[1]-c2[1])**2.+(n1.coordinates[2]-c2[2])**2.)
             if(dist<0.1*self.carac):  
              found=True
              nodes2elem.remove(j)
              break
           if(not found):
             print "The mesh periodicity is required "
             return False
         return True  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def VerifRefCoord(self,s1,s2,a):
         nodes1=a.allSets[s1].referencePoints
         nodes2=a.allSets[s2].referencePoints
         nodes2elem=range(len(nodes2))

         for n1 in nodes1:
           found=False
           for j in nodes2elem:
             n2coord=a.getCoordinates(nodes2[j])
             n1coord=a.getCoordinates(n1)
             c2=(n2coord[0]-self.Vx,n2coord[1]-self.Vy,n2coord[2]-self.Vz)
             dist=sqrt((n1coord[0]-c2[0])**2+(n1coord[1]-c2[1])**2+(n1coord[2]-c2[2])**2)
             if(dist<0.1*self.carac):  
              found=True
              nodes2elem.remove(j)
              break
           if(not found):
             print "The mesh periodicity is required"
             return False
         return True       
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def CheckInput(self):
       if(self.RP1==None): return False
       elif(self.set1==None): return False
       elif(self.set2==None): return False
       return True  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setRefpoint1(self,RP1): 
       self.RP1=RP1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setRefpoint2(self,RP2): 
       self.RP2=RP2
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setRefpoint3(self,RP3): 
       self.RP3=RP3
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setGroup1(self,set1):
          self.set1=set1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setGroup2(self,set2): 
          self.set2=set2
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def setGroup3(self,set3): 
       self.set3=set3
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getCurrentModel(self):
      vpName = session.currentViewportName
      modelName = session.sessionState[vpName]['modelName']
      return mdb.models[modelName]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getCurrentViewport(self):
      vpName = session.currentViewportName
      return session.viewports[vpName]

