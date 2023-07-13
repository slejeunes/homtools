from abaqus import *
from abaqusConstants import *
from math import ceil, fabs, atan, degrees, pi
from numpy import *
import __main__
import random

# VARIABLES GLOBALES
# NOMBRE MAX D ITERATION SUR LES TIRAGES ALEATOIRES
ITEMAX=20
# CONCENTRATION DE CHARGE EN VOLUME
C = 0.2
# RAYON MOYEN DES CHARGES
Rb = 0.05 

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
import os
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getCurrentViewport():

    vpName = session.currentViewportName
    return session.viewports[vpName]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getCurrentModel():

    vpName = session.currentViewportName
    modelName = session.sessionState[vpName]['modelName']
    return mdb.models[modelName]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getNewModel():

    vpName = session.currentViewportName
    modelName = session.sessionState[vpName]['modelName']
    if(len(mdb.models[modelName].parts)==0):
      return mdb.models[modelName]    
    modelName=modelName+"-1"  
    return mdb.Model(modelName)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mean(m,axis=0):
    m = asarray(m)
    return add.reduce(m,axis)/float(m.shape[axis])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def std(m,axis=0):
    x = asarray(m)
    n = float(x.shape[axis])
    mx = asarray(mean(x,axis))
    if axis < 0:
        axis = len(x.shape) + axis
    mx.shape = mx.shape[:axis] + (1,) + mx.shape[axis:]
    x = x - mx
    return sqrt(add.reduce(x*x,axis)/(n-1.0))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def dist(p1,p2):
   return sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_point(p,spheres,DMAX,ep):
  for s in spheres:
        d=dist(s[0],p)-ep
        if(d<s[1] or d-s[1]<2*ep): 
            return False
  return True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test_mask(p,r,DMAX,EP):
   d=dist(r,p)
   if(r[0]-d<DMAX/4+EP/5 or d>r[0]): 
       return False
   return True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def SphereCell(name, mask, width,height,DMAX,NPOINT,EP,REMOVEB=True):
   BOUNDSIZE=(width,height)
   myModel = getCurrentModel()
   sa = myModel.ConstrainedSketch(name='Sketch A', sheetSize=200.0)
   sr = myModel.ConstrainedSketch(name='Sketch Rect', sheetSize=200.0)
   if(mask==1):
     sr.rectangle(point1=(0.0, 0.0), point2=(BOUNDSIZE[0], BOUNDSIZE[1]))
     geomcarac=[BOUNDSIZE[0],BOUNDSIZE[1]]
   elif(mask==2):
     sr.CircleByCenterPerimeter(center=(width*0.5, width*0.5), point1=( width*0.5, width))
     geomcarac=[width*0.5,width*0.5]
   elif(mask==3):
     sr.EllipseByCenterPerimeter(center=(width*0.5,height*0.5), axisPoint1=(width,height*0.5),axisPoint2=( width*0.5, height))
     geomcarac=[width*0.5,height*0.5]
   bound = [(0., 0.), (0., BOUNDSIZE[1]), (BOUNDSIZE[0], BOUNDSIZE[1]), (BOUNDSIZE[0], 0.)]

   # tirage aleatoire avec controle de la distance entre germe
   spheres=[]
   p1=(random.uniform(0+geomcarac[0]*0.25,BOUNDSIZE[0]-geomcarac[0]*0.25),random.uniform(0+geomcarac[0]*0.25,BOUNDSIZE[1]-geomcarac[0]*0.25)) 
   cercle1=(p1,DMAX/2-EP/2)
   spheres.append(cercle1)
   for i in range(0,NPOINT-1):
      p=(random.uniform(0,BOUNDSIZE[0]),random.uniform(0,BOUNDSIZE[1])) 
      iter=0
      test=test_point(p,spheres,DMAX,EP) 
      if(REMOVEB and mask==2): 
          test=test and test_mask(p,geomcarac,DMAX,EP)
      while((not test) and iter<900):             
         p=(random.uniform(0,BOUNDSIZE[0]),random.uniform(0,BOUNDSIZE[1])) 
         test=test_point(p,spheres,DMAX,EP)
         if(REMOVEB and mask==2): 
          test=test and test_mask(p,geomcarac,DMAX,EP)
         iter+=1
      if(iter<900 and test):
         distances={}
         for s in spheres:
            distances[dist(s[0],p)-s[1]]=s
         dists=distances.keys()
         dists.sort()
         r=dists[0]-EP
         if(r>0.5*DMAX):
            r=0.5*DMAX-0.5*EP
         if(REMOVEB and mask==2): 
            if(geomcarac[0]-dist(p,geomcarac)<r):
               r=geomcarac[0]-dist(p,geomcarac)-1.5*EP
         cercle=(p,r)
         spheres.append(cercle)
                        
   # trace du diagramme
   for s in spheres:
      sa.CircleByCenterPerimeter(center=(s[0][0], s[0][1]), point1=( s[0][0]+s[1], s[0][1]))
 

   myModel.convertAllSketches()
   p = myModel.Part(name=name, dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   p.BaseShell(sketch=sr)
   sa.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=sa)

   vp = getCurrentViewport()
   vp.setValues(displayedObject=p)

   a = myModel.rootAssembly
   a.Instance(name=name, part=p, dependent=ON)
   i=0
   vol_charge=0
   rayons=[]
   surfspec=[]
   for s in spheres:
      r1=a.ReferencePoint(point=(s[0][0], s[0][1],0))
      p.DatumPointByCoordinate(coords=(s[0][0], s[0][1], 0.0))
      f1 = a.instances[name].faces.findAt(((s[0][0], s[0][1],0),))
      vol_charge+=f1[0].getSize()
      region2=(f1,)
      ref = a.referencePoints
      refPoints1=(ref[r1.id], )
      region1=regionToolset.Region(referencePoints=refPoints1)
      nom='Constraint'+str(i+1)
      i+=1
      rayons.append(s[1])
      surfspec.append(2/s[1])
      myModel.RigidBody(name=nom, refPointRegion=region1,bodyRegion=region2)

   tf = a.instances[name].faces
   mstot=a.getMassProperties(regions=tf[0:len(tf)])
   print 'Taux volumique de charge :',(vol_charge/mstot['area'])
   print 'nombre de grain:',len(spheres) 
   print 'Rayon moyen',mean(array(rayons))
   print 'Surf Spec moyen',mean(array(surfspec))
   print 'Ecart Type Surf Spec',std(array(surfspec))



