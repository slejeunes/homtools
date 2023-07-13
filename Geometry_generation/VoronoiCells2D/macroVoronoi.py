#-----------------------------------------------------------------------
#     Plugin pour la generation d'un schema Voronoi avec epaisseur
#     entre les cellules
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#     Rev: 14/06/10  
#-----------------------------------------------------------------------
from abaqus import *
from abaqusConstants import *
#from Numeric import *
from math import ceil, fabs, atan, degrees, pi, asin, cos, sin, radians
import __main__
from voronoi import *
import random

# Iterations max sur les tirages aleatoires
ITEMAX=20


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
# Construction d'un schema de voronoi aleatoire avec controle de la distance
# mini entre grain (le schemas peut-etre periodique ou non...)
def constrainedVoronoi(bound,npoint,DMIN,PERIODIC,radius,center):
   nbound=bound
   if(PERIODIC):
      nbound = [(-bound[2][0],-bound[2][1]), (-bound[2][0],2.*bound[2][1]), (2.*bound[2][0],2.*bound[2][1]), (2.*bound[2][0],-bound[2][1])]
   points = []
   voronoi = bounded_voronoi(nbound, points)
   ite = 0
   for i in range(1,npoint+1):
     if(len(bound)<=4):
        newpoint = (random.uniform(0,bound[2][0]), random.uniform(0,bound[2][1]))
     else:
        r=random.uniform(0,radius)
        angle= random.uniform(0,radians(360))
        newpoint = (center[0]+r*cos(angle),center[1]+r*sin(angle))
     j = 0
     if i>1:
       ite = 0
       while ((j < len(points)) and (ite < ITEMAX)):
        d = dist(points[j],newpoint)
        j = j+1
        if (d < DMIN):
           if(len(bound)<=4):
              newpoint = (random.uniform(0,bound[2][0]), random.uniform(0,bound[2][1]))
           else:
              r=random.uniform(0,radius)
              angle= random.uniform(0,radians(360))
              newpoint = (center[0]+r*cos(angle),center[1]+r*sin(angle))
           j = 0
           ite = ite + 1
     if (ite>=ITEMAX): 
        npoint=npoint-1
     else:
       if(PERIODIC):
        for ix in (-1,1,0):
         for iy in (1,-1,0):   
           pointadd=(newpoint[0]+float(ix)*bound[2][0],newpoint[1]+float(iy)*bound[2][1])
           points.append(pointadd)
           update_diagram(voronoi, pointadd, nbound)
           voronoi = bounded_voronoi(nbound, points)
       else:
        points.append(newpoint)
        update_diagram(voronoi, newpoint, nbound)
        voronoi = bounded_voronoi(nbound, points)     
   return points,voronoi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def point_on_bound(pt,bound):
   for i in range(0,len(bound)-1):
      seg1=(bound[i],bound[i+1])
      a1=slope(pt,bound[i+1])
      a2=slope(bound[i],pt)
      if(a1==a2 and point_on_segment(pt,seg1)):
         return True
   seg1=(bound[len(bound)-1],bound[0])
   a1=slope(pt,bound[len(bound)-1])
   a2=slope(bound[0],pt)
   if(a1==a2 and point_on_segment(pt,seg1)):
      return True
   return False      
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ajoute une epaisseur au diagramme de Voronoi en supposant
# que toute les cellules ont le meme vecteur normal !
def addThickness(voronoi,bound,points,EP,s):
   refpoint=[]
   cellpoint=[]
   for i in range(0,len(points)):   
     poly = voronoi[points[i]] 
     lg=[]
     for j in range(0,len(poly)-1):
       if (point_on_bound(poly[j],bound) and point_on_bound(poly[j+1],bound)):
        scalar=0.0
       else:
        scalar=1.0
       pointmil=midpoint(poly[j],poly[j+1])
       vect=(poly[j+1][0]-poly[j][0], poly[j+1][1]-poly[j][1])
       vect=normalize(vect)
       ml1=slope(poly[j],poly[j+1]) 
       ptl=(pointmil[0]+scalar*(EP/2.0)*vect[1],pointmil[1]-scalar*(EP/2.0)*vect[0] )     
       lig=(ptl,ml1)
       lg.append(lig)

     if(point_on_bound(poly[len(poly)-1],bound) and point_on_bound(poly[0],bound)):
      scalar=0.0
     else:
      scalar=1.0
     pointmil=midpoint(poly[len(poly)-1],poly[0])
     vect=(poly[len(poly)-1][0]-poly[0][0], poly[len(poly)-1][1]-poly[0][1]) 
     vect=normalize(vect)
     ml1=slope(poly[len(poly)-1],poly[0]) 
     ptl=(pointmil[0]-scalar*(EP/2.0)*vect[1], pointmil[1]+scalar*(EP/2.0)*vect[0])
     lig=(ptl,ml1)
     lg.append(lig)
   
     newpointj=[]
     for j in range(0,len(lg)-1):
       ap=line_intersect(lg[j], lg[j+1])
       if(ap!=None): newpointj.append(ap)
     ap=line_intersect(lg[len(lg)-1], lg[0])
     if(ap!=None): newpointj.append(ap)
     newpointj.append(ap)
     refflag=True 
     for j in range(0,len(newpointj)-1):
#         lig=(newpointj[j],newpointj[j+1])
#         if(lineInBound(lig,bound)):
#            refflag=True 
#         if(pointsInBound(newpointj[j],bound)):
#           refflag=True
         s.Line(newpointj[j],newpointj[j+1])
#     if(pointsInBound(newpointj[len(newpointj)-1],bound)):
#         refflag=True
#     lig=(newpointj[len(newpointj)-1],newpointj[0])
#     if(lineInBound(lig,bound)):
#       refflag=True 
     s.Line(newpointj[len(newpointj)-1],newpointj[0])    
     if(refflag): 
       refpoint.append(points[i]) 
       cellpoint.append(newpointj)
   return refpoint,cellpoint

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def drawVoronoi(voronoi,bound,points,s):
   for i in range(0,len(points)):   
     poly = voronoi[points[i]] 
     if(len(poly)>0):
       for j in range(0,len(poly)-1):
         s.Line(poly[j+1],poly[j])
       s.Line(poly[len(poly)-1],poly[0])  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pointsInBound(point,bound):
  if((point[0]>=bound[0][0]) and (point[0]<=bound[2][0]) and (point[1]>=bound[0][1]) and (point[1]<=bound[1][1])):
    return True
  else:
    return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def lineInBound(line,bound):
   point1=line[0]
   point2=line[1]  
   if(pointsInBound(point1,bound)): return True
   if(pointsInBound(point2,bound)): return True
   ml1=slope(point1,point2) 
   l1=(point1, ml1)
   for i in (0,3):
      if(i<3):
        seg=(bound[i],bound[i+1])
      else:
        seg=(bound[i],bound[0]) 
      pt3=segment_intersect_line(seg,l1)
      if(pt3!= None):
        if point_on_segment(pt3, line):
           return True
   return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def removePoints(voronoi,points,thepoint,eps):
  for i in range(0,len(points)):   
    poly = voronoi[points[i]]    
    npoly=len(poly)
    j=0
    while (j<npoly):
      d=dist(poly[j],thepoint)
      if(d<eps):
         voronoi[points[i]][j]=thepoint
      j=j+1
   # suppression des redondances
  for i in range(0,len(points)):   
    poly = voronoi[points[i]]    
    npoly=len(poly)
    j=0
    while (j<npoly):
      d=dist(poly[j],thepoint)
      if(j==len(poly)-1): 
        d=dist(poly[j],poly[0])
      else: 
         d=dist(poly[j+1],poly[j])
      if(d<1.e-10):  
         voronoi[points[i]].remove(poly[j])    
         npoly=npoly-1
      j=j+1
  return voronoi   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def verifVoronoi(voronoi,points,eps):
   for i in range(0,len(points)):   
     poly = voronoi[points[i]]    
     npoly=len(poly)
     j=0
     while (j<npoly):
       if(j==npoly-1): 
         d=dist(poly[j],poly[0])
       else: 
         d=dist(poly[j+1],poly[j])
       if(d<eps):
         thepoint=poly[j]
         voronoi[points[i]].remove(poly[j])   
         voronoi=removePoints(voronoi,points,thepoint,eps)   
         poly = voronoi[points[i]]      
         npoly=len(poly)
       j=j+1
   return voronoi 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def pointsinlist(pt,liste):
  for i in range(0,len(liste)):
      if(dist(pt,liste[i])<1.e-4): return True
  return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# recherche le segment partagant le noeud prev dans une liste
def findnextseg(prev,seglist):
   for i in range(0,len(seglist)):
      seg=seglist[i]
      if(dist(prev,seg[0])<1.e-4):
        next=seg[1]   
        return True,next,i
   for i in range(0,len(seglist)):
      seg=seglist[i]
      if(dist(prev,seg[1])<1.e-4):
        next=seg[0]   
        return True,next,i
   return False,None,0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def verifMask(voronoi,points,mask,geomcarac):
   if(mask==1):   # cas du carre rien a faire par defaut
     return voronoi  
   elif(mask==2): # cas du cercle
     r=geomcarac[0]
     center=[r,r]
     # suppression des cellules qui ont un point ou plusieurs points sur la frontiere
     i=0
     numcell=[]
     npoint=len(points)
     while (i<npoint):
       int_point=0  
       poly=voronoi[points[i]]
       for j in range(0,len(poly)):       
          d=dist(poly[j],center)
          if(d>=r or fabs(d-r)<=1.e-4): int_point=int_point+1
       if(int_point>=1): 
         numcell.append(i)
       i=i+1
     pt=len(numcell)
     for i in range(0,len(numcell)):  
       del voronoi[points[numcell[pt-1-i]]] 
       del points[numcell[pt-1-i]]

     # cherche les segments qui forment le contour i.e les 
     # segments qui n'ont pas leur deux noeuds partages avec une autre cellule
     seglist=[] 
     for i in range(0,len(points)): 
       polyi=voronoi[points[i]]
       for j in range(0,len(polyi)-1):
         seg1=(polyi[j+1],polyi[j])
         cont=0
         for k in range(0,len(points)): 
           if(i!=k):
             polyk=voronoi[points[k]]
             if(pointsinlist(seg1[0],polyk) and  pointsinlist(seg1[1],polyk)): cont=cont+1
         if(cont<1): seglist.append(seg1)
       seg1=(polyi[len(polyi)-1],polyi[0])
       cont=0
       for k in range(0,len(points)): 
         if(i!=k):
           polyk=voronoi[points[k]]
           if(pointsinlist(seg1[0],polyk) and  pointsinlist(seg1[1],polyk)): cont=cont+1
       if(cont<1): seglist.append(seg1)

     # on range les segments dans l'ordre pour faire une nouvelle cellule de voronoi
     if(len(seglist)>2):
       extcell=[]
       seg=seglist[0]
       extcell.append(seg[0]) 
       extcell.append(seg[1])
       del seglist[0]
       prev=seg[1]
       next=seg[0]
       i=0
       while True:
         test,next,i=findnextseg(prev,seglist)
         if(test==False or len(seglist)<=1): break  
         extcell.append(next)
         prev=next   
         del seglist[i]

       # on ajoute le nouveau polygone 
       center=(r,r)
       points.append(center)
       voronoi[center]=extcell
        
     return voronoi
   elif(mask==3): # cas de l'ellipse (a ecrire ...)
     return voronoi           

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# permet de recuperer un point appartenant a l'interphase...
def getInteriorpoint(voronoi,points,BOUNDSIZE,mask):
   getpt=False
   center=(BOUNDSIZE[0]*0.5,BOUNDSIZE[1]*0.5) 
   r=BOUNDSIZE[0]*0.5
   intpt=center
   for i in range(0,len(points)): 
     poly=voronoi[points[i]]
     if(getpt): break
     for pt in poly:
        if(dist(center,pt)<r):
           intpt=pt 
           getpt=True
           break
   return intpt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generation d'un modele complet avec choix du type de geometrie de base
# rectangle, cercle, ellipse de l'epaisseur et du nombre de cellules

def VoronoiCell(name, mask, width,height,DMIN,NPOINT,EP,PERIODIC=False):
   BOUNDSIZE=(width,height)
   myModel = getCurrentModel()
   s = myModel.ConstrainedSketch(name='Sketch A', sheetSize=200.0)
   s0 = myModel.ConstrainedSketch(name='Sketch Ini', sheetSize=200.0)
   sr = myModel.ConstrainedSketch(name='Sketch Rect', sheetSize=200.0)
   if(mask==1):
     sr.rectangle(point1=(0.0, 0.0), point2=(BOUNDSIZE[0], BOUNDSIZE[1]))
     geomcarac=[BOUNDSIZE[0],BOUNDSIZE[1]]
   elif(mask==2):
     sr.CircleByCenterPerimeter(center=(width*0.5, width*0.5), point1=( width*0.5, width))
     geomcarac=[width*0.5]
   elif(mask==3):
     sr.EllipseByCenterPerimeter(center=(width*0.5,height*0.5), axisPoint1=(width,height*0.5),axisPoint2=( width*0.5, height))
     geomcarac=[width*0.5,height*0.5]
   bound = [(0., 0.), (0., BOUNDSIZE[1]), (BOUNDSIZE[0], BOUNDSIZE[1]), (BOUNDSIZE[0], 0.)]

   # tirage aleatoire avec controle de la distance entre germe
   center=(BOUNDSIZE[0]/2.,BOUNDSIZE[0]/2.)
   points,voronoi = constrainedVoronoi(bound,NPOINT,DMIN,PERIODIC,BOUNDSIZE[1]/2,center) 
   # supression des cellules qui n'ont pas leurs germes dans le masque   
   voronoi = verifMask(voronoi,points,mask,geomcarac)  
   # supression des sommets proches
   voronoi = verifVoronoi(voronoi,points,1.5*EP)

   # trace du diagramme
   drawVoronoi(voronoi,bound,points,s0) 
   myModel.convertAllSketches()
   p = myModel.Part(name="Initial Voronoi", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   p.BaseShell(sketch=sr)
   s.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=s0)
   pinit = p
   # recuperation d'un sommet pour la selection de la gomme 
   intpt=getInteriorpoint(voronoi,points,BOUNDSIZE,mask)
   # ajout d'une epaisseur sur le motif
   refpoint,cellpoint = addThickness(voronoi,bound,points,EP,s)

   myModel.convertAllSketches()
   p = myModel.Part(name=name, dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   p.BaseShell(sketch=sr)
   s.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
   for i in range(0,len(refpoint)):   
     xd=refpoint[i][0]
     yd=refpoint[i][1]
     p.DatumPointByCoordinate(coords=(xd, yd, 0.0))
   vp = getCurrentViewport()
   vp.setValues(displayedObject=p)

   a = myModel.rootAssembly
   a.Instance(name=name, part=p, dependent=ON)
   dat = a.instances[name].datums

   vp.assemblyDisplay.setValues(interactions=ON,constraints=ON, connectors=ON, engineeringFeatures=ON)
   vp.assemblyDisplay.setValues(renderStyle=WIREFRAME)
   InclusionFace=[]
#   for i in range(0,len(dat)):
#     r1=a.ReferencePoint(point=(refpoint[i][0],refpoint[i][1],0))
#     lpoints=[]
#     for j in range(0,len(cellpoint[i])):
#         apoints=(cellpoint[i][j][0],cellpoint[i][j][1],0.0)
#         lpoints.append(apoints)
#     center=[0,0]
#     for pts in voronoi[refpoint[i]]:
#          center[0]=pts[0]+center[0]
#          center[1]=pts[1]+center[1]
#     if(len(voronoi[refpoint[i]])>0):
#       center[0]=center[0]/len(voronoi[refpoint[i]])
#       center[1]=center[1]/len(voronoi[refpoint[i]])
#     f1 = a.instances[name].faces.findAt(((center[0],center[1],0),))
#     region2=(f1,)
#     ref = a.referencePoints
#     refPoints1=(ref[r1.id], )
#     region1=regionToolset.Region(referencePoints=refPoints1)
#     nom='Constraint'+str(i+1)
#     myModel.RigidBody(name=nom, refPointRegion=region1,bodyRegion=region2)

#   gommeface = a.instances[name].faces.findAt(((intpt[0],intpt[1],0),))
#   gommeregion=(gommeface,)
#   a.Set(faces=[gommeface], name='Gomme')
#   ms=a.getMassProperties(regions=gommeregion)
#   tf = a.instances[name].faces
#   mstot=a.getMassProperties(regions=tf[0:len(tf)])
   
#   for fc in tf:
#      if(fc.index!=gommeface[0].index):  
#        InclusionFace.append(tf[fc.index:fc.index+1]) 
#   a.Set(faces=tuple(InclusionFace), name='Inclusion')


 #  if(mask==1):
 #     perim=-(2.*BOUNDSIZE[1]+2.*BOUNDSIZE[0])
 #  elif(mask==2):
 #     perim=-2.*pi*(width*0.5) 
 #  elif(mask==3):
 #     perim=-pi*(3.*0.5*(width+height)-sqrt(0.25*(3.*width+height)*(width+3.*height)))
    
#   for i in range(0,len(p.edges)):
#      perim=perim+p.edges[i].getSize(False)   

#   print 'Taux volumique de charge :',1-(ms['area']/mstot['area'])
#   print 'Surface specifique d echange :',perim
#   print 'nombre de grain:',len(points) 
#  makedidDiagram(pinit,BOUNDSIZE,a,name)
   


def angleceil(x): return (ceil(x[0]),x[1])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def makedidDiagram(part,BOUNDSIZE,assembly=None,name=None):
   f=open('did.txt','w')
   L=[]
   for i in range(0,len(part.edges)):
      siz=part.edges[i].getSize(False)   
      pts=part.edges[i].getVertices()
      coord1=part.vertices[pts[0]].pointOn[0]
      coord2=part.vertices[pts[1]].pointOn[0]
      angle=-200         
      # si on est vertical et sur un bord
      if(coord1[0]==coord2[0]) and (coord2[0]==0 or coord2[0]==BOUNDSIZE[0]):
        angle=-1000
      # si on est vertical et pas sur un bord
      if(coord1[0]==coord2[0]) and not (coord2[0]==0 or coord2[0]==BOUNDSIZE[0]):
        angle=90
      # si on est horizontal et sur un bord 
      if(coord1[1]==coord2[1]) and (coord2[1]==0 or coord2[1]==BOUNDSIZE[1]):
        angle=-1000
      # si on est horizontal et pas sur un bord 
      if(coord1[1]==coord2[1]) and not (coord2[1]==0 or coord2[1]==BOUNDSIZE[1]):
        angle=0
      if(angle==-200):     
         angle=degrees(atan((coord2[1]-coord1[1])/(coord2[0]-coord1[0])))
      if(angle!=-1000):
         if(angle<0): angle=angle+180;
         elem=(angle,siz)
         L.append(elem)


   L=map(angleceil,L)
   L.sort() 

   for i in range(len(L)-1,0,-1):
      anglei,sizi=L[i]
      anglej,sizj=L[i-1]
      if(anglei==anglej):
        sizi=sizj+sizi
        L[i-1]=(anglei,sizi)  
        L.pop(i) 
   for i in range(0,len(L)):
      angle,siz=L[i]
      f.write(str(angle)+'   '+str(siz)+"\n")

   f.close()

   if(assembly!=None):
      f=open('vgz.txt','w')
      f.write("No grain      Surface\n 0.0        0.0\n")
      tf=assembly.instances[name].faces
      mstot=assembly.getMassProperties(regions=tf[0:len(tf)])
      vol=[]
      for i in range(0,len(tf)-1):    
         ms=assembly.getMassProperties(regions=tf[i:i+1])
         vol.append((round(ms['area'],3),1))

      vol.sort()   
      for i in range(len(vol)-1,0,-1):  
         szi,nbi=vol[i]
         szj,nbj=vol[i-1]
         if(fabs(szi-szj)<0.0001):
            nbi=nbi+nbj
            vol[i-1]=(szi,nbi)
            vol.pop(i)
      for i in range(0,len(vol)):
          szi,nbi=vol[i]
          f.write(str(nbi)+'  '+str(szi)+"\n")
      f.write('0      '+str(mstot['area'])+"\n")
      f.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def InterpolEp(listEp,FRATIO):
  if len(listEp)<2:
    x1,y1=listEp[0]
    EP=y1*(1.+1.5*(x1-FRATIO)/FRATIO)
    return EP
  if len(listEp)==2:
    # interpolation lineaire
    x1,y1=listEp[0]
    x2,y2=listEp[1]
    a=(y2-y1)/(x2-x1)      
    b=y1-a*x1
    EP=a*FRATIO+b
    return EP
  else:
    # on s'arrete au quadratique
    x1,y1=listEp[len(listEp)-3]
    x2,y2=listEp[len(listEp)-2]
    x3,y3=listEp[len(listEp)-1]    
    L1=(FRATIO-x2)*(FRATIO-x3)/((x1-x2)*(x1-x3))
    L2=(FRATIO-x1)*(FRATIO-x3)/((x2-x1)*(x2-x3))
    L3=(FRATIO-x1)*(FRATIO-x2)/((x3-x1)*(x3-x2))
    EP=L1*y1+L2*y2+L3*y3
    return EP

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generation d'un modele complet avec choix du type de geometrie de base
# rectangle, cercle, ellipse du taux de charge et du nombre initiale cellules

def VoronoiCell2(name, mask,width,height,DMIN,NPOINT,FRATIO,filename="",PERIODIC=False):
   BOUNDSIZE=(width,height)
   myModel = getNewModel()
   s = myModel.ConstrainedSketch(name='Sketch A', sheetSize=200.0)
   s2 = myModel.ConstrainedSketch(name='Sketch B', sheetSize=200.0)
   sr = myModel.ConstrainedSketch(name='Sketch Rect', sheetSize=200.0)
   if(mask==1):
     sr.rectangle(point1=(0.0, 0.0), point2=(BOUNDSIZE[0], BOUNDSIZE[1]))
     geomcarac=[BOUNDSIZE[0],BOUNDSIZE[1]]
   elif(mask==2):
     sr.CircleByCenterPerimeter(center=(width*0.5, width*0.5), point1=( width*0.5, width))
     geomcarac=[width*0.5]
   elif(mask==3):
     sr.EllipseByCenterPerimeter(center=(width*0.5,height*0.5), axisPoint1=(width,height*0.5),axisPoint2=( width*0.5, height))
     geomcarac=[width*0.5,height*0.5]
   bound = [(0., 0.), (0., BOUNDSIZE[1]), (BOUNDSIZE[0], BOUNDSIZE[1]), (BOUNDSIZE[0], 0.)]
   if (filename==""):
     # tirage aleatoire avec controle de la distance entre germe
     center=(BOUNDSIZE[0]/2.,BOUNDSIZE[0]/2.)
     points,voronoi = constrainedVoronoi(bound,NPOINT,DMIN,PERIODIC,BOUNDSIZE[1]/2,center) 
     # supression des cellules qui n'ont pas leurs germes dans le masque   
     voronoi = verifMask(voronoi,points,mask,geomcarac)  
   else:
     points=[]
     prevmodel=mdb.openAuxMdb(filename)  
     if(not prevmodel):
        print "Impossible d'ouvrir ", filename
        return
     oldname=mdb.getAuxMdbModelNames()
     mdb.copyAuxMdbModel(oldname[0],"OldModel")
     prevmodel=mdb.models['OldModel']
     for k, pa in prevmodel.parts.items():      
        if(len(pa.datums)>1):
          for k, da in pa.datums.items():
            apoint=(da.pointOn[0],da.pointOn[1])
            points.append(apoint)
          break 
     voronoi = bounded_voronoi(bound, points)

   listEp=[]
   points_ini=points
   voronoi_ini=voronoi
   drawVoronoi(voronoi,bound,points,s) 
 
   myModel.convertAllSketches()
   p = myModel.Part(name="Initial Voronoi", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   pinit = p
   p.BaseShell(sketch=sr)
   s.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)

   perim=-(2.*BOUNDSIZE[1]+2.*BOUNDSIZE[0])
   for i in range(0,len(p.edges)):
      perim=perim+p.edges[i].getSize(False)   
   EP=(1-FRATIO)/perim
   voronoi = verifVoronoi(voronoi,points,0.90*EP)

   # recuperation d'un sommet pour la selection de la gomme 
   intpt=getInteriorpoint(voronoi,points,BOUNDSIZE,mask)
   refpoint,cellpoint = addThickness(voronoi,bound,points,EP,s2)
   myModel.convertAllSketches()
   p = myModel.Part(name=name+str(0), dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   p.BaseShell(sketch=sr)
   s2.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=s2)

   a = myModel.rootAssembly
   a.Instance(name=name+str(0), part=p, dependent=ON)
   tf = a.instances[name+str(0)].faces
   gommeface = a.instances[name+str(0)].faces.findAt(((intpt[0],intpt[1],0),))
   gommeregion=(gommeface,)
   a.Set(faces=[gommeface], name='Gomme')
   ms=a.getMassProperties(regions=gommeregion)
   mstot=a.getMassProperties(regions=tf[0:len(tf)])
   cr=1.-(ms['area']/mstot['area'])

   itera=0
   listEp.append((cr,EP))

   while (abs(cr-FRATIO)/FRATIO>0.005 and itera<8):
     a.suppressFeatures(featureNames=(name+str(itera),))
     itera=itera+1
     print "iteration ",itera
     voronoi=voronoi_ini
     points=points_ini
     EP=InterpolEp(listEp,FRATIO)  
     s2 = myModel.ConstrainedSketch(name='Sketch B'+str(itera), sheetSize=200.0)
     voronoi = verifVoronoi(voronoi,points,0.90*EP)
     refpoint,cellpoint = addThickness(voronoi,bound,points,EP,s2)    
     myModel.convertAllSketches()
     p = myModel.Part(name=name+str(itera), dimensionality=TWO_D_PLANAR,  type=DEFORMABLE_BODY)
     p.BaseShell(sketch=sr)
     s2.unsetPrimaryObject()
     f = p.faces
     pickedFaces = f[0:1]   #f.getSequenceFromMask(mask=('[#1 ]', ), )
     p.PartitionFaceBySketch(faces=pickedFaces, sketch=s2)
     a = myModel.rootAssembly
     a.Instance(name=name+str(itera), part=p, dependent=ON)

     tf = a.instances[name+str(itera)].faces
     gommeface = a.instances[name+str(itera)].faces.findAt(((intpt[0],intpt[1],0),))
     gommeregion=(gommeface,)
     a.Set(faces=[gommeface], name='Gomme')
     ms=a.getMassProperties(regions=gommeregion)
     mstot=a.getMassProperties(regions=tf[0:len(tf)])
     cr=1.-(ms['area']/mstot['area'])
     print 'Taux volumique de charge :',cr
     listEp.append((cr,EP))

   for i in range(0,len(refpoint)):   
     xd=refpoint[i][0]
     yd=refpoint[i][1]
     p.DatumPointByCoordinate(coords=(xd, yd, 0.0))
   vp = getCurrentViewport()
   vp.setValues(displayedObject=p)

   dat = a.instances[name+str(itera)].datums

   vp.assemblyDisplay.setValues(interactions=ON,constraints=ON, connectors=ON, engineeringFeatures=ON)
   vp.assemblyDisplay.setValues(renderStyle=WIREFRAME)
   InclusionFace=[]
   for i in range(0,len(dat)):
     r1=a.ReferencePoint(point=(refpoint[i][0],refpoint[i][1],0))
     lpoints=[]
     for j in range(0,len(cellpoint[i])):
         apoints=(cellpoint[i][j][0],cellpoint[i][j][1],0.0)
         lpoints.append(apoints)
     center=[0,0]
     for pts in voronoi[refpoint[i]]:
          center[0]=pts[0]+center[0]
          center[1]=pts[1]+center[1]
     if(len(voronoi[refpoint[i]])>0):
       center[0]=center[0]/len(voronoi[refpoint[i]])
       center[1]=center[1]/len(voronoi[refpoint[i]])
     f1 = a.instances[name+str(itera)].faces.findAt(((center[0],center[1],0),))
     region2=(f1,)
     ref = a.referencePoints
     refPoints1=(ref[r1.id], )
     region1=regionToolset.Region(referencePoints=refPoints1)
     nom='Constraint'+str(i+1)
     myModel.RigidBody(name=nom, refPointRegion=region1,bodyRegion=region2)

   gommeface = a.instances[name+str(itera)].faces.findAt(((intpt[0],intpt[1],0),))
   gommeregion=(gommeface,)
   a.Set(faces=[gommeface], name='Gomme')
   ms=a.getMassProperties(regions=gommeregion)
   tf = a.instances[name+str(itera)].faces
   mstot=a.getMassProperties(regions=tf[0:len(tf)])
   
   for fc in tf:
      if(fc.index!=gommeface[0].index):  
        InclusionFace.append(tf[fc.index:fc.index+1]) 
   a.Set(faces=tuple(InclusionFace), name='Inclusion')

   if(mask==1):
      perim=-(2.*BOUNDSIZE[1]+2.*BOUNDSIZE[0])
   elif(mask==2):
      perim=-2.*pi*(width*0.5) 
   elif(mask==3):
      perim=-pi*(3.*0.5*(width+height)-sqrt(0.25*(3.*width+height)*(width+3.*height)))
    
   for i in range(0,len(p.edges)):
      perim=perim+p.edges[i].getSize(False)   

   print 'Taux volumique de charge :',1-(ms['area']/mstot['area'])
   print 'Surface specifique d echange :',perim
   print 'nombre de grain:',len(points) 
   makedidDiagram(pinit,BOUNDSIZE,a,name+str(itera))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def voronoicellarea(poly,a,faces):
   #calcul du centre du polygone
   center=[0,0]
   for pt in poly:
     center[0]=center[0]+pt[0]
     center[1]=center[1]+pt[1]
   center[0]=center[0]/len(poly)
   center[1]=center[1]/len(poly)
   cellface=faces.findAt(((center[0],center[1],0),))
   cellregion=(cellface,)
   ms=a.getMassProperties(regions=cellregion)
   return ms['area']
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def voronoicellface(poly,faces):
   #calcul du centre du polygone
   center=[0,0]
   for pt in poly:
     center[0]=center[0]+pt[0]
     center[1]=center[1]+pt[1]
   center[0]=center[0]/len(poly)
   center[1]=center[1]/len(poly)
   cellface=faces.findAt(((center[0],center[1],0),))
   return cellface

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generation d'un modele complet avec choix du type de geometrie de base
# rectangle, cercle, ellipse du taux de charge, du nombre initiale cellules,
# du rayon moyen et de l'epaisseur : tirage aletoire de gomme occluse...

def VoronoiCell3(name, mask,width,height,RMEAN,EP,NPOINT,FRATIO):
   BOUNDSIZE=(width,height)
   myModel = getNewModel()
   s0 = myModel.ConstrainedSketch(name='Sketch ini', sheetSize=200.0)
   s = myModel.ConstrainedSketch(name='Sketch without thickness', sheetSize=200.0)
   s2 = myModel.ConstrainedSketch(name='Sketch with thickness', sheetSize=200.0)
   sr = myModel.ConstrainedSketch(name='Sketch Mask', sheetSize=200.0)
   bound=[]
   if(mask==1):
     sr.rectangle(point1=(0.0, 0.0), point2=(BOUNDSIZE[0], BOUNDSIZE[1]))
     geomcarac=[BOUNDSIZE[0],BOUNDSIZE[1]]
     bound = [(0., 0.), (0., BOUNDSIZE[1]), (BOUNDSIZE[0], BOUNDSIZE[1]), (BOUNDSIZE[0], 0.)]
   elif(mask==2):
     ndisk=50
     for i in range(0,ndisk):
        angle=radians(i*(360./ndisk))
        center=(BOUNDSIZE[0]/2.,BOUNDSIZE[1]/2.)  
        r=BOUNDSIZE[0]/2.
        apoint=(center[0]-r*cos(angle),center[1]+r*sin(angle))
        bound.append(apoint)    
     sr.CircleByCenterPerimeter(center=(width*0.5, width*0.5), point1=( width*0.5, width))
     geomcarac=[width*0.5]
   elif(mask==3):
     bound = [(0., 0.), (0., BOUNDSIZE[1]), (BOUNDSIZE[0], BOUNDSIZE[1]), (BOUNDSIZE[0], 0.)]
     sr.EllipseByCenterPerimeter(center=(width*0.5,height*0.5), axisPoint1=(width,height*0.5),axisPoint2=( width*0.5, height))
     geomcarac=[width*0.5,height*0.5]

   # tirage aleatoire avec controle de la distance entre germe
   DMIN=2.*RMEAN+EP
   center=(BOUNDSIZE[0]/2.,BOUNDSIZE[1]/2.)
   points,voronoi = constrainedVoronoi(bound,NPOINT,DMIN,False,BOUNDSIZE[1]/2.,center) 

   drawVoronoi(voronoi,bound,points,s0) 
   # supression des cellules qui n'ont pas leurs germes dans le masque+construction du contour   
   voronoi = verifMask(voronoi,points,mask,geomcarac)  
   drawVoronoi(voronoi,bound,points,s) 
 
   myModel.convertAllSketches()
   p = myModel.Part(name="Initial Voronoi", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   pinit = p
   p.BaseShell(sketch=sr)
   s.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)

   voronoi = verifVoronoi(voronoi,points,0.9*EP)
   # recuperation d'un sommet pour la selection de la gomme 
   intpt=getInteriorpoint(voronoi,points,BOUNDSIZE,mask)
   refpoint,cellpoint = addThickness(voronoi,bound,points,EP,s2)
   del refpoint[len(refpoint)-1]

   myModel.convertAllSketches()
   p = myModel.Part(name="Voronoi", dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)
   p.BaseShell(sketch=sr)
   s2.unsetPrimaryObject()
   f = p.faces
   pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
   p.PartitionFaceBySketch(faces=pickedFaces, sketch=s2)

   # calcul du taux de charge initial
   a = myModel.rootAssembly
   a.Instance(name="Voronoi", part=p, dependent=ON)
   tf = a.instances["Voronoi"].faces
   gommeface = a.instances["Voronoi"].faces.findAt(((intpt[0],intpt[1],0),))
   gommeregion=(gommeface,)
   ms=a.getMassProperties(regions=gommeregion)
   mstot=a.getMassProperties(regions=tf[0:len(tf)])
   cr=0
   for i in range(0,len(points)-1):
     polyi=voronoi[points[i]]  
     cr=cr+voronoicellarea(polyi,a,tf)
   cr=cr/mstot['area']

   rubbermatrix=[]
   pointsnum=(range(0,len(refpoint)-1))

   if(cr<=FRATIO):
       print 'Valeur cible non atteinte: C=',cr, '\nIl faut augementer le nombre de grain ou diminuer la valeur cible'
   else:
     itera=0
     newarea=0 
     while (abs(cr-FRATIO)/FRATIO>0.005 and itera<NPOINT):
       if(len(pointsnum)<2): break
       aparticule=random.choice(pointsnum)  
       itera=itera+1
       polyi=voronoi[refpoint[aparticule]]  
       newratio=cr-(voronoicellarea(polyi,a,tf)/mstot['area'])
       print 'iteration ',itera,' C=',newratio
       if(newratio>=FRATIO):
          rubbermatrix.append(aparticule)
          pointsnum.remove(aparticule)
          cr=newratio 
   particuleslist={}
   k=0
   radius=[]
   for i in range(0,len(refpoint)): 
     if(len(rubbermatrix)<=0 or i not in rubbermatrix):
       xd=refpoint[i][0]
       yd=refpoint[i][1]
       particuleslist[k]=i  
       p.DatumPointByCoordinate(coords=(xd, yd, 0.0))
       k=k+1
       polyi=voronoi[refpoint[i]]
       radius.append(sqrt(voronoicellarea(polyi,a,tf))/pi)

   dat = a.instances["Voronoi"].datums

   InclusionFace=[]


   for i in range(0,len(dat)):
     r1=a.ReferencePoint(point=(refpoint[particuleslist[i]][0],refpoint[particuleslist[i]][1],0))
     lpoints=[]
     for j in range(0,len(cellpoint[particuleslist[i]])):
         apoints=(cellpoint[particuleslist[i]][j][0],cellpoint[particuleslist[i]][j][1],0.0)
         lpoints.append(apoints)
     center=[0,0]
     for pts in voronoi[refpoint[particuleslist[i]]]:
          center[0]=pts[0]+center[0]
          center[1]=pts[1]+center[1]
     if(len(voronoi[refpoint[i]])>0):
       center[0]=center[0]/len(voronoi[refpoint[particuleslist[i]]])
       center[1]=center[1]/len(voronoi[refpoint[particuleslist[i]]])
     f1 = tf.findAt(((center[0],center[1],0),))
     fp = p.faces.findAt(((center[0],center[1],0),))
     region2=(f1,)
     ref = a.referencePoints
     refPoints1=(ref[r1.id], )
     region1=regionToolset.Region(referencePoints=refPoints1)
     InclusionFace.append(fp) 
     nom='Constraint'+str(i+1)
     myModel.RigidBody(name=nom, refPointRegion=region1,bodyRegion=region2)



   rubberocclused=[]
   tf = p.faces
   for i in rubbermatrix:
      polyi=voronoi[refpoint[i]]  
      f1=voronoicellface(polyi,tf)
      rubberocclused.append(f1)

   gommeface = p.faces.findAt(((intpt[0],intpt[1],0),))
   gommeregion=(gommeface,)
   if(len(gommeface)>0):
     p.Set(faces=[gommeface], name='Interphace')
   if(len(InclusionFace)>0):
     p.Set(faces=InclusionFace, name='Inclusion')
   if(len(rubberocclused)>0):
     p.Set(faces=rubberocclused, name='Occluded Rubber')

   rubberface= p.faces.findAt(((BOUNDSIZE[0]/2.,BOUNDSIZE[1]-EP/10,0),))
   if(len(rubberface)>0):
     p.Set(faces=[rubberface], name='Rubber')   

   print 'Taux volumique de charge :',cr
   print 'nombre de grain:',len(dat)
   print 'Rayon moyen',mean(array(radius))
   print 'Rayon moyen/ep',mean(array(radius))/float(EP)
   print 'Ecart Type Rayon moyen/ep',std(array(radius))/float(EP)
 #  makedidDiagram(pinit,BOUNDSIZE,a,"Voronoi")



