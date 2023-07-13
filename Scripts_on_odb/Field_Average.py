#-----------------------------------------------------------------------
#     Script to compute the avarage of a field on a domain
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#     
#-----------------------------------------------------------------------
#from abaqus import *
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import __main__

# to be defined by the user
odbname='Job.odb' 
fieldname='U'
elsetname=' ALL ELEMENTS'
stepname='Step-1'
framenumber=-1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getinterpts(ty):

    gpts=()
    if(ty[0:4]=='CPS4' or ty[0:4]=='CPE4'):
        gpts=((-0.577350269189626,-0.577350269189626),
              (0.577350269189626,-0.577350269189626),
              (0.577350269189626,0.577350269189626),
              (-0.577350269189626,0.577350269189626))

    elif(ty[0:4]=='CPS6' or ty[0:4]=='CPE6'):
        gpts=((1./6.,1./6.),(2./3.,1./6.),(1./6.,2./3.))

    elif(ty[0:4]=='CPS3' or ty[0:4]=='CPE3'):
        gpts=((0.577350269189626,0.577350269189626),)

    elif(ty[0:4]=='CPS8' or ty[0:4]=='CPE8'):
        gpts=((-0.774596669241483,-0.774596669241483),
             (0.0,-0.774596669241483),
             (0.774596669241483,-0.774596669241483),
             (-0.774596669241483,0.0),
             (0.0,0.0),
             (0.774596669241483,0.0),
             (-0.774596669241483,0.774596669241483),
             (0.0,0.774596669241483),
             (0.774596669241483,0.774596669241483))
  

    elif(ty[0:4]=='C3D4'):
          gpts=((0.577350269189626,0.577350269189626,0.577350269189626),)
        
    elif(ty[0:5]=='C3D10'):
        gpts=((1./6.,1./6.,1./6.),
              (2./3.,1./6.,1./6.),
              (1./6.,2./3.,1./6.),
              (1./6.,1./6.,2./3.))

    elif(ty[0:4]=='C3D8'):
        gpts=((-0.577350269189626,-0.577350269189626,-0.577350269189626),
              (0.577350269189626,-0.577350269189626,-0.577350269189626),
              (0.577350269189626,0.577350269189626,-0.577350269189626),
              (-0.577350269189626,0.577350269189626,-0.577350269189626),
              (-0.577350269189626,-0.577350269189626,0.577350269189626),
              (0.577350269189626,-0.577350269189626,0.577350269189626),
              (0.577350269189626,0.577350269189626,0.577350269189626),
              (-0.577350269189626,0.577350269189626,0.577350269189626))

    elif(ty[0:5]=='C3D20'):
        gpts=((-0.774596669241483,-0.774596669241483,-0.774596669241483),
             (0.0,-0.774596669241483,-0.774596669241483),
             (0.774596669241483,-0.774596669241483,-0.774596669241483),
             (-0.774596669241483,0.0,-0.774596669241483),
             (0.0,0.0,-0.774596669241483),
             (0.774596669241483,0.0,-0.774596669241483),
             (-0.774596669241483,0.774596669241483,-0.774596669241483),
             (0.0,0.774596669241483,-0.774596669241483),
             (0.774596669241483,0.774596669241483,-0.774596669241483),
             (-0.774596669241483,-0.774596669241483,0.0),
             (0.0,-0.774596669241483,0.0),
             (0.774596669241483,-0.774596669241483,0.0),
             (-0.774596669241483,0.0,0.0),
             (0.0,0.0,0.0),
             (0.774596669241483,0.0,0.0),
             (-0.774596669241483,0.774596669241483,0.0),
             (0.0,0.774596669241483,0.0),
             (0.774596669241483,0.774596669241483,0.0),
             (-0.774596669241483,-0.774596669241483,0.774596669241483),
             (0.0,-0.774596669241483,0.774596669241483),
             (0.774596669241483,-0.774596669241483,0.774596669241483),
             (-0.774596669241483,0.0,0.774596669241483),
             (0.0,0.0,0.774596669241483),
             (0.774596669241483,0.0,0.774596669241483),
             (-0.774596669241483,0.774596669241483,0.774596669241483),
             (0.0,0.774596669241483,0.774596669241483),
             (0.774596669241483,0.774596669241483,0.774596669241483))


    elif(ty[0:4]=='C3D6'):
        gpts=((0.577350269189626,0.577350269189626,-0.577350269189626),
              (0.577350269189626,0.577350269189626,0.577350269189626))


    elif(ty[0:5]=='C3D15'):
        gpts=((1./6.,1./6.,-0.577350269189626),
              (2./3.,1./6.,-0.577350269189626),
              (1./6.,2./3.,-0.577350269189626),
              (1./6.,1./6.,0),
              (2./3.,1./6.,0),
              (1./6.,2./3.,0),
              (1./6.,1./6.,0.577350269189626),
              (2./3.,1./6.,0.577350269189626),
              (1./6.,2./3.,0.577350269189626))

    return gpts

 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getinterfunct(ty,chi):

    Ni=()
    if(ty[0:4]=='CPS4' or ty[0:4]=='CPE4'):
        Ni=((1.0-chi[0]-chi[1]+chi[1]*chi[0])*0.25,
           (1.0+chi[0]-chi[1]-chi[1]*chi[0])*0.25,
           (1.0+chi[0]+chi[1]+chi[1]*chi[0])*0.25,
           (1.0-chi[0]+chi[1]-chi[1]*chi[0])*0.25)

    elif(ty[0:4]=='CPS6' or ty[0:4]=='CPE6'):
        u=1.0-chi[0]-chi[1]
        Ni=(u*(2.0*u-1.0),
               chi[0]*(2.0*chi[0]-1.0),
               chi[1]*(2.0*chi[1]-1.0),
               4.0*chi[0]*u,               
               4.0*chi[0]*chi[1],
               4.0*chi[1]*u)

    elif(ty[0:4]=='CPS3'or ty[0:4]=='CPE3'):
        Ni=(1.0-chi[0]-chi[1],
               chi[0],
               chi[1])



    elif(ty[0:4]=='CPS8' or ty[0:4]=='CPE8'):
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


    elif(ty[0:4]=='C3D4'):
        Ni=((1-chi[0]-chi[1]-chi[2]),
              chi[0],
              chi[1],
              chi[2])

    elif(ty[0:4]=='C3D8'):
        x=1+chi[0]
        y=1+chi[1]
        z=1+chi[2]
        r=1-chi[0]
        s=1-chi[1]
        t=1-chi[2]
        c=1./8
        Ni=(c*r*s*t,
             c*x*s*t,
             c*x*y*t,
             c*r*y*t,
             c*r*s*z,
             c*x*s*z,
             c*x*y*z,
             c*r*y*z)

    elif(ty[0:4]=='C3D6'):
        l=1-chi[0]-chi[1]
        a=(1-chi[2])/2
        b=(1-chi[2])/2
        Ni=(l*a,
            chi[0]*a,
            chi[1]*a,
            l*b,
            chi[0]*b,
            chi[1]*b)

    elif(ty[0:5]=='C3D10'):
        l=1-chi[0]-chi[1]-chi[2]
        Ni=(-l*(1-2*l),
             -chi[0]*(1-2*chi[0]),
             -chi[1]*(1-2*chi[1]),
             -chi[2]*(1-2*chi[2]),
             4*l*chi[0],
             4*chi[0]*chi[1],
             4*l*chi[1],
             4*l*chi[2],
             4*chi[0]*chi[2],
             4*chi[1]*chi[2])


    elif(ty[0:5]=='C3D15'):
        l=1-chi[0]-chi[1]
        x=-chi[2]*(1-chi[2])/2
        y=chi[2]*(1+chi[2])/2
        Ni=(-l*(1-2*l)*x,
            -chi[0]*(1-2*chi[0])*x,
            -chi[1]*(1-2*chi[1])*x,
            -l*(1-2*l)*y,
            -chi[0]*(1-2*chi[0])*y,
            -chi[1]*(1-2*chi[1])*y,
            4*l*chi[0]*x,
            4*chi[1]*chi[0]*x,
            4*l*chi[1]*x,
            4*l*chi[0]*y,
            4*chi[1]*chi[0]*y,
            4*l*chi[1]*y,
            l*(1-chi[2]*chi[2]),
            chi[0]*(1-chi[2]*chi[2]),
            chi[1]*(1-chi[2]*chi[2]))


    elif(ty[0:5]=='C3D20'):
         s=chi[0]
         t=chi[1]
         u=chi[2]
         s0=1.+s
         s1=1.-s
         s2=1.-s*s
         t0=1.+t
         t1=1.-t
         t2=1.-t*t
         u0=1.+u
         u1=1.-u
         u2=1.-u*u
         Ni=(-0.1250*s1*t1*u1*(2.00+s+t+u),
               -0.1250*s0*t1*u1*(2.00-s+t+u),
               -0.1250*s0*t0*u1*(2.00-s-t+u),
               -0.1250*s1*t0*u1*(2.00+s-t+u),
               -0.1250*s1*t1*u0*(2.00+s+t-u),
               -0.1250*s0*t1*u0*(2.00-s+t-u),
               -0.1250*s0*t0*u0*(2.00-s-t-u),
               -0.1250*s1*t0*u0*(2.00+s-t-u),
                0.250*s2*t1*u1,
                0.250*s0*t2*u1,
                0.250*s2*t0*u1,
                0.250*s1*t2*u1,
                0.250*s2*t1*u0,
                0.250*s0*t2*u0,
                0.250*s2*t0*u0,
                0.250*s1*t2*u0,
                0.250*s1*t1*u2,
                0.250*s0*t1*u2,
                0.250*s0*t0*u2,
                0.250*s1*t0*u2)
#        Ni=((1-chi[0])*(1-chi[1])*(1-chi[2])*(-2-chi[0]-chi[1]-chi[2])/8.0,
#            (1+chi[0])*(1-chi[1])*(1-chi[2])*(-2+chi[0]-chi[1]-chi[2])/8.0,
#            (1+chi[0])*(1+chi[1])*(1-chi[2])*(-2+chi[0]+chi[1]-chi[2])/8.0,
#            (1-chi[0])*(1+chi[1])*(1-chi[2])*(-2-chi[0]+chi[1]-chi[2])/8.0,
#            (1-chi[0])*(1-chi[1])*(1+chi[2])*(-2-chi[0]-chi[1]+chi[2])/8.0,
#            (1+chi[0])*(1-chi[1])*(1+chi[2])*(-2+chi[0]-chi[1]+chi[2])/8.0,
#            (1+chi[0])*(1+chi[1])*(1+chi[2])*(-2+chi[0]+chi[1]+chi[2])/8.0,
#            (1-chi[0])*(1+chi[1])*(1+chi[2])*(-2-chi[0]+chi[1]+chi[2])/8.0,
#            (1-chi[0]*chi[0])*(1-chi[1])*(1-chi[2])/4.0,
#            (1+chi[0])*(1-chi[1]*chi[1])*(1-chi[2])/4.0,
#            (1-chi[0]*chi[0])*(1+chi[1])*(1-chi[2])/4.0,
#            (1-chi[0])*(1-chi[1]*chi[1])*(1-chi[2])/4.0,
#            (1-chi[0]*chi[0])*(1-chi[1])*(1+chi[2])/4.0,
#            (1+chi[0])*(1-chi[1]*chi[1])*(1+chi[2])/4.0,
#            (1-chi[0]*chi[0])*(1+chi[1])*(1+chi[2])/4.0,
#            (1-chi[0])*(1-chi[1]*chi[1])*(1+chi[2])/4.0,
#            (1-chi[0])*(1-chi[1])*(1-chi[2]*chi[2])/4.0,
#            (1+chi[0])*(1-chi[1])*(1-chi[2]*chi[2])/4.0,
#            (1+chi[0])*(1+chi[1])*(1-chi[2]*chi[2])/4.0,
#            (1-chi[0])*(1+chi[1])*(1-chi[2]*chi[2])/4.0)

    return Ni

 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




odb = openOdb(path=odbname)

lastFrame = odb.steps[stepname].frames[framenumber]

#definition d'une region
phase=odb.rootAssembly.elementSets[elsetname]
varfield = lastFrame.fieldOutputs[fieldname]
field = varfield.getSubset(position= NODAL)
Values = field.values

#creation d'un dictionnaire ayant pour chaque nodeLabel, la position associe dans le tableau Values
place={}
g=0
for i in range(len(Values)):
    if Values[i].instance!= None:
        place[Values[i].nodeLabel]=g
    g= g+1

l=len(Values[0].data)

jacobien= lastFrame.fieldOutputs['IVOL']
jac= jacobien.getSubset(position= INTEGRATION_POINT, region=phase)
det = jac.values

el=odb.rootAssembly.elementSets[elsetname].elements[0]

z=[0]*l
a = [0.]*l
b= 0

for t in det:
    b = b +t.data

integer=[0,0,0]
val=[0.]*len(el)

#pour chaque element
for i in range(len(el)):
    # on lance la fonction getvalues avec le type de lelement
    ipts=getinterpts(el[i].type)
    nodes=el[i].connectivity
   
    # on fait la boucle sur le nombre de compososantes de la variables
    for k in range(l):
        f=[0.]*l
        #on fait une boucle sur le nombre de pts dintegration pour avoir les fonctions Ni
        for g in ipts:
            if g!= None:
                t=0
                s=[0.]*len(nodes)
                N=getinterfunct(el[i].type, g)
                #on fait la somme en bouclant sur le nombre de noeuds dans l element
                for j in range(len(nodes)):
                    s[t]= s[t]+ N[j]*Values[place[nodes[j]]].data[k]
                f[k]= f[k]+s[t]*det[z[k]].data
                z[k]=z[k]+1
                if (t=='len(ipts)-1'):
                    t=0
                else:
                    t=t+1

        a[k]=a[k] +f[k]

# On imprime la somme des vitesses.
for j in range (l):
    a[j]=a[j]/b  
    print fieldname,j,'=', a[j]
print 'volume =', b

