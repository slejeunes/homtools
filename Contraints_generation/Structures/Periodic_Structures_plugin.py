#-----------------------------------------------------------------------
#     Plugin for periodic condition on Structures
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 08/04/11
#
#-----------------------------------------------------------------------
from abaqusGui import *
from shellicon import *
from beamicon import *
from kernelAccess import mdb, session
import i18n
########################################################################
# Class definition
########################################################################
class PeriodicStructuresForm(AFXForm):
    
    ID_WARNING=AFXForm.ID_LAST
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):

        # Construct the base class.
        #
        AFXForm.__init__(self, owner)
                
        # Command
        #
        self.cmd = AFXGuiCommand(self, 'PeriodicStruct', 'periodicBoundaryStruct')

        self.pvectorX = AFXFloatKeyword(self.cmd, 'Vx', TRUE,1.0)
        self.pvectorY = AFXFloatKeyword(self.cmd, 'Vy', TRUE,0.0)
        self.pvectorZ = AFXFloatKeyword(self.cmd, 'Vz', TRUE,0.0)
        self.is_plate = AFXBoolKeyword(self.cmd, 'is_plate', AFXBoolKeyword.TRUE_FALSE, True,True)
        self.owner=owner
        self.ind=-1
        self.is_plate.setValue(True)
        FXMAPFUNC(self, SEL_COMMAND, self.ID_WARNING,PeriodicStructuresForm.onCmdWarning)
        self.ok=False
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):
        self.ind=-1
        self.db1 = PeriodicDBStruct(self,self.is_plate.getValue())
        return self.db1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getNextDialog(self,previous):
      if previous == self.db1:
          self.ind=self.ind+1
          self.db2 = PeriodicDBStruct2(self)
          return self.db2
      else:
          return None
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def doCustomChecks(self):
        if(self.ok==False):
          sendCommand('periodicBoundaryStruct.showMidPlane()')
          showAFXWarningDialog( self.getCurrentDialog(),'Is the Mid Plane Correct?',AFXDialog.YES | AFXDialog.NO,self,self.ID_WARNING)
          return False
        else:
          return True
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdWarning(self, sender, sel, ptr):
       if sender.getPressedButtonId() == AFXDialog.ID_CLICKED_YES:
         self.ok=True
       if(self.is_plate.getValue()):
           sendCommand('periodicBoundaryStruct.removeMidPlane()')
       else:
           sendCommand('periodicBoundaryStruct.removeRefAxis()')
       return 1
     #  elif sender.getPressedButtonId() == AFXDialog.ID_CLICKED_NO:
      #   self.show()



class PeriodicStructuresForm_Beam(PeriodicStructuresForm):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):
       PeriodicStructuresForm.__init__(self,owner) 
       self.is_plate.setValue(False)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def doCustomChecks(self):
        if(self.ok==False):
          sendCommand('periodicBoundaryStruct.showRefAxis()')
          showAFXWarningDialog( self.getCurrentDialog(),'Is the Ref Axis Correct?',AFXDialog.YES | AFXDialog.NO,self,self.ID_WARNING)
          return False
        else:
          return True


class PeriodicDBStruct2(AFXDataDialog):
    [
        ID_S1,
        ID_S2
    ] = range(AFXDataDialog.ID_LAST, AFXDataDialog.ID_LAST+2)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):
        AFXDataDialog.__init__(self, form, 'Periodic conditions for Shell', self.OK|self.CANCEL)
        vf = FXVerticalFrame(self,  LAYOUT_FILL_Y|LAYOUT_FILL_X)
        gb2 = FXGroupBox(vf, 'Boundary Sets', LAYOUT_FILL_X|FRAME_GROOVE)
        gb3 = FXGroupBox(vf, 'Periodicity Vector (in global CSYS)', LAYOUT_FILL_X|FRAME_GROOVE)
        # boundary set
        hf2 = FXMatrix(gb2, 2,opts=MATRIX_BY_COLUMNS|LAYOUT_FILL_X|LAYOUT_FILL_Y)        
        FXLabel(hf2,'Select the first face',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        FXButton(hf2, 'Edit...',None,self,self.ID_S1,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
        FXLabel(hf2, 'Select the second face',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        FXButton(hf2, 'Edit...',None,self,self.ID_S2,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
 
        AFXTextField(gb3, 5, 'X:                         ', form.pvectorX, opts=JUSTIFY_LEFT|LAYOUT_FILL_X|AFXTEXTFIELD_FLOAT)
        AFXTextField(gb3, 5, 'Y:                         ', form.pvectorY, opts=JUSTIFY_LEFT|LAYOUT_FILL_X|AFXTEXTFIELD_FLOAT)
        AFXTextField(gb3, 5, 'Z:                         ', form.pvectorZ, opts=JUSTIFY_LEFT|LAYOUT_FILL_X|AFXTEXTFIELD_FLOAT)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_S1,PeriodicDBStruct2.onCmdS1)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_S2,PeriodicDBStruct2.onCmdS2)

        self.form=form
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdS1(self,sender,sel,ptr):
         self.step2=getS1Procedure(self.form)
         self.step2.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdS2(self,sender,sel,ptr):
         self.step2=getS2Procedure(self.form)
         self.step2.activate()
         return 1

class PeriodicDBStruct(AFXDataDialog):
    [
        ID_M1,
        ID_M2,
       # ID_PL1,
        ID_CSYS
    ] = range(AFXDataDialog.ID_LAST, AFXDataDialog.ID_LAST+3)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form,is_plate):
        if is_plate:
           AFXDataDialog.__init__(self, form, 'Periodic conditions for Plates', self.CONTINUE|self.CANCEL)
        else:
           AFXDataDialog.__init__(self, form, 'Periodic conditions for Beam', self.CONTINUE|self.CANCEL)
        AFXNote(self, 'First Step for the definition of equations on the boundary\n')

        vf = FXVerticalFrame(self,  LAYOUT_FILL_Y|LAYOUT_FILL_X|LAYOUT_CENTER_X)
        if is_plate:
          icon = FXXPMIcon(getAFXApp(), shellicon)
          FXLabel(vf, '', icon, pl=0, pr=0, pt=0, pb=0,opts=LAYOUT_CENTER_X)
        else:
          icon = FXXPMIcon(getAFXApp(), beamicon)
          FXLabel(vf, '', icon, pl=0, pr=0, pt=0, pb=0,opts=LAYOUT_CENTER_X)
        gb0 = FXGroupBox(vf, 'Geometry', LAYOUT_FILL_X|FRAME_GROOVE)                      
        gb1 = FXGroupBox(vf, 'Macro Nodes', LAYOUT_FILL_X|FRAME_GROOVE)

        # geometry
        hf0 = FXMatrix(gb0, 2,opts=MATRIX_BY_COLUMNS|LAYOUT_FILL_X|LAYOUT_FILL_Y)        
        #self.label1=FXLabel(hf0, 'Select a reference plane',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        #FXButton(hf0, 'Edit...',None,self,self.ID_PL1,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
        self.label2=FXLabel(hf0, 'Select local coordinate system (z will define the normal to mid plane):',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        if not is_plate:
           self.label2.setText('Select local coordinate system (z will define the beam neutral axis):')
        FXButton(hf0, 'Edit...',None,self,self.ID_CSYS,opts=LAYOUT_RIGHT|BUTTON_NORMAL)

        # refpoints
        hf1 = FXMatrix(gb1, 2,opts=MATRIX_BY_COLUMNS|LAYOUT_FILL_X|LAYOUT_FILL_Y)     
        label1 = FXLabel(hf1, 'Select the first RefPoint (E11,E22,2E12)',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)  
        FXButton(hf1, 'Edit...',None,self,self.ID_M1,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
        label2=FXLabel(hf1, 'Select the second RefPoint (K11,K22,2K12)',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        FXButton(hf1, 'Edit...',None,self,self.ID_M2,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
        if not is_plate:        
          label1.setText('Select the first RefPoint (E,not used,not used)')
          label2.setText('Select the second RefPoint (K1,K2,KT)')
 

        FXMAPFUNC(self,SEL_COMMAND,self.ID_M1,PeriodicDBStruct.onCmdM1)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_M2,PeriodicDBStruct.onCmdM2)
       # FXMAPFUNC(self,SEL_COMMAND,self.ID_PL1,PeriodicDBStruct.onCmdPL1)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_CSYS,PeriodicDBStruct.onCmdCSYS)

        self.form=form


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdM1(self,sender,sel,ptr):
         self.step2=getRP1Procedure(self.form.owner)
         self.step2.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdM2(self,sender,sel,ptr):
         self.step2=getRP2Procedure(self.form.owner)
         self.step2.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdPL1(self,sender,sel,ptr):
         self.step2=getPL1Procedure(self.form.owner)
         self.step2.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdCSYS(self,sender,sel,ptr):
         self.step2=getCSYSProcedure(self.form.owner)
         self.step2.activate()
         return 1

class getRP1Procedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        self.cmd = AFXGuiCommand(self, 'setRefpoint1', 'periodicBoundaryStruct')
        self.setM1 = AFXObjectKeyword(self.cmd, 'RP1', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setM1, i18n.tr('Select a Reference Point'), REFERENCE_POINTS, ONE)
        return self.step1


class getRP2Procedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        self.cmd = AFXGuiCommand(self, 'setRefpoint2', 'periodicBoundaryStruct')
        self.setM2 = AFXObjectKeyword(self.cmd, 'RP2', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setM2, i18n.tr('Select a Reference Point'), REFERENCE_POINTS, ONE)
        return self.step1


class getPL1Procedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        self.cmd = AFXGuiCommand(self, 'setShellRefplane', 'periodicBoundaryStruct')
        self.setPL1 = AFXObjectKeyword(self.cmd, 'PL1', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setPL1, i18n.tr('Select a Datum Plane'), DATUM_PLANES, ONE)
        return self.step1

class getCSYSProcedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        self.cmd = AFXGuiCommand(self, 'setCSYS', 'periodicBoundaryStruct',True)
        self.setCSYS = AFXObjectKeyword(self.cmd, 'CSYS', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setCSYS, i18n.tr('Select a local Cordinate System'), DATUM_CSYS, ONE) 
        return self.step1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  def getNextStep(self,previous):
      #  sendCommand('periodicBoundaryStruct.showMidPlane()')


class getS1Procedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner.owner)
        self.cmd = AFXGuiCommand(self, 'setGroup1', 'periodicBoundaryStruct')
        self.setS1 = AFXObjectKeyword(self.cmd, 'set1', True)
        self.ind = AFXIntKeyword(self.cmd, 'ind', True, owner.ind)
        self.owner=owner
        self.step1=None
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.ind.setValue(self.owner.ind)
        self.step1=AFXPickStep(self, self.setS1, i18n.tr('Select Entities'), FACES|EDGES|REFERENCE_POINTS|NODES|VERTICES, MANY,1,TUPLE)
        return self.step1

class getS2Procedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner.owner)
        self.cmd = AFXGuiCommand(self, 'setGroup2', 'periodicBoundaryStruct')
        self.setS2 = AFXObjectKeyword(self.cmd, 'set2', True)
        self.ind = AFXIntKeyword(self.cmd, 'ind', True, owner.ind)
        self.owner=owner
        self.step1=None
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setS2, i18n.tr('Select Entities'), FACES|EDGES|REFERENCE_POINTS|NODES|VERTICES, MANY,1,TUPLE)
        self.ind.setValue(self.owner.ind)
        return self.step1


import os
absPath = os.path.abspath(__file__)
absDir  = os.path.dirname(absPath)
helpUrl = os.path.join(absDir, 'www.lma.cnrs-mrs.fr')

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()

# Register a GUI plug-in in the Plug-ins menu.
#
toolset.registerGuiMenuButton(
    object=PeriodicStructuresForm(toolset), buttonText='Homtools|Periodic Bounday Conditions for Plates',
    kernelInitString='import periodicBoundaryStruct; periodicBoundaryStruct=periodicBoundaryStruct.PeriodicBoundaryStruct()',
    version='1.0', author='S. Lejeunes & S. Bourgeois (LMA-CNRS UPR7051)',
    applicableModules = ['Interaction'],
    description='A simple Gui to define periodic Boundary Conditions for Structures'
                "This plug-in's files may be copied from " + absDir,
    helpUrl=helpUrl
)

toolset.registerGuiMenuButton(
    object=PeriodicStructuresForm_Beam(toolset), buttonText='Homtools|Periodic Bounday Conditions for Beams',
    kernelInitString='import periodicBoundaryStruct; periodicBoundaryStruct=periodicBoundaryStruct.PeriodicBoundaryStruct()',
    version='1.0', author='S. Lejeunes & S. Bourgeois (LMA-CNRS UPR7051)',
    applicableModules = ['Interaction'],
    description='A simple Gui to define periodic Boundary Conditions for Structures'
                "This plug-in's files may be copied from " + absDir,
    helpUrl=helpUrl
)

