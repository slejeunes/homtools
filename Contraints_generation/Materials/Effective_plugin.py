#-----------------------------------------------------------------------
#     Plugin for linear deformation on the boundary
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#     
#-----------------------------------------------------------------------
from abaqusGui import *
from kernelAccess import mdb, session
import i18n
########################################################################
# Class definition
########################################################################
class EffectiveForm(AFXForm):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):

        # Construct the base class.
        # 
        AFXForm.__init__(self, owner)
                
        # Command
        #
        self.cmd = AFXGuiCommand(self, 'Effective', 'homogeneBoundary')
        self.is_smallstrain=AFXBoolKeyword(self.cmd, 'is_smallstrain', AFXBoolKeyword.TRUE_FALSE, True,True)
        self.dim=AFXIntKeyword(self.cmd, 'dim',3)
        self.owner=owner
        self.setS1=None
        self.setM3=None
        self.setM2=None
        self.setM1=None
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):
        vpName = session.currentViewportName
        modelName = session.sessionState[vpName]['modelName']
        m=mdb.models[modelName]
        a = m.rootAssembly
        self.dim.setValue(3)
        if(a.getMassProperties()['volume']==None): 
           self.dim.setValue(2)

        db = EffectiveDB(self,self.dim)
        return db

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class EffectiveDB(AFXDataDialog):
    [
        ID_M1,
        ID_M2,
        ID_M3,
        ID_S1,
    ] = range(AFXDataDialog.ID_LAST, AFXDataDialog.ID_LAST+4)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form,dim):
        sendCommand("homogeneBoundary.__init__()")
        form.setS1=None

        self.step2=None
        AFXDataDialog.__init__(self, form, 'Linear displacements on the boundary', self.OK|self.CANCEL)

                      
        vf = FXVerticalFrame(self, LAYOUT_FILL_Y|LAYOUT_FILL_X)
        self.comboBox = AFXComboBox(vf, 0, 2, 'Formulation:      ',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_Y)
        self.comboBox.appendItem('Small strain')
        self.comboBox.appendItem('Finite strain')
        form.is_smallstrain.setValue(True)

        gb = FXGroupBox(vf, 'Macro Nodes', LAYOUT_FILL_X|FRAME_GROOVE)
        gb1 = FXGroupBox(vf, 'Boundary Sets', LAYOUT_FILL_X|FRAME_GROOVE)
        hf = FXMatrix(gb, 2,opts=MATRIX_BY_COLUMNS|LAYOUT_FILL_X|LAYOUT_FILL_Y)        

        self.dim=dim.getValue()
        self.label1=FXLabel(hf, 'Select the first RefPoint (E11,E22)',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        FXButton(hf, 'Edit...',None,self,self.ID_M1,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
        self.label2=FXLabel(hf, 'Select the second RefPoint (GAMMA12,Unused)',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        FXButton(hf, 'Edit...',None,self,self.ID_M2,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
        if(self.dim==3): 
          self.label1.setText('Select the first RefPoint (E11,E22,E33)')
          self.label2.setText('Select the second RefPoint (GAMMA12=2xE12,GAMMA13=2xE13,GAMMA23=2xE23)')
          self.label3=FXLabel(hf, 'Select the Third RefPoint (Unused)',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
          self.label3.disable()
          self.button3=FXButton(hf, 'Edit...',None,self,self.ID_M3,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
          self.button3.disable()

        hf2 = FXMatrix(gb1, 2,opts=MATRIX_BY_COLUMNS|LAYOUT_FILL_X|LAYOUT_FILL_Y)        
        label4=FXLabel(hf2,'Select the Edges',opts=JUSTIFY_LEFT|LAYOUT_FILL_X|LAYOUT_FILL_COLUMN)
        if(self.dim==3): 
           label4.setText('Select the Faces')
        FXButton(hf2, 'Edit...',None,self,self.ID_S1,opts=LAYOUT_RIGHT|BUTTON_NORMAL)
   

        FXMAPFUNC(self,SEL_COMMAND,self.ID_M1,EffectiveDB.onCmdM1)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_M2,EffectiveDB.onCmdM2)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_M3,EffectiveDB.onCmdM3)
        FXMAPFUNC(self,SEL_COMMAND,self.ID_S1,EffectiveDB.onCmdS1)
        self.form=form

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def processUpdates(self):
        if(self.comboBox.getCurrentItem()==0):
            if(self.dim==2):
              self.label1.setText('Select the first RefPoint (E11,E22)')
              self.label2.setText('Select the second RefPoint (GAMMA12=2xE12,Unused)')
            else:
              self.label1.setText('Select the first RefPoint (E11,E22,E33)')
              self.label2.setText('Select the second RefPoint (GAMMA12=2xE12,GAMMA13=2xE13,GAMMA23=2xE23)')
              self.label3.setText('Select the third RefPoint (Unused)')
              self.label3.disable()
              self.button3.disable()
            self.form.is_smallstrain.setValue(True)
        else:
            if(self.dim==2):
              self.label1.setText('Select the first RefPoint (F11-1,F12)')
              self.label2.setText('Select the second RefPoint (F21,F22-1)')
            else:
              self.label1.setText('Select the first RefPoint (F11-1,F12,F13)')
              self.label2.setText('Select the second RefPoint (F21,F22-1,F23)')
              self.label3.setText('Select the third RefPoint (F31,F32,F33-1)')
              self.label3.enable()
              self.button3.enable()
            self.form.is_smallstrain.setValue(False)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdM1(self,sender,sel,ptr):
         self.step1=getRP1rocedure(self.form.owner)
         self.step1.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdM2(self,sender,sel,ptr):
         self.step2=getRP2rocedure(self.form.owner)
         self.step2.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdM3(self,sender,sel,ptr):
         self.step3=getRP3rocedure(self.form.owner)
         self.step3.activate()
         return 1
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdS1(self,sender,sel,ptr):
         self.step2=getS1Procedure(self.form.owner,self.dim)
         self.step2.activate()
         return 1


class getRP1rocedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        cmd = AFXGuiCommand(self,'setRefpoint1', 'homogeneBoundary')
        self.setM1 = AFXObjectKeyword(cmd, 'RP1', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setM1, i18n.tr('Select a Reference Point'), REFERENCE_POINTS, ONE)
        return self.step1


class getRP2rocedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        cmd = AFXGuiCommand(self,'setRefpoint2', 'homogeneBoundary')
        self.setM2 = AFXObjectKeyword(cmd, 'RP2', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setM2, i18n.tr('Select a Reference Point'), REFERENCE_POINTS, ONE)
        return self.step1

class getRP3rocedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner):
        AFXProcedure.__init__(self, owner=owner)
        cmd = AFXGuiCommand(self,'setRefpoint3', 'homogeneBoundary')
        self.setM3 = AFXObjectKeyword(cmd, 'RP3', TRUE)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        self.step1=AFXPickStep(self, self.setM3, i18n.tr('Select a Reference Point'), REFERENCE_POINTS, ONE)
        return self.step1

class getS1Procedure(AFXProcedure):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self,  owner,dim):
        self.dim=dim
        AFXProcedure.__init__(self, owner=owner)
        self.cmd = AFXGuiCommand(self, 'setGroup1', 'homogeneBoundary')
        self.setS1 = AFXObjectKeyword(self.cmd, 'set1', TRUE)
        self.step1=None
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstStep(self):
        if(self.dim==2):
          self.step1=AFXPickStep(self, self.setS1, i18n.tr('Select Edges'), EDGES|NODES|VERTICES, MANY,1,TUPLE)
        else:
          self.step1=AFXPickStep(self, self.setS1, i18n.tr('Select Faces'), FACES|EDGES|NODES|VERTICES, MANY,1,TUPLE)
        return self.step1



import os
absPath = os.path.abspath(__file__)
absDir  = os.path.dirname(absPath)
helpUrl = os.path.join(absDir, 'www.lma.cnrs-mrs.fr')

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()

# Register a GUI plug-in in the Plug-ins menu.
#
toolset.registerGuiMenuButton(
    object=EffectiveForm(toolset), buttonText='Homtools|Kinematic Uniform Boundary Condition',
    kernelInitString='import homogeneBoundary; homogeneBoundary=homogeneBoundary.HomogeneBoundary()',
    version='1.0', author='S. Lejeunes & S. Bourgeois (LMA-CNRS UPR7051)',
    applicableModules = ['Interaction'],
    description='A simple Gui to define linear deformation conditions on boundary'
                "This plug-in's files may be copied from " + absDir,
    helpUrl=helpUrl
)



