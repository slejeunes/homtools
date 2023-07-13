#-----------------------------------------------------------------------
#     Plugin for generating circular cells (2D)
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#     
#-----------------------------------------------------------------------

from abaqusGui import *
from sphereDB import sphereDB


###########################################################################
# Class definition
###########################################################################

class SphereForm(AFXForm):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):

        # Construct the base class.
        #
        AFXForm.__init__(self, owner)
                
        # Command
        #
        self.cmd = AFXGuiCommand(self, 'SphereCell', 'macroSphere')
        
        self.nameKw = AFXStringKeyword(self.cmd, 'name', TRUE,"Sphere")
        self.maskKw = AFXIntKeyword(self.cmd, 'mask',TRUE,1,FALSE)
        self.widthKw = AFXFloatKeyword(self.cmd, 'width', TRUE,1.0)
        self.heightKw = AFXFloatKeyword(self.cmd, 'height', TRUE,1.0)
        self.dmaxKw = AFXFloatKeyword(self.cmd, 'DMAX', TRUE,0.1)
        self.npointKw = AFXFloatKeyword(self.cmd, 'NPOINT', TRUE,10)
        self.epKw = AFXFloatKeyword(self.cmd, 'EP', TRUE,0.03)
        self.removeBKw = AFXBoolKeyword(self.cmd, 'REMOVEB', AFXBoolKeyword.ON_OFF, FALSE, FALSE)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):

        return sphereDB(self)


import os
absPath = os.path.abspath(__file__)
absDir  = os.path.dirname(absPath)
helpUrl = os.path.join(absDir, 'www.lma.cnrs-mrs.fr')

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()

# Register a GUI plug-in in the Plug-ins menu.
#
toolset.registerGuiMenuButton(
    object=SphereForm(toolset), buttonText='Homtools|Spheric Inclusions (2D)',
    kernelInitString='import macroSphere',
    version='1.0', author='S. Lejeunes & S. Bourgeois (LMA-CNRS UPR7051)',
    applicableModules = ['Part'],
    description='A simple Gui to generate Spheric inclusions geometry '
                "This plug-in's files may be copied from " + absDir,
    helpUrl=helpUrl
)
    


