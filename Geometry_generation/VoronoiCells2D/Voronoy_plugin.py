#-----------------------------------------------------------------------
#     Plugin for generating voronoy cells (2D)
#-----------------------------------------------------------------------
#     Authors: Stephane Lejeunes,  Stephane Bourgeois
#     Institute: LMA UPR7051, CNRS
#     Date: 24/02/10
#     
#-----------------------------------------------------------------------

from abaqusGui import *
from voronoiDB import *


###########################################################################
# Class definition
###########################################################################

class VoronoiForm(AFXForm):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):

        # Construct the base class.
        #
        AFXForm.__init__(self, owner)
                
        # Command
        #
        self.cmd = AFXGuiCommand(self, 'VoronoiCell', 'macroVoronoi')
        
        self.nameKw = AFXStringKeyword(self.cmd, 'name', TRUE,"Voronoi")
        self.maskKw = AFXIntKeyword(self.cmd, 'mask',TRUE,1,FALSE)
        self.widthKw = AFXFloatKeyword(self.cmd, 'width', TRUE,1.0)
        self.heightKw = AFXFloatKeyword(self.cmd, 'height', TRUE,1.0)
        self.dminKw = AFXFloatKeyword(self.cmd, 'DMIN', TRUE,0.1)
        self.npointKw = AFXFloatKeyword(self.cmd, 'NPOINT', TRUE,10)
        self.epKw = AFXFloatKeyword(self.cmd, 'EP', TRUE,0.03)
        self.periodicKw = AFXBoolKeyword(self.cmd, 'PERIODIC', AFXBoolKeyword.ON_OFF, FALSE, FALSE)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):

        return VoronoiDB(self)


import os
absPath = os.path.abspath(__file__)
absDir  = os.path.dirname(absPath)
helpUrl = os.path.join(absDir, 'www.lma.cnrs-mrs.fr')

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()

# Register a GUI plug-in in the Plug-ins menu.
#
toolset.registerGuiMenuButton(
    object=VoronoiForm(toolset), buttonText='Homtools|Voronoi cells (2D)',
    kernelInitString='import macroVoronoi',
    version='1.0', author='S. Lejeunes & S. Bourgeois (LMA-CNRS UPR7051)',
    applicableModules = ['Part'],
    description='A simple Gui to generate a Voronoi geometry '
                "This plug-in's files may be copied from " + absDir,
    helpUrl=helpUrl
)
    



