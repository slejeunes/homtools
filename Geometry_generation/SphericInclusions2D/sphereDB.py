from abaqusGui import *
from sphericIcon import *

###########################################################################
# Class definition
###########################################################################

class sphereDB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #
        AFXDataDialog.__init__(self, form, 'Create Spheric inclusions', self.OK|self.CANCEL)
            
        AFXNote(self, 'This GUI is used for the definition of a part with random spheric inclusion.')
                      
        hf = FXHorizontalFrame(self, LAYOUT_FILL_X, 0,0,0,0, 0,0,0,0)
        
        gb = FXGroupBox(hf, 'Parameters', LAYOUT_FILL_Y|FRAME_GROOVE)
        va = AFXVerticalAligner(gb)
        AFXTextField(va, 12, 'Name:', form.nameKw, 0)
        combo = AFXComboBox(va, 8, 3, 'Mask:', form.maskKw, 0)
        combo.appendItem('Rectangle', 1)
        combo.appendItem('Circle', 2)
        combo.appendItem('Ellipse', 3)
        AFXTextField(va, 12, 'Width (w):', form.widthKw, 0)
        AFXTextField(va, 12, 'Height (h):', form.heightKw, 0)
        AFXTextField(va, 12, 'Max distance (dmax):', form.dmaxKw, 0)
        AFXTextField(va, 12, 'Number of cells:', form.npointKw, 0)
        AFXTextField(va, 12, 'Thickness (ep):', form.epKw, 0)
        FXCheckButton(va,  'Remove sphere on the bound:',form.removeBKw,0)

        gb = FXGroupBox(hf, 'Spheric Geom', LAYOUT_FILL_Y|FRAME_GROOVE)
        icon = FXXPMIcon(getAFXApp(), spheric_xpm)
        FXLabel(gb, '', icon, pl=0, pr=0, pt=0, pb=0)


