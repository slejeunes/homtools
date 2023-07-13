from abaqusGui import *
from voronoiIcon import *

###########################################################################
# Class definition
###########################################################################

class VoronoiDB(AFXDataDialog):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #
        AFXDataDialog.__init__(self, form, 'Create Voronoi', self.OK|self.CANCEL)
            
        AFXNote(self, 'This GUI is used for the definition of a Voronoi part.')
                      
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
        AFXTextField(va, 12, 'Min distance (dmin):', form.dminKw, 0)
        AFXTextField(va, 12, 'Number of cells:', form.npointKw, 0)
        AFXTextField(va, 12, 'Thickness (ep):', form.epKw, 0)
        FXCheckButton(va,  'Is periodic:',form.periodicKw,0)

        gb = FXGroupBox(hf, 'Voronoi Geom', LAYOUT_FILL_Y|FRAME_GROOVE)
        icon = FXXPMIcon(getAFXApp(), voronoi_xpm)
        FXLabel(gb, '', icon, pl=0, pr=0, pt=0, pb=0)


class VoronoiDB2(AFXDataDialog):
    [
        ID_1
    ] = range(AFXDataDialog.ID_LAST, AFXDataDialog.ID_LAST+1) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #
        self.db=None
        self.form=form
        AFXDataDialog.__init__(self, form, 'Create Voronoi with constraint', self.OK|self.CANCEL)
            
        AFXNote(self, 'This GUI is used for the definition of a Voronoi part.')
                      
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
        AFXTextField(va, 12, 'Min distance (dmin):', form.dminKw, 0)
        AFXTextField(va, 12, 'Number of cells:', form.npointKw, 0)
        AFXTextField(va, 12, 'Filler ratio (in volume):', form.frationKw, 0)
 
        FXCheckButton(va,  'Is periodic:',form.periodicKw,0)

        gb = FXGroupBox(hf, 'Voronoi Geom', LAYOUT_FILL_Y|FRAME_GROOVE)
        icon = FXXPMIcon(getAFXApp(), voronoi_xpm)
        FXLabel(gb, '', icon, pl=0, pr=0, pt=0, pb=0)

        gb2 = FXGroupBox(self, 'Use previous diagram', LAYOUT_FILL_X|FRAME_GROOVE)
        hf2 = FXHorizontalFrame(gb2, LAYOUT_FILL_X, 0,0,0,0, 0,0,0,0)
        AFXTextField(hf2, 32, 'Select a file:', form.fileNameKw, 0)
        FXButton(hf2,'Browse...',None,self,self.ID_1)

        FXMAPFUNC(self,SEL_COMMAND,self.ID_1,VoronoiDB2.onCmdFile)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def onCmdFile(self,sender,sel,ptr):
        patterns= 'Cae Database (*.cae)'
        if not self.db:
          self.db= AFXFileSelectorDialog(self,'Select an input file',self.form.fileNameKw,None,
                  AFXSELECTFILE_EXISTING,patterns)
          self.db.create()

        self.db.showModal()
        return 1 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class VoronoiDB3(AFXDataDialog):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, form):

        # Construct the base class.
        #
        self.db=None
        self.form=form
        AFXDataDialog.__init__(self, form, 'Create Voronoi Aggregate', self.OK|self.CANCEL)
            
        AFXNote(self, 'This GUI is used for the definition of a Voronoi Aggregate.')
                      
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
        AFXTextField(va, 12, 'Particule Radius (r):', form.rmeanKw, 0)
        AFXTextField(va, 12, 'Initial Number of cells:', form.npointKw, 0)
        AFXTextField(va, 12, 'Thickness (ep):', form.epKw, 0)
        AFXTextField(va, 12, 'Filler ratio (in volume):', form.frationKw, 0)
 
        gb = FXGroupBox(hf, 'Voronoi Geom', LAYOUT_FILL_Y|FRAME_GROOVE)
        icon = FXXPMIcon(getAFXApp(), voronoi_xpm)
        FXLabel(gb, '', icon, pl=0, pr=0, pt=0, pb=0)


