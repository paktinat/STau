from array import *

class sms():

    def __init__(self, modelname):
        if modelname.find("TChipChimSTauSnu") != -1: self.TChipChimSTauSnu()


    def TChipChimSTauSnu(self):
        # model name                   
        self.modelname = "TChipChimSTauSnu"
        # decay chain
        self.label= "#scale[0.8]{pp #rightarrow #tilde{#chi^{+}} #tilde{#chi^{-}},#tilde{#chi} #rightarrow #tilde{#tau} #rightarrow #tau #chi^{0}_{1}}";
        # scan range to plot
        self.Xmin = 100
        self.Xmax = 500
        self.Ymin = 0
        self.Ymax = 500
        # produce sparticle
        self.sParticle = "m_{#tilde{#chi}^{#pm}_{1}} (GeV)"
        # LSP
        self.LSP = "m_{#tilde{#chi}^{0}_{1}} (GeV)"        
        # diagonal position: mLSP = mgluino - 2mtop 
        self.diagX = array('d',[200,100,500])
        self.diagY = array('d',[0,100,500])        



