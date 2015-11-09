from ROOT import TCanvas, TH1, TFile, TPad, TEfficiency, TAxis, kBlack, kBlue , kRed,TF1

from optparse import OptionParser
import fnmatch
import os


class EffCanvas:
    def __init__(self , fileName):
        self.File = TFile.Open(fileName)
        #self.File.ls()
        self.Canvas = self.File.Get("muTau_log_comp_overlayc_ratio")

        if self.Canvas == None :
            self.Canvas = self.File.Get("muTau_log_comp_overlay")

        self.Pad = self.Canvas.GetListOfPrimitives().At(0)

        for i in range( 0 , self.Pad.GetListOfPrimitives().GetSize() ):
            if self.Pad.GetListOfPrimitives().At( i ).GetName() == "h2_copy" :
                self.Histo = self.Pad.GetListOfPrimitives().At( i )


class Efficiency:
    def __init__(self, Trigger , Sample , Dir , Initial="HighPtTau"):

        if Sample.find("Endcap") > -1:
            self.Color = kRed
        elif Sample.find("Barrel") > -1:
            self.Color = kBlue
        else:
            self.Color = kBlack

        if Sample.find("Loose") > -1:
            if Sample.find("Endcap") > -1 :
                self.Color = kBlue
            if Sample.find("Barrel") > -1 :
                self.Color = kRed    
            
        filenameformat = Initial + "{TrgName}Trg{Sample}{AllOrPass}*.root"

        PassFile_ = filenameformat.format( TrgName=Trigger , Sample=Sample , AllOrPass="Pass" )
        PassFile = ""
        AllFile_ = filenameformat.format( TrgName=Trigger , Sample=Sample , AllOrPass="All"  )
        AllFile = ""
        self.Title = "Efficiency of {TrgName} in {Sample} sample".format( TrgName=Trigger , Sample=Sample )
        
        for file in os.listdir( Dir  ):
            if fnmatch.fnmatch(file , PassFile_ ):
                PassFile = Dir + "/" + file
            if fnmatch.fnmatch(file , AllFile_ ):
                AllFile = Dir + "/" + file

        self.Pass = EffCanvas( PassFile )
        self.All  = EffCanvas( AllFile  )

    
#     def __init__(self , filePass , fileAll ):
#         self.Pass = EffCanvas( filePass )
#         self.All  = EffCanvas( fileAll  )

    def GetEfficiency( self , newNBinGroups , bins = None ):
        if bins == None:
            self.LastPass = self.Pass.Histo.Rebin( newNBinGroups , "newPassHisto" )
            self.LastAll = self.All.Histo.Rebin( newNBinGroups , "newAllHisto" )
        else :
            self.LastPass = self.Pass.Histo.Rebin( newNBinGroups , "newPassHisto" , bins )
            self.LastAll = self.All.Histo.Rebin( newNBinGroups , "newAllHisto" , bins )

        self.LastEff = TEfficiency( self.LastPass , self.LastAll )
        self.LastEff.SetLineColor( self.Color )
        self.LastEff.SetLineWidth( 3 )
        self.LastPass.Divide( self.LastAll )
        self.LastPass.SetLineColor( self.Color )
        self.LastPass.SetLineWidth( 3 )
        self.LastPass.GetYaxis().SetRangeUser( -.2 , 1.2 )

     
effMuTau = Efficiency( "muTau" , "MuEG" , "/home/hbakhshi/TrgEff" )
effBMuTau = Efficiency( "muTau" , "MuEGBarrel" , "/home/hbakhshi/TrgEff" )
effEMuTau = Efficiency( "muTau" , "MuEGEndcap" , "/home/hbakhshi/TrgEff" )

from array import array

newBins = array('d' , [20,30,40, 50 , 60, 80 , 100,120,140,200 , 450 ] )
nNewBins = 10
effBMuTau.GetEfficiency(nNewBins , newBins)

newBins = array('d' , [20,30,40,  50,  100,200 , 450 ] )
nNewBins = 6
effEMuTau.GetEfficiency(nNewBins , newBins)

newBins = array('d' , [0,50,140,450] )
nNewBins = 3
effMuTau.GetEfficiency(nNewBins , newBins)

c = TCanvas()

effMuTau.LastEff.SetTitle( effMuTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
effBMuTau.LastEff.SetTitle( effBMuTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
effEMuTau.LastEff.SetTitle( effEMuTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")


effBMuTau.LastEff.Draw("AP")
effEMuTau.LastEff.Draw("P SAME")
effMuTau.LastEff.Draw("P SAME")

#exit()
#eff.LastPass.Draw()

effEleTau = Efficiency( "eleTau" , "MuEG" , "/home/hbakhshi/TrgEff" )
effBEleTau = Efficiency( "eleTau" , "MuEGBarrel" , "/home/hbakhshi/TrgEff" )
effEEleTau = Efficiency( "eleTau" , "MuEGEndcap" , "/home/hbakhshi/TrgEff" )

effEleTau.GetEfficiency(4)
effBEleTau.GetEfficiency(4)
effEEleTau.GetEfficiency(4)

c2 = TCanvas()

effEleTau.LastEff.SetTitle( effEleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
effBEleTau.LastEff.SetTitle( effBEleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
effEEleTau.LastEff.SetTitle( effEEleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")

effEleTau.LastEff.Draw("AP")
effBEleTau.LastEff.Draw("P SAME")
effEEleTau.LastEff.Draw("P SAME")



effDoubleTau = Efficiency( "DoubleTau" , "MuEG" , "/home/hbakhshi/TrgEff" )
effDoubleTau.GetEfficiency(4)
c3 = TCanvas()

effDoubleTau.LastEff.SetTitle( effDoubleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
effDoubleTau.LastEff.Draw("AP")

SMeffMuTau = Efficiency( "muTau" , "SingleMu" , "/home/hbakhshi/TrgEff" , "TauEff" )
SMeffBMuTau = Efficiency( "muTau" , "SingleMuBarrel" , "/home/hbakhshi/TrgEff" , "TauEff"  )
SMeffEMuTau = Efficiency( "muTau" , "SingleMuEndcap" , "/home/hbakhshi/TrgEff" , "TauEff"  )

SMeffMuTau.GetEfficiency(10)
SMeffBMuTau.GetEfficiency(10)
SMeffEMuTau.GetEfficiency(10)

c4 = TCanvas()

SMeffMuTau.LastEff.SetTitle( SMeffMuTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
SMeffBMuTau.LastEff.SetTitle( SMeffBMuTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
SMeffEMuTau.LastEff.SetTitle( SMeffEMuTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")

SMeffMuTau.LastEff.Draw("AP")
SMeffBMuTau.LastEff.Draw("P SAME")
SMeffEMuTau.LastEff.Draw("P SAME")

SMeffEleTau = Efficiency( "eleTau" , "SingleEle" , "/home/hbakhshi/TrgEff" , "TauEff")
SMeffBEleTau = Efficiency( "eleTau" , "SingleEleBarrel" , "/home/hbakhshi/TrgEff" , "TauEff" )
SMeffEEleTau = Efficiency( "eleTau" , "SingleEleEndcap" , "/home/hbakhshi/TrgEff" , "TauEff" )

SMeffEleTau.GetEfficiency(4)
SMeffBEleTau.GetEfficiency(4)
SMeffEEleTau.GetEfficiency(4)

c5 = TCanvas()

SMeffEleTau.LastEff.SetTitle( SMeffEleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
SMeffBEleTau.LastEff.SetTitle( SMeffBEleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
SMeffEEleTau.LastEff.SetTitle( SMeffEEleTau.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")

SMeffEleTau.LastEff.Draw("AP")
SMeffBEleTau.LastEff.Draw("P SAME")
SMeffEEleTau.LastEff.Draw("P SAME")



SMeffMuTauL = Efficiency( "muTau" , "LooseSingleMu" , "/home/hbakhshi/TrgEff" , "TauEff" )
SMeffBMuTauL = Efficiency( "muTau" , "LooseSingleMuBarrel" , "/home/hbakhshi/TrgEff" , "TauEff"  )
SMeffEMuTauL = Efficiency( "muTau" , "LooseSingleMuEndcap" , "/home/hbakhshi/TrgEff" , "TauEff"  )

SMeffMuTauL.GetEfficiency(50)
SMeffBMuTauL.GetEfficiency(50)
SMeffEMuTauL.GetEfficiency(50)

c6 = TCanvas()

SMeffMuTauL.LastEff.SetTitle( SMeffMuTauL.Title + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
SMeffBMuTauL.LastEff.SetTitle( "#mu#tau trigger efficiency in (#mu + #tau) selected events in SingleMu dataset (Barrel)" + ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")
SMeffEMuTauL.LastEff.SetTitle( "#mu#tau trigger efficiency in (#mu + #tau) selected events in SingleMu dataset (Endcap)"+ ";p_{T}^{#tau};#epsilon" ) #options.Title + ";p_{T}^{#tau};#epsilon")

SMeffMuTauL.LastEff.Draw("AP")
SMeffBMuTauL.LastEff.Draw("P SAME")
SMeffEMuTauL.LastEff.Draw("P SAME")

EfficiencyCurves = TFile.Open("EfficiencyCurves.root" )

f1 = EfficiencyCurves.Get("Eff_ETauTrg_Tau_Data_2012_Barrel")
f2 = EfficiencyCurves.Get("Eff_ETauTrg_Tau_Data_2012_EndCap")

f1.Draw("SAME")
f2.Draw("SAME")

c7 = TCanvas()
SMeffBMuTauL.LastEff.Draw("AP")
SMeffBMuTau.LastEff.Draw("P SAME")
f1.Draw("SAME")


c8 = TCanvas()
SMeffEMuTauL.LastEff.Draw("AP")
SMeffEMuTau.LastEff.Draw("P SAME")
f2.Draw("SAME")

