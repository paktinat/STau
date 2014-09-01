#!/usr/bin/env python

from ROOT import TFile, TH1D, TH2D, TTree, TGraphAsymmErrors, TCanvas, gROOT, TMultiGraph, TLatex, gStyle, gDirectory, TLegend
import os
import numpy as n
import math
import array
import datetime

LogFile = open( 'LogFile' , 'w' ) 

class Cuts:
    def __init__(self ):
        self.File = None
        self.fTotalNumbers = TFile('/afs/cern.ch/work/h/hbakhshi/STau/CMSSW_6_1_1/src/HiggsAnalysis/all_Histos.root')
        self.hNEvents = self.fTotalNumbers.Get('h_SMSEvents')

        PDFCTEQ66 = {100:[5823.40, 0.0 , +3.4 , -.6 , -3.2],
                     125:[2434.10, 0.0, 3.6 , -.6 , -3.5 ],
                     150:[1194.60, 0.3, +3.9 , -.5 , -3.8],
                     175:[649.58, 0.3 , 4.2, -.5 , -4. ] ,
                     200:[379.24, 0.4, 4.5 , -.4 , -4.4] ,
                     225:[233.41 , 0.5, 5.0 , -.3 , -4.4],
                     250:[149.86 , 0.3 ,5.1 , -.4 , -4.8],
                     275:[99.27, 0.1 , 5.5 , -.4 , -5. ],
                     300:[67.51 , 0.0 , 5.9 , -.2 , -5.1],
                     325:[46.97 , 0.1 , 6.1 , -.2 , -5.5],
                     350:[33.28 , 0.0 , 6.4 , -.2 , -5.6],
                     375:[23.95 , 0.0 , 7.0 , -.1 , -5.7],
                     400:[17.51 , 0.0 , 6.8 , -.3 , -6.3],
                     425:[12.93 , 0.0 , 7.5 , -.3 , -6.1],
                     450:[9.66, 0.0 , 7.5 , -.5 , -6.7],
                     475:[7.28 , 0.1 , 7.8, -1 , -6.8],
                     500:[5.53 , 0.0 , 8.1 , -.9 , -7.0]
                     }

        self.gXSections = TGraphAsymmErrors( len(PDFCTEQ66) )
        step = 25
        cc = 0
        for mass in range(100 , 500 , step):
            self.gXSections.SetPoint( cc , mass , PDFCTEQ66[mass][0] )
        
            errh = math.hypot( PDFCTEQ66[mass][1] , PDFCTEQ66[mass][2] )
            errl = math.hypot( PDFCTEQ66[mass][3] , PDFCTEQ66[mass][4] )

            self.gXSections.SetPointError( cc , errl , errh , step/2. , step/2. )
            
            cc += 1

        self.hWeights = TH2D("hWeights" , "Weights" , 20 , 100 , 500 , 25 , 0 , 500 )

        self.hMedian = TH2D("hMedian" , "Median" , 20 , 100 , 500 , 25 , 0 , 500 )
        self.hMedianCurve = TH2D("hMedianCurve" , "Expected Exclusion" , 20 , 100 , 500 , 25 , 0 , 500 )

        self.hp1Sigma = TH2D("hp1Sigma" , "+1#sigma" , 20 , 100 , 500 , 25 , 0 , 500 )
        self.hp2Sigma = TH2D("hp2Sigma" , "+2#sigma" , 20 , 100 , 500 , 25 , 0 , 500 )
        self.hm1Sigma = TH2D("hm1Sigma" , "-1#sigma" , 20 , 100 , 500 , 25 , 0 , 500 )
        self.hm2Sigma = TH2D("hm2Sigma" , "-2#sigma" , 20 , 100 , 500 , 25 , 0 , 500 )
        self.hObserved = TH2D("hObserved" , "Observed" , 20 , 100 , 500 , 25 , 0 , 500 )

        self.hEmpty  = TH2D("hEmpty" , "Empty" , 20 , 100 , 500 , 25 , 0 , 500 )

        self.FixTH2D( self.hMedianCurve , 3,2,1)
        self.FixTH2D( self.hm1Sigma , 2,4,2)
        self.FixTH2D( self.hm2Sigma , 2,4,3)
        self.FixTH2D( self.hp1Sigma , 2,4,2)
        self.FixTH2D( self.hp2Sigma , 2,4,3)
        self.FixTH2D( self.hObserved , 3,1,1)

        self.hTimeMonitor = TH2D("hTimeMonitor" , "Time" , 20 , 100 , 500 , 25 , 0 , 500 )

        

    def getVarName(self , var):
        varNames = { 0:'METModMPZMod',
                     1:'METModPPZMod',
                     2:'TauPt',
                     3:'EleTauPt',
                     4:'MET',
                     5:'EleMT',
                     6:'MT2'
                     }
        return varNames[var]

    class BinResults:
        def __init__(self , index , hBKG , hSignal , hSigNorm ):
            self.Index = index
            self.BKG = hBKG
            self.DY = self.BKG.GetBinContent(1)
            self.WJets = self.BKG.GetBinContent(2)
            self.Top = self.BKG.GetBinContent(3)
            self.QCD = self.BKG.GetBinContent(4)
            self.Data = self.BKG.GetBinContent(5)
            self.TotalBkg= self.DY + self.WJets + self.Top + self.QCD
            
            self.Signal = hSignal
            self.hSignalNorm = hSigNorm

            self.hNSignals = TH2D("hNSignals_%d" % (index) , "NSignals Bin#%d" % (index) , 20 , 100 , 500 , 25 , 0 , 500 )


        def BinName(self):
            return 'bin%d' % (self.Index)

        def GetIndices(self):
            return '\t'.join( s for s in [self.BinName()]*5 )

        def getSystematics(self, systematics):
            return '\t'.join( str(ss) for ss in systematics )
            
        def getRates(self , MGlu, MLSP , weight):
            mgluBin = self.Signal.GetXaxis().FindBin( MGlu )
            mlspBin = self.Signal.GetYaxis().FindBin( MLSP )
            currentSignal = self.Signal.GetBinContent( mgluBin , mlspBin )

            mgluBin = self.hNSignals.GetXaxis().FindBin( MGlu )
            mlspBin = self.hNSignals.GetYaxis().FindBin(MLSP)

            signals = weight*currentSignal

            self.hNSignals.SetBinContent( mgluBin, mlspBin , signals )

            signals = self.hSignalNorm.GetBinContent( mgluBin , mlspBin )

            return '%(sig)f\t%(dy)f\t%(w)f\t%(top)f\t%(qcd)f' % {'sig':signals,
                                                                 'dy':self.DY,
                                                                 'w':self.WJets,
                                                                 'top':self.Top,
                                                                 'qcd':self.QCD
                                                                 }

    def makecardfile( self , MGlu, MLSP , weight , systematics ) :
        filename = 'cardfile'
        cardfile = open( filename , 'w' )
        print >> cardfile, ('imax\t%d\tnumber of channels' % (self.NBins) ).expandtabs()
        print >> cardfile, 'jmax\t4\tnumber of backgrounds'
        print >> cardfile, 'kmax\t1\tnumber of nuisance parameters (sources of systematic uncertainties)'
        print >> cardfile, "----------------------------------------------"  
        print >> cardfile, ('bin\t' + '\t'.join( a.BinName() for a in self.BinsInfo  ) ).expandtabs()
        print >> cardfile, ('observation\t' + '\t'.join( str(a.TotalBkg) for a in self.BinsInfo) ).expandtabs()
        print >> cardfile, "----------------------------------------------"  
        print >> cardfile, ('bin\t' + '\t'.join( a.GetIndices() for a in self.BinsInfo) ).expandtabs()
        print >> cardfile, ('process\t' + '\t'.join( 'SUSY\tDY\tW\tTOP\tQCD' for a in self.BinsInfo) ).expandtabs()
        print >> cardfile, ('process\t' + '\t'.join( '0\t1\t2\t3\t4' for a in self.BinsInfo) ).expandtabs()
        print >> cardfile, ('rate\t' + '\t'.join( a.getRates(MGlu, MLSP , weight) for a in self.BinsInfo) ).expandtabs()
        print >> cardfile, "----------------------------------------------"  
        print >> cardfile, ('syst\tlnN\t' + '\t'.join( a.getSystematics(systematics) for a in self.BinsInfo) ).expandtabs()
        cardfile.close()
        return filename



    def run_cuts( self , cuts , binVariable , bins , name ):
        #if self.File:
        #    self.File.Close()

        self.Name  = name

        METmPZ , METpPZ , TauPt , zPt , MET , EleMT , MT2 = cuts

        command =  ' '.join( str(i) for i in ["DrawHists" , -100 , name + ".root" ] + cuts + [ binVariable ] + bins ) + " >& " + name + "_out"
        print command
        self.OutStreamDrawHist = os.system(command )

        self.File = TFile.Open( name + '.root'  )

        self.CutsStr = ''
        for nCut in range( 0, len(cuts) ):
            if cuts[nCut] == -1 or cuts[nCut] == -999 :
                continue
            if not self.CutsStr == '' :
                self.CutsStr += " And "
            if binVariable == nCut :
                self.CutsStr += '%s > %f' % ( self.getVarName(nCut) , bins[0] )
            else:
                self.CutsStr += '%s > %f' % ( self.getVarName(nCut) , cuts[nCut] )

        if binVariable == -1:
            self.NBins = 1
            self.BinsStr = "No Binning"
        else:
            self.BinVar = self.getVarName( binVariable )
            self.NBins = len(bins)
            self.BinsStr = 'No binning'
            if not self.NBins == 1:
                self.BinsStr = 'bins :' + self.BinVar
                allBins = bins 
                self.BinsStr += ' | '.join( '>' + str(bbb) for bbb in allBins )

        self.BinsInfo = []

        rangeOfBins = range(1, self.NBins+1)
        if self.NBins == 1:
            rangeOfBins = [0]
        for binid in rangeOfBins:
            print >> LogFile, binid
            self.File.ls()
            hhBKG = self.File.Get("hBKG_%d"%(binid)).Clone("hAllBKG___%d"%(binid))
            hhSignal = self.File.Get("hSignal_%d"%(binid)).Clone("hMGluMLSP___%d"%(binid))
            hhsignorm = self.File.Get("hSignal_%d_Normalized"%(binid)).Clone("hMGluMLSP___%d_Normalized"%(binid))
            bbinin = self.BinResults(binid , hhBKG , hhSignal , hhsignorm )
            
            self.BinsInfo.append( bbinin )

    def FixTH2D(self, h1 , lw , lc , ls , clone=False):
        h = h1
        if clone:
            h = h1.Clone( "Cloned" + h1.GetName() )

        h.SetLineWidth(lw)
        h.SetLineColor(lc)
        h.SetLineStyle(ls)

        nbnx = h.GetNbinsX()
        nbny = h.GetNbinsY()
        for i in range(1 , nbnx+1):
            for j in range(1 , nbny+1):
                z = h.GetBinContent(i, j)
                if z == 0:
                    h.SetBinContent(i, j, 10000)
        h.SetMinimum(0.01)
        h.SetMaximum(100)
        h.SetContour(1)
        h.SetContourLevel(0, 1)
        
        return h

    def GetContour(self, histo):
        ctemp = TCanvas('temp' , 'tmp canvas' , 800 , 600 )
        ctemp.cd()
        
        histo.Draw('CONT Z LIST')
        ctemp.Update()
        conts = gROOT.GetListOfSpecials().FindObject("contours")

        ret = TMultiGraph()

        if conts.GetSize() == 0: 
            ctemp.Close()
            return ret
        if conts.At(0).GetSize() == 0:
            ctemp.Close()
            return ret

        for ggid in range(0,conts.At(0).GetSize()) :
            graph = conts.At(0).At(ggid).Clone( histo.GetName() + "_Con1Graph_%d" % (ggid) )
            graph.SetTitle( histo.GetTitle() )

            graph.SetLineWidth( histo.GetLineWidth() )
            graph.SetLineColor( histo.GetLineColor() )
            graph.SetLineStyle( histo.GetLineStyle() )

            graph.Write()
            ret.Add(graph)

        ctemp.Close()
        return ret

    def Write(self):
        self.hWeights.Write()
        for m in self.BinsInfo:
            m.hNSignals.Write()
        self.hMedian.Write()
        self.hp1Sigma.Write()
        self.hp2Sigma.Write()
        self.hm1Sigma.Write()
        self.hm2Sigma.Write()
        self.hObserved.Write()
        self.hTimeMonitor.Write()

        gStyle.SetOptTitle(0)

        cc = TCanvas('All' , 'All' , 800 , 600 )

        mg = TMultiGraph()

        self.hMedian.SetStats(False)

        self.hMedian.Draw('colz')
        mg.Add( self.GetContour( self.hMedianCurve ) )
        mg.Add( self.GetContour( self.hm1Sigma) )
        mg.Add( self.GetContour( self.hm2Sigma ))
        mg.Add( self.GetContour( self.hp1Sigma ))
        mg.Add( self.GetContour( self.hp2Sigma ))
        mg.Add( self.GetContour( self.hObserved))
        cc.cd()
        mg.Draw('C')
        cc.Update()

        leg = cc.BuildLegend()
        for lentry in leg.GetListOfPrimitives():
            lentry.SetOption('l')

        l = TLatex()
        l.SetTextSize(0.03)
        l.SetNDC()
        l.DrawLatex( .03 , .97 , self.CutsStr ).Draw()
        l.DrawLatex( .03 , .94 , self.BinsStr ).Draw()


        cc.SaveAs( gDirectory.GetName() + ".C" )

        #cc.PaintTextNDC( 
        #cc.PaintTextNDC(
        cc.Write()

    def getFakeValue(self, XSection , MGlu , MLSP ):
        if MGlu < 300:
            if MGlu - MLSP > 100 :
                if MLSP >50: 
                    return .5

        return 5.

    def upperlimit( self , MGlu , MLSP , systematics , fake=False):
        mgluBin = self.hNEvents.GetXaxis().FindBin( MGlu )
        mlspBin = self.hNEvents.GetYaxis().FindBin( MLSP )
        TotalNumberOfSignals = self.hNEvents.GetBinContent( mgluBin , mlspBin )
        
        if TotalNumberOfSignals < 100:
            return

        XSection = self.gXSections.Eval( MGlu )
        Lumi = 19.6
        self.WeightSignal = Lumi*XSection / TotalNumberOfSignals

        #print( '%d -- %d , W:%f' % (mgluBin , mlspBin , self.WeightSignal) )

        cardfilename = self.makecardfile( MGlu , MLSP , self.WeightSignal , systematics  )

        mgluBin = self.hWeights.GetXaxis().FindBin( MGlu )
        mlspBin = self.hWeights.GetYaxis().FindBin(MLSP)

        self.hWeights.SetBinContent( mgluBin, mlspBin , self.WeightSignal )

        before = datetime.datetime.now()
        if not fake:
            self.CombineStream = os.system( ' '.join( str(i) for i in ['combine' , '-M' , 'Asymptotic' , cardfilename ] ))
        after = datetime.datetime.now()
        
        delta = after - before
        self.hTimeMonitor.SetBinContent( mgluBin , mlspBin , delta.total_seconds() )


        if fake:
            value = self.getFakeValue( XSection , MGlu , MLSP )

            self.hMedian.SetBinContent( mgluBin, mlspBin , value )
            self.hMedianCurve.SetBinContent( mgluBin, mlspBin , value )
            self.hp2Sigma.SetBinContent( mgluBin, mlspBin , value+2*value )
            self.hp1Sigma.SetBinContent( mgluBin, mlspBin , value+1*value )
            self.hm2Sigma.SetBinContent( mgluBin, mlspBin , value-.2*value )
            self.hm1Sigma.SetBinContent( mgluBin, mlspBin , value-.1*value )
            self.hObserved.SetBinContent( mgluBin, mlspBin , value+.2*value )
        else:
            file = TFile.Open( 'higgsCombineTest.Asymptotic.mH120.root' )
            limit = file.Get('limit')

            blimits = array.array('d', [0.0])
            quantileExpected = array.array( 'f' , [0.0])

            limit.SetBranchAddress( 'limit' , blimits )
            limit.SetBranchAddress( 'quantileExpected' , quantileExpected )
            
            for i in range( 0, limit.GetEntries() ):
                limit.GetEntry( i )
                q = quantileExpected[0]
                l = blimits[0]
                if q > 0.9 :
                    self.hp2Sigma.SetBinContent( mgluBin, mlspBin , l )
                elif q > 0.8 :
                    self.hp1Sigma.SetBinContent( mgluBin, mlspBin , l )
                elif q > 0.4 :
                    self.hMedian.SetBinContent( mgluBin, mlspBin , l )
                    self.hMedianCurve.SetBinContent( mgluBin, mlspBin , l )
                elif q > 0.1 :
                    self.hm1Sigma.SetBinContent( mgluBin, mlspBin , l )
                elif q > 0.0 : 
                    self.hm2Sigma.SetBinContent( mgluBin, mlspBin , l)
                elif q > -2.0 : 
                    self.hObserved.SetBinContent( mgluBin, mlspBin , l)
                #print ( '%d : %f' % (i , l) )

            os.system('/bin/rm -f roostats-*')
            file.Close()
        
            #os.system('/bin/rm -f higgsCombineTest.Asymptotic.mH120.root')
            #os.system('/bin/rm -f %s' % (cardfilename) )

if __name__ == "__main__":
    fout = TFile('ExclusionOutput.root', 'recreate')


    for i in range(0,1):
        m = Cuts()

        #METmPZ , METpPZ , TauPt , zPt , MET , EleMT , MT2 
        try:
            if i==0:
                m.run_cuts( [-50,175,-1,-1,-1,-1,-1] , 6 , [50,70,90,110,130] , 'SyncCuts' )
            elif i==1:
                m.run_cuts( [-999,-1,40,-1,-1,-1,20] , -1 , [50] , 'MET_TauPT' )
            elif i==2:
                m.run_cuts( [-999 , -1 , 50 , -1 , -1 , -1 , 10] , 4 , [  40 , 60 , 80] , 'METBin_TauPt'  )
            elif i==3:
                continue
                m.run_cuts( [-50 , 150 , -1 , -1 , -1 , -1 , 50] , -1 , [] , 'MT2_MET_PZ' )
            elif i==4:
                m.run_cuts( [-50 , -1 , -1 , -1 , -1 , -1 , 50] , 1 , [180 , 200 , 220 , 250] , 'MT2_Met_Pz_Binned' )
            elif i==5:
                m.run_cuts( [-999 , -1 , -1 , -1 , -1 , -1 , -1] , 6 , [80 , 100 , 120 ], 'MT2Binned' )
            elif i==6:
                m.run_cuts( [-50 , 175 , -1 , -1 , -1 , -1 , -1] , 6 , [ 80, 100 , 120 ], 'MT280Binned_Pz_MET' )
            elif i==7:
                m.run_cuts( [-50 , 175 , -1 , -1 , -1 , -1 , -1] , 6 , [ 50 , 80 , 100 , 120 ], 'MT250Binned_Pz_MET' )

        except:
            print >> LogFile, m.Name + " Error " 
            continue

        print >> LogFile , m.Name 

        # for mchargino in range(100 , 500 , 20):
        #     LogFile.write( "\n%d" % (mchargino) )
        #     for mlsp in range( 0 , 500 , 20 ):
        #         LogFile.write("|%d" % (mlsp) )
        #         LogFile.flush()
        #         m.upperlimit(mchargino , mlsp , [ 1.1 ]*5 )
        
        m.upperlimit( 200 , 0 , [1.1]*5 )
                

        fout.mkdir( m.Name  ).cd()
        m.Write()

    fout.Close()
