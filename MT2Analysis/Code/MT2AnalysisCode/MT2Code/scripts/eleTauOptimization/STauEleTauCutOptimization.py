#!/usr/bin/env python
# Standard python import
import os
import sys    # exit
import time   # time accounting
from optparse import OptionParser

import itertools #combinatiories
import math
from array import array
import xml.etree.ElementTree as etree 

# Import ROOT classes
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut, gDirectory , TH1F , TFormula, std

# Import TMVA classes from ROOT
from ROOT import TMVA

os.environ['DISPLAY'] = '0'
os.environ['TMVASYS'] = os.environ['ROOTSYS']

# Default settings for command line arguments
DEFAULT_INFNAME  = "AfterBVetoMET_Trees.root"
DEFAULT_TREENAME = "bVeto"


VARIABLES = {"ElePt":['F'] , "EleMT":['F'] ,  "MT2":['F'] , "MCT":['F'] , "MT2Imb":['F'] ,"TauPt":['F'] , "DPt":['F'] , "MET":['F'] ,"EleTauPt":['F'] ,"JPTModMZPTMod":['F'] , 'METModPPZMod':['F'] , 'METModMPZMod':['F'] , 'ModMETmPZ':['F'] }

Category1 = [ 'EleTauPt' , 'METModPPZMod' , 'ModMETmPZ' ]
Category2 = [ 'MCT' , 'MT2' , 'MT2Imb' ]


SUSYCats = {"50-100":"SUSYCategory >= 1  && SUSYCategory < 2" ,  "100-150":"SUSYCategory >= 2  && SUSYCategory < 3" , 
            "150-200":"SUSYCategory >= 3  && SUSYCategory < 4" , "200-250":"SUSYCategory >= 4  && SUSYCategory < 5" , 
            "250-300":"SUSYCategory >= 5  && SUSYCategory < 6" , "300-350":"SUSYCategory >= 6  && SUSYCategory < 7" , 
            "350-400":"SUSYCategory >= 7  && SUSYCategory < 8" , "400-450":"SUSYCategory >= 8  && SUSYCategory < 9" , 
            "450-500":"SUSYCategory >= 9  && SUSYCategory < 10" }

class MethodInfo :

    def __init__(self , titDir , name , title):
          self.methodName = name
          self.methodTitle = title

          hname = "MVA_" + self.methodTitle

          self.sig = titDir.Get( hname + "_S" )
          self.bgd = titDir.Get( hname + "_B" )
          self.origSigE = titDir.Get( hname + "_effS" )
          self.origBgdE = titDir.Get( hname + "_effB" )
          self.sigE = None
          self.bgdE = None
          self.purS = None
          self.sSig = None
          self.effpurS = None
          self.maxSignificance = 0.0
          self.maxSignificanceErr = 0.0

          self.maxbin = -1
          self.fNSignal = 1000
          self.fNBackground = 1000

          self.SetResultHists()
          self.UpdateSignificanceHists('S/sqrt(B)')

    def SetResultHists(self):
  
        pname    = "purS_"         + self.methodTitle;
        epname   = "effpurS_"      + self.methodTitle;
        ssigname = "significance_" + self.methodTitle;

        self.sigE = self.origSigE.Clone("sigEffi");
        self.bgdE = self.origBgdE.Clone("bgdEffi");
      
        nbins = self.sigE.GetNbinsX();
        low = self.sigE.GetBinLowEdge(1);
        high = self.sigE.GetBinLowEdge(nbins+1);
        self.purS = TH1F(pname, pname, nbins, low, high);
        self.sSig  = TH1F(ssigname, ssigname, nbins, low, high);
        self.effpurS = TH1F(epname, epname, nbins, low, high);        

        
        self.sigE.SetFillStyle( 0 );
        self.bgdE.SetFillStyle( 0 );
        self.sSig.SetFillStyle( 0 );
        self.sigE.SetLineWidth( 3 );
        self.bgdE.SetLineWidth( 3 );
        self.sSig.SetLineWidth( 3 );

        self.purS.SetFillStyle( 0 );
        self.purS.SetLineWidth( 2 );
        self.purS.SetLineStyle( 5 );
        self.effpurS.SetFillStyle( 0 );
        self.effpurS.SetLineWidth( 2 );
        self.effpurS.SetLineStyle( 6 );
  
    def UpdateSignificanceHists(self, significance_formula , doprint = False):
        f = TFormula("sigf", significance_formula.replace("S", "x").replace("B", 'y') )
        cname = "Classifier"
        maxLenTitle = len(cname)

        maxSig    = -1.;
        maxSigErr = -1.;
        for i in range( 1 , self.origSigE.GetNbinsX()+1 )  :
            eS = self.origSigE.GetBinContent( i )
            S = eS * self.fNSignal
            B = self.origBgdE.GetBinContent( i ) * self.fNBackground
            purity = 0.0
            if not S+B ==0 :
                purity = S/(S+B)
            self.purS.SetBinContent( i, purity)
         
            sig = f.Eval(S,B);
            if sig > maxSig :
                maxSig = sig;
                if significance_formula == "S/sqrt(B)" and not B == 0 and not S == 0:
                    maxSigErr = sig * math.sqrt( 1./S + 1./(2.*B))
                                
            self.sSig.SetBinContent( i, sig )
            self.effpurS.SetBinContent( i, eS*self.purS.GetBinContent( i ) )
    
      
        self.maxSignificance    = self.sSig.GetMaximum()
        self.maxbin = self.sSig.GetMaximumBin();
        self.maxSignificanceErr = maxSigErr if maxSigErr > 0 else 0
        self.sSig.Scale(1/self.maxSignificance);
        self.optimalcut = self.sSig.GetXaxis().GetBinCenter( self.maxbin )
        self.NSig = self.origSigE.GetBinContent( self.maxbin )*self.fNSignal
        self.NBkg = self.origBgdE.GetBinContent( self.maxbin )*self.fNBackground
        self.EffSig = self.origSigE.GetBinContent( self.maxbin ) 
        self.EffBkg = self.origBgdE.GetBinContent( self.maxbin )
        
        if doprint:
            str = "%(name)s   (  #signal, #backgr.)  Optimal-cut  %(sig_formula)s      NSig      NBkg   EffSig   EffBkg" % {'name':cname , 'sig_formula':significance_formula }
            print str

            str = "%(name)s   (  %(signal)d, %(bkg)d )  %(optimalcut)f  %(sig)f+- %(sigerr)f      %(NSig)f       %(NBkg)f   %(EffSig)f   %(EffBkg)f" % {'name':self.methodName , 'signal':self.fNSignal , 'bkg':self.fNBackground , 'optimalcut':self.optimalcut , 'sig':self.maxSignificance , 'sigerr':self.maxSignificanceErr , 'NSig':self.NSig , 'NBkg':self.NBkg , 'EffSig':self.EffSig , 'EffBkg':self.EffBkg }
            print str

    def ReadCuts( self, xmlfilename):
        tree = etree.parse(xmlfilename)
        root = tree.getroot()
        variables = []
        self.Variables = {}
        vars = root.find( 'Variables' )
        for var in vars:
            varname = var.attrib[ 'Title' ]
            variables.append( varname )
        weights = root.find( "Weights")
        for bins in weights:
            effs = float(bins.attrib['effS'])
            if round(self.EffSig,2) == round( effs , 2) :
                cuts = bins[0].attrib
                i = 0
                for varname in variables:
                    self.Variables[ varname ] = []
                    self.Variables[ varname ].append( round( float( cuts[ 'cutMin_%d' % (i) ] ) , 2 ) )
                    self.Variables[ varname ].append( round( float( cuts[ 'cutMax_%d' % (i) ] ) , 2 ) )
                    i+=1

    def Print( self ):
        ret = []
        ret.append( '%.2f' % ( 100* self.EffSig ) )
        ret.append( '%.2f' % ( 100* self.EffBkg ) )
        ret.append( '%.2f' % ( self.maxSignificance ) )
        ret.append( '%.2f' % ( self.maxSignificanceErr ) )
        for var in self.Variables:
            lowercut = self.Variables[var][0]
            uppercut = self.Variables[var][1]
            if lowercut < -10000000:
                ret.append( '%(var)s < %(val).2f' % {'var':var , 'val':uppercut} )
            elif uppercut > 10000000:
                ret.append( '%(val).2f < %(var)s' % {'var':var , 'val':lowercut} )
            else:
                ret.append( '%(val).2f < %(var)s < %(val2).2f' % {'var':var , 'val':lowercut , 'val2':uppercut} )

        return ','.join( [str(item) for item in ret] )

if __name__ == "__main__":
    infname     = DEFAULT_INFNAME
    treeName    = DEFAULT_TREENAME
    verbose     = False

    parser = OptionParser()
    parser.add_option("-n", "--ncombs", dest="ncombs" , type='int' , default=-1)    
    parser.add_option("-s", "--signal", dest="signal" , type='int')
    parser.add_option("-d", "--dir", dest="dir" , default='.' )
    (options, args) = parser.parse_args() 

    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    gROOT.Macro       ( "TMVAlogon.C" )    
    gROOT.LoadMacro   ( "TMVAGui.C" )

    input = TFile.Open( infname )
    TheTree = input.Get( treeName )

    catIndex = -1
    for (SUSCatName,SUSCatCut) in SUSYCats.iteritems():
        print SUSCatName + " : "
        catIndex += 1
        if catIndex == options.signal:
            print "ok"
        else:
            print "skipped"
            continue

        for COMB_L in range(2,5):
            if options.ncombs > -1:
                if not COMB_L == options.ncombs:
                    continue
                
            outtxtfile = open( 'csvfiles/%s_%d.csv' % (SUSCatName, COMB_L) , 'w' )

            ALLComs = itertools.combinations( VARIABLES , COMB_L )
            print str(COMB_L) + ":"

            ccc = 0
            for comb in ALLComs:
                ccc +=1
                Name = "TMVA_" + SUSCatName
                for cut in comb:
                    Name += "_" + cut

                print "Combination #%d : %s " % (ccc , Name)

                nCategory1 = 0
                nCategory2 = 0
                for cut in comb:
                    if cut in Category1:
                        nCategory1 += 1
                    if cut in Category2:
                        nCategory2 += 1
                if (nCategory1 > 1) or (nCategory2 > 2):
                    print "two members from the same category, skipped"
                    continue


                outfname = options.dir + "/" + Name + ".root"
                outputFile = TFile( outfname, 'RECREATE' )

                factory = TMVA.Factory( Name, outputFile, 
                                        "!V:Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" )

                # Set verbosity
                factory.SetVerbose( verbose )

                #it crashes when these variables are being studied
                forbidden_vars = [ 'METModPPZMod' , 'JPTModMZPTMod' , 'MET' , 'METModMPZMod' ]
                forbidden_vars2 = [ 'EleTauPt' , 'MET' , 'METModMPZMod' ]
                n_forbidden_vars = 0
                n_forbidden_vars2 = 0

                selection_cut = ''
                for var in VARIABLES:
                    if var in comb:
                        if not selection_cut == '':
                            selection_cut += " && "
                        factory.AddVariable( var, VARIABLES[ var ][0] )
                        selection_cut += var + " > -1000000"
                        if var in forbidden_vars:
                            n_forbidden_vars += 1
                        if var in forbidden_vars2:
                            n_forbidden_vars2 += 1
                    else:
                        factory.AddSpectator( var )


                if n_forbidden_vars > 2 or n_forbidden_vars2 == 3:
                    continue

                SUSYSignalCut = TCut( SUSCatCut )
                BKGCut = TCut( "SUSYCategory < 0 && SUSYCategory > -10" )

                factory.SetInputTrees( TheTree , SUSYSignalCut  , BKGCut  )
                factory.SetSignalWeightExpression("W")
                factory.SetBackgroundWeightExpression("W")

                mycutSig = TCut( selection_cut ) 
                mycutBkg = TCut( selection_cut )
                try:
                    err = ['err']*4
                    err += selection_cut.split("&&")
                    print ','.join( err )
                    factory.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" )
                except:
                    err = ['err']*4
                    err += selection_cut.split("&&")
                    print >> outtxtfile, ','.join( err )
                    continue


                factory.BookMethod( TMVA.Types.kCuts, "Cuts",
                                   "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" )

                # Train MVAs
                factory.TrainAllMethods()

                # Test MVAs
                factory.TestAllMethods()

                # Evaluate MVAs
                factory.EvaluateAllMethods()    

                # Save the output.
                outputFile.Close()

                fin = TFile( outfname , 'read') 
                fin.cd('Method_Cuts/Cuts/')
                mInfo = MethodInfo( gDirectory , "Cuts" , "Cuts" )
                mInfo.ReadCuts('weights/' + Name + '_Cuts.weights.xml')
                print >> outtxtfile, mInfo.Print()
                fin.Close()
            outtxtfile.close()
    input.Close()
