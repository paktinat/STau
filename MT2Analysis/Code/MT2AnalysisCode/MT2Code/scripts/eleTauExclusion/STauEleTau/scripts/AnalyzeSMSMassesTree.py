from ROOT import TFile, TH1D, TH2D, TTree, TGraphAsymmErrors, TCanvas, gROOT, TMultiGraph, TLatex, gStyle, gDirectory, TLegend
import array

class allvars :
    def __init__(self , filename):
        self.File = TFile.Open( filename )
        self.Tree = self.File.Get("demo/tree")
        tree = self.Tree

        self.M22_1 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000022_1' , self.M22_1 )
        self.M22_2 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000022_2' , self.M22_2 )
        self.M22_avg = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000022' , self.M22_avg )
        self.M22_Best = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000022_Best' , self.M22_Best )
        self.M22_lv = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000022_lv' , self.M22_lv )

        self.M24_1 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000024_1' , self.M24_1 )
        self.M24_2 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000024_2' , self.M24_2 )
        self.M24_avg = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000024' , self.M24_avg )
        self.M24_Best = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000024_Best' , self.M24_Best )
        self.M24_lv = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000024_lv' , self.M24_lv )
        self.M24_15 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000024_15' , self.M24_15 )

        self.M15_1 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000015_1' , self.M15_1 )
        self.M15_2 = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000015_2' , self.M15_2 )
        self.M15_avg = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000015' , self.M15_avg )
        self.M15_lv = array.array( 'f' , [0.0])
        tree.SetBranchAddress( 'M1000015_lv' , self.M15_lv )

        self.MLSPComment = array.array( 'i' , [0])
        tree.SetBranchAddress( 'mlsp_comment' , self.MLSPComment )
        self.MCHComment = array.array( 'i' , [0])
        tree.SetBranchAddress( 'mch_comment' , self.MCHComment )

        self.LastMLSPComment = array.array( 'i' , [0])
        tree.SetBranchAddress( 'lastmlsp_comment' , self.LastMLSPComment )
        self.LastMCHComment = array.array( 'i' , [0])
        tree.SetBranchAddress( 'lastmch_comment' , self.LastMCHComment )


    def roundit( self,  a , b):
        c = a/b
        cr = round(c)
        return cr*b

    def diff2round( self , a , b):
        return abs( a - self.roundit( a , b ) )

    def getBetterLSP(self):
        if self.MLSPComment[0] == -1:
            diff1 = self.diff2round( self.M22_1[0] , 25.0 )
            diff2 = self.diff2round( self.M22_2[0] , 25.0 )
            if diff1 < diff2 :
                return self.M22_1[0]
            else:
                return self.M22_2[0]
        else:
            return self.MLSPComment[0]

    def getBetterChi(self):
        if self.MCHComment[0] == 85:
            diff1 = self.diff2round( self.M24_1[0] , 20.0 )
            diff2 = self.diff2round( self.M24_2[0] , 20.0 )
            if diff1 < diff2 :
                return self.M24_1[0]
            else:
                return self.M24_2[0]
        else:
            return self.MCHComment[0]

    def getBetterSTauMass(self):
        chi1 = 2*self.M15_1[0] - self.getBetterLSP() 
        chi2 = 2*self.M15_2[0] - self.getBetterLSP() 
        diff1 = self.diff2round( chi1 , 20.0 )
        diff2 = self.diff2round( chi2 , 20.0 )
        if diff1 < diff2 :
            return self.M15_1, chi1
        else:
            return self.M15_2, chi2

        
    def loop1(self):

        hlsp_12 = TH2D("hlsp_12" , "LSPMass;LSP1;LSP2" , 21 , -25 , 500 ,21 , -25 , 500 )
        hlsp_best1 = TH2D("hlsp_best1" , "LSPMass;BestLSP;LSP1" , 21 , -25 , 500 ,21 , -25 , 500 )
        hlsp_bestcomment = TH2D("hlsp_1comment" , "LSPMass;LSP1;FromComment" , 21 , -25 , 500 ,21 , -25 , 500 )
        hlsp_bestlastcomment = TH2D("hlsp_bestlastcomment" , "LSPMass;BestLSP;FromLastComment" , 21 , -25 , 500 ,21 , -25 , 500 )

        hch_12 = TH2D("hch_12" , "#chi mass;#chi_{1};#chi_{2}" , 22 , 80 , 520 ,22 , 80 , 520 )
        hch_best1 = TH2D("hch_best1" , "#chi mass;#chi_{best};#chi_{1}" , 22 , 80 , 520 ,22 , 80 , 520 )
        hch_bestcomment = TH2D("hch_bestcomment" , "#chi mass;#chi_{best};#chi_{comment}" , 22 , 80 , 520 ,22 , 80 , 520 )
        hch_commentfromstau = TH2D("hch_commentfromstau" , "#chi mass;#chi_{comment};#chi_{from stau}" , 22 , 80 , 520 ,22 , 80 , 520 )
        hch_commentavg = TH2D("hch_commentavg" , "#chi mass;#chi_{comment};#chi_{avg}" , 22 , 80 , 520 ,22 , 80 , 520 )

        h_chibestlspbest = TH2D("h_chibestlspbest" , ";#chi_{best};lsp_{best}" , 22 , 80 , 520 ,21 , -25 , 500 )
        h_chicommentlspcomment = TH2D("h_chicommentlspcomment" , ";#chi_{comment};lsp_{comment}" , 22 , 80 , 520 ,21 , -25 , 500 )
        h_chilastcommentlsplastcomment = TH2D("h_chilastcommentlsplastcomment" , ";#chi_{lastcomment};lsp_{lastcomment}" , 22 , 80 , 520 ,21 , -25 , 500 )

        h_betterchibetterlsp =  TH2D("h_betterchibetterlsp" , ";#chi_{better};lsp_{better}" , 22 , 70 , 510 ,21 , -35 , 490 )
        h_chifromstau_betterlsp =  TH2D("h_chifromstau_betterlsp" , ";#chi_{from stau};lsp_{better}" , 22 , 70 , 510 ,21 , -35 , 490 )

        

        for i in range(0 , self.Tree.GetEntries() ):
            self.Tree.GetEntry( i )

            if self.MCHComment[0] == -1 :
                self.MCHComment[0] = 85

            stau,chifromstau = self.getBetterSTauMass()

            hlsp_12.Fill( self.M22_1[0] , self.M22_2[0] )
            hlsp_best1.Fill( self.M22_Best[0] , self.M22_1[0] )
            hlsp_bestcomment.Fill( self.M22_1[0] , self.MLSPComment[0] )
            hlsp_bestlastcomment.Fill( self.M22_Best[0] , self.LastMLSPComment[0] )
            
            hch_12.Fill( self.M24_1[0] , self.M24_2[0] )
            hch_best1.Fill( self.M24_Best[0] , self.M24_1[0] )
            hch_bestcomment.Fill( self.M24_Best[0] , self.MCHComment[0] )
            hch_commentfromstau.Fill( self.MCHComment[0] , chifromstau )
            hch_commentavg.Fill( self.MCHComment[0] , (self.M24_1[0] + self.M24_2[0]) / 2.0 )

            h_chibestlspbest.Fill( self.M24_Best[0] , self.M22_Best[0] )
            h_chicommentlspcomment.Fill( self.MCHComment[0] ,  self.MLSPComment[0] )
            h_chilastcommentlsplastcomment.Fill( self.M24_2[0] , self.M22_1[0] )

            h_betterchibetterlsp.Fill( self.getBetterChi() , self.getBetterLSP() )
            h_chifromstau_betterlsp.Fill( chifromstau , self.getBetterLSP() )

        fout = TFile.Open('fout.root', 'recreate')

        hlsp_bestlastcomment.Write()
        hlsp_bestcomment.Write()
        hlsp_best1.Write()
        hlsp_12.Write()

        hch_12.Write()
        hch_commentfromstau.Write()
        hch_bestcomment.Write()
        hch_best1.Write()
        hch_commentavg.Write()


        h_chibestlspbest.Write()
        h_chicommentlspcomment.Write()
        h_chilastcommentlsplastcomment.Write()

        h_betterchibetterlsp.Write()
        h_chifromstau_betterlsp.Write()

        fout.Close()

a = allvars("~/Desktop/histo.root")
a.loop1()
a.File.Close()
