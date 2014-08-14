from ROOT import TFile, TH1D, TH2D, TTree, TGraphAsymmErrors, TCanvas, gROOT, TMultiGraph, TLatex, gStyle, gDirectory, TLegend
import os
import numpy as n
import math
import array
import datetime


gROOT.ProcessLine( "struct combination_info { double limit ; float qe ; } ;" ) 


file = TFile.Open( 'higgsCombineTest.Asymptotic.mH120.root' )
limit = file.Get('limit')

limit.Scan('*')

#info = combination_info()

blimits = array.array('d', [0.0])
#quantileExpected = n.zeros( 1 , dtype=float )
quantileExpected = array.array( 'f' , [0.0])

#limit.SetBranchAddress( 'limit' , AddressOf( info , 'limit') )
limit.SetBranchAddress( 'limit' , blimits )
limit.SetBranchAddress( 'quantileExpected' , quantileExpected )
            
for i in range( 0, limit.GetEntries() ):
    limit.GetEntry( i )

    print ( '%d : %f - %f' % (i , blimits[0] , quantileExpected[0]) )

file.Close()
