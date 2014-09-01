#!/usr/bin/env python
# Standard python import
from optparse import OptionParser
import os

import array

# Import ROOT classes
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut, gDirectory , TH1F , TFormula, std, TPad , TCanvas, TGraph, TLegend
from ROOT import kYellow , kBlue , kOrange , kCyan , kRed , kSpring , kMagenta, kTeal , kViolet , kGreen 

SUSYNAMES = [ '00_50' , '50_100' , '100_150' , '150_200' , '200_250' , '250_300' , '300_350' , '350_400' , '400_450' , '450_500' ]
COLOR =  [ kYellow-2 , kBlue , kOrange , kCyan , kRed , kSpring , kMagenta, kTeal , kViolet , kGreen ]

if __name__ == "__main__":
    #to disable x session during the run
    os.environ['DISPLAY'] = '0'

    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="dir" )
    parser.add_option("-f", "--file", dest="file" )
    (options, args) = parser.parse_args() 

    method = 0

    file = TFile( options.file , 'read' )
    directory = file.GetDirectory(options.dir)

    contents = directory.GetListOfKeys()
    for dirid in range(0 , contents.GetSize()):
        dirname = contents[dirid].GetName()
        if dirname == 'MT2' :
            dir = directory.GetDirectory( dirname )
            print dir.GetName()
            for lower in range( 1 , 2):
                canvas_name = '%(varname)s_%(lower)d_%(method)d_%(cut)s' % {'varname':dirname , 'lower':lower , 'method':method , 'cut':options.dir}
                c1 = TCanvas( canvas_name , canvas_name , 800 , 800 )

                mmm = 0
                for susyname in SUSYNAMES:
                    name = '%(varname)s_%(lower)d_%(method)d_SUSY_%(susy)s_%(cut)s' % {'varname':dirname , 'lower':lower , 'method':method , 'susy':susyname , 'cut':options.dir}
                    histo = dir.Get( name )
                    #c1.Clear()
                    #histo.Draw()
                    #c1.Print( canvas_name + '.gif+100' )
                    
                    x = array.array( 'd' , [0.0])
                    y = array.array( 'd' , [0.0])

                    newx = []
                    newy = []
                    maximum = -1.0
                    for point in range(0 , histo.GetN() ):
                        histo.GetPoint( point , x , y )
                        if x[0] < 300 :
                            newx.append(x[0])
                            newy.append(y[0])
                            if y[0] > maximum :
                                maximum = y[0]
                        else:
                            histo.SetPoint( point , x[0] , 0.0 )
                    for point in range(0 , len(newx) ):
                        histo.SetPoint( point , newx[point] , 5.0*newy[point]/maximum )
                    option = "AL"
                    c1.cd()
                    histo.SetTitle( '%s < #Delta M < %s' % tuple( susyname.split("_") ) )
                    print histo.GetTitle()
                    histo.SetLineColor( COLOR[mmm] )
                    histo.SetLineWidth( 3 )
                    histo.Draw( option )
                    mmm += 1
                    histo.GetXaxis().SetLimits( 0 , 300 )
                    histo.GetXaxis().SetTitle( "MT2 Cut")
                    histo.GetYaxis().SetTitle( "Normalized Significance" )
                    c1.SaveAs( susyname + ".gif" )
