#!/usr/bin/env python
# Standard python import
from optparse import OptionParser
import os

# Import ROOT classes
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut, gDirectory , TH1F , TFormula, std, TPad , TCanvas

SUSYNAMES = [ '00_50' , '50_100' , '100_150' , '150_200' , '200_250' , '250_300' , '300_350' , '350_400' , '400_450' , '450_500' ]

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
        if not dirname == 'EventLists' :
            dir = directory.GetDirectory( dirname )
            print dir.GetName()
            for lower in range( 0 , 2):
                canvas_name = '%(varname)s_%(lower)d_%(method)d_%(cut)s' % {'varname':dirname , 'lower':lower , 'method':method , 'cut':options.dir}
                c1 = TCanvas( canvas_name , canvas_name , 800 , 800 )
                for susyname in SUSYNAMES:
                    name = '%(varname)s_%(lower)d_%(method)d_SUSY_%(susy)s_%(cut)s' % {'varname':dirname , 'lower':lower , 'method':method , 'susy':susyname , 'cut':options.dir}
                    histo = dir.Get( name )
                    c1.Clear()
                    histo.Draw()
                    c1.Print( canvas_name + '.gif+100' )
