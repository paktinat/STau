#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import os
#from dbs.apis.dbsClient import DbsApi
from eostools import *
import getpass

from ROOT import TFile, TTree, gSystem, gROOT

def root_logon():
    print "1"
    gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreFWLite.so");
    print "1"
    gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreUtilities.so");
    print "1"
    gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsCommon.so");
    print "1"
    gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsFWLite.so");

    print "1"
    gSystem.Load("libPhysics");

    print "1"
    gSystem.Load("shlib/libDiLeptonAnalysis.so");



if __name__ == "__main__":
    parser = OptionParser(usage='%prog --dir=</specify/eos/dir>')
    parser.add_option("-d", "--dir", dest="directory", help="Directory", metavar="/specify/eos/directory")
    parser.add_option("-n", type="int", dest="nfiles", default=1)
    (options, args) = parser.parse_args()

    if not (options.directory):
        parser.print_help()
        parser.error('Mandatory option is --dir')
        
        
    path = "/eos/cms/store/user/" + getpass.getuser() + options.directory
    eosFiles = listFiles( path )

    outeospathbasename = os.path.basename( path ) + "_MT2Trees" 
    outeospath = os.path.normpath( path + "/../" + outeospathbasename + "/" )
    print path
    #print eosFiles

    OutFileName = options.directory.split("/")[1]

    InputFileList = {}
    if os.path.isfile( OutFileName ):
        with open(OutFileName) as f:
            for line in f:
                index_,filename,nevents_,ntotalevents_ = line.split(",")
                index = int( index_)
                nevents = int(nevents_)
                ntotalevents = int(ntotalevents_)
                InputFileList[ filename ] = [ index , nevents , ntotalevents ] 
    else:    
        f1 = open(OutFileName , 'w')
        N = 0
        NEntriesTotal = 0
        for fileName in eosFiles:
            fileName = "root://eoscms/" + fileName
            file = TFile.Open(fileName)
            if file:
                t = file.Get("Events")
                N = N + 1
                NEntriesTotal = NEntriesTotal + t.GetEntries()
                print >> f1 ,"%(N)d,%(filename)s,%(NEntries)d,%(NTotal)d" % {"N":N , "filename":fileName , "NEntries":t.GetEntries(), "NTotal":NEntriesTotal}
                InputFileList[ fileName ] = [N,t.GetEntries() , NEntriesTotal]
                file.Close()
            else:
                N = N+1
                print str(N) + " " + fileName + " is not found"



    nJobs = (len(InputFileList)/options.nfiles) + 1
        

    submitcommands = open("submitall" , 'w')
    queue = "8nh"
    infileites = InputFileList.iterkeys()
    for jobindex in range( 0 , nJobs):
        output = "MT2tree_%(jobindex)d.root" % {"jobindex":jobindex}
        command = "bsub -q %(queue)s -J job%(jobindex)d runMt2Analyzer.sh %(output)s %(outeospath)s " % {"queue":queue , "jobindex":jobindex , "output":output , "outeospath":outeospath}
        countor = 0
        for fileid in infileites:
            countor+=1
            command += fileid + " "
            if( countor == options.nfiles ):
                break

        print >> submitcommands, command
