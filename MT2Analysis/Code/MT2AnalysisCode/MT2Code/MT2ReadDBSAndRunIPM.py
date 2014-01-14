#!/usr/bin/python

### WATCH OUT: NEEDS TO BE LAUNCHED AS python runManager.py configuration.cfg and NOT ./runManager configuration.cfg

# System libraries to interact with the shell
import sys
import os
from os import popen
import commands
import time

from ROOT import TFile, TTree, gSystem, gROOT

def root_logon():
  gSystem.SetIncludePath(" -I./include/ -I../TESCO/include/");
  #gSystem->Load("/public/V_CMSSW64/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/slc5_amd64_gcc462/libFWCoreFWLite.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreFWLite.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreUtilities.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsCommon.so");
  gSystem.Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsFWLite.so");

  gSystem.Load("libPhysics");
  gSystem.Load("shlib/libDiLeptonAnalysis.so");
  gROOT.ProcessLine(".x SetStyle_PRD.C");

def showMessage(message):
  print "[runManager] " + message


def process(dataset):
  dpmPath = "rfio:///dpm/particles.ipm.ac.ir/home/cms/"
  dbsUrl = 'http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet'
  list = [] # Initialize file list

  # Check if cache file exists and if we should use it
  cacheFile = '.'+dataset.replace('/','_').replace('?','_')
  allmyfiles=[]
  if os.path.isfile(cacheFile):
      print "Reading from cached file"
      f = open(cacheFile)
      showMessage('Reading files pertaining to '+dataset+'from cache file '+cacheFile)
      allmyfiles = f.readlines()
      f.close()
  else:
      # If not: rebuild the list
      showMessage("Going to fetch list of files pertaining to "+str(dataset))
      command = 'dbs search --query="find file where dataset='+dataset+'" --noheader --url='+dbsUrl
      print command
      theList = os.popen(command)
      list = theList.readlines()
      for li in list:
              allmyfiles.append(li)
      f = open(cacheFile,'w')
      f.write(''.join(allmyfiles))
  
  numberOfFiles = len(allmyfiles)
  if(numberOfFiles == 0):
    showMessage("No files found in "+str(task[0]))
    return "Error"

  correctList = [];
  for fileName in allmyfiles:
    auxiliar = dpmPath + fileName
    auxiliar = auxiliar[0:len(auxiliar)-1]
    correctList.append(auxiliar)



  root_logon()

  
  N = 0
  NEntriesTotal = 0
  for fileName in correctList:
    file = TFile.Open(fileName)
    t = file.Get("Events")
    N = N + 1
    NEntriesTotal = NEntriesTotal + t.GetEntries()
    print str(N) + ": " + fileName + ", NEntries:" + str(t.GetEntries()) + ", NTotal = " + str(NEntriesTotal)
    file.Close()

#process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012C-24Aug2012-v1-1c5ff4400f073c7d50e40dc88ea10d23/USER")    
process("/MultiJet/paktinat-V03-09-01_MultiJet-Run2012B-13Jul2012-v1-61aef9e3bbf721ce333708258d830dbb/USER")    
