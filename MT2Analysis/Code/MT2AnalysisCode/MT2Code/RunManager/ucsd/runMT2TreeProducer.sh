#! /bin/bash 
#-vx

#sample command
# "/DoubleMu/hbakhshi-V03-09-01_RunA-EleTau-9bc3a6bb0ab58621fcde9a392f5df1a3/USER" hbakhshi 8 1 0 data TauFR 22def0c673eac3687b5be93a0b60dbee775adb8d /tmp/hbakhshi/

dataset=$1
usernametoreplace=$2
nfiles=$3
iteration=$4

processid=$5
sampletype=$6
cutset=$7

commitid=$8

outputdir=$9

debugmode=0

WORKINGDIR=`pwd`

export CMS_PATH=/code/osgcode/cmssoft/cms
VERSION_SLC=$(lsb_release -r -s)
if (( `expr $VERSION_SLC \< 6.0` )); then
    export SCRAM_ARCH=slc5_amd64_gcc462
else
    export SCRAM_ARCH=slc6_amd64_gcc472
fi
echo $SCRAM_ARCH
. ${CMS_PATH}/cmsset_default.sh

export GLITE_VERSION=gLite-3.2.11-1 
. /code/osgcode/ucsdt2/$GLITE_VERSION/etc/profile.d/grid-env.sh 

export LCG_GFAL_INFOSYS=lcg-bdii.cern.ch:2170 
export GLOBUS_TCP_PORT_RANGE=20000,25000 
. /code/osgcode/ucsdt2/Crab/etc/crab.sh

scramv1 project CMSSW CMSSW_5_3_14
cd CMSSW_5_3_14/src/
eval `scramv1 runtime -sh`

vomsfile=`voms-proxy-info | grep path | awk -F ":" '{print $2}' |  tr -d ' '`
#this variable is needed by dbsApi for authentication
export X509_USER_PROXY=${vomsfile}
echo $X509_USER_PROXY

git clone "https://github.com/paktinat/STau" 
cd STau/
export STauWD=`pwd`
git checkout $commitid
cd MT2Analysis/Code/MT2AnalysisCode/MT2Code
make

cd RunManager/ucsd/
ListOfFilesToRunOn=`python getlistoffiles.py -n $nfiles -i $iteration -p $WORKINGDIR/Files`

cd ../..


otherarguments=""

if [ "$sampletype" = "data" ]; then
    otherarguments=" -j Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt "
elif [ "$sampletype" = "scan" ]; then
    otherarguments=" -s MC2012 -p Cert_8TeV_13Jul_06Aug_24Aug_11Dec_ReReco_and_Rest_PromptReco_Collisions12_PU_60bins_true.root -P MC2012PU_S10_60bins.root -b BEfficiencies_ge2jets_CSVM_HT750andMET30_or_HT450andMET200_MinDPhi4_newTTbar.root -u TauEfficiencies_LooseIsoLooseEleTightMuo_fromTauPOGplots.root -f "
elif [ "$sampletype" = "mc" ]; then
    otherarguments=" -s MC2012 -p Cert_8TeV_13Jul_06Aug_24Aug_11Dec_ReReco_and_Rest_PromptReco_Collisions12_PU_60bins_true.root -P MC2012PU_S10_60bins.root -b BEfficiencies_ge2jets_CSVM_HT750andMET30_or_HT450andMET200_MinDPhi4_newTTbar.root -u TauEfficiencies_LooseIsoLooseEleTightMuo_fromTauPOGplots.root "
fi

echo $otherarguments


set -vx

if [ "$debugmode" == "1" ]; then

    COUNTER=0
    for file in $ListOfFilesToRunOn
      do
      echo ./RunMT2Analyzer -d . -i $processid -t $sampletype -m $cutset $otherarguments  -e -E -c -o MT2tree_$COUNTER.root "$file"
      let COUNTER=COUNTER+1
    done
    
    echo hadd ./MT2tree_$iteration.root MT2tree_*.root
  
    echo lcg-cp -v -D srmv2 file://`pwd`/MT2tree_$iteration.root $outputdir/MT2tree_$iteration.root
    echo rm -rf ./MT2tree_$iteration.root
else
    COUNTER=0
    for file in $ListOfFilesToRunOn
      do
      while [ ! -f ./IN.root ]
	do
	xrdcp "$file" ./IN.root
      done
      ./RunMT2Analyzer -d . -i $processid -t $sampletype -m $cutset $otherarguments  -e -E -c -o MT2treeS_$COUNTER.root ./IN.root
      rm -rf ./IN.root
      let COUNTER=COUNTER+1
    done
    
    hadd ./MT2tree_$iteration.root MT2treeS_*.root
    
    lcg-cp -v -D srmv2 file://`pwd`/MT2tree_$iteration.root $outputdir/MT2tree_$iteration.root
    rm -rf ./MT2tree*.root
fi
