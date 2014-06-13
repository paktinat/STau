#! /bin/bash -vx

#sample command
# "/DoubleMu/hbakhshi-V03-09-01_RunA-EleTau-9bc3a6bb0ab58621fcde9a392f5df1a3/USER" hbakhshi 8 1 /tmp/hbakhshi/ 0 data TauFR 

dataset=$1
usernametoreplace=$2
nfiles=$3
iteration=$4

outputdir=$5

processid=$6
sampletype=$7
cutset=$8

debugmode=1

export CMS_PATH=/code/osgcode/cmssoft/cms
export SCRAM_ARCH=slc5_amd64_gcc462
. ${CMS_PATH}/cmsset_default.sh

export GLITE_VERSION=gLite-3.2.11-1 
. /code/osgcode/ucsdt2/$GLITE_VERSION/etc/profile.d/grid-env.sh 

export LCG_GFAL_INFOSYS=lcg-bdii.cern.ch:2170 
export GLOBUS_TCP_PORT_RANGE=20000,25000 
. /code/osgcode/ucsdt2/Crab/etc/crab.sh

scramv1 project CMSSW CMSSW_5_3_14
cd CMSSW_5_3_14/src/
eval `scramv1 runtime -sh`

git clone "https://github.com/paktinat/STau" 
cd STau/
export STauWD=`pwd`
cd MT2Analysis/Code/MT2AnalysisCode/MT2Code
make

cd RunManager/ucsd/
ListOfFilesToRunOn=`python getlistoffiles.py -u $usernametoreplace  -d $dataset  -n $nfiles -i $iteration`

cd ../..
if [ "$debugmode" == "1" ]; then
    echo "./RunMT2Analyzer -d . -i $processid -t $sampletype -m $cutset  -e -E -c -o MT2tree_$iteration.root $ListOfFilesToRunOn"
    cp ./MT2tree_$iteration.root $outputdir
    rm -rf ./MT2tree_$iteration.root
else
    echo "./RunMT2Analyzer -d . -i $processid -t $sampletype -m $cutset  -e -E -c -o MT2tree_$iteration.root $ListOfFilesToRunOn"
    echo "cp ./MT2tree_$iteration.root $outputdir"
    echo "rm -rf ./MT2tree_$iteration.root"
fi



##python getlistoffiles.py -d "/DoubleMu/hbakhshi-V03-09-01_RunA-EleTau-9bc3a6bb0ab58621fcde9a392f5df1a3/USER" -u hbakhshi -n 8 -i 1 -c