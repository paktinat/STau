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
ListOfFilesToRunOn=`python getlistoffiles.py -u $usernametoreplace  -d $dataset  -n $nfiles -i $iteration -x $X509_USER_PROXY -c`

cd ../..
if [ "$debugmode" == "1" ]; then
    echo "./RunMT2Analyzer -d . -i $processid -t $sampletype -m $cutset  -e -E -c -o MT2tree_$iteration.root $ListOfFilesToRunOn"
    echo "mkdir $outputdir"
    echo "cp ./MT2tree_$iteration.root $outputdir"
    echo "rm -rf ./MT2tree_$iteration.root"
else
    ./RunMT2Analyzer -d . -i $processid -t $sampletype -m $cutset  -e -E -c -o MT2tree_$iteration.root $ListOfFilesToRunOn
    mkdir $outputdir
    cp ./MT2tree_$iteration.root $outputdir
    rm -rf ./MT2tree_$iteration.root
fi



##python getlistoffiles.py -d "/DoubleMu/hbakhshi-V03-09-01_RunA-EleTau-9bc3a6bb0ab58621fcde9a392f5df1a3/USER" -u hbakhshi -n 8 -i 1 -c
