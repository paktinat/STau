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

cp /cvmfs/cms.cern.ch/$SCRAM_ARCH/external/expat/2.0.1*/lib/libexpat.so ./libexpat.so.0

export LD_LIBRARY_PATH=${PWD}:$(scram tool tag expat LIBDIR):$LD_LIBRARY_PATH

export PATH=$(scram tool tag expat BINDIR):$PATH

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

cp $(scram tool tag expat LIBDIR)/libexpat.so ./libexpat.so.0
export LD_LIBRARY_PATH=${PWD}:$(scram tool tag expat LIBDIR):$LD_LIBRARY_PATH

otherarguments=""

if [ "$sampletype" = "data" ]; then
    otherarguments=" -j Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt "
elif [ "$sampletype" = "scan" ]; then
    otherarguments=" -s MC2012 -p Cert_8TeV_13Jul_06Aug_24Aug_11Dec_ReReco_and_Rest_PromptReco_Collisions12_PU_60bins_true.root -P MC2012PU_S10_60bins.root -b BEfficiencies_ge2jets_CSVM_HT750andMET30_or_HT450andMET200_MinDPhi4_newTTbar.root -u TauEfficiencies_LooseIsoLooseEleTightMuo_fromTauPOGplots.root -f "
elif [ "$sampletype" = "mc" ]; then
    otherarguments=" -s MC2012 -p Cert_8TeV_13Jul_06Aug_24Aug_11Dec_ReReco_and_Rest_PromptReco_Collisions12_PU_60bins_true.root -P MC2012PU_S10_60bins.root -b BEfficiencies_ge2jets_CSVM_HT750andMET30_or_HT450andMET200_MinDPhi4_newTTbar.root -u TauEfficiencies_LooseIsoLooseEleTightMuo_fromTauPOGplots.root "
fi

echo $otherarguments


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
    set -vx
    COUNTER=0
    for file in $ListOfFilesToRunOn
    do
	
	COUNTER2=0
	while [ ! -f ./IN.root ]
	do
	    if [ $COUNTER2 -gt 20 ]; then
		break
	    fi

	    if [[ ${file:0:1} == "/" ]]; then
		cp "$file" ./IN.root
	    else
		xrdcp "$file" ./IN.root
	    fi
	    let COUNTER2=COUNTER2+1
	    sleep 3
	done

	if [ ! -f ./IN.root ]; then
	    echo "The $file can't be copied after 20 tries and skipped"
	else
	    NEventsIn=`python RunManager/ucsd/getnevents.py`
	    NEventsInRun=500
	    NJobs=`expr $NEventsIn / $NEventsInRun`
	    let NJobs=NJobs+1
	    FIRSTEVENT=0
	    for n in $(seq 1 $NJobs)
	    do
		./RunMT2Analyzer -n $NEventsInRun -a $FIRSTEVENT -d . -i $processid -t $sampletype -m $cutset $otherarguments  -e -E -c -o MT2treeS_${COUNTER}_$n.root ./IN.root
		let FIRSTEVENT=FIRSTEVENT+$NEventsInRun
	    done
	    rm -rf ./IN.root
	fi
	let COUNTER=COUNTER+1
    done
    
    hadd ./MT2tree_$iteration.root MT2treeS_*.root
    
    HADOOPPATH=`echo "$outputdir" | cut -d '=' -f 2`
    COUNTER2=0
    while [ ! -f  $HADOOPPATH/MT2tree_$iteration.root ]
    do
        if [ $COUNTER2 -gt 20 ]; then
	    break
	fi
	lcg-cp -v -D srmv2 file://`pwd`/MT2tree_$iteration.root $outputdir/MT2tree_$iteration.root
	let COUNTER2=COUNTER2+1
	echo ${COUNTER2}th Try
	sleep 10
    done
    
    rm -rf ./MT2tree*.root

    if [ ! -f  $HADOOPPATH/MT2tree_$iteration.root ]; then
	echo "The file was not copied to /hadoop after 20 tries"
	exit 1
    fi
fi
