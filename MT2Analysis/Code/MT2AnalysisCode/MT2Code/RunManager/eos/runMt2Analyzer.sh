#! /bin/bash -vx

output=$1
shift

eosoutdir=$1
shift

input=$@

cd ~/work/TauFR/CMSSW_5_3_14_patch1/src/
eval `scramv1 runtime -sh`
cd ~/work/TauFR/STau/MT2Analysis/Code/MT2AnalysisCode/MT2Code/

eos1="/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select"

./RunMT2Analyzer -d $WORKDIR -i 0  -t data -m TauFR  -e -E -c -o $output $input
 
$eos1 mkdir -p $eosoutdir
$eos1 cp $WORKDIR/$output $eosoutdir/
rm -rf $WORKDIR/$output


