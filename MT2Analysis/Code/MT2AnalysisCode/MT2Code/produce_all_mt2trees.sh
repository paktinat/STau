#! /bin/bash

for sample in 500_300_Tauola_1514 500_300_Tauolapp_1000 800_0_Tauola_1336 800_0_Tauolapp_1257
do
  echo $sample
    ./RunMT2Analyzer -d . -i 4 -o MassTree_$sample.root -t mc -m MC -e -E /dataLOCAL/MT2Tau/RootFiles/TChipChimSTauSnu_Test/NTupleProducer_CMSSW_5_3_7_$sample.root > $sample.out
done
