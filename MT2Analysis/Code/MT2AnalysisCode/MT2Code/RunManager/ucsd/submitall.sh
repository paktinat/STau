#! /bin/bash


export GLITE_VERSION=gLite-3.2.11-1 
. /code/osgcode/ucsdt2/$GLITE_VERSION/etc/profile.d/grid-env.sh 

export LCG_GFAL_INFOSYS=lcg-bdii.cern.ch:2170 
export GLOBUS_TCP_PORT_RANGE=20000,25000 
. /code/osgcode/ucsdt2/Crab/etc/crab.sh


username=`whoami`
vomsfile=`voms-proxy-info | grep path | awk -F ":" '{print $2}' |  tr -d ' '`
vomslefthours=$(voms-proxy-info | grep timeleft | awk -F ":" '{print $2}')

#echo "the current proxy file $vomsfile is valid for $vomslefthours hours, user:$username"

if [[ ! -f $vomsfile ]]; then
    echo "You first neet to get a proxy"
    voms-proxy-init -voms cms
elif (( `expr $vomslefthours \< 3` )); then 
    echo "You need to first get a proxy"
    voms-proxy-init -voms cms
else
    echo "the current proxy is valid and will be used"
fi

vomsfile=`voms-proxy-info | grep path | awk -F ":" '{print $2}'`
vomslefthours=$(voms-proxy-info | grep timeleft | awk -F ":" '{print $2}')

echo "the current proxy file $vomsfile is valid for $vomslefthours hours, user:$username"


echo looping over $1
while read p; do
    p=${p##*( )}
   if [[ ${p:0:1} == '#' ]]; then
       echo -e "skipped commented line \n \n"
   elif [  ${#p} -eq 0 ]; then
       echo -e "skipped empty line \n \n"
   else
       splitted=($(echo $p | tr " " "\n"))
       jobname=${splitted[0]}
       jobname=($(echo $jobname | tr -d ' '))
       len=${#jobname}
       len=$(expr $len + 1)
       arguments=${p:$len}
       
       datasetname=${splitted[1]}
       nfilesinds=$(python getnfiles.py -d $datasetname)
       
       nfilesperjob=${splitted[3]}

       njobs=$(expr $nfilesinds \/ $nfilesperjob)
       njobs=$(expr $njobs + 1)

       usernametoreplace=${splitted[2]}
       commitid=${splitted[8]}
       outdir=$(python createoutdir.py -d $datasetname -a $commitid)

     
       arguments="$arguments $outdir"

       echo "creating jobs for $jobname"
       echo -e "\t DS : $datasetname"
       echo -e "\t NFilesPerJob : $nfilesperjob"
       echo -e "\t NFilesInDS : $nfilesinds"
       echo -e "\t NJobs : $njobs"
       echo -e "\t OutDir : $outdir"
       
       sed -e "s|USERNAME|$username|g" -e "s|PROXYFILENAME|$vomsfile|g" -e "s|JOBNAME|$jobname|g" -e "s|ARGS|$arguments|g" condorTemplate -e "s|NJOBS|$njobs|g"  > condor_$jobname

       echo -e "\t condor_$jobname file is created"

       if [[ ! -d $jobname ]]; then
	   mkdir $jobname
	   echo -e "\t $jobname directory is created"
       else
	   echo -e "\t the current $jobname directory will be used"
       fi

       echo -e "\t Producing the list of files"
       python getlistoffiles.py -u $usernametoreplace -d $datasetname -c -x -j $jobname

       if [[ $2 = "submit" ]]; then
	   condor_submit condor_$jobname
	   echo -e "\t Submitted"
       fi
   fi
done < $1
