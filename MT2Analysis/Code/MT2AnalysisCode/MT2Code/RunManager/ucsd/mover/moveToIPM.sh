#! /bin/bash                                                    
 

export PREFIXLS="srm"

export SOURCEADD="://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN="
export DESTADD="://se1.particles.ipm.ac.ir:8446/srm/managerv2?SFN="

export lfn_pfnsrc="hadoop/cms"
export lfn_pfndest="dpm/particles.ipm.ac.ir/home/cms"


#OldUserNameToReplace=$2
#if [[ $OldUserNameToReplace == "" ]]; then
#    OldUserNameToReplace=`whoami`
#fi

#CurrentUserName=`whoami`

for DSName in $(cat $1) ; do
    src=$DSName
#    dest=$( sed -e "s|$OldUserNameToReplace|$CurrentUserName|g" <<< $src )
    dest=$( sed -e "s|$lfn_pfnsrc|$lfn_pfndest|g" <<< $src )
 
    #echo $OldUserNameToReplace $CurrentUserName $src $dest
    
    lcg-cp -v -D srmv2  "$PREFIXLS$SOURCEADD/$src" "$PREFIXLS$DESTADD/$dest"
done

