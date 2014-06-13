#! /bin/bash                                                    
 
export PREFIXCP="gsiftp"
export PREFIXLS="srm"

export SOURCEADD="://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN="
export DESTADD="://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN="

export lfn_pfnsrc="pnfs/lcg.cscs.ch/cms/trivcat"
export lfn_pfndest="hadoop/cms"

OldUserNameToReplace=$2
if [[ $OldUserNameToReplace == "" ]]; then
    OldUserNameToReplace=`whoami`
fi
CurrentUserName=`whoami`

for DSName in $(cat $1) ; do
    src=$DSName
    dest=$( sed -e "s|$OldUserNameToReplace|$CurrentUserName|g" <<< $src )

    #echo $OldUserNameToReplace $CurrentUserName $src $dest
    
    lcg-cp -v -D srmv2  "$PREFIXLS$SOURCEADD/$lfn_pfnsrc/$src" "$PREFIXLS$DESTADD/$lfn_pfndest/$dest"
done

