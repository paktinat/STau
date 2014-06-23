from optparse import OptionParser
import logging
import sys
import os
from dbs.apis.dbsClient import DbsApi

import getpass

if __name__ == "__main__":
    parser = OptionParser(usage='%prog --dataset=</specify/dataset/name> --username=<username> --appendix=<appendix>')
    parser.add_option("-d", "--dataset", dest="ds", help="dataset", metavar="/specify/dataset/name")
    parser.add_option("-a", "--appendix", dest="appendix" )
    (options, args) = parser.parse_args()    

    hadoopbasename = "/hadoop/cms/"

    url = "https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter"
    api = DbsApi(url)

    files_ = api.listFiles(dataset=options.ds)

    if len( files_ ) == 0:
        url = "https://cmsweb.cern.ch/dbs/prod/phys02/DBSReader"
        api = DbsApi( url )
        files_ = api.listFiles(dataset=options.ds)

    if len( files_ ) == 0:
        print "Dataset wasn't found"

    producerusername = ''
    sss = files_[0]['logical_file_name'].split('/')
    for i in range(0,len(sss)):
        if sss[i] == 'user':
            producerusername = sss[i+1]


    username = getpass.getuser()

    ntupledir = os.path.dirname( files_[0]['logical_file_name'] )
    ntupledir = hadoopbasename + ntupledir.replace( producerusername , username , 1 )

    ntupledir += "_MT2tree_" + options.appendix
    while os.path.isdir( ntupledir ) and len( os.listdir( ntupledir ) ) > 0 :
        ntupledir += "_"
        

    if not os.path.isdir( ntupledir ):
        os.mkdir( ntupledir )

    ntupledir = "srm://bsrm-1.t2.ucsd.edu:8443/srm/v2/server?SFN=" + ntupledir
    print ntupledir
