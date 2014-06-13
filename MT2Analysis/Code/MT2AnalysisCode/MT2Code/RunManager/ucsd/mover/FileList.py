#!/usr/bin/env python
"""
_setDatasetStatus_

Give the dataset path, the new status and DBS Instance url (writer), it will
set the new status. If the children option is used, the status of all its
children will be changed as well

"""
from optparse import OptionParser
import logging
import sys
import os
from dbs.apis.dbsClient import DbsApi

if __name__ == "__main__":
    parser = OptionParser(usage='%prog --query=</specify/dataset/query>')
    parser.add_option("-q", "--query", dest="query", help="Query", metavar="/specify/dataset/query")
    (options, args) = parser.parse_args()

    if not (options.query):
        parser.print_help()
        parser.error('Mandatory option is --query')

    url = "https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter"

    api = DbsApi(url)

    datasets = api.listDatasets(dataset=options.query)

    for dataset in datasets:
        files_ = api.listFiles(dataset=dataset['dataset'])
        files = []
        for fff in files_:
            files.append( fff['logical_file_name'] )
        sss = dataset['dataset'].split("/")

        oldfiles = []
        outfile = None
        totallen = len(files)

        outfilename = sss[1]

        while os.path.isfile(outfilename):
            with open(outfilename) as f:
                oldfiles += f.read().splitlines()
            outfilename+="_"

        files = list( set(files) - set(oldfiles) )
        print "%(ds)s : %(len)d (%(totallen)d)" % {"ds":outfilename, "len":len(files) , "totallen":totallen}

        if(len(files)==0):
            continue
        outfile = open( outfilename , 'w')
        for f in files:
            print >>outfile, f

