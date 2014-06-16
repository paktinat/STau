from optparse import OptionParser
import logging
import sys
import os
from dbs.apis.dbsClient import DbsApi

def GetFileIndex( filename ):
    fn = os.path.basename( filename )
    return int( fn.split("_")[3] )
     

if __name__ == "__main__":
    parser = OptionParser(usage='%prog --dataset=</specify/dataset/name> --username=<username>')
    parser.add_option("-d", "--dataset", dest="ds", help="dataset", metavar="/specify/dataset/name")
    parser.add_option("-u", "--username", dest="username", help="the username in ucsd se")
    parser.add_option("-n", "--nfiles", dest="nfiles", help="number of files to return" , type="int")
    parser.add_option("-i", "--iteration", dest="iteration", help="the number of job" , type="int")
    parser.add_option("-c", "--check", dest="check", help="check if the file exists in /hadoop/" , action="store_true" ,  default=False)    
        
    (options, args) = parser.parse_args()


    hadoopbasename = "/hadoop/cms/"

    url = "https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter"
    api = DbsApi(url)

    files_ = api.listFiles(dataset=options.ds)

    producerusername = ''
    sss = files_[0]['logical_file_name'].split('/')
    for i in range(0,len(sss)):
        if sss[i] == 'user':
            producerusername = sss[i+1]


    if(not options.username ):
        options.username = producerusername
        

    files = []
    for fff in files_:
        newfilename = hadoopbasename + fff['logical_file_name'].replace( producerusername , options.username , 1 )
        if options.check :
            if os.path.isfile( newfilename ):
                files.append( newfilename  )
            else:
                print newfilename + " doesn't exist"
        else:
            files.append(newfilename)

    files.sort( key=GetFileIndex )
    sss = options.ds.split("/")
    outfilename = sss[1]
    print outfilename

    printvar = "" 
    for i in range( options.iteration * options.nfiles , options.iteration * options.nfiles + options.nfiles ):
        if i < len(files):
            printvar += files[i] + " "

    print printvar
