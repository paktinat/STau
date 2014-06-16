from optparse import OptionParser
from dbs.apis.dbsClient import DbsApi
     
if __name__ == "__main__":
    parser = OptionParser(usage='%prog --dataset=</specify/dataset/name> --username=<username>')
    parser.add_option("-d", "--dataset", dest="ds", help="dataset", metavar="/specify/dataset/name")
        
    (options, args) = parser.parse_args()

    url = "https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter"
    api = DbsApi(url)

    files_ = api.listFiles(dataset=options.ds)

    if len( files_ ) == 0:
        url = "https://cmsweb.cern.ch/dbs/prod/phys02/DBSReader"
        api = DbsApi( url )
        files_ = api.listFiles(dataset=options.ds)


    print len(files_)
