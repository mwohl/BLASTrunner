import argparse
import re
import time

import backoff
import requests

BLAST_QUERY_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

@backoff.on_predicate(backoff.fibo, lambda status: status == "WAITING", max_value=12)
def _check_status(RID):
    blast_params = {}
    blast_params["CMD"] = "Get"
    blast_params["FORMAT_OBJECT"] = "SearchInfo"
    blast_params["RID"] = RID

    response = requests.post(BLAST_QUERY_URL, params=blast_params)
    response_text = response.text

    status_block = re.search("Status=(.*)\n", response_text)
    if status_block:
        status = status_block.group(1)
        print("Status = " + status)
    
    return status

#Using RID from search query, send request to retrieve results
def _fetch_results(RID):
    results = {}

    blast_params = {}
    blast_params["CMD"] = "Get"
    blast_params["FORMAT_OBJECT"] = "Alignment"
    blast_params["RID"] = RID

    return results

#Using input fasta_file, send search request to web BLAST
#If search returns results, fetch them and insert into SQLite db
def query_blast(fasta_file):
    blast_params = {}
    blast_params["CMD"] = "Put"
    blast_params["DATABASE"] = "nr"
    blast_params["PROGRAM"] = "blastn"


    with open(fasta_file, "r") as seq:
        blast_params["QUERY"] = seq.read()
        print(blast_params["QUERY"])

    response = requests.post(BLAST_QUERY_URL, params=blast_params)
    response_text = response.text
    
    rid_block = re.search("RID = (.*)\n", response_text)
    if rid_block:
        RID = rid_block.group(1)
        print("RID = " + RID)
    
    rtoe_block = re.search("RTOE = (.*)\n", response_text)
    if rtoe_block:
        RTOE = int(rtoe_block.group(1))
        print("RTOE = " + str(RTOE) + " seconds")

    if RTOE:
        print("Sleeping for {} seconds while awaiting results".format(RTOE))
        time.sleep(RTOE)
    
    if RID:
        status = _check_status(RID)
        print(status)
        import pdb;pdb.set_trace()
        # results = _fetch_results(RID)
    
    #if results:
        # insert results into SQLite db

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="Path to the fasta file to be used to query NCBI BLAST")
    #parser.add_argument("-d", "--blast_database", default="nr", help="Which BLAST database to query, default is nr")
    args = parser.parse_args()
    if args.input_file:
        #print("Querying NCBI BLAST {} database using fasta file {}".format(args.blast_database, args.input_file))
        print("Performing web blastp query against nr database with fasta file {}".format(args.input_file))
        query_blast(args.input_file)