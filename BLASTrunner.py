import argparse
import re
import sys
import time
from xml.etree import ElementTree

import backoff
import requests

BLAST_QUERY_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

@backoff.on_predicate(backoff.fibo, lambda status: status == "WAITING", max_value=60)
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
def fetch_results(RID):
    results = {}

    blast_params = {}
    blast_params["CMD"] = "Get"
    blast_params["FORMAT_OBJECT"] = "Alignment"
    blast_params["FORMAT_TYPE"] = "XML"
    blast_params["RID"] = RID

    #response = requests.post(BLAST_QUERY_URL, params=blast_params)
    #response_text = response.text
    #root = ElementTree.fromstring(response_text)

    root = ElementTree.parse("results.xml").getroot()

    queries = []
    hits = []
    hsps = []

    for BlastOutput_iteration in root.findall("BlastOutput_iterations"):
        for iteration in BlastOutput_iteration.findall("Iteration"):
            query_data = {}
            query_data["queryID"] = iteration.find("Iteration_query-ID").text
            query_data["queryDef"] = iteration.find("Iteration_query-def").text
            query_data["queryLength"] = int(iteration.find("Iteration_query-len").text)
            queries.append(query_data)

            for hit in iteration.findall("Iteration_hits/Hit"):
                hit_data = {}
                hit_data["hitID"] = hit.find("Hit_id").text
                hit_data["hitDef"] = hit.find("Hit_def").text
                hit_data["accession"] = hit.find("Hit_accession").text
                hits.append(hit_data)

                for hsp in hit.findall("Hit_hsps/Hsp"):
                    hsp_data = {}
                    hsp_data["alignmentLength"] = int(hsp.find("Hsp_align-len").text)
                    hsp_data["bitScore"] = float(hsp.find("Hsp_bit-score").text)
                    hsp_data["eValue"] = float(hsp.find("Hsp_evalue").text)
                    hsp_data["gaps"] = int(hsp.find("Hsp_gaps").text)
                    hsp_data["percentID"] = 100*(hsp_data["alignmentLength"]/query_data["queryLength"])
                    hsps.append(hsp_data)
            
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
    if not RID:
        print("Something went wrong. Please try search again.")
        sys.exit(1)
    
    rtoe_block = re.search("RTOE = (.*)\n", response_text)
    if rtoe_block:
        RTOE = int(rtoe_block.group(1))
        print("RTOE = " + str(RTOE) + " seconds")

    if RTOE:
        print("Sleeping for {} seconds while awaiting results...".format(RTOE))
        time.sleep(RTOE)
    
    print("Checking status of web blast search {}".format(RID))
    status = _check_status(RID)
    print(status)

    if status == "FAILED":
        print("Web blast search {} failed.".format(RID))
        print("Report error at https://support.nlm.nih.gov/support/create-case/")
        sys.exit(1)
    if status == "UNKNOWN":
        print("Web blast search {} has expired; try re-running a new search.".format(RID))
        sys.exit(1)
    if status == "READY":
        print("Search complete; Retrieving results...")
        results = _fetch_results(RID)
    
    #if results:
        # insert results into SQLite db

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="Path to the fasta file to be used to query NCBI BLAST")
    parser.add_argument("-f", "--fetch_results", help="Fetch results for search with input RID")
    #parser.add_argument("-d", "--blast_database", default="nr", help="Which BLAST database to query, default is nr")
    args = parser.parse_args()
    if args.input_file:
        #print("Querying NCBI BLAST {} database using fasta file {}".format(args.blast_database, args.input_file))
        print("Performing web blastp query against nr database with fasta file {}".format(args.input_file))
        query_blast(args.input_file)
    if args.fetch_results:
        print("Fetching results for RID {}".format(args.fetch_results))
        fetch_results(args.fetch_results)
    