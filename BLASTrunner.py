import argparse
import requests
import re

BLAST_QUERY_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

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
    RID_BLOCK = re.search("RID = (.*)\n", response_text)
    RID = RID_BLOCK.group(1)
    print("RID = " + RID)
    RTOE_BLOCK = re.search("RTOE = (.*)\n", response_text)
    RTOE = RTOE_BLOCK.group(1)
    print("RTOE = " + RTOE + " seconds") 
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="Path to the fasta file to be used to query NCBI BLAST")
    parser.add_argument("-d", "--blast_database", default="nr", help="Which BLAST database to query, default is nr")
    args = parser.parse_args()
    if args.input_file:
        #print("Querying NCBI BLAST {} database using fasta file {}".format(args.blast_database, args.input_file))
        print("Performing web blastp query against nr database with fasta file {}".format(args.input_file))
        query_blast(args.input_file)