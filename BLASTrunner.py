import argparse
import re
import sqlite3
import sys
import time
from xml.etree import ElementTree

import backoff
import requests


# Set global variables
BLAST_QUERY_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

CREATE_QUERIES_TABLE = (
    "CREATE TABLE IF NOT EXISTS queries "
    "(queryID TEXT PRIMARY KEY, queryDef TEXT, queryLength INTEGER)"
)
CREATE_HITS_TABLE = (
    "CREATE TABLE IF NOT EXISTS hits "
    "(hitID TEXT PRIMARY KEY, hitDef TEXT, accession TEXT, queryID TEXT, "
    "FOREIGN KEY (queryID) REFERENCES queries (queryID))"
)
CREATE_HSPS_TABLE = (
    "CREATE TABLE IF NOT EXISTS hsps "
    "(hspID INTEGER PRIMARY KEY AUTOINCREMENT, alignLength INTEGER, "
    "bitScore REAL, eValue REAL, gaps INTEGER, percentID REAL, hitID TEXT, "
    "FOREIGN KEY (hitID) REFERENCES hits (hitID))"
)

CREATE_STATEMENTS = [CREATE_QUERIES_TABLE, CREATE_HITS_TABLE, CREATE_HSPS_TABLE]

INSERTS = {
    "queries": "INSERT INTO queries(queryID, queryDef, queryLength) values (?,?,?)",
    "hits": "INSERT INTO hits(hitID, hitDef, accession, queryID) values (?,?,?,?)",
    "hsps": "INSERT INTO hsps(alignLength, bitScore, eValue, gaps, percentID, hitID) "
    "values (?,?,?,?,?,?)",
}


@backoff.on_predicate(backoff.fibo, lambda status: status == "WAITING", max_value=60)
def _check_status(RID):
    """Check the status of a query submitted to web BLAST using query's RID. Backoff/retry
    while status equals 'WAITING'.

    Parameters
        RID (str): the RID of the query for which to perform status check
    
    Returns
        status (str): query status reported by BLAST ("WAITING", "FAILED", "UNKNOWN", or "READY")
    """
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


def _fetch_results(RID):
    """Use RID from search query to retrieve results from web BLAST
    and convert response text to XML ElementTree object format.

    Parameters:
        RID (str): RID to identify search query from which to retrieve results

    Returns:
        root (obj of class xml.etree.ElementTree): ElementTree object representing full XML data
    """
    blast_params = {}
    blast_params["CMD"] = "Get"
    blast_params["FORMAT_OBJECT"] = "Alignment"
    blast_params["FORMAT_TYPE"] = "XML"
    blast_params["RID"] = RID

    # TODO uncomment real request code here and remove static response from file
    # response = requests.post(BLAST_QUERY_URL, params=blast_params)
    # response_text = response.text
    # root = ElementTree.fromstring(response_text)
    root = ElementTree.parse("results.xml").getroot()

    return root


def _parse_xml_results(root):
    """Parse the XML results fetched from BLAST into data structures for insertion into database

    Parameters:
        root (obj of class xml.etree.ElementTree): ElementTree object representing full XML data
    
    Returns:
        queries (list): a list of tuples containing information about each query performed
        hits (list): a list of tuples containing information about each hit returned from BLAST
        hsps (list): a list of tuples containing information about each hsp returned from BLAST
    """
    queries = []
    hits = []
    hsps = []

    for BlastOutput_iteration in root.findall("BlastOutput_iterations"):
        for iteration in BlastOutput_iteration.findall("Iteration"):
            # build a tuple consisting of (queryID, queryDef, queryLength)
            query_data = (
                iteration.find("Iteration_query-ID").text,
                iteration.find("Iteration_query-def").text,
                int(iteration.find("Iteration_query-len").text),
            )
            queries.append(query_data)

            for hit in iteration.findall("Iteration_hits/Hit"):
                # build a tuple consisting of (hitID, hitDef, accession, queryID)
                hit_data = (
                    hit.find("Hit_id").text,
                    hit.find("Hit_def").text,
                    hit.find("Hit_accession").text,
                    iteration.find("Iteration_query-ID").text,
                )
                hits.append(hit_data)

                for hsp in hit.findall("Hit_hsps/Hsp"):
                    # build a tuple consisting of:
                    # (alignmentLength, bitScore, eValue, gaps, percentID, hitID)
                    hsp_data = (
                        int(hsp.find("Hsp_align-len").text),
                        float(hsp.find("Hsp_bit-score").text),
                        float(hsp.find("Hsp_evalue").text),
                        int(hsp.find("Hsp_gaps").text),
                        100
                        * (
                            int(hsp.find("Hsp_align-len").text)
                            / int(iteration.find("Iteration_query-len").text)
                        ),
                        hit.find("Hit_id").text,
                    )
                    hsps.append(hsp_data)

    return queries, hits, hsps


def _initialize_database():
    """
    Create a SQLite database named blastresults.db containing queries, hits, and hsps tables

    Parameters:
        (None)
    
    Returns:
        bool: True on success, or prints an error and exits the program on Exception
    """
    try:
        conn = sqlite3.connect("blastresults.db")

        for create in CREATE_STATEMENTS:
            conn.execute(create)

        conn.commit()
        conn.close()

    except Exception:
        print("An error occurred when trying to initialize the SQLite database.")
        sys.exit(1)

    return True


def _load_results_into_database(result_data, db_table):
    """
    Load a dataset into the specified table in SQLite database

    Parameters
        result_data (list): A list of tuples, each containing a result row to insert into database
        db_table (str): name of database table into which result_data should be inserted

    Returns
        bool: True on success, or prints an error and exits the program on Exception
        
    """
    try:
        conn = sqlite3.connect("blastresults.db")

        conn.executemany(INSERTS[db_table], result_data)

        conn.commit()
        conn.close()

    except Exception:
        print("An error occurred when trying to insert data into the {} table".format(db_table))
        sys.exit(1)

    return True


# Using input fasta_file, send search request to web BLAST
# If search returns results, fetch them and insert into SQLite db
def run_blast(fasta_file):
    """TODO add docstring here"""
    # blast_params = {}
    # blast_params["CMD"] = "Put"
    # blast_params["DATABASE"] = "nr"
    # blast_params["PROGRAM"] = "blastn"

    # with open(fasta_file, "r") as seq:
    #     blast_params["QUERY"] = seq.read()
    #     print(blast_params["QUERY"])

    # response = requests.post(BLAST_QUERY_URL, params=blast_params)
    # response_text = response.text

    # rid_block = re.search("RID = (.*)\n", response_text)
    # if rid_block:
    #     RID = rid_block.group(1)
    #     print("RID = " + RID)
    # if not RID:
    #     print("Something went wrong. Please try search again.")
    #     sys.exit(1)

    # rtoe_block = re.search("RTOE = (.*)\n", response_text)
    # if rtoe_block:
    #     RTOE = int(rtoe_block.group(1))
    #     print("RTOE = " + str(RTOE) + " seconds")

    # if RTOE:
    #     print("Sleeping for {} seconds while awaiting results...".format(RTOE))
    #     time.sleep(RTOE)

    # print("Checking status of web blast search {}".format(RID))
    # status = _check_status(RID)
    # print(status)

    # if status == "FAILED":
    #     print("Web blast search {} failed.".format(RID))
    #     print("Report error at https://support.nlm.nih.gov/support/create-case/")
    #     sys.exit(1)
    # if status == "UNKNOWN":
    #     print("Web blast search {} has expired; try re-running a new search.".format(RID))
    #     sys.exit(1)
    # if status == "READY":
    #     print("Search complete; Retrieving results...")
    RID = "123ABC"
    root = _fetch_results(RID)
    queries, hits, hsps = _parse_xml_results(root)

    # Initialize the SQLite database
    if _initialize_database():
        print("initialized the database")

    # Load BLAST results into SQLite database
    if _load_results_into_database(queries, "queries"):
        print("loaded queries into database")

    if _load_results_into_database(hits, "hits"):
        print("loaded hits into database")

    if _load_results_into_database(hsps, "hsps"):
        print("loaded hsps into database")

    print("Successfully loaded BLAST results into SQLite database.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input_file", help="Path to the fasta file to be used to query NCBI BLAST"
    )

    # TODO remove extra arg used for development
    parser.add_argument("-f", "--fetch_results")

    args = parser.parse_args()
    if args.input_file:
        print(
            "Performing web blastp query against nr database with fasta file {}".format(
                args.input_file
            )
        )
        run_blast(args.input_file)

    # TODO remove extra arg used for development
    if args.fetch_results:
        print("Fetching results for RID {}".format(args.fetch_results))
        _fetch_results(args.fetch_results)
