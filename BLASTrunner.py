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


def _submit_query(fasta_file):
    """Build query using input fasta file and submit to web BLAST to run blastn against nr database

    Parameters
        fasta_file (str): fasta file to use for querying web BLAST

    Returns
        response_text (str): response text from the query request
    """
    blast_params = {}
    blast_params["CMD"] = "Put"
    blast_params["DATABASE"] = "nr"
    blast_params["PROGRAM"] = "blastn"

    with open(fasta_file, "r") as seq:
        blast_params["QUERY"] = seq.read()

    response = requests.post(BLAST_QUERY_URL, params=blast_params)
    response_text = response.text
    print("Query submitted to web BLAST:")

    return response_text


def _parse_RID_RTOE(response_text):
    """Parse out the Request ID (RID) and RTOE from the response received from web BLAST upon query submission

    Parameters:
        response_text (str): response text from the query request

    Returns:
        RID (str): RID of the query submitted to web BLAST
        RTOE (int): estimated time until search is completed, in seconds
    """
    rid_block = re.search("RID = (.*)\n", response_text)
    if rid_block:
        RID = rid_block.group(1)
        print("RID = " + RID)
    else:
        RID = ""

    rtoe_block = re.search("RTOE = (.*)\n", response_text)
    if rtoe_block:
        RTOE = int(rtoe_block.group(1))
        print("RTOE = " + str(RTOE) + " seconds")
    else:
        RTOE = 0
    return RID, RTOE


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
        print(status + "...")

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

    response = requests.post(BLAST_QUERY_URL, params=blast_params)
    response_text = response.text
    root = ElementTree.fromstring(response_text)
    # root = ElementTree.parse("results.xml").getroot()

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


def _initialize_database(db_name):
    """
    Create a SQLite database containing queries, hits, and hsps tables

    Parameters:
        db_name (str): Name of output SQLite database
    
    Returns:
        bool: True on success, or prints an error and exits the program on Exception
    """
    try:
        conn = sqlite3.connect(db_name)

        for create in CREATE_STATEMENTS:
            conn.execute(create)

        conn.commit()
        conn.close()

    except Exception:
        print("An error occurred when trying to initialize the SQLite database.")
        sys.exit(1)

    return True


def _load_results_into_database(db_name, result_data, db_table):
    """
    Load a dataset into the specified table in SQLite database

    Parameters
        db_name (str): Name of output SQLite database
        result_data (list): A list of tuples, each containing a result row to insert into database
        db_table (str): name of database table into which result_data should be inserted

    Returns
        bool: True on success, or prints an error and exits the program on Exception
        
    """
    try:
        conn = sqlite3.connect(db_name)

        conn.executemany(INSERTS[db_table], result_data)

        conn.commit()
        conn.close()

    except Exception:
        print("An error occurred when trying to insert data into the {} table".format(db_table))
        sys.exit(1)

    return True


# Using input fasta_file, send search request to web BLAST
# If search returns results, fetch them and insert into SQLite db
def run_blast(fasta_file, output_db_name):
    """Procedure for BLASTrunner
        - queries web BLAST with input fasta file
        - fetches results in XML format when ready
        - parses XML results
        - initializes SQLite database
        - inserts results (queries, hits, and hsps) into SQLite database

    Parameters
        fasta_file (str): fasta file to use for querying web BLAST
        output_db_name (str): name for local results database
    
    Returns
        None
    """
    response_text = _submit_query(fasta_file)

    RID, RTOE = _parse_RID_RTOE(response_text)
    if not RID:
        print("Something went wrong. Please try search again.")
        sys.exit(1)
    if RTOE:
        print("Sleeping for {} seconds while awaiting results...".format(RTOE))
        time.sleep(RTOE)

    print("Checking status of web BLAST search: RID {}".format(RID))
    status = _check_status(RID)

    if status == "FAILED":
        print("Web BLAST search {} failed.".format(RID))
        print("Report error at https://support.nlm.nih.gov/support/create-case/")
        sys.exit(1)
    if status == "UNKNOWN":
        print("Web BLAST search {} has expired; try re-running a new search.".format(RID))
        sys.exit(1)
    if status == "READY":
        print("Retrieving results...")

    root = _fetch_results(RID)
    queries, hits, hsps = _parse_xml_results(root)

    if _initialize_database(output_db_name):
        print("Initialized SQLite database")

    if _load_results_into_database(output_db_name, queries, "queries"):
        print("Loaded query data into database")

    if _load_results_into_database(output_db_name, hits, "hits"):
        print("Loaded hit data into database")

    if _load_results_into_database(output_db_name, hsps, "hsps"):
        print("Loaded hsp data into database")

    print("Successfully loaded BLAST results into SQLite database!")
    print("See README for help with querying local results database")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Path to the fasta file to be used to query NCBI BLAST")
    parser.add_argument(
        "-o", "--output_db_name", default="blastresults.db", help="name for local results database"
    )

    args = parser.parse_args()
    if args.input_file:
        print(
            "Performing web BLAST blastn query against nr database with fasta file {}".format(
                args.input_file
            )
        )
        run_blast(args.input_file, args.output_db_name)
