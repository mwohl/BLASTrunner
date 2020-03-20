# BLASTrunner
**Perform a web BLAST query and store results in a SQLite database**

BLASTrunner provides a simple mechanism for running a blastn query against the BLAST nr database using the NCBI-BLAST Common URL API.  Using a fasta file as input, BLASTrunner composes the query, makes the API call to submit it, checks status, and retrieves the XML result.  It then parses the XML and inserts the results into an SQLite database.

## Setup/Installation

After cloning the repository, setup consists of a few simple steps:

1. Create a python3 virtual environment in which to run the program

    ` python3 -m venv .venv `

2. Activate the virtual environment

    ` source .venv/bin/activate `

3. Install BLASTrunner's dependencies using pip

    ` pip install -r requirements.txt `

## Usage

BLASTrunner requires a fasta file, which can contain one or more DNA sequences, as input.  For example:

    python BLASTrunner.py /path/to/myseq.fasta

It optionally accepts a name to use for the SQLite output database.  For example:

    python BLASTrunner.py /path/to/myseq.fasta -o nrblast20200320.db

If no database name is provided, the SQLite database will be named **blastresults.db** by default.

## Output

The output from BLASTrunner is a SQLite database consisting of three tables:
| queries | info about queries submitted to web BLAST |
| ----------- | ----------- |
| queryID | TEXT |
| queryDef | TEXT |
| queryLength | INTEGER|

| hits | hits returned from query |
| ----------- | ----------- |
| hitID | TEXT |
| hitDef | TEXT |
| accession | TEXT |
| queryID | TEXT |

| hsps | hsps found per hit |
| ----------- | ----------- |
| hspID | INTEGER |
| alignLength | INTEGER |
| bitScore | REAL |
| eValue | REAL |
| gaps | INTEGER |
| percentID | REAL |
| hitID | TEXT |

## Querying Results Database

To access the SQLite results database via command line, invoke sqlite3 and provide the name of the database.  For example:

    sqlite3 blastresults.db

#### Example Queries
To see BLAST queries and get their IDs:

    SELECT * FROM queries;

To see hits returned for a particular BLAST query of interest:

    SELECT * FROM hits WHERE queryID = "Query_14919";

To see hsps associated with a particular hit:

    SELECT * FROM hsps WHERE hitID = "gi|1772680595|gb|CP045560.1|";

To see all hsps with a percent ID of 100.0:

    SELECT * FROM hsps WHERE percentID = 100.0;

To narrow that down to hsps with percent ID of 100.0 that came from a specific query:

    SELECT * FROM hsps p JOIN hits h on p.hitID = h.hitID JOIN queries q on h.queryID = q.queryID WHERE q.queryID = "Query_14919" AND p.percentID = 100.0;

Etc.

Happy querying!
