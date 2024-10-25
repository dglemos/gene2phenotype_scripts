import os
import re
import sys
import datetime
import argparse
import MySQLdb

attrib_mapping = [
    "gain_of_function_mp",
    "loss_of_function_mp",
    "dominant_negative_mp",
]

KEY = "Badonyi_score"

def get_details_from_file(file):
    """
        Extracts lines from a file, skipping the header.
        It stores the lines after the header in the variable `skip_header`.

        Args:
            file : The path to the file to be read.
    """
    with open(file, "r") as opened_file:
        lines = opened_file.readlines()
        skip_header = lines[1:]
    
    return skip_header

def get_locus_id_from_g2p_db(list_lines, host, port, db, password, user):
    """
        Retrieves locus IDs from a gene-to-phenotype (G2P) database and appends them to the input list.

        Args:
            list_lines : A list of strings or lists, where each element contains a gene symbol as the first item.
            host       : The database hostname.
            port       : The port number to connect to the database.
            db         : The name of the database.
            password   : The password.
            user       : The username.

        Returns:
            The modified `list_lines` with locus IDs appended to each line.
            If no locus ID is found for a gene symbol, `None` is appended.

        Example:
            list_lines = ["BRCA1 info", "TP53 info"]
            updated_lines = get_locus_id_from_g2p_db(list_lines, "localhost", 3306, "g2p_db", "password", "user")
    """

    get_locus_source = """ SELECT id from locus where name = %s """
    get_locus_attrib_source = """SELECT locus_id from locus_attrib where value = %s """

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    for i, line in enumerate(list_lines):
        if isinstance(line, str):
            line = line.split()
            gene_symbol = line[0]
            cursor.execute(get_locus_source, (gene_symbol,))
            locus_id = cursor.fetchone()
    

            if locus_id:
                line.append(locus_id[0])
            elif not locus_id:
                cursor.execute(get_locus_attrib_source, (gene_symbol,))
                locus_attrib_id = cursor.fetchone()
                line.append(locus_attrib_id[0] if locus_attrib_id else None)
            else:
                line.append(None)
            
            list_lines[i] = line
        
    cursor.close()
    database.close()

    return list_lines

def get_attrib_id(host, port, db, password, user, attrib):
    """
        Retrieve the ID corresponding to a specific attribute value.

        Args:
            host     : The hostname or IP address of the MySQL server.
            port     : The port number on which the MySQL server.
            db       : The name of the database to connect to.
            password : The password
            user     : The username
            attrib   : The attribute value to search for in the 'attrib' table.

        Returns:
            The ID associated with the given attribute value from the 'attrib' table.

        Raises:
            ValueError: Source ID not found in the database
    """
    get_attrib_id_query = """SELECT id from attrib where value = %s"""

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    cursor.execute(get_attrib_id_query, (attrib,))

    attrib_id = cursor.fetchone()

    if attrib_id is None:
        raise ValueError("Source ID not found in the database.")

    return attrib_id[0]


def insert_details_into_meta(host, port, db, password, user, attrib):
    """
        Inserts details into the 'meta' table of the given database.

        Args:
            host     : the hostname or IP address of the database server
            port     : the port number used to connect to the database server
            db       : the name of the database
            password : the password
            user     : the username
    """

    source_id = get_source_details(host, port, db, password, user)
    description = "Baydoni & Marsh probabilities"
    
    insert_into_meta_query = """ INSERT into meta(`key`, date_update, description, version, source_id, is_public)
                                 VALUES (%s, %s, %s, %s, %s, %s) """

    #current_datetime = datetime.now()
    meta_key = KEY + "_" + attrib

    #formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S.%f')
    version = 1
    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()
    
    cursor.execute(insert_into_meta_query, (meta_key, datetime.datetime.now(), description, version, source_id, 0))
    database.commit()
    cursor.close()
    database.close()


def get_source_details(host, port, db, password, user):
    """
        Retrieves the source ID for 'Marsh Mechanism probabilities' from the 'source' table.

        Args:
            host     : the hostname or IP address of the MySQL database server.
            port     : the port number used to connect to the MySQL database server.
            db       : the name of the database
            password : sthe password
            user     : the username

        Returns:
            source_id : the ID of the source with the name 'Marsh Mechanism probabilities'.

        Raises:
            ValueError : Source ID not found in the database
    """

    get_source_query = """ SELECT id from source where name = 'Marsh Mechanism probabilities'"""

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    cursor.execute(get_source_query)
    source_id = cursor.fetchone()
    
    if source_id is None:
        raise ValueError("Source ID not found in the database.")
 
    source_id = source_id[0]

    cursor.close()
    database.close()

    return source_id

def insert_into_gene_stats(list_lines, host, port, db, password, user, attrib):

    """
        Inserts gene-related data into the `gene_stats` table.

        Parameters:
        ----------
        list_lines : a list containing gene-related data, where each sublist represents a row of values with the following elements:
            [0]: gene_symbol (str)
            [1]: Not used in the query
            [2]: score (float)
            [3]: Not used in the query
            [4]: gene_id (int)
            [5]: locus_identifier_id (int or None)

        host     : the hostname or IP address of the MySQL server.
        port     : the port number on which the MySQL server.
        db       : the name of the database to connect to.
        password : the password
        user     : the username
        attrib   : the attribute value

        Example usage:
        --------------
        list_lines = [
            ['gene1', 'x', 0.85, 'x', 101, 201],
            ['gene2', 'y', 0.67, 'y', 102, None],
            ['gene3', 'z', 0.92, 'z', 103, 202]
        ]

        insert_into_gene_stats(list_lines, 'localhost', 3306, 'genedb', 'password123', 'user123', 'attrib_key')
    """

    insert_into_gene_stats_query = """ INSERT into gene_stats (gene_symbol, gene_id, score, source_id, description_attrib_id)
                                       VALUES (%s, %s, %s, %s, %s)"""

    source_id = get_source_details(host, port, db, password, user)

    attrib_value = get_attrib_id(host,port,db,password,user,attrib)

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    for line in list_lines:
        if line[4] is not None:
            cursor.execute(insert_into_gene_stats_query, (line[0], line[4], line[2], source_id, attrib_value ))

    database.commit()
    cursor.close()
    database.close()


def main():
    """
        Main function to import gene probability data into the gene_stats table of a G2P database.

        Command-line Arguments:
        -----------------------
        --host : str (required)
            The hostname or IP address of the database server where the data will be imported.

        -p, --port : int (required)
            The port number on which the database server is listening.

        -d, --database : str (required)
            The name of the G2P database into which the data will be imported.

        -pwd, --password : str (required)
            The password for the database user.

        -u, --user : str (required)
            The username for connecting to the G2P database.

        -f, --file : str (required)
            Path to the file containing the gene probability data.

        -a, --attrib : str (required)
            The attribute type in the file being inserted. Allowed types include:
            'gain_of_function_mp', 'loss_of_function_mp', and 'dominant_negative_mp'.
            This argument will be used to map the `attrib` value to the corresponding description ID.

        Example usage:
        --------------
        python script.py --host localhost --port 3306 --database g2p_db --password secret_pwd --user db_user --file data.csv --attrib gain_of_function_mp
    """

    parser = argparse.ArgumentParser(description="This script is used to import the probabilities from a file and imports it to the gene_stats table in the G2P DB")
    parser.add_argument("--host", required=True, help="Host of the database were the data is to be imported")
    parser.add_argument("-p", "--port", required=True, help="Port information of the database to be imported")
    parser.add_argument("-d", "--database", required=True, help="G2P Database to import the information into")
    parser.add_argument("-pwd", "--password", required=True, help="Password for the G2P database information")
    parser.add_argument("-u", "--user", required=True, help="User information for the G2P database")
    parser.add_argument("-f", "--file", required=True, help="File containing the information for the score")
    parser.add_argument('-a', "--attrib", required=True, help="The attrib type in the file you are trying to insert. Allowed types are gain_of_function_mp, loss_of_function_mp, dominant_negative_mp ")

    args = parser.parse_args()

    host = args.host
    port = args.port
    db = args.database
    pwd = args.password
    user = args.user
    file = args.file
    port = int(port)
    attrib = args.attrib

    if attrib not in attrib_mapping:
        print("The type applied is not permitted")
        sys.exit()

    if file:
        print("Getting details from file")
        file_lines = get_details_from_file(file)
        print("Getting locus id from the G2P DB")
        get_locus_id_from_g2p_db(file_lines, host, port, db, pwd, user)
        print("Inserting into gene stats")
        insert_into_gene_stats(file_lines, host, port, db, pwd, user, attrib)
        print("Inserting into meta table")
        insert_details_into_meta(host, port, db, pwd, user, attrib)
        print("File has been loaded")
        sys.exit


if __name__ == '__main__':
    main()