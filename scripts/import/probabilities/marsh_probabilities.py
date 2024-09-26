import os
import re
import sys
from datetime import datetime
import argparse
import MySQLdb

attrib_mapping = {
    "gain_of_function_mp" : 63,
    "loss_of_function_mp" : 64,
    "dominant_negative_mp" : 65
}

key = "Badonyi_score"
def get_details_from_file(file):
    """
        Extracts lines from a file, skipping the header.

        This function reads the contents of a file, skipping the first line
        (assumed to be a header) and returning the remaining lines.

        Args:
            file (str): The path to the file to be read.

        Returns:
            None: This function does not return a value, but it stores the
            lines after the header in the variable `skip_header`.

        Example:
            get_details_from_file("data.txt")
    """
    with open(file, "r") as opened_file:
        lines = opened_file.readlines()
        skip_header = lines[1:]

    
    return skip_header

def get_locus_id_from_g2p_db(list_lines, host, port, db, password, user):
    """
        Retrieves locus IDs from a gene-to-phenotype (G2P) database and appends them to the input list.

        This function connects to a MySQL database and executes a query to fetch the locus ID
        for each gene symbol (first item in each line) from the 'locus' table. It appends the
        fetched locus ID to the corresponding line in `list_lines`. If a locus ID is not found,
        it appends `None`.

        Args:
            list_lines (list of str): A list of strings or lists, where each element contains
                                  a gene symbol as the first item.
            host (str): The database host address.
            port (int): The port number to connect to the database.
            db (str): The name of the database.
            password (str): The password for authenticating with the database.
            user (str): The username for authenticating with the database.

        Returns:
            list: The modified `list_lines` with locus IDs appended to each line.
              If no locus ID is found for a gene symbol, `None` is appended.

        Example:
            list_lines = ["BRCA1 info", "TP53 info"]
            updated_lines = get_locus_id_from_g2p_db(list_lines, "localhost", 3306, "g2p_db", "password", "user")
    """

    get_locus_source = """ SELECT id from locus where name = %s
"""
    get_locus_attrib_source = """SELECT locus_id from locus_attrib where value = %s
"""
    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    for i, line in enumerate(list_lines):
        if isinstance(line, str):
            line = line.split()
            gene_symbol = line[0]
            cursor.execute(get_locus_source, (gene_symbol,))
            cursor.execute(get_locus_attrib_source, (gene_symbol,))
            locus_id = cursor.fetchone()
            locus_attrib_id = cursor.fetchone()
    

            if locus_id:
                line.append(locus_id[0])
            elif locus_attrib_id:
                line.append(locus_attrib_id[0])
            else:
                line.append(None)
            
            list_lines[i] = line
        
    cursor.close()
    database.close()


    return list_lines


def insert_details_into_meta(host, port, db, password, user):
    """
        Inserts details into the 'meta' table of the given database.

        This function establishes a connection to a MySQL database and inserts a record into the 'meta' table
        with the specified key, current timestamp, description, version, and source ID. The source ID is 
        retrieved via the `get_source_details` function.

        Parameters:
        ----------
        host : str
            The hostname or IP address of the MySQL database server.
        
        port : int
            The port number used to connect to the MySQL database server.
        
        db : str
            The name of the MySQL database where the 'meta' table resides.
        
        password : str
            The password to authenticate the MySQL user.
        
        user : str
            The username to authenticate with the MySQL database server.

        Other Variables:
        ----------
        source_id : int
            The ID obtained from the source details function used to populate the 'source_id' column.
        
        description : str
            A brief description of the entry (default: "Baydoni & Marsh probabilities").
        
        formatted_datetime : str
            The current date and time in the format 'YYYY-MM-DD HH:MM:SS.microseconds'.
        
        version : int
            The version of the entry being inserted (default: 1).
        
        insert_into_meta_query : str
            The SQL query used to insert data into the 'meta' table.

        Raises:
        ----------
        MySQLdb.Error
            If any error occurs during the database connection or query execution.
        
        Example:
        ----------
        insert_details_into_meta('localhost', 3306, 'mydatabase', 'password123', 'root')

        This inserts a record into the 'meta' table of the 'mydatabase' database using the current timestamp, 
        a static description, and source details retrieved from the database.

    """
    
    source_id = get_source_details(host, port, db, password, user)
    description = "Baydoni & Marsh probabilities"
    
    insert_into_meta_query = """ INSERT into meta(key, date_update, description, version, source_id) VALUES (%s, %s, %s, %s, %s)

"""
    current_datetime = datetime.now()


    formatted_datetime = current_datetime.strftime('%Y-%m-%d %H:%M:%S.%f')
    version = 1
    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()
    
    cursor.execute(insert_into_meta_query, (key, formatted_datetime, description, version, source_id ))
    cursor.close()
    database.close()

    

def get_source_details(host, port, db, password, user):
    """
        Retrieves the source ID from the 'source' table in the database.

        This function connects to a MySQL database and retrieves the `id` of the source where the `name` 
        matches 'Marsh Mechanism probabilities'. It returns the corresponding `source_id`.

        Parameters:
        ----------
        host : str
            The hostname or IP address of the MySQL database server.
        
        port : int
            The port number used to connect to the MySQL database server.
        
        db : str
            The name of the MySQL database containing the 'source' table.
        
        password : str
            The password to authenticate the MySQL user.
        
        user : str
            The username to authenticate with the MySQL database server.

        Returns:
        ----------
        source_id : int
            The ID of the source with the name 'Marsh Mechanism probabilities'.

        Raises:
        ----------
        MySQLdb.Error
            If any error occurs during the database connection or query execution.

        Example:
        ----------
        source_id = get_source_details('localhost', 3306, 'mydatabase', 'password123', 'root')

        This retrieves the `id` from the 'source' table where `name` equals 'Marsh Mechanism probabilities'
        and returns the corresponding `source_id`.
    """

    get_source_query = """ SELECT id from source where name = 'Marsh Mechanism probabilities'"""

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    cursor.execute(get_source_query)
    source_id = cursor.fetchone()
    source_id = source_id[0]

    cursor.close()
    database.close()

    return source_id

def insert_into_gene_stats(list_lines, host, port, db, password, user, attrib):

    """
        Inserts gene-related data into the `gene_stats` table of a MySQL database.

        Parameters:
        ----------
        list_lines : list of lists
            A list containing gene-related data, where each sublist represents a row of values with the following elements:
            [0]: gene_symbol (str)
            [1]: Not used in the query
            [2]: score (float)
            [3]: Not used in the query
            [4]: gene_id (int)
            [5]: locus_identifier_id (int or None)
        
        host : str
            The hostname or IP address of the MySQL server.

        port : int
            The port number on which the MySQL server is listening.

        db : str
            The name of the database to connect to.

        password : str
            The password for the database user.

        user : str
            The username for authenticating with the MySQL server.

        attrib : str
            An attribute key that will be mapped to a specific description ID from a predefined `attrib_mapping` dictionary.
            This is used as the `description_id` in the query.

        Functionality:
        -------------
        - Connects to the MySQL database using the provided connection details.
        - Retrieves the `source_id` corresponding to the 'Marsh Mechanism probabilities' source from the `source` table.
        - Iterates over each sublist in `list_lines`, and for each line where the `locus_identifier_id` (index 5) is not None,
        it inserts the gene data into the `gene_stats` table using the `insert_into_gene_stats_query`.
        - Closes the database connection after all the insert operations are completed.

        Example usage:
        --------------
        list_lines = [
            ['gene1', 'x', 0.85, 'x', 101, 201],
            ['gene2', 'y', 0.67, 'y', 102, None],
            ['gene3', 'z', 0.92, 'z', 103, 202]
        ]
        
        insert_into_gene_stats(list_lines, 'localhost', 3306, 'genedb', 'password123', 'user123', 'attrib_key')
        
        Notes:
        -----
        The function assumes that the table and column names, as well as the structure of the input list, align with the schema.
    """

    insert_into_gene_stats_query = """ INSERT into gene_stats (gene_symbol, gene_id, score, source_id, description_id) VALUES (%s, %s, %s, %s, %s, %s)
"""
    source_id = get_source_details(host, port, db, password, user)

    attrib_value = attrib_mapping.get(attrib)

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()


    for line in list_lines:
        if line[5] is not None:
            cursor.execute(insert_into_gene_stats_query, (line[0], line[4], line[2], source_id, attrib_value ))

    cursor.close()
    database.close()
   


def main():
    """
        Main function to import gene probability data into the gene_stats table of a G2P database.

        This script parses command-line arguments to obtain the necessary parameters for connecting to a 
        MySQL database, retrieves data from a specified file, processes it, and then inserts the data into 
        the `gene_stats` table.

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
            Path to the file containing the gene probability data. The file can be a single file or multiple 
            files separated by commas.

        -a, --attrib : str (required)
            The attribute type in the file being inserted. Allowed types include:
            'gain_of_function_mp', 'loss_of_function_mp', and 'dominant_negative_mp'.
            This argument will be used to map the `attrib` value to the corresponding description ID.

        Functionality:
        -------------
        1. Parses the command-line arguments to obtain the necessary database connection details, file path, and attribute type.
        2. Checks if the `attrib` value is valid by confirming that it exists in the `attrib_mapping` dictionary.
        If the attribute type is not valid, the script exits with an error message.
        3. Reads the file and retrieves its contents using `get_details_from_file()`.
        4. Obtains the necessary `locus_identifier_id` values from the G2P database using the `get_locus_id_from_g2p_db()` function.
        5. Retrieves the HGNC (HUGO Gene Nomenclature Committee) locus identifier from the G2P database using the `get_hgnc_id_from_g2p_db()` function.
        6. Finally, the processed data is inserted into the `gene_stats` table using the `insert_into_gene_stats()` function.

        Example usage:
        --------------
        python script.py --host localhost --port 3306 --database g2p_db --password secret_pwd --user db_user --file data.csv --attrib gain_of_function_mp

        Notes:
        -----
        - Ensure that the database and tables are correctly configured before running the script.
        - The file provided should follow the format expected by the script for proper processing.
        - This script requires the `attrib_mapping` dictionary to be predefined, which maps attribute types to description IDs.
    """

    parser = argparse.ArgumentParser(description="This script is used to import the probabilities from a file and imports it to the gene_stats table in the G2P DB")
    parser.add_argument("--host", required=True, help="Host of the database were the data is to be imported")
    parser.add_argument("-p", "--port", required=True, help="Port information of the database to be imported")
    parser.add_argument("-d", "--database", required=True, help="G2P Database to import the information into")
    parser.add_argument("-pwd", "--password", required=True, help="Paaword for the G2P database information")
    parser.add_argument("-u", "--user", required=True, help="User information for the G2P database")
    parser.add_argument("-f", "--file", required=True, help="File containing the information for the score, can either be one file or files seperated by ,")
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
        print("Getting HGNC locus identifier from G2P DB")
        get_hgnc_id_from_g2p_db(file_lines, host, port, db, pwd, user)
        print("Inserting into gene stats")
        insert_into_gene_stats(file_lines, host, port, db, pwd, user, attrib)



if __name__ == '__main__':
    main()