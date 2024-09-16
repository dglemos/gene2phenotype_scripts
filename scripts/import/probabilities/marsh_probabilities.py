import os
import re 
import argparse
import MySQLdb 


def get_details_from_file(file):
    with open(file, "r") as opened_file:
        lines = opened_file.readlines()
        skip_header = lines[1:]

    
    return skip_header

def get_locus_id_from_g2p_db(list_lines, host, port, db, password, user):

    get_locus_source = """ SELECT id from locus where name = %s
"""

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
            else:
                line.append(None)
            
            list_lines[i] = line
        
    cursor.close()
    database.close()


    return list_lines


def get_hgnc_id_from_g2p_db(list_lines, host, port, db, password, user):
    
    get_locus_identifier = """ SELECT id from locus_identifier where locus_id = %s """

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    for i,line in enumerate(list_lines):
        locus_id = line[4]
        if locus_id is not None:
            cursor.execute(get_locus_identifier, (locus_id,))
            locus_identifer_id = cursor.fetchone()

            if locus_identifer_id:
                line.append(locus_identifer_id[0])
            else: 
                line.append(None)
        else:
            line.append(None)
            
            list_lines[i] = line

    cursor.close()
    database.close()

    return list_lines


def insert_into_gene_stats(list_lines, host, port, db, password, user):
    insert_into_gene_stats_query = """ INSERT into gene_stats (gene_id, gene_symbol, hgnc_id, score, source_id) VALUES (%s, %s, %s, %s, %s)
"""
    get_source_query = """ SELECT id from source where name = 'Marsh Mechanism probabilities'"""

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    cursor.execute(get_source_query)
    source_id = cursor.fetchone()
    source_id = source_id[0]


    for line in list_lines:
        if line[5] is not None:
            cursor.execute(insert_into_gene_stats_query, (line[4], line[0], line[5], line[2], source_id))

    cursor.close()
    database.close()
   


def main():
    parser = argparse.ArgumentParser(description="This script is used to import the probabilities from a file and imports it to the gene_stats table in the G2P DB")
    parser.add_argument("--host", required=True, help="Host of the database were the data is to be imported")
    parser.add_argument("-p", "--port", required=True, help="Port information of the database to be imported")
    parser.add_argument("-d", "--database", required=True, help="G2P Database to import the information into")
    parser.add_argument("-pwd", "--password", required=True, help="Paaword for the G2P database information")
    parser.add_argument("-u", "--user", required=True, help="User information for the G2P database")
    parser.add_argument("-f", "--file", required=True, help="File containing the information for the score, can either be one file or files seperated by ,")

    args = parser.parse_args()

    host = args.host
    port = args.port
    db = args.database
    pwd = args.password
    user = args.user
    file = args.file
    port = int(port)

    if file: 
        print("Getting details from file")
        file_lines = get_details_from_file(file)
        print("Getting locus id from the G2P DB")
        get_locus_id_from_g2p_db(file_lines, host, port, db, pwd, user)
        print("Getting HGNC locus identifier from G2P DB")
        get_hgnc_id_from_g2p_db(file_lines, host, port, db, pwd, user)
        print("Inserting into gene stats")
        insert_into_gene_stats(file_lines, host, port, db, pwd, user)



if __name__ == '__main__':
    main()