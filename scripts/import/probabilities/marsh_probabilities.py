import os
import re 
import argparse
import MySQLdb 


def get_details_from_file(file):
    with open(file, "r") as opened_file:
        lines = opened_file.readlines()
        skip_header = lines[1:]

    
    return skip_header

def get_locus_id_from_g2p_db(list_lines,host,port,db,password,user):

    get_locus_source = """ SELECT locus_id from locus where name = %s
"""

    database = MySQLdb.connect(host=host,port=port,user=user,passwd=password,db=db)
    cursor = database.cursor()

    for line in list_lines:
        gene_symbol = line[0]
        cursor.execute(get_locus_source, (gene_symbol,))
        locus_id = cursor.fecthone()

        if locus_id:
            line.append(locus_id)
        else:
            line.append(None)
        
    cursor.close()
    database.close()


    return list_lines


def get_hgnc_id_from_g2p_db(list_lines,host,port,db,password,user)





def main():
    parser = argparse.ArgumentParser(description="This script is used to import the probabilities from a file and imports it to the gene_stats table in the G2P DB")
    parser.add_argument("-h", "--host", required=True, help="Host of the database were the data is to be imported")
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