import requests
from requests.adapters import HTTPAdapter, Retry
import re
import MySQLdb
import datetime
import argparse

# Uniprot data fetch URL
url = 'https://rest.uniprot.org/uniprotkb/search?query=(organism_id:9606)&(reviewed:true)&fields=accession,cc_function,xref_mim,xref_hgnc,gene_primary&size=500' # To be changed

# Configuration to fetch Uniprot data
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)

# Function to fetch Uniprot data
def fetch_all_data():
    total_items = []
    for batch, total in get_batch(url):
        current_batch_json = batch.json()
        for item in current_batch_json["results"]:
            current_item = {}
            current_item["gene_symbol"] = item["genes"][0]["geneName"]["value"]
            current_item["accession"] = item["primaryAccession"]
            current_item["protein_function"] = item["comments"][0]["texts"][0]["value"]
            current_item["HGNC"] = item["uniProtKBCrossReferences"][0]["id"]
            current_item["MIM"] = item["uniProtKBCrossReferences"][1]["id"]
            total_items.append(current_item)
    print("Uniprot data successfully fetched via Uniprot API.")
    return total_items

# Function to insert Uniprot data to database
def insert_uniprot_data(db_host, db_port, db_name, db_user, db_password, total_items):
    sql_source = f""" SELECT id, name FROM source WHERE name = 'UniProt' or name = 'HGNC'
                 """
    sql_identifier = f""" SELECT identifier, locus_id FROM locus_identifier WHERE source_id = %s
                 """
    insert_sql = f""" INSERT INTO uniprot_annotation(uniprot_accession, gene_id, hgnc, gene_symbol, mim, protein_function, source_id)
                  VALUES(%s, %s, %s, %s, %s, %s, %s)
              """
    sql_meta = """ INSERT INTO meta(`key`, date_update, is_public, description, source_id, version)
               VALUES(%s,%s,%s,%s,%s,%s)
           """

    db = MySQLdb.connect(host=db_host, port=db_port, user=db_user, passwd=db_password, db=db_name)
    cursor = db.cursor()
    # Fetch source ids
    source_ids = {}
    cursor.execute(sql_source)
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            source_ids[row[1]] = row[0]
    # Fetch locus identifiers
    identifier_to_locus_id_map = {}
    cursor.execute(sql_identifier, [source_ids['HGNC']])
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            identifier_to_locus_id_map[row[0]] = row[1]
    # Insert Uniprot data
    insert_count = 0
    for item in total_items:
        if item["HGNC"] in identifier_to_locus_id_map:
            cursor.execute(insert_sql, [item["accession"],
                                  identifier_to_locus_id_map[item["HGNC"]],
                                  int(item["HGNC"][5:]), # To be changed once column datatype changes
                                  item["gene_symbol"],
                                  item["MIM"],
                                  item["protein_function"],
                                  source_ids['UniProt']])
            insert_count+=1
    # Insert import info into meta table
    cursor.execute(sql_meta, ['import_uniprot',
                              datetime.datetime.now(),
                              0,
                              'Import Uniprot data',
                              source_ids['UniProt'],
                              'ensembl_111'])
    db.commit()
    db.close()
    print("Uniprot data successfully inserted into G2P database.")
    print(f'Total Uniprot entries fetched: {len(total_items)}')
    print(f'Total Uniprot entries inserted: {insert_count}')
    print("Note: Only Uniprot data entries with existing Gene information in the database will be inserted.")

def main():
    # Fetch Database configuration from User
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--g2p_host", required=True, help="G2P database host")
    parser.add_argument("--g2p_port", required=True, help="G2P database host port")
    parser.add_argument("--g2p_database", required=True, help="G2P database name")
    parser.add_argument("--g2p_user", required=True, help="G2P database Username")
    parser.add_argument("--g2p_password", default='', help="G2P database Password (default: '')")

    args = parser.parse_args()

    g2p_db_host = args.g2p_host
    g2p_db_port = int(args.g2p_port)
    g2p_db_name = args.g2p_database
    g2p_user = args.g2p_user
    g2p_password = args.g2p_password

    total_items = fetch_all_data()
    insert_uniprot_data(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, total_items)

if __name__ == '__main__':
    main()
