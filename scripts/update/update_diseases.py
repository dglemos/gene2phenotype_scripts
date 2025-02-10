#!/usr/bin/env python3

import os.path
import sys
import argparse
import MySQLdb
import requests
import configparser

lgd_disease_to_update = []

def dump_data(db_host, db_port, db_name, user, password):
    gene_records = {} # key = gene symbol; value = list of records
    diseases = {} # key = disease name; value = list of record ids

    sql = """   SELECT l.name, g2p.stable_id, d.name, d.id, a.value, a1.value, lgd.id
                FROM locus_genotype_disease lgd
                LEFT JOIN locus l ON l.id = lgd.locus_id
                LEFT JOIN g2p_stableid g2p ON g2p.id = lgd.stable_id
                LEFT JOIN disease d ON d.id = lgd.disease_id
                LEFT JOIN attrib a ON a.id = lgd.genotype_id
                LEFT JOIN attrib a1 ON a1.id = lgd.confidence_id
                ORDER BY l.name
          """
    
    sql_disease = """
                        SELECT d.name, d.id, lgd.id
                        FROM disease d
                        LEFT JOIN locus_genotype_disease lgd ON lgd.disease_id = d.id
                        ORDER BY d.name
                  """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        if row[0] not in gene_records:
            gene_records[row[0]] = [{
                "stable_id": row[1],
                "disease_name": row[2],
                "disease_id": row[3],
                "genotype": row[4],
                "confidence": row[5],
                "record_id": row[6]
            }]
        else:
            gene_records[row[0]].append({
                "stable_id": row[1],
                "disease_name": row[2],
                "disease_id": row[3],
                "genotype": row[4],
                "confidence": row[5],
                "record_id": row[6]
            })

    cursor.execute(sql_disease)
    data_disease = cursor.fetchall()
    for row in data_disease:
        if row[0] not in diseases:
            diseases[row[0]] = {
                "disease_id": row[1],
                "records": [row[2]]
            }
        else:
            diseases[row[0]]["records"].append(row[2])

    db.close()

    return gene_records, diseases

def read_file(file, gene_records, diseases):
    file_output = "gene_disease_not_updated.txt"
    file_not_updated = "diseases_not_updated.txt"
    diseases_to_update = []

    with open(file, "r") as fh, open(file_output, "w") as wr, open(file_not_updated, "w") as wr_diseases:
        wr.write("gene symbol\tdisease name\tdisease name formatted\tdiseases found in G2P linked to gene\n")

        # Header:
        # gene symbol, gene mim, disease name, disease_name_formatted, Updated , disease mim, confidence category,
        # allelic requirement, mutation consequence, phenotypes, organ specificity list, pmids, panel, prev symbols,
        # hgnc id, gene disease pair entry date, cross cutting modifier, mutation consequence flag, confidence value,
        # flag, comments, variant consequence, disease ontology
        for line in fh:
            if not line.startswith("gene symbol"):
                data = line.split("\t")
                gene_symbol = data[0]
                current_disease = data[1] # before searching for the disease we have to check if 'gene-related' is used
                new_disease = data[2]
                is_updated = data[3]

                if is_updated == "Yes":
                    continue

                # Check if the record can be found in G2P
                try:
                    db_data = gene_records[gene_symbol]
                except:
                    print(f"WARNING: {gene_symbol} not found in G2P")
                else:
                    # print(f"\n{gene_symbol}; {current_disease} (New: {new_disease})")
                    to_update = 0
                    list_disease = []
                    for record in db_data:
                        if record["disease_name"].endswith(current_disease) and record["disease_name"] != new_disease:
                            to_update = 1
                            # print("-> to update (record from db):", record)
                            # Check if new disease name is already in the db
                            if new_disease in diseases:
                                # If the new disease name already exists then add it to replace existing name
                                # Before running the updating, the API checks if records are not the same
                                print(f"Update LGD records: replace disease_id {record['disease_id']} by {diseases[new_disease]['disease_id']}")
                                lgd_disease_to_update.append(
                                    {
                                        "disease_id": record["disease_id"],
                                        "new_disease_id": diseases[new_disease]["disease_id"]
                                    }
                                )
                            elif gene_symbol not in new_disease:
                                # This should not happen
                                # Action: print this warning and not update
                                wr_diseases.write("Gene not found in new disease name\t"+line)
                            else:
                                # Create list of diseases to update
                                diseases_to_update.append(
                                    {
                                        "id": record["disease_id"],
                                        "name": new_disease
                                    }
                                )

                        else:
                            list_disease.append(record["disease_name"])

                    # Print to file genes that won't have disease updates
                    if not to_update:
                        wr.write(f"{gene_symbol}\t{current_disease}\t{new_disease}\t{list_disease}\n")

    return diseases_to_update

def update_diseases(diseases_to_update, config_api):
    """
        Method to update the disease names.
        This method calls the gene2phenotype api and updates the database defined in gene2phenotype_project/config.ini
    """
    disease_url = "update/diseases/"
    lgd_disease_url = "lgd_disease_updates/"
    login_url = "login/"

    data = {
        "username": config_api["username"],
        "password": config_api["password"]
    }

    response = requests.post(config_api["api_url"] + login_url, json=data)
    if response.status_code == 200:
        try:
            response_update = requests.post(config_api["api_url"] + disease_url, json=diseases_to_update, cookies=response.cookies)
            if response_update.status_code == 200:
                response_json = response_update.json()
                print("Diseases updated successfully:", response_json)
                # Some diseases were probably not updated: if they are already in the db
                # For those we can update the disease_id in lgd to point to the existing disease
                if "errors" in response_json:
                    for error in response_json["errors"]:
                        print(f"Update LGD records: replace disease_id {error['id']} by {error['existing_id']}")
                        lgd_disease_to_update.append(
                            {
                                "disease_id": error['id'],
                                "new_disease_id": error['existing_id']
                            }
                        )
                    response_update_lgd = requests.post(config_api["api_url"] + lgd_disease_url, json=lgd_disease_to_update, cookies=response.cookies)
                    if response_update_lgd.status_code == 200:
                        print("LGD records updated successfully:", response_update_lgd.json())
                    else:
                        print("Failed to update LGD records:", response_update_lgd.status_code, response_update_lgd.json())
            else:
                print("Failed to update diseases:", response_update.status_code, response_update.json())
        except Exception as e:
            print("Error:", e)
    else:
        print("Error: cannot login into G2P")


def main():
    """
        Params:
                --config : Config file name containing the database and API connection info (mandatory)
                        File format is the following: 
                            [database]
                            host = <>
                            port = <>
                            user = <>
                            password = <>
                            name = <>

                            [api]
                            api_url = <>
                            username = <>
                            password = <>

                --file : Tab delimited file with all diseases to be updated (mandatory)
                    File format is the following:
                        gene symbol\tdisease name\tdisease name formatted\tUpdated
    """

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", required=True, help="Config file")
    parser.add_argument("--file", required=True, help="Tab delimited file with all diseases to be updated")
    args = parser.parse_args()

    file = args.file
    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config['database']['host']
    db_port = config['database']['port']
    db_name = config['database']['name']
    user = config['database']['user']
    password = config['database']['password']

    print("Dump data from G2P...")
    gene_records, diseases = dump_data(db_host, int(db_port), db_name, user, password)
    print("Dump data from G2P... done\n")

    if os.path.isfile(file):
        expected_columns = ['gene symbol', 'disease name', 'disease name formatted', 'Updated']

        with open(file, "r", encoding="utf-8") as fh:
            header = fh.readline().strip().split("\t")
            if header[:4] != expected_columns:
                sys.exit(f"Error: File format is incorrect. Found: {header[:4]}; Expected: {expected_columns}")

        print("Parsing diseases to update...")
        diseases_to_update = read_file(file, gene_records, diseases)
        print("Parsing diseases to update... done\n")

        print("Updating disease names...")
        update_diseases(diseases_to_update, config['api'])
        print("Updating disease names... done\n")

    else:
        print(f"Input file is invalid '{file}'")


if __name__ == '__main__':
    main()