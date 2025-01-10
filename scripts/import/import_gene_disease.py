#!/usr/bin/env python3

import os
import re
import argparse
import datetime
import MySQLdb
import xml.etree.ElementTree as ET
import warnings

"""
    Script      : import_gene_disease.py

    Description : Imports gene-disease associations from OMIM (Ensembl) and Mondo
                  The script can be run in two modes:
                    - import (default) : Runs a full import of the data. This mode should only be run on a newly created db.
                    - update : Updates the data if necessary. This mode should be used in the release cycle.

    Options
                --host : Ensembl server where core db is stored
                --port : Ensembl server's port
                --user : Ensembl server's username
                --database : Ensembl core db
                --g2p_host : G2P server
                --g2p_port : G2P server's port
                --g2p_database : G2P database
                --g2p_user : G2P server's username
                --g2p_password : G2P server's password
                --omim : option to import/update OMIM gene-disease data (default: 1)
                --mondo : option to import/update Mondo gene-disease data (default: 1)
                --mondo_file : Mondo owl file with gene-disease data (default: '')
                --update : option to run the update mode (default: 0)
"""

debug = 0

"""
    Retrieve meta info from G2P db
    Retrieve the latest OMIM and Mondo updates info
"""
def get_g2p_meta(db_host, db_port, db_name, user, password):
    meta_info = {}

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    sql = """ SELECT s.name, m.version, m.date_update
              FROM meta m
              LEFT JOIN source s ON s.id = m.source_id
              WHERE m.`key` = 'import_gene_disease'
              ORDER BY m.date_update desc
          """

    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        source_name = row[0]
        version = row[1]
        date_update = row[2]

        if source_name not in meta_info:
            meta_info[source_name] = { 
                "data_version": version,
                "date_update": date_update
            }
        elif date_update > meta_info[source_name]["date_update"]:
            meta_info[source_name] = { 
                "data_version": version,
                "date_update": date_update
            }

    return meta_info

"""
    Retrieve OMIM gene disease data from Ensembl core db
"""
def get_mim_gene_diseases(db_host, db_port, db_name, user, password):
    gene_diseases = {}

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    sql = """ SELECT a.value, g.stable_id, x.dbprimary_acc, x.description 
              FROM external_db e, xref x, object_xref o, gene g, gene_attrib a 
              WHERE e.external_db_id = x.external_db_id AND x.xref_id = o.xref_id 
              AND o.ensembl_id = g.gene_id AND e.db_name = 'MIM_MORBID' AND a.gene_id = g.gene_id 
              AND a.attrib_type_id = 4
          """

    sql_variation = """ SELECT pf.object_id, po.accession FROM phenotype_feature pf
                        LEFT JOIN phenotype p ON p.phenotype_id = pf.phenotype_id
                        LEFT JOIN phenotype_ontology_accession po ON po.phenotype_id = pf.phenotype_id
                        LEFT JOIN source s ON s.source_id = pf.source_id
                        WHERE pf.type = 'Gene' AND pf.object_id like 'ENSG%' AND s.name = 'MIM morbid' AND
                        po.accession like 'MONDO%'
                    """

    cursor.execute(sql)
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            if row[0] not in gene_diseases.keys():
                gene_diseases[row[0]] = [{ 'stable_id':row[1],
                                           'mim_id':row[2],
                                           'disease':row[3] }]
            else:
                gene_diseases[row[0]].append({ 'stable_id':row[1],
                                               'mim_id':row[2],
                                               'disease':row[3] })

    db.close()
    return gene_diseases


"""
    Insert OMIM gene disease data into G2P database
"""
def insert_gene_diseases(db_host, db_port, db_name, user, password, gene_diseases, version):
    sql_gene = f""" SELECT l.id, i.identifier FROM locus l
                    LEFT JOIN locus_identifier i on i.locus_id = l.id
                    LEFT JOIN source s on s.id = i.source_id
                    WHERE l.name = %s AND s.name = 'Ensembl'
                """

    sql_source = f""" SELECT id, name FROM source WHERE name = 'OMIM' OR name = 'Mondo' or name = 'Ensembl'
                 """

    sql_insert = f""" INSERT IGNORE INTO gene_disease(gene_id, disease, identifier, source_id)
                      VALUES(%s, %s, %s, %s)
                  """

    sql_meta = """ INSERT INTO meta(`key`, date_update, is_public, description, source_id, version)
                   VALUES(%s,%s,%s,%s,%s,%s)
               """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Fetch source id
    source_ids = {}
    cursor.execute(sql_source)
    data = cursor.fetchall()
    if len(data) != 0:
        for row in data:
            source_ids[row[1]] = row[0]

    for gd, gd_info in gene_diseases.items():
        cursor.execute(sql_gene, [gd])
        data = cursor.fetchall()
        gene_id = None
        gene_stable_id = None
        if len(data) != 0:
            for row in data:
                gene_id = row[0]
                gene_stable_id = row[1]
        if gene_id is not None:
            for info in gd_info:
                disease = info['disease'].split(';')[0]
                if info['stable_id'] == gene_stable_id:
                    cursor.execute(sql_insert, [gene_id, disease, info['mim_id'], source_ids['OMIM']])

    # Insert import info into meta
    cursor.execute(sql_meta, ['import_gene_disease',
                              datetime.datetime.now(),
                              0,
                              'Import OMIM gene disease associations from Ensembl',
                              source_ids['Ensembl'],
                              version])

    db.commit()
    db.close()


"""
    Retrieve the gene-disease association from the Mondo file in owl format
    Returns a dict "mondo_id": { "disease_name", "hgnc_id" }
"""
def get_mondo_gene_diseases(file):
    with open(file, "r") as fh:
        results = {}
        mondo_ids = []
        mondo_diseases = {}
        start_class = 0

        for event, elem in ET.iterparse(fh, events=("start", "end")):
            is_element = None
            disease_name = None
            disease_synonym = None
            hgnc_id = None

            if event == "start" and elem.tag == "{http://www.w3.org/2002/07/owl#}Class":
                # Some entries have a class inside a class
                # We want to skip those - only the first class is relevant
                if start_class == 1 and debug == 1:
                    print("Start class when another class still open")

                if start_class == 0:
                    start_class = 1
                    
                    if debug == 1:
                        ET.dump(elem)
                    
                    for i in elem.iter():
                        # Get mondo ID
                        if i.tag == "{http://www.geneontology.org/formats/oboInOwl#}id" and i.text is not None:
                            is_element = i.text

                        # Get disease name
                        if i.tag == "{http://www.w3.org/2000/01/rdf-schema#}label" and i.text is not None:
                            disease_name = i.text

                        # Get disease synonym
                        # This is necessary if for some reason it can't get the name from the label
                        if i.tag == "{http://www.geneontology.org/formats/oboInOwl#}hasExactSynonym" and i.text is not None:
                            disease_synonym = i.text

                        # Get HGNC ID from the class
                        # This is probably not necessary
                        if(i.tag == "{http://www.w3.org/2002/07/owl#}someValuesFrom"
                        and i.attrib is not None and "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource" in i.attrib.keys()
                        and "hgnc" in i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"]):
                            hgnc_id = i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].split("/")[-1]

                        if is_element is not None and disease_name is not None and not disease_name.startswith("obsolete"):
                            mondo_ids.append(is_element)
                        
                        # This is a special case
                        elif is_element is not None and disease_name is None and disease_synonym is not None:
                            mondo_ids.append(is_element)
                            disease_name = disease_synonym

            # Class is over
            if event == "end" and elem.tag == "{http://www.w3.org/2002/07/owl#}Class":
                start_class = 0

            # Save the disease name - this is necessary because the HGNC is outside of the Class
            if is_element is not None and disease_name is not None and is_element not in mondo_diseases:
                mondo_diseases[is_element] = disease_name

            # Get HGNC ID from outside the elem "Class"
            # At this point the is_element is None - we don't have the Mondo ID anymore
            # but we know that the Mondo ID is the previous ID analysed
            if event == "end" and elem.tag == "{http://www.w3.org/2002/07/owl#}Restriction" and is_element is None and len(mondo_ids) >= 1:
                is_element = mondo_ids[-1]

                for i in elem.iter():
                    if(i.tag == "{http://www.w3.org/2002/07/owl#}someValuesFrom" 
                       and "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource" in i.attrib.keys()
                       and "hgnc" in i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"]):
                        hgnc_id = i.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].split("/")[-1]
                        disease_name = mondo_diseases[is_element]
 
            if is_element is not None and disease_name is not None and hgnc_id is not None and is_element not in results:
                results[is_element] = { "disease": disease_name, "hgnc_id": hgnc_id }

    fh.close()

    return results

"""
    Insert Mondo gene disease data into G2P database
    This method runs a full import
"""
def insert_mondo_gene_diseases(db_host, db_port, db_name, user, password, gene_diseases, mondo_version):
    sql_gene = """ SELECT l.id FROM locus l
                    LEFT JOIN locus_identifier i on i.locus_id = l.id
                    LEFT JOIN source s on s.id = i.source_id
                    WHERE i.identifier = %s AND s.name = 'HGNC'
                """

    sql_source = """ SELECT id, name FROM source WHERE name = 'Mondo'
                 """

    sql_insert = """ INSERT IGNORE INTO gene_disease(gene_id, disease, identifier, source_id)
                      VALUES(%s, %s, %s, %s)
                  """

    sql_meta = """ INSERT INTO meta(`key`, date_update, is_public, description, source_id, version)
                   VALUES(%s,%s,%s,%s,%s,%s)
               """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Fetch source id
    source_id = None
    cursor.execute(sql_source)
    data = cursor.fetchall()
    for row in data:
        source_id = row[0]

    for gd, gd_info in gene_diseases.items():
        cursor.execute(sql_gene, [f"HGNC:{gd}"])
        data = cursor.fetchall()
        gene_id = None
        for row in data:
            gene_id = row[0]
        if gene_id is not None:
            for info in gd_info:
                cursor.execute(sql_insert, [gene_id, info['disease'], info['mondo_id'], source_id])

    # Insert import info into meta
    cursor.execute(sql_meta, ['import_gene_disease',
                              datetime.datetime.now(),
                              0,
                              'Import Mondo gene disease associations',
                              source_id,
                              mondo_version])

    db.commit()
    db.close()

"""
    Updates the Mondo gene-disease associations
    Affected tables: gene_disease, meta

    'gene_diseases' format example:
        MONDO:1040051 : {'disease': 'IMPDH1-related retinopathy', 'hgnc_id': '6052'}
"""
def update_mondo_gene_diseases(db_host, db_port, db_name, user, password, gene_diseases, mondo_version):
    g2p_current_data = {}

    sql_query = """ SELECT gd.disease, gd.identifier, li.identifier
                    FROM gene_disease gd
                    LEFT JOIN locus l on l.id = gd.gene_id
                    LEFT JOIN locus_identifier li on li.locus_id = l.id
                    LEFT JOIN source s ON s.id = li.source_id
                    LEFT JOIN source s_gd ON s_gd.id = gd.source_id
                    WHERE s.name = 'HGNC' and s_gd.name = 'Mondo'
                """

    sql_ins_gene_disease = """  """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    cursor.execute(sql_query)
    data = cursor.fetchall()
    for row in data:
        mondo_id = row[1]

        if mondo_id not in g2p_current_data:
            g2p_current_data[mondo_id] = {
                "disease": row[0],
                "hgnc_id": re.sub("HGNC:", "", row[2])
            }
        else:
            warnings.warn(f"Warning: Duplicated Mondo ID found in G2P '{mondo_id}'")

    db.close()

    # Iterate through the new mondo ids extracted from the owl file and
    # compare them with the existing ids in G2P
    for mondo, data in gene_diseases.items():
        print("\n", mondo, ":", data)
        if mondo in g2p_current_data:
            print("Found in g2p:", mondo, "; from g2p:", g2p_current_data[mondo])
            if data["disease"] != g2p_current_data[mondo]["disease"]:
                print("> to update")
        else:
            print("New:", mondo)

    # TODO: when we don't use the owl format, check which ids are obsolet, i.e. found in g2p but not in mondo

"""
    Fetch the Mondo data version from the owl file
"""
def fetch_mondo_version(file):
    version = None

    with open(file, "r") as fh:
        for event, elem in ET.iterparse(fh, events=("start", "end")):
            if event == "start" and elem.tag == "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}RDF":
                for i in elem.iter():
                    if i.tag == "{http://www.w3.org/2002/07/owl#}versionIRI":
                        resource = i.get("{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource")
                        version = re.search("[0-9]+\\-[0-9]+\\-[0-9]+", resource)
                        return version.group()

    return version


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--host", required=False, help="Ensembl core database host")
    parser.add_argument("--port", required=False, help="Ensembl core host port")
    parser.add_argument("--database", required=False, help="Ensembl core database name")
    parser.add_argument("--user", required=False, help="Username")
    parser.add_argument("--password", default='', help="Password (default: '')")
    parser.add_argument("--g2p_host", required=True, help="G2P database host")
    parser.add_argument("--g2p_port", required=True, help="G2P host port")
    parser.add_argument("--g2p_database", required=True, help="G2P database name")
    parser.add_argument("--g2p_user", required=True, help="Username")
    parser.add_argument("--g2p_password", default='', help="Password (default: '')")
    parser.add_argument("--omim", default=1, help="Import OMIM gene-disease data")
    parser.add_argument("--mondo", default=1, help="Import MONDO gene-disease data")
    parser.add_argument("--mondo_file", default='', help="MONDO owl file")
    parser.add_argument("--update", default=0, help="Run data update. By default, the script imports all data from scratch")

    args = parser.parse_args()

    db_host = args.host
    db_port = args.port
    db_name = args.database
    user = args.user
    password = args.password
    g2p_db_host = args.g2p_host
    g2p_db_port = int(args.g2p_port)
    g2p_db_name = args.g2p_database
    g2p_user = args.g2p_user
    g2p_password = args.g2p_password
    import_omim = int(args.omim)
    import_mondo = int(args.mondo)
    mondo_file = args.mondo_file
    update_mode = int(args.update)

    # Set the type of run
    run_import = 0
    if update_mode == 0:
        run_import = 1

    # Fetch gene-disease import info from G2P
    g2p_meta_info = get_g2p_meta(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password)

    ### Import or update OMIM ###
    if import_omim == 1:
        # Get version from Ensembl db name
        version = re.search("[0-9]+", db_name)

        # Compare G2P OMIM last update with Ensembl (OMIM) data version
        if "Ensembl" in g2p_meta_info:
            if run_import == 1:
                print("Cannot run import: G2P already has OMIM gene-disease associations. Please run the script in update mode (--update 1)")
                run_import = 0 # set flag to 0 to make sure we don't insert the data again

            if g2p_meta_info["Ensembl"]["data_version"] >= version.group() and update_mode == 1:
                print("Omim gene-disease data is the latest. Skipping update.")
                update_mode = 0 # set flag to 0 because data in G2P is already updated

        # Only fetch OMIM data from Ensembl if import/update can be done
        if update_mode == 1 or run_import == 1:
            # Retrive OMIM gene-disease from Ensembl core db
            print("Getting OMIM gene-disease associations...")
            gene_diseases = get_mim_gene_diseases(db_host, int(db_port), db_name, user, password)
            print("Getting OMIM gene-disease associations... done")

        # Run the import
        if run_import == 1:
            print("> Running OMIM full import")
            # Insert OMIM gene-disease data into G2P db
            print("Inserting OMIM gene-disease associations into G2P...")
            insert_gene_diseases(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, gene_diseases, version.group())
            print("Inserting OMIM gene-disease associations into G2P... done")

        # Run the update
        elif update_mode == 1:
            print("> Running OMIM update")
            # Update OMIM gene-disease data in G2P db
            # TODO

    ### Import or update Mondo ###
    if import_mondo == 1:
        """
            Input file: it is necessary to perform some cleaning before the import.
            The Mondo file should only have the following type of entries:
                <!-- http://purl.obolibrary.org/obo/MONDO_0011584 -->
            
            The owl file can be downloaded from https://mondo.monarchinitiative.org/pages/download/
        """

        # Fetch Mondo version from file
        mondo_version = fetch_mondo_version(mondo_file)
        mondo_version_obj = datetime.datetime.strptime(mondo_version, "%Y-%m-%d")

        # Compare G2P OMIM last update with Ensembl (OMIM) data version
        if "Mondo" in g2p_meta_info:
            if run_import == 1:
                print("Cannot run import: G2P already has Mondo gene-disease associations. Please run the script in update mode (--update 1)")
                run_import = 0 # set flag to 0 to make sure we don't insert the data again

            g2p_mondo_version_obj = datetime.datetime.strptime(g2p_meta_info["Mondo"]["data_version"], "%Y-%m-%d")
            if g2p_mondo_version_obj >= mondo_version_obj and update_mode == 1:
                print("Mondo gene-disease data is the latest. Skipping update.")
                update_mode = 0 # set flag to 0 because data in G2P is already updated

        # Only fetch OMIM data from Ensembl if import/update can be done
        if update_mode == 1 or run_import == 1:
            # Retrive Mondo gene-disease from Mondo owl file
            print("Getting Mondo gene-disease associations...")
            mondo_gene_diseases = get_mondo_gene_diseases(mondo_file)
            print("Getting Mondo gene-disease associations... done")

        # Run the import
        if run_import == 1:
            print("> Running Mondo full import")

            # Format data for the full import
            new_mondo_data = {}
            for mondo, data in mondo_gene_diseases.items():
                if data["hgnc_id"] not in new_mondo_data:
                    new_mondo_data[data["hgnc_id"]] = [{ "mondo_id": mondo, "disease": data["disease"] }]
                else:
                    new_mondo_data[data["hgnc_id"]].append({ "mondo_id": mondo, "disease": data["disease"] })

            # Insert Mondo gene-disease data into G2P db
            print("Inserting Mondo gene-disease associations into G2P...")
            insert_mondo_gene_diseases(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, new_mondo_data, mondo_version)
            print("Inserting Mondo gene-disease associations into G2P... done")

        # Run the update
        elif update_mode == 1:
            print("> Running Mondo update")

            # Update Mondo gene-disease data in G2P db
            update_mondo_gene_diseases(g2p_db_host, g2p_db_port, g2p_db_name, g2p_user, g2p_password, mondo_gene_diseases, mondo_version)

if __name__ == '__main__':
    main()
