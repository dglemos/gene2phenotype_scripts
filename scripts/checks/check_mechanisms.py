#!/usr/bin/env python3

import configparser
import argparse
import MySQLdb

def dump_old_data(db_host, db_port, db_name, user, password):
    records = {}
    attribs = {}

    sql = """   select gfd.genomic_feature_disease_id, gf.gene_symbol, d.name, gfd.allelic_requirement_attrib, gfd.mutation_consequence_attrib, gfd.variant_consequence_attrib from genomic_feature_disease gfd
                left join genomic_feature gf on gfd.genomic_feature_id = gf.genomic_feature_id
                left join disease d on gfd.disease_id = d.disease_id
                where gfd.variant_consequence_attrib like "%,%"
          """
    
    sql_attribs = """
                        SELECT attrib_id, value
                        FROM attrib
                  """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    # Get attribs list
    cursor.execute(sql_attribs)
    attribs_data = cursor.fetchall()
    for row in attribs_data:
        attribs[row[0]] = row[1]

    # Get records
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        allelic_requirements = row[3].split(",")
        mutation_consequences = row[4].split(",")
        variant_consequences = row[5].split(",")

        allelic_requirement_list = []
        for ar in allelic_requirements:
            allelic_requirement_list.append(attribs[int(ar)])
        
        mutation_consequence_list = []
        for mc in mutation_consequences:
            mutation_consequence_list.append(attribs[int(mc)])

        variant_consequence_list = []
        for vc in variant_consequences:
            consequence = attribs[int(vc)]
            if consequence == "gain_of_function_variant" or consequence == "loss_of_function_variant":
                variant_consequence_list.append(consequence)

        if variant_consequence_list:
            records[row[0]] = {
                "gene_symbol": row[1],
                "disease_name": row[2],
                "allelic_requirement": allelic_requirement_list,
                "mutation_consequence": mutation_consequence_list,
                "variant_consequence": variant_consequence_list
            }

    db.close()

    return records

def dump_data(db_host, db_port, db_name, user, password):
    records = {}

    sql = """   select l.name, d.name, a1.value, m.value, g2p.stable_id, g.old_g2p_id from locus_genotype_disease lgd
                left join locus l on l.id = lgd.locus_id
                left join disease d on d.id = lgd.disease_id
                left join attrib a1 on a1.id = lgd.genotype_id
                left join cv_molecular_mechanism m on m.id = lgd.mechanism_id
                left join gencc_submission g on g.g2p_stable_id = lgd.stable_id
                left join g2p_stableid g2p on g2p.id = lgd.stable_id
          """

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()
    # Get records
    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        records[row[4]] = {
            "gene_symbol": row[0],
            "disease_name": row[1],
            "allelic_requirement": row[2],
            "mechanism": row[3],
            "old_g2p_id": row[5]
        }

    db.close()

    return records

def compare_records(old_records, records):
    for old_record_id in old_records:
        for new_record_id in records:
            old_genotype = old_records[old_record_id]['allelic_requirement'][0]
            if old_genotype == "monoallelic_X_het":
                old_genotype = "monoallelic_X_heterozygous"
            elif old_genotype == "monoallelic_X_hem":
                old_genotype = "monoallelic_X_hemizygous"

            if old_record_id == records[new_record_id]["old_g2p_id"]:
                print(f"-> found: {old_record_id} ({old_records[old_record_id]['variant_consequence']}), {new_record_id} ({records[new_record_id]['mechanism']})")
            elif (old_records[old_record_id]["gene_symbol"] == records[new_record_id]["gene_symbol"] 
                  and old_records[old_record_id]['disease_name'].lower() == records[new_record_id]['disease_name'].lower()
                  and old_genotype == records[new_record_id]['allelic_requirement']):
                print(f"-> checking: old_id {old_record_id} ({old_records[old_record_id]['variant_consequence']}) new_id {new_record_id} ({records[new_record_id]['mechanism']})")


def main():
    parser = argparse.ArgumentParser(description="Script to compare mechanisms between new and old system")
    parser.add_argument("--config", default='', help="Config file with details to dbs")
    args = parser.parse_args()

    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config['database']['host']
    db_port = config['database']['port']
    db_name = config['database']['name']
    user = config['database']['user']
    password = config['database']['password']
    old_db_host = config['old_database']['host']
    old_db_port = config['old_database']['port']
    old_db_name = config['old_database']['name']
    old_user = config['old_database']['user']
    old_password = config['old_database']['password']

    old_records = dump_old_data(old_db_host, int(old_db_port), old_db_name, old_user, old_password)

    # for r in old_records:
    #     print("->", r, ":", old_records[r])

    records = dump_data(db_host, int(db_port), db_name, user, password)

    compare_records(old_records, records)

if __name__ == '__main__':
    main()