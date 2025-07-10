#!/usr/bin/env python3

import sys
import re
import argparse
import MySQLdb
import configparser
from pathlib import Path

"""
    Script to dump the table publication_families to a file.
    The file is going to be used to populate the new schema.
    Update done in release v7.
"""


def fetch_data_from_db(db_host, db_port, db_name, user, password):
    """
    Fetch the contents of the table publication families.
    """
    families_data = {}

    db = MySQLdb.connect(host=db_host, port=db_port, user=user, passwd=password, db=db_name)
    cursor = db.cursor()

    sql = """
        select pf.families, pf.affected_individuals, pf.ancestries, pf.is_deleted, a.value,
        p.id, p.pmid, g.stable_id, lgd.id
        from publication_families pf
        left join publication p on p.id = pf.publication_id
        left join attrib a on a.id = pf.consanguinity_id
        left join lgd_publication lp on lp.publication_id = p.id
        left join locus_genotype_disease lgd on lgd.id = lp.lgd_id
        left join g2p_stableid g on g.id = lgd.stable_id
          """

    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        pmid = row[6]

        if pmid not in families_data:
            families_data[pmid] = [{ 
                "families": row[0],
                "affected_individuals": row[1],
                "ancestries": row[2].replace("\t", ", "),
                "is_deleted": row[3],
                "consanguinity": row[4],
                "publication_id": row[5],
                "g2p_id": row[7],
                "lgd_id": row[8]
            }]
        else:
            families_data[pmid].append({ 
                "families": row[0],
                "affected_individuals": row[1],
                "ancestries": row[2].replace("\t", ", "),
                "is_deleted": row[3],
                "consanguinity": row[4],
                "publication_id": row[5],
                "g2p_id": row[7],
                "lgd_id": row[8]
            })

    return families_data


def generate_file(families_data):
    output_file = "publication_families_release_v7.txt"

    with open(output_file, "w") as wr:
        wr.write("g2p id\tlgd_id\tpmid\tpublication_id\tnumber of families\taffected individuals\tancestries\tconsanguinity\tis_deleted\n")
    
        for pmid in families_data:
            for data in families_data[pmid]:
                wr.write(
                    str(data["g2p_id"])+"\t"+
                    str(data["lgd_id"])+"\t"+
                    str(pmid)+"\t"+
                    str(data["publication_id"])+"\t"+
                    str(data["families"])+"\t"+
                    str(data["affected_individuals"])+"\t"+
                    data["ancestries"]+"\t"+
                    data["consanguinity"]+"\t"+
                    str(data["is_deleted"])+"\n"
                )


def main():
    parser = argparse.ArgumentParser(description="Dump the table publication_families to a file")
    parser.add_argument("--config", required=True, help="Config file with details to the G2P database")

    args = parser.parse_args()

    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config['g2p_database']['g2p_host']
    db_port = int(config['g2p_database']['g2p_port'])
    db_name = config['g2p_database']['g2p_database']
    user = config['g2p_database']['g2p_user']
    password = config['g2p_database']['g2p_password']

    families_data = fetch_data_from_db(db_host, db_port, db_name, user, password)
    generate_file(families_data)

if __name__ == '__main__':
    main()
