#!/usr/bin/env python3

import sys
import re
import argparse
import MySQLdb
import configparser
from pathlib import Path

"""
    Script to dump the table publication_comment to a file.
    The file is going to be used to populate the new schema.
    Update done in release v8.
"""


def fetch_data_from_db(db_host, db_port, db_name, user, password):
    """
    Fetch the contents of the table publication_comment.
    """
    comments_data = {}

    db = MySQLdb.connect(
        host=db_host, port=db_port, user=user, passwd=password, db=db_name
    )
    cursor = db.cursor()

    sql = """
        select pc.comment, pc.is_public, pc.is_deleted, pc.date, p.id, p.pmid, g.stable_id,
        lgd.id, pc.user_id, u.username
        from publication_comment pc
        left join publication p on p.id = pc.publication_id
        left join lgd_publication lp on lp.publication_id = p.id
        left join locus_genotype_disease lgd on lgd.id = lp.lgd_id
        left join g2p_stableid g on g.id = lgd.stable_id
        left join user u on u.id = pc.user_id
        """

    cursor.execute(sql)
    data = cursor.fetchall()
    for row in data:
        pmid = row[5]

        if pmid not in comments_data:
            comments_data[pmid] = [
                {
                    "comment": row[0],
                    "is_public": row[1],
                    "is_deleted": row[2],
                    "date": row[3],
                    "publication_id": row[4],
                    "g2p_id": row[6],
                    "lgd_id": row[7],
                    "user_id": row[8],
                    "username": row[9],
                }
            ]
        else:
            comments_data[pmid].append(
                {
                    "comment": row[0],
                    "is_public": row[1],
                    "is_deleted": row[2],
                    "date": row[3],
                    "publication_id": row[4],
                    "g2p_id": row[6],
                    "lgd_id": row[7],
                    "user_id": row[8],
                    "username": row[9],
                }
            )

    return comments_data


def generate_file(comments_data):
    output_file = "publication_comments_release_v8.txt"

    with open(output_file, "w") as wr:
        wr.write(
            "g2p id\tlgd_id\tpmid\tpublication_id\tcomment\tuser_id\tusername\tdate\tis_deleted\n"
        )

        for pmid in comments_data:
            for data in comments_data[pmid]:
                wr.write(
                    str(data["g2p_id"])
                    + "\t"
                    + str(data["lgd_id"])
                    + "\t"
                    + str(pmid)
                    + "\t"
                    + str(data["publication_id"])
                    + "\t"
                    + data["comment"]
                    + "\t"
                    + str(data["user_id"])
                    + "\t"
                    + data["username"]
                    + "\t"
                    + str(data["date"])
                    + "\t"
                    + str(data["is_deleted"])
                    + "\n"
                )


def main():
    parser = argparse.ArgumentParser(
        description="Dump the table publication_comment to a file"
    )
    parser.add_argument(
        "--config", required=True, help="Config file with details to the G2P database"
    )

    args = parser.parse_args()

    config_file = args.config

    # Load the config file
    config = configparser.ConfigParser()
    config.read(config_file)

    db_host = config["g2p_database"]["g2p_host"]
    db_port = int(config["g2p_database"]["g2p_port"])
    db_name = config["g2p_database"]["g2p_database"]
    user = config["g2p_database"]["g2p_user"]
    password = config["g2p_database"]["g2p_password"]

    comments_data = fetch_data_from_db(db_host, db_port, db_name, user, password)
    generate_file(comments_data)


if __name__ == "__main__":
    main()
