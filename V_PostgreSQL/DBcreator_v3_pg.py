#import psycopg2
import subprocess


creation_dict={'ReferenceSequence':
'CREATE TABLE IF NOT EXISTS ReferenceSequence (    name       VARCHAR PRIMARY KEY,    length     INTEGER);',
'Metadata':
'CREATE TABLE IF NOT EXISTS Metadata (    id              INTEGER PRIMARY KEY,    main_tag        VARCHAR,    tag_key         VARCHAR,tag_value       VARCHAR,sub_tag_key     VARCHAR,    sub_tag_value   VARCHAR);',
'ReadGroup':
'CREATE TABLE IF NOT EXISTS ReadGroup (    group_id   INTEGER PRIMARY KEY,    sample     VARCHAR,    platform   VARCHAR);',
'Read':
'CREATE TABLE IF NOT EXISTS Read (    id         INTEGER PRIMARY KEY,    qname      VARCHAR,    flag       INTEGER,    rname      VARCHAR  ,    pos        INTEGER,    mapq       INTEGER,    cigar      VARCHAR,    rnext      VARCHAR,    pnext      INTEGER,    tlen       INTEGER,    seq        VARCHAR,    qual       VARCHAR,    RG         VARCHAR  );',
'Tags':
'CREATE TABLE IF NOT EXISTS Tags (    id              VARCHAR PRIMARY KEY,    tags_values     VARCHAR,    read_id         INTEGER );'}


def main():

    for table,query in creation_dict.items():

        status=subprocess.run(["PGPASSWORD=pysam8192 psql -U postgres -h 127.0.0.1  -d pysamdb -c \""+query+"\""],shell=True, capture_output=True)
        print(status)
        print('Table '+table+' successfully created')


if __name__ == "__main__":
    main()