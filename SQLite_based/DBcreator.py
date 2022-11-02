
from sqlite3 import connect

creation_dict={'ReferenceSequence':
'''
CREATE TABLE IF NOT EXISTS ReferenceSequence (
    name       STRING PRIMARY KEY,
    length     INTEGER
);
''',
'Metadata':
'''
CREATE TABLE IF NOT EXISTS Metadata (
    id              INTEGER PRIMARY KEY,
    main_tag        STRING,
    tag_key         STRING,
    tag_value       STRING,
    sub_tag_key     STRING,
    sub_tag_value   STRING
);
''',
'ReadGroup':
'''
CREATE TABLE IF NOT EXISTS ReadGroup (
    group_id   INTEGER PRIMARY KEY,
    sample     STRING,
    platform   STRING
);
''',
'Read':
'''
CREATE TABLE IF NOT EXISTS Read (
    id         INTEGER PRIMARY KEY,
    qname      STRING,
    flag       INTEGER,
    rname      STRING  REFERENCES ReferenceSequence (name),
    pos        INTEGER,
    mapq       INTEGER,
    cigar      STRING,
    rnext      STRING,
    pnext      INTEGER,
    tlen       INTEGER,
    seq        STRING,
    qual       STRING,
    RG         STRING  REFERENCES ReadGroup (group_id)
);
''',
'Tags':
'''
CREATE TABLE IF NOT EXISTS Tags (
    id              STRING PRIMARY KEY,
    tags_values     STRING,
    read_id         INTEGER REFERENCES Read (id) 
);
'''}


def main():

    connection=connect('PysamDB.db')
    db_cursor=connection.cursor()
    for table,query in creation_dict.items():
        db_cursor.execute(query)

    connection.commit()
    db_cursor.close()

    connection.close()


if __name__ == "__main__":
    main()