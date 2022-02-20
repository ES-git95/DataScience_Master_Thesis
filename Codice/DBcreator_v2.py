
from sqlite3 import connect

creation_dict={'ReferenceSequence':
'''
CREATE TABLE IF NOT EXISTS ReferenceSequence (
    id         INTEGER PRIMARY KEY,
    name       STRING,
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
    rname      STRING  REFERENCES ReferenceSequence (id),
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
    id              INTEGER PRIMARY KEY,
    tags_values     STRING,
    read_id         INTEGER REFERENCES Read (id) 
);
'''}


def DBConnection(DBName,operation,connection=None):
    if operation=='OPEN':
        connection=connect(DBName)
        print('Connected to '+DBName+'\n')
    if operation=='CLOSE':
        try:
            connection.close()
            print("The connection to "+DBName+" has close \n")
        except:
            connection=None
            print("ERROR: The connection to "+ DBName+" has not been open yet, or the 'connection' parameter is missing \n")
    return connection


def main():

    n='n' #this is usefur if you input n without '' from the terminal
    y='y'

    while True:
        choice=str(input('This program will create the table of PysamDB.db database schema in SQLite, do you want to continue? [y/n] : ')).lower()
        if choice not in ['y','n']:
            print("Please, enter 'y' to continue, 'n' to exit the program. Other inputs are not valid. \n")
        else:
            break

    print('\n')
    if choice=='y':

        connection=DBConnection('PysamDB.db','OPEN')
        db_cursor=connection.cursor()
        for table,query in creation_dict.items():
            db_cursor.execute(query)
            print('Table '+table+' successfully created')

        connection.commit()
        db_cursor.close()
        print('\n') 
        DBConnection('PysamDB.db','CLOSE',connection)
        print('Congratulations! Database PysamDB.db schema is now populated with your tables. \n')

    else:
        print('PysamDB.db database will not be created. See you! \n')


if __name__ == "__main__":
    main()