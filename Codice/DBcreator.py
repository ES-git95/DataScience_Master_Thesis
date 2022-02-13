
from DBfunctions import *


creation_dict={'SamFile':
'''
CREATE TABLE IF NOT EXISTS SamFile (
    id       INTEGER PRIMARY KEY AUTOINCREMENT,
    filename STRING,
    index_filename STRING,
    description STRING
);
''',
'RS_Dictionary':
'''
CREATE TABLE IF NOT EXISTS RS_Dictionary (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    name       STRING,
    length     INTEGER,
    samfile_id INTEGER REFERENCES SamFile (id) 
);
''',
'Metadata':
'''
CREATE TABLE IF NOT EXISTS Metadata (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    main_tag   STRING,
    tag        STRING,
    value      STRING,
    samfile_id INTEGER REFERENCES SamFile (id) 
);
''',
'ReadGroup':
'''
CREATE TABLE IF NOT EXISTS ReadGroup (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    group_id   INTEGER,
    sample     STRING,
    platform   STRING,
    samfile_id INTEGER REFERENCES SamFile (id) 
);
''',
'Read':
'''
CREATE TABLE IF NOT EXISTS Read (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    qname      STRING,
    flag       INTEGER,
    rname      STRING  REFERENCES RS_Dictionary (id),
    pos        INTEGER,
    mapq       INTEGER,
    cigar      STRING,
    rnext      STRING,
    pnext      INTEGER,
    tlen       INTEGER,
    seq        STRING,
    qual       STRING,
    group_id   INTEGER REFERENCES ReadGroup (id),
    samfile_id INTEGER REFERENCES SamFile (id) 
);
''',
'Tags':
'''
CREATE TABLE IF NOT EXISTS Tags (
    id      INTEGER PRIMARY KEY AUTOINCREMENT,
    tags_values STRING,
    read_id INTEGER REFERENCES Read (id) 
);
'''}



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