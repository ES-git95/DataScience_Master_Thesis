
from sqlite3 import connect

def main():

    connection=connect('PysamDB.db')
    db_cursor=connection.cursor()
    
    query='CREATE INDEX Read_idx ON Read(rname, pos);'

    db_cursor.execute(query)

    connection.commit()
    db_cursor.close()

    connection.close()


if __name__ == "__main__":
    main()
