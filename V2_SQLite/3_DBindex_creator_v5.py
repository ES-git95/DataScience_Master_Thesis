
from sqlite3 import connect

def main():

    connection=connect('PysamDB.db')
    db_cursor=connection.cursor()
    
    query1='CREATE INDEX Read_idx ON Read(rname, pos);'
    query2='CREATE INDEX InsDel_idx ON InsDel(read_id);'

    db_cursor.execute(query1)
    db_cursor.execute(query2)

    connection.commit()
    db_cursor.close()

    connection.close()


if __name__ == "__main__":
    main()