from sqlite3 import connect

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


def GetMaxIdFromTable(samfile_name,table_name,operation,connection):

    if operation=='MAX':
        SelectQuery="SELECT max(id) FROM '"+table_name+"'"

    if operation=='SAMFILE_ID':
        SelectQuery="SELECT max(id) FROM '"+table_name+"' WHERE filename='"+samfile_name+"'"
    
    cursor=connection.cursor()

    try:
        cursor.execute(SelectQuery)
        id=cursor.fetchone()[0]
        if id is None:
            id=0
    except:
        id=0

    cursor.close()

    return id