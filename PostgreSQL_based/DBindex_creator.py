import psycopg2

def main():

    query='CREATE INDEX idx_chr_pos ON read(rname,pos);'

    conn = psycopg2.connect(database="pysamdb", user='postgres', password='pysam8192', host='127.0.0.1', port= '5432')
    conn.autocommit = True
    cursor = conn.cursor()
    cursor.execute(query)
    conn.close()


if __name__ == "__main__":
    main()