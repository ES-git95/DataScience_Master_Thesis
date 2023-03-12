#import psycopg2
import subprocess

def main():

    query='CREATE INDEX idx_chr_pos ON read(rname,pos);'
    status=subprocess.run(["PGPASSWORD=pysam8192 psql -U postgres -h 127.0.0.1  -d pysamdb -c \""+query+"\""],shell=True, capture_output=True)
    print(status)


if __name__ == "__main__":
    main()