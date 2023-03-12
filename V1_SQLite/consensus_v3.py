import pandas as pd
import sqlite3
import re
import sys, getopt
import time

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"r:p",["region=","perfomance="])
    except getopt.GetoptError: 
        print("mpileup.py -r <chr region> (ex -r chr:1-100)")
        sys.exit() 

    for opt, arg in opts: 
        performance=False
        if opt in ('-r','-region'):
            arg_list=arg.split(':')
            Chr_name=arg_list[0]
            Starting_pos=arg_list[1].split('-')[0]
            Ending_pos=arg_list[1].split('-')[1]
        if opt in ('-p','-performance'):
            performance=True

    begin=time.time()
    con = sqlite3.connect('PysamDB.db')
    cur = con.cursor()
    v_range=str(cur.execute("select length(seq) from Read limit 1;").fetchall()[0][0])

    query="select mapq, flag, seq,pos,cigar,qual from Read where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-"+v_range+"+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)

    db_stop=time.time()
    db_duration=(db_stop-begin)

    pos_list=list(sorted(set(df.pos)))

    df_before_seq=pd.DataFrame()
    Starting_pos=int(Starting_pos)
    Ending_pos=int(Ending_pos)
    freq_base=''
    for position in range(Starting_pos,Ending_pos+1):

        if not not pos_list:
            if position == pos_list[0]:
                df_idx=pos_list.pop(0)
                df_sub=df[(df.pos==df_idx)]
                df_seq_tolist=pd.DataFrame(df_sub.seq.apply(list).tolist()).set_index([df_sub.pos])
                df_qual_tolist=pd.DataFrame(df_sub.qual.apply(list).tolist()).set_index([df_sub.pos])
                df_seq_final=df_seq_tolist 

                df_seq_t=df_seq_final.T.agg(lambda x: ''.join(x).replace(' ',''), axis=1)
                new_idxs=list(range(position,position+len(df_seq_t)))
                df_seq_t.index=new_idxs

                if df_before_seq.empty==False:
                    df_before_seq=pd.concat([df_before_seq, df_seq_t], axis=1, join="outer").fillna('').agg(''.join, axis=1).sort_index()
                else:
                    df_before_seq=df_seq_t

        if df_before_seq.empty==False:
            seq_final=list(df_before_seq.pop(position))
            if seq_final != []:
                freq_base=freq_base+max(set(seq_final), key = seq_final.count)
                if freq_base == '*':
                    freq_base='N'
            else:
                freq_base='N'
        else:
            freq_base=freq_base+'N'

    stop=time.time()
    duration=(stop-begin)

    if performance==False:
        print(Chr_name+'\n'+freq_base)
    else:
        print(db_duration,duration)

if __name__ == "__main__":
    main(sys.argv[1:])