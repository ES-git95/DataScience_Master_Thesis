import pandas as pd
import sqlite3
import sys, getopt
import numpy as np
import itertools as itr
from collections import Counter

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"r:",["region="])
    except getopt.GetoptError:
        print("coverage.py -r <chr region> (ex -r chr:1-100)")
        sys.exit()

    for opt, arg in opts:
        if opt in ('-r','-region'):
            arg_list=arg.split(':')
            Chr_name=arg_list[0]
            Starting_pos=arg_list[1].split('-')[0]
            Ending_pos=arg_list[1].split('-')[1]

    con = sqlite3.connect('PysamDB.db')
    query="select Read.id id,mapq,qual,qual_insdel,pos,seq,seq_insdel from Read left join InsDel on Read.id=InsDel.read_id where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-50+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)  

    idx_list=df[df['seq_insdel'].notnull()].index.tolist()
    df.update(pd.DataFrame({'seq':list(df.loc[idx_list,'seq_insdel']),'qual':list(df.loc[idx_list,'qual_insdel'])},index=idx_list))

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
                df_bool=~(df_qual_tolist=='-')
                df_seq_final=df_seq_tolist.where(df_bool,' ') #--> rivedere questione dimensionalità

                df_seq_t=df_seq_final.T.agg(lambda x: ''.join(x).replace(' ',''), axis=1)
                new_idxs=list(range(position,position+len(df_seq_t)))
                df_seq_t.index=new_idxs

                if df_before_seq.empty==False:
                    df_before_seq=pd.concat([df_before_seq, df_seq_t], axis=1, join="outer").fillna('').agg(''.join, axis=1).sort_index()
                else:
                    df_before_seq=df_seq_t

        if df_before_seq.empty==False:
            seq_final=list(df_before_seq.pop(position))
            freq_base=freq_base+max(set(seq_final), key = seq_final.count)
        else:
            freq_base=freq_base+'N'

            
    print(Chr_name+'\n'+freq_base)


if __name__ == "__main__":
    main(sys.argv[1:])