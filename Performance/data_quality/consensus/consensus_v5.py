import pandas as pd
import sqlite3
import sys, getopt
import time

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"r:s:p",["region=","show=","perfomance="])
    except getopt.GetoptError:
        print("consensus.py -r <chr region> -show_dels <yes,no>")
        sys.exit()

    for opt, arg in opts:
        show_ins=False
        performance=False
        if opt in ('-r','-region'):
            arg_list=arg.split(':')
            Chr_name=arg_list[0]
            Starting_pos=arg_list[1].split('-')[0]
            Ending_pos=arg_list[1].split('-')[1]
        if opt in ('-s',"-show"):
            if arg=='no':
                show_ins=False
            if arg=='yes':
                show_ins=True
        if opt in ('-p','-performance'):
            performance=True

    begin=time.time()
    con = sqlite3.connect('PysamDB.db')
    cur = con.cursor()
    v_range=str(cur.execute("select length(seq) from Read limit 1;").fetchall()[0][0])

    query="select Read.id id,mapq,qual,qual_insdel,pos,seq,seq_insdel from Read left join InsDel on Read.id=InsDel.read_id where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-"+v_range+"+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)  

    db_stop=time.time()
    db_duration=(db_stop-begin)

    if show_ins==True: 
        idx_list=df[df['seq_insdel'].notnull()].index.tolist()
        df.update(pd.DataFrame({'seq':list(df.loc[idx_list,'seq_insdel']),'qual':list(df.loc[idx_list,'qual_insdel'])},index=idx_list))

    pos_list=list(sorted(set(df.pos)))

    df_before_seq=pd.DataFrame()
    Starting_pos=int(Starting_pos)
    Ending_pos=int(Ending_pos)
    freq_base=''
    long_base=''
    for position in range(Starting_pos-int(v_range)+1,Ending_pos+1):


        base=''
        if not not pos_list:
            if position == pos_list[0]:
                df_idx=pos_list.pop(0)
                df_sub=df[(df.pos==df_idx)]
                df_seq_tolist=pd.DataFrame(df_sub.seq.apply(list).tolist()).set_index([df_sub.pos])
                df_qual_tolist=pd.DataFrame(df_sub.qual.apply(list).tolist()).set_index([df_sub.pos])
                df_bool=~(df_qual_tolist=='-')
                df_seq_final=df_seq_tolist.where(df_bool,' ')

                df_seq_t=df_seq_final.T.agg(lambda x: ''.join(x).replace(' ',''), axis=1)
                new_idxs=list(range(position,position+len(df_seq_t)))
                df_seq_t.index=new_idxs

                if df_before_seq.empty==False:
                    df_before_seq=pd.concat([df_before_seq, df_seq_t], axis=1, join="outer").fillna('').agg(''.join, axis=1).sort_index()
                else:
                    df_before_seq=df_seq_t

        if (df_before_seq.empty==False) and (position>=Starting_pos) and (position in df_before_seq.index):
            seq_final=list(df_before_seq.pop(position))
            if seq_final != []:
                base=max(set(seq_final), key = seq_final.count) 
                if (base == '*') and (show_ins==False):
                    base='N' 
                if (base == '*') and (show_ins==True):
                    base='*' 
                freq_base=freq_base+base
            if seq_final == [] and (position>=Starting_pos):
                base='N'
                freq_base=freq_base+base

        if (df_before_seq.empty==True) and (position>=Starting_pos):
            base='N'
            freq_base=freq_base+base
        
        if base!='':
            long_base=long_base+base
            if len(long_base)==50 and performance==False: 
                print(long_base)
                long_base=''

    stop=time.time()
    duration=(stop-begin)

    if performance==True:       
        print(db_duration,duration)


if __name__ == "__main__":
    main(sys.argv[1:])