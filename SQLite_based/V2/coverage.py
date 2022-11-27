import pandas as pd
import sqlite3
import sys, getopt
import numpy as np
import itertools as itr
from collections import Counter
import re

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
    query="select Read.id id,mapq,qual,qual_insdel,pos,seq,seq_insdel,pos_ixd,attachment from Read left join InsDel on Read.id=InsDel.read_id where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-50+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)  

    idx_list=df[df['seq_insdel'].notnull()].index.tolist()
    df.update(pd.DataFrame({'seq':list(df.loc[idx_list,'seq_insdel']),'qual':list(df.loc[idx_list,'qual_insdel'])},index=idx_list))

    read_pos=df.pos

    #1. estrazione dei read.id --> numreads
    read_count=len(read_pos)

    #2. estrazione read_qual e media --> meanmapq
    read_qual_avg=np.mean(df['mapq'])

    #4. estrazione tutte qual, conversione ASCII > media --> meanbaseq
    qual_vals=[]
    qual_exclude=[]
    counter_tot=0
    for base_quals in df[['pos','qual']].iterrows():
        position=int(base_quals[1]['pos'])
        for counter,base_qual in enumerate(base_quals[1]['qual']):
            if base_qual != '-':
                counter_tot+=1
                if base_qual != '<':
                    qual_vals.append(ord(base_qual)-33)
            else:
                qual_exclude.append(str(position+counter))
    mean_quals=sum(qual_vals)/counter_tot

    #3. iterazione su ciascun read.id+49> unione lista> Counter --> covbases
    pos_list=[list(range(i,i+50)) for i in read_pos]
    pos_list_join=list(filter(lambda x : x >= int(Starting_pos) and x <= int(Ending_pos),itr.chain.from_iterable(pos_list)))


    sub_df=df[['pos_ixd','attachment']].dropna()
    sub_df_del=list(sub_df[sub_df['attachment'].str.contains("-")]['pos_ixd'])+qual_exclude

    pos_to_exclude=sub_df_del
    for counter,x in enumerate(sub_df_del):
        if ',' in x:
            pos_list_rm=pos_to_exclude.pop(counter)
            multi_pos=pos_list_rm.split(',')
            for pos in multi_pos:
                pos_list_join.remove(int(pos))
        else:
            single_pos=pos_to_exclude.pop(counter)
            pos_list_join.remove(int(single_pos))

    num_bases=len(set(pos_list_join))    

    #coverage
    coverage=num_bases/(int(Ending_pos)-int(Starting_pos))*100

    pos_list_counter=Counter(pos_list_join)
    mean_bases_pos=sum(pos_list_counter.values())/(int(Ending_pos)-int(Starting_pos)) #--> meandepth
    #--> questo si può aggiustare con l'estrazione degli indici relativi alle cancellazioni e rimozione dalla prima lista prima del counter

    print(Chr_name,Starting_pos,Ending_pos,read_count,num_bases,coverage,mean_bases_pos,mean_quals,read_qual_avg)


if __name__ == "__main__":
    main(sys.argv[1:])