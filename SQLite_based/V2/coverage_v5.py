import pandas as pd
import sqlite3
import sys, getopt
import numpy as np
import itertools as itr
from collections import Counter
import re
import warnings
import copy
import time

warnings.filterwarnings("ignore", category=FutureWarning) 


def main(argv):

    try:
        opts, args = getopt.getopt(argv,"r:p",["region=","perfomance="])
    except getopt.GetoptError:
        print("coverage.py -r <chr region> (ex -r chr:1-100)")
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
    con = sqlite3.connect('PysamDB.db')#V5/
    cur = con.cursor()
    v_range=str(cur.execute("select length(seq) from Read limit 1;").fetchall()[0][0])

    query="select Read.id id,mapq,qual,qual_insdel,pos,seq,seq_insdel,pos_ixd,attachment from Read left join InsDel on Read.id=InsDel.read_id where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-"+v_range+"+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)  

    db_stop=time.time()
    db_duration=(db_stop-begin)

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
        if position >= int(Starting_pos):
            for counter,base_qual in enumerate(base_quals[1]['qual']):
                counter_tot+=1
                if base_qual != '-':
                    #counter_tot+=1
                    #if base_qual != '<':
                    qual_vals.append(ord(base_qual)-33)
                        
                else:
                    qual_exclude.append(str(position+counter))
    mean_quals=sum(qual_vals)/counter_tot

    #3. iterazione su ciascun read.id+49> unione lista> Counter --> covbases
    pos_list=[list(range(i,i+int(v_range))) for i in read_pos]
    pos_list_join=list(filter(lambda x : x >= int(Starting_pos) and x <= int(Ending_pos),itr.chain.from_iterable(pos_list)))
    pos_list_join_1=copy.copy(pos_list_join)

    sub_df=df[['pos_ixd','attachment']].dropna().reset_index(drop=True)
    for row in sub_df[sub_df['attachment'].str.contains(",")].iterrows():
        for idx,pos in enumerate(row[1]['pos_ixd'].split(',')):
            attach=row[1]['attachment'].split(',')
            sub_df = sub_df.append({'pos_ixd':pos, 'attachment':attach[idx]}, ignore_index=True)
        sub_df=sub_df.drop(row[0]).reset_index(drop=True)


    sub_df_insertions=pd.concat([sub_df[sub_df['attachment'].str.contains("\+")]['attachment'].str.extract('(\d+)'),sub_df[sub_df['attachment'].str.contains("\+")]['pos_ixd']],axis=1)
    sub_df_ins=[]
    for i in sub_df_insertions.iterrows():
        for element in range(int(i[1]['pos_ixd'])+1,int(i[1]['pos_ixd'])+int(i[1][0])+1):
            sub_df_ins.append(str(element))
    
    #sub_df_del=list(sub_df[sub_df['attachment'].str.contains("-")]['pos_ixd'])+qual_exclude
    sub_df_deletions=pd.concat([sub_df[sub_df['attachment'].str.contains("-")]['attachment'].str.extract('(\d+)'),sub_df[sub_df['attachment'].str.contains("-")]['pos_ixd']],axis=1)
    sub_df_del=[]
    for i in sub_df_deletions.iterrows():
        for element in range(int(i[1]['pos_ixd'])+1,int(i[1]['pos_ixd'])+int(i[1][0])+1):
            sub_df_del.append(str(element))
    sub_df_del=sub_df_del#+qual_exclude
    sub_df_del_qual=sub_df_del+qual_exclude
    #sub_df_ins=list(sub_df[sub_df['attachment'].str.contains("\+")]['pos_ixd'])
    #final_count=0
    if sub_df_del != []:
        for x in sub_df_del:
            #if ',' in x:
            #    pos_list_rm=pos_to_exclude.pop(counter)
            #    multi_pos=pos_list_rm.split(',')
            #    for pos in multi_pos:
            #        pos_list_join.remove(int(pos))
            #else:
            pos_list_join.remove(int(x))
        for x in sub_df_del_qual:
            pos_list_join_1.remove(int(x))
            

    
    if sub_df_ins != []:
        for x in sub_df_ins:
            #if ',' in x:
            #    pos_list_rm=pos_to_add.pop(counter)
            #    multi_pos=pos_list_rm.split(',')
            #    for pos in multi_pos:
            #        pos_list_join.append(int(pos))
            #else:
            pos_list_join.append(int(x))
            pos_list_join_1.append(int(x))

    num_bases=len(set(pos_list_join))    

    #coverage
    coverage=num_bases/(int(Ending_pos)-int(Starting_pos)+1)*100

    #meandepth
    pos_list_counter=Counter(pos_list_join_1)
    mean_bases_pos=(sum(pos_list_counter.values()))/(int(Ending_pos)-int(Starting_pos)+1) 
    #--> questo si può aggiustare con l'estrazione degli indici relativi alle cancellazioni e rimozione dalla prima lista prima del counter

    stop=time.time()
    duration=(stop-begin)

    if performance==False:    
        print(Chr_name,Starting_pos,Ending_pos,read_count,num_bases,coverage,mean_bases_pos,mean_quals,read_qual_avg)
    else:
        print(db_duration,duration)


if __name__ == "__main__":
    main(sys.argv[1:])