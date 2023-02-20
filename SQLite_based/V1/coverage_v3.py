import pandas as pd
import sqlite3
import re
import sys, getopt
import numpy as np
import itertools as itr
from collections import Counter
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

    #df['new_pos']=df['pos']

    for row in df.iterrows():
        cigar=row[1]['cigar'] 
        digits=list(map(int,re.findall('\d+',cigar)))
        letters=re.findall("[a-zA-Z]+", cigar)
        seq_insdel=None

        if 'D' in letters:
            cigar_extended=''.join([a*b for a,b in zip(digits,letters)])
            cigar_ext_idx=enumerate(cigar_extended)
            del_idxs= [idx for idx, op in cigar_ext_idx if op == 'D']
            counter=0

            for i,del_idx in enumerate(del_idxs):
                seq_insdel=df.loc[row[0],'seq'][:del_idx+counter] + '*' + df.loc[row[0],'seq'][del_idx+counter:]
                qual_insdel=df.loc[row[0],'qual'][:del_idx+counter+1] + df.loc[row[0],'qual'][del_idx+counter:]
                counter+=1

            #df.loc[row[0],'new_pos']=df.loc[row[0],'new_pos']+len(del_idxs)


        if 'I' in letters:
            cigar_extended=''.join([a*b for a,b in zip(digits,letters)])
            cigar_ext_idx=enumerate(cigar_extended)
            ins_idxs= [idx for idx, op in cigar_ext_idx if op == 'I']
            counter=0
            if seq_insdel is None:
                seq_insdel=df.loc[row[0],'seq']
                qual_insdel=df.loc[row[0],'qual']
            else:
                pass

            for i,ins_idx in enumerate(ins_idxs):
                df.loc[row[0],'seq']=seq_insdel[:ins_idx+counter] + df.loc[row[0],'seq'][ins_idx+1+counter:] 
                df.loc[row[0],'qual']=qual_insdel[:ins_idx+counter] + df.loc[row[0],'qual'][ins_idx+1+counter:] 
                counter-=1

        if ('I' in letters or 'D' in letters):
            df.loc[row[0],'seq']=seq_insdel
            df.loc[row[0],'qual']=qual_insdel

    read_pos=df.pos

    #1. estrazione dei read.id --> numreads
    read_count=len(read_pos)

    #2. estrazione read_qual e media --> meanmapq
    read_qual_avg=np.mean(df['mapq'])

    #4. estrazione tutte qual, conversione ASCII > media --> meanbaseq
    qual_vals=[]
    qual_exclude=[]
    seq_exclude=[]
    counter_tot=0
    for base_quals in df[['pos','seq','qual']].iterrows():
        position=int(base_quals[1]['pos'])
        for counter in range(0,len(base_quals[1]['qual'])):
            base_qual=base_quals[1]['qual'][counter]
            base_seq=base_quals[1]['seq'][counter]
            if base_qual != '-':
                counter_tot+=1
                if base_qual != '<':
                    qual_vals.append(ord(base_qual)-33)
                else:
                    qual_exclude.append(str(position+counter))
            if base_seq == '*':
               seq_exclude.append(position+counter) 

    mean_quals=sum(qual_vals)/counter_tot 


    #3. iterazione su ciascun read.id+49> unione lista> Counter --> covbases
    pos_list=[list(range(i,i+int(v_range))) for i in read_pos]
    pos_list_join=list(filter(lambda x : x >= int(Starting_pos) and x <= int(Ending_pos),itr.chain.from_iterable(pos_list)))


    pos_to_exclude=seq_exclude+qual_exclude

    sub_df_del=pos_to_exclude
    for counter,x in enumerate(sub_df_del):
        single_pos=pos_to_exclude.pop(counter)
        pos_list_join.remove(int(single_pos))

    num_bases=len(set(pos_list_join))    

    #coverage
    coverage=num_bases/(int(Ending_pos)-int(Starting_pos))*100

    pos_list_counter=Counter(pos_list_join)
    mean_bases_pos=sum(pos_list_counter.values())/(int(Ending_pos)-int(Starting_pos)) #--> meandepth
    #--> questo si può aggiustare con l'estrazione degli indici relativi alle cancellazioni e rimozione dalla prima lista prima del counter

    stop=time.time()
    duration=(stop-begin)

    if performance==False:
        print(Chr_name,Starting_pos,Ending_pos,read_count,num_bases,coverage,mean_bases_pos,mean_quals,read_qual_avg)
    else:
        print(db_duration,duration)

if __name__ == "__main__":
    main(sys.argv[1:])