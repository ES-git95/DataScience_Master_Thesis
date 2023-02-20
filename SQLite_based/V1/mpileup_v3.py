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

    df['new_pos']=df['pos']

    del_dic={}
    ins_dic={}

    for row in df.iterrows():
        cigar=row[1]['cigar']
        digits=list(map(int,re.findall('\d+',cigar)))
        letters=re.findall("[a-zA-Z]+", cigar)

        if 'D' in letters:
            cigar_extended=''.join([a*b for a,b in zip(digits,letters)])
            cigar_ext_idx=enumerate(cigar_extended)
            del_idxs= [idx for idx, op in cigar_ext_idx if op == 'D']
            del_num=[a for a,b in zip(digits,letters) if b=='D']
            counter=0

            for i,del_idx in enumerate(del_idxs):
                df.loc[row[0],'seq']=df.loc[row[0],'seq'][:del_idx+counter] + '*' + df.loc[row[0],'seq'][del_idx+counter:]
                df.loc[row[0],'qual']=df.loc[row[0],'qual'][:del_idx+counter+1] + df.loc[row[0],'qual'][del_idx+counter:]
                ins_count=cigar_extended[:del_idx+counter].count('I')
                dic_index=row[1]['pos']+del_idx+counter-1-ins_count
                if (i==0) | (del_idxs[i-1] != del_idxs[i]-1):
                    if  (dic_index not in del_dic):
                        if (row[1]['flag']==16):
                            del_dic[dic_index]=['-'+str(del_num.pop(0))+'n']
                        else:
                            del_dic[dic_index]=['-'+str(del_num.pop(0))+'N'] 
                    else:
                        if (row[1]['flag']==16):
                            del_dic[dic_index].append('-'+str(del_num.pop(0))+'n')
                        else:
                            del_dic[dic_index].append('-'+str(del_num.pop(0))+'N') 
                    counter+=1
                else:
                    continue 

            df.loc[row[0],'new_pos']=df.loc[row[0],'new_pos']+len(del_idxs)


        if 'I' in letters:
            cigar_extended=''.join([a*b for a,b in zip(digits,letters)])
            cigar_ext_idx=enumerate(cigar_extended)
            ins_idxs= [idx for idx, op in cigar_ext_idx if op == 'I']
            ins_num=[a for a,b in zip(digits,letters) if b=='I']
            counter=0

            for i,ins_idx in enumerate(ins_idxs):
                df.loc[row[0],'seq']=df.loc[row[0],'seq'][:ins_idx+counter] + df.loc[row[0],'seq'][ins_idx+1+counter:] 
                df.loc[row[0],'qual']=df.loc[row[0],'qual'][:ins_idx+counter] + df.loc[row[0],'qual'][ins_idx+1+counter:] 
                dic_index=row[1]['pos']+ins_idx+counter-1
                if (i==0) | (ins_idxs[i-1] != ins_idxs[i]-1):
                    num=ins_num.pop(0)
                    if  (dic_index not in ins_dic):
                        if (row[1]['flag']==16):
                            ins_dic[dic_index]=['+'+str(num)+(row[1]['seq'][ins_idx+counter:ins_idx+counter+num]).lower()]
                        else:
                            ins_dic[dic_index]=['+'+str(num)+(row[1]['seq'][ins_idx+counter:ins_idx+counter+num]).upper()] 
                    else:
                        if (row[1]['flag']==16):
                            ins_dic[dic_index].append('+'+str(num)+(row[1]['seq'][ins_idx+counter:ins_idx+counter+num]).lower())
                        else:
                            ins_dic[dic_index].append('+'+str(num)+(row[1]['seq'][ins_idx+counter:ins_idx+counter+num]).upper()) 
                    counter-=1
                else:
                    continue   


    idxs=list(df.pos)
    idxs_new=list(df.new_pos)

    for position in range(int(Starting_pos),int(Ending_pos)+1): 
        ixds_filtered=sorted(set([a for a,b in zip(idxs,idxs_new) if (b>=position-int(v_range)+1 and a<=position)]))
        seq_final=[]
        qual_final=[]
        for idx in ixds_filtered:
            pos=position-idx
            df_seq=df[df['pos']==idx]
            if df_seq.empty==True:
                continue
            for row in df_seq.iterrows():
                initial=''
                if pos==0:
                    initial='^'+chr(row[1]['mapq']+33)
                last=''
                if pos==len(row[1]['seq'])-1:
                    last='$'

                attach=''
                try:
                    base_qual=row[1]['qual'][pos]
                    base=row[1]['seq'][pos]+attach+last
                    if base_qual=='-':
                        continue

                    if row[1]['flag']==16:
                        seq_final.append(initial+base.lower())
                    else:
                        seq_final.append(initial+base.upper())
                    
                    qual_final.append(base_qual)
                except:
                    pass

        seq_final=''.join(seq_final)
        qual_final=''.join(qual_final)

        if (position in del_dic.keys()): #(row[1]['flg_del']==1) & 
            attach=''.join(del_dic[position])
            seq_final=seq_final+attach

        if (position in ins_dic.keys()) : #(row[1]['flg_ins']==1) & 
            attach=''.join(ins_dic[position])
            seq_final=seq_final+attach

        length=len(qual_final)

        if (length != 0) and (performance==False):
            print(Chr_name,position,'N',length,seq_final,qual_final)
    
    stop=time.time()
    duration=(stop-begin)

    if performance==True:
        print(db_duration,duration)


if __name__ == "__main__":
    main(sys.argv[1:])