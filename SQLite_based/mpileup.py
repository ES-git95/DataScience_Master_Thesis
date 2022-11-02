import pandas as pd
import sqlite3
import re
import sys, getopt

def main(argv):

    Chr_name='chrX'
    Starting_pos='12245035'
    Ending_pos='12245084'

    try:
        opts, args = getopt.getopt(argv,"r:",["region="])
    except getopt.GetoptError:
        print("mpileup.py -r <chr region> (ex -r chr:1-100)")
        sys.exit()

    for opt, arg in opts:
        if opt in ('-r','-region'):
            arg_list=arg.split(':')
            Chr_name=arg_list[0]
            Starting_pos=arg_list[1].split('-')[0]
            Ending_pos=arg_list[1].split('-')[1]

    con = sqlite3.connect('PysamDB.db')
    query="select mapq, flag, seq,pos,cigar,qual from Read where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-50+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)

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
                df.loc[row[0],'seq']=row[1]['seq'][:del_idx+counter] + '*' + row[1]['seq'][del_idx+counter:]
                df.loc[row[0],'qual']=row[1]['qual'][:del_idx+counter+1] + row[1]['qual'][del_idx+counter:]
                dic_index=row[1]['pos']+del_idx+counter-1
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
                df.loc[row[0],'seq']=row[1]['seq'][:ins_idx+counter] + row[1]['seq'][ins_idx+1+counter:]
                df.loc[row[0],'qual']=row[1]['qual'][:ins_idx+counter] + row[1]['qual'][ins_idx+1+counter:]
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

    for position in range(int(Starting_pos),int(Ending_pos)+1): #range(30000547,30000548):#
        #ixds_filtered=sorted(set(filter(lambda idx: (idx>=position-49 and idx<=position), idxs)))
        ixds_filtered=sorted(set([a for a,b in zip(idxs,idxs_new) if (b>=position-49 and a<=position)]))
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

        if length != 0:
            print(Chr_name,position,'N',length,seq_final,qual_final)


if __name__ == "__main__":
    main(sys.argv[1:])