import pandas as pd
import sqlite3
import sys, getopt
import ast
import re

def main(argv):

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
    query="select pos,pos_shifted,mapq,seq,seq_insdel,qual,qual_insdel,flag,pos_ixd,attachment from Read left join InsDel on Read.id=InsDel.read_id where Read.rname='"+Chr_name+"' and Read.pos between ("+Starting_pos+"-50+1) and "+Ending_pos+";"

    df = pd.read_sql(query, con)  

    idx_list=df[df['seq_insdel'].notnull()].index.tolist()
    df.update(pd.DataFrame({'seq':list(df.loc[idx_list,'seq_insdel']),'qual':list(df.loc[idx_list,'qual_insdel'])},index=idx_list))

    df_attachments=df.loc[idx_list,['pos_ixd','attachment','flag']]
    try:
        df_attachments['pos_ixd']=df_attachments['pos_ixd'].astype('int64').astype('str')
    except ValueError:
        pass

    attach_idx_set=set(','.join(list(df_attachments['pos_ixd'])).split(','))

    pos_list=list(sorted(set(df.pos)))

    df_before_seq=pd.DataFrame()
    Starting_pos=int(Starting_pos)
    Ending_pos=int(Ending_pos)
    for position in range(Starting_pos,Ending_pos+1):

        if not not pos_list:
            if position == pos_list[0]:
                df_idx=pos_list.pop(0)
                df_sub=df[(df.pos==df_idx)]
                df_seq_tolist=pd.DataFrame(df_sub[['seq','flag']].apply(lambda x: list(x['seq'].lower()) if x['flag']==16 else list(x['seq'].upper()),axis=1).tolist()).set_index([df_sub.pos])
                df_qual_tolist=pd.DataFrame(df_sub.qual.apply(list).tolist()).set_index([df_sub.pos])
                df_bool=~(df_qual_tolist=='-')
                df_seq_final=df_seq_tolist.where(df_bool,' ') #--> rivedere questione dimensionalità
                df_qual_final=df_qual_tolist.where(df_bool,' ')

                df_seq_t=df_seq_final.T.agg(lambda x: ''.join(x).replace(' ',''), axis=1)
                new_idxs=list(range(position,position+len(df_seq_t)))
                df_seq_t.index=new_idxs
                df_seq_t[position]=('^'+chr(df_sub.iloc[0]['mapq']+33)).join(' '+df_seq_t[position])[1:]
                df_seq_t[position+len(df_seq_t)-1]='$'.join(df_seq_t[position+len(df_seq_t)-1]+' ')[:-1]
                df_qual_t=df_qual_final.T.agg(lambda x: ''.join(x).replace(' ',''), axis=1)
                df_qual_t.index=new_idxs

                if df_before_seq.empty==False:
                    df_before_seq=pd.concat([df_before_seq, df_seq_t], axis=1, join="outer").fillna('').agg(''.join, axis=1).sort_index()
                    df_before_qual=pd.concat([df_before_qual, df_qual_t], axis=1, join="outer").fillna('').agg(''.join, axis=1).sort_index()
                else:
                    df_before_seq=df_seq_t
                    df_before_qual=df_qual_t

        if df_before_seq.empty==False:
            attachment=''
            if str(position) in attach_idx_set:
                df_sub_attach=df_attachments[df_attachments['pos_ixd'].str.contains(str(position))]
                for attach in df_sub_attach.iterrows():
                    attach_idx=attach[1]['pos_ixd'].split(',').index(str(position))
                    if attach[1]['flag']==16:
                        attachment+=attach[1]['attachment'].split(',')[attach_idx].lower()
                    else:
                        attachment+=attach[1]['attachment'].split(',')[attach_idx].upper()

            seq_final=df_before_seq.pop(position)
            qual_final=df_before_qual.pop(position)
            length=len(qual_final)
            if seq_final.strip() != '':
                print(Chr_name,position,'N',length,seq_final+attachment,qual_final)


if __name__ == "__main__":
    main(sys.argv[1:])