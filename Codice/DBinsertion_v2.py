#importazione pacchetto
import pysam
from tqdm import tqdm
import pandas as pd
from tqdm import tqdm
import csv
from pathlib import Path
import subprocess
import os
import warnings


warnings.filterwarnings("ignore")

def CreateInsertionDF(dataframe,tag_dict,tag):
    if tag!='CO':
        df=dataframe.set_index(tag_dict[tag]).stack().reset_index()
        df.insert(0,'id',tag_dict[tag])
        df.insert(0,'tag',tag)
    else:
        df=pd.DataFrame([tag,None,None,dataframe[tag].values[0][0]]).T
        
    df=df.set_axis(['main_tag','tag_key','tag_value','sub_tag_key','sub_tag_value'], axis=1, inplace=False)

    return df

def MetadataFilesCreation(samfile,tag_dict,df_metadata,files_dict):

    for tag in tag_dict.keys():

        try:
            df=pd.DataFrame(samfile.header[tag])
        except ValueError:
            df=pd.DataFrame(samfile.header[tag],index=[0])
        except KeyError:
            continue

        if tag=='SQ':

            df_RS=df[['SN','LN']].reset_index()
            df_RS.to_csv('insert_RS.csv', header=None, index=False)
            RS_dict=df_RS[['index','SN']].set_index('SN').to_dict()['index']
            files_dict['ReferenceSequence']='insert_RS.csv'
            df=df.drop(['LN'], axis=1)

        if tag=='RG':

            df_RG=df[['ID','PL','SM']].reset_index()
            df_RG.to_csv('insert_RG.csv', header=None, index=False)
            files_dict['ReadGroup']='insert_RG.csv'
            df=df.drop(['PL','SM'], axis=1)

        if (len(df.columns) == 1) and (tag != 'CO') :
            continue
        else:
            df_metadata=df_metadata.append(CreateInsertionDF(df,tag_dict,tag))

    return df_metadata,RS_dict,files_dict


def InsertReadTagsRows(sam_file,key,writer_reads,writer_tags,read_id_last,counter_id,rname_id):

    read_id=read_id_last

    for read in sam_file.fetch(key): #remembre that * is not included
        d=read.to_dict()
        d.update({'ref_name':rname_id})

        tags=d.pop('tags')

        if tags != []:
            counter_id += 1
            writer_tags.writerow([counter_id,str(tags),read_id])

        writer_reads.writerow([read_id]+list(d.values())+[None]) 

        read_id += 1

    return read_id,counter_id


def main():

    n='n'
    y='y'

    while True:
        start=input('This program will populate PysamDB database tables from your input. Do you want to continue? [y|n]: ')
        if start not in ['y','n']:
            print("Please, enter 'y' to continue, 'n' to exit the program. Other inputs are not valid. \n")
        else:
            break

    print('\n')
    #inserire un comando per cambiare la directory da cui estrarre i file
    input_file=input('Please enter the single .bam file name: \n')
    index_file=input('Please enter the corrisponding .bam.bai index file: \n')

    samfile=pysam.AlignmentFile(input_file,"rb",filepath_index=index_file)

    tag_dict={'HD':'VN','SQ':'SN','RG':'ID','PG':'ID'}
    files_dict={}

    df_metadata=pd.DataFrame(columns=['main_tag','tag_key','tag_value','sub_tag_key','sub_tag_value'])
    df_metadata,RS_dict,files_dict=MetadataFilesCreation(samfile,tag_dict,df_metadata,files_dict)
    df_metadata=df_metadata.reset_index(drop=True).reset_index().rename(columns={'index':'id'})
    df_metadata.to_csv('insert_meta.csv', header=None, index=False)
    files_dict['Metadata']='insert_meta.csv'

    print("\nThe Metadata files are created, ready for insertion in the DB.\nNow the Read file will be created:")

    Index_stat=samfile.get_index_statistics()
    pbar = tqdm(total=len(Index_stat))

    with open('insert_read.csv','w',encoding='UTF8') as read_file,open('insert_tags.csv','w',encoding='UTF8') as tags_file:
        writer_reads=csv.writer(read_file)
        writer_tags=csv.writer(tags_file)
        read_id_last=0
        counter_id=0
        for index in Index_stat:  
            if index != '*':
                rname_id=str(RS_dict[index[0]])
            else:
                rname_id= None #questo per adesso 
            read_id_last,counter_id=InsertReadTagsRows(samfile,index[0],writer_reads,writer_tags,read_id_last,counter_id,rname_id)
            pbar.update(n=1)

        files_dict['Read']='insert_read.csv'
        files_dict['Tags']='insert_tags.csv'

    print("\n\nThe Read and Tags files are created, concerning "+str(read_id_last)+" reads, ready for insertion in the DB")

    for table,csv_file in files_dict.items():
        result = subprocess.run(['sqlite3',str(Path('PysamDB.db').resolve()),'-cmd','.mode csv','.import ' + str(Path(csv_file).resolve()).replace('\\','\\\\')+' '+str(table)],
                        capture_output=True)
        #os.remove("csv_file")
        print(str(result))

    print('\n Congratulations! All the records of '+input_file+' file has been inserted in the database.')


if __name__ == "__main__":
    main()