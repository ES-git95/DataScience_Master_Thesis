#importazione delle librerie
import pysam
import pandas as pd
import csv
from pathlib import Path
import subprocess
import os
import warnings
import multiprocessing as mp
import string
import time
import sys, getopt

warnings.filterwarnings("ignore")

def CreateInsertionDF(dataframe,tag_dict,tag):
    '''Funzione di creazione del Dataframe pandas relativo ai metadati'''
    if tag!='CO':
        df=dataframe.set_index(tag_dict[tag]).stack().reset_index()
        df.insert(0,'id',tag_dict[tag])
        df.insert(0,'tag',tag)
    else:
        df=pd.DataFrame([tag,None,None,dataframe[tag].values[0][0]]).T
        
    df=df.set_axis(['main_tag','tag_key','tag_value','sub_tag_key','sub_tag_value'], axis=1, inplace=False)

    return df

def MetadataFilesCreation(samfile,tag_dict,df_metadata,files_dict):
    '''Funzione che crea i file csv per il popolamento delle tabelle ReferenceSequence e ReadGroup 
    e lancia la funzione CreateInsertionDF per il popolamento del dataframe dei metadati'''
    for tag in tag_dict.keys():

        try:
            df=pd.DataFrame(samfile.header[tag])
        except ValueError:
            df=pd.DataFrame(samfile.header[tag],index=[0])
        except KeyError:
            continue

        if tag=='SQ':

            df_RS=df[['SN','LN']]
            df_RS.to_csv('insert_RS.csv', header=None, index=False)
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

    return df_metadata,files_dict


def InsertReadTagsRows(sam_file,start,end,counter): 
    '''Funzione che popola  i file csv relativi alle tabelle Read e Tags'''
    letter=string.ascii_uppercase[counter]
    counter_id=0
    with open('samfile_read.csv','a',encoding='UTF8') as read_file,open('samfile_tags.csv','a',encoding='UTF8') as tags_file: #
        writer_reads=csv.writer(read_file)
        writer_tags=csv.writer(tags_file)
        for read_id,read in enumerate(sam_file.fetch(until_eof=True)):

            if (read_id>=start) & (read_id<end):
                d=read.to_dict()
    
                tags=d.pop('tags')
    
                if tags != []:
                    writer_tags.writerow([letter+str(counter_id),str(tags),read_id])
                    counter_id+=1

                writer_reads.writerow([read_id]+list(d.values())+[None])                 


def processor(csv_file,table):
    '''Processo di inserimento file nel database, instanziato per ciascuna tabella dello schema di database'''
    status=subprocess.run(["PGPASSWORD=pysam8192 psql -U postgres -h 127.0.0.1  -d pysamdb -c \"\copy "+ str(table) +" FROM '"+str(Path(csv_file).resolve()).replace('\\','\\\\')+"' delimiter ',' csv header\""],
                        shell=True, capture_output=True)
    if performance==False:
        print(status)
    os.remove(csv_file)


def main(argv):

    #si legge l'argomento in input che consiste nel nome del file
    try:
        opts, args = getopt.getopt(argv,"f:p",["file=","performance="])
    except getopt.GetoptError:
        print("DBinsertion.py -f <file name>")
        sys.exit()

    for opt, arg in opts:
        global performance
        performance=False
        if opt in ('-f','-file'):
            input_file=arg
        if opt in ('-p','-performance'):
            performance=True

    #viene estratto il numero di thread della macchina corrente
    tot_cpus=mp.cpu_count()
    if tot_cpus < 6: used_cpus=tot_cpus
    else: used_cpus=6

    files_number=used_cpus+1

    #vengono creati i file per il popolamento delle tabelle Read e Tag, e i rispettivi dataframe pandas
    files_dict={}
    files_dict['Read']='samfile_read.csv'
    files_dict['Tags']='samfile_tags.csv'
    pd.DataFrame(list()).to_csv('samfile_read.csv',header=False)
    pd.DataFrame(list()).to_csv('samfile_tags.csv',header=False)

    #si duplica il file di imput per il numero di processi che sfrutteranno il multithread
    save = pysam.set_verbosity(0) 
    input_dict={}
    for i in range(0,files_number):
        input_dict["samfile_{0}".format(i+1)]=pysam.AlignmentFile(input_file,"rb")
    pysam.set_verbosity(save)
    
    samfile=input_dict['samfile_1']

    begin=time.time()
    tag_dict={'HD':'VN','SQ':'SN','RG':'ID','PG':'ID'}

    #si crea il dataframe per Metadata
    #si lancia la funzione MetadataFilesCreation che popola il file relativo ai metadati
    df_metadata=pd.DataFrame(columns=['main_tag','tag_key','tag_value','sub_tag_key','sub_tag_value'])
    df_metadata,files_dict=MetadataFilesCreation(samfile,tag_dict,df_metadata,files_dict)
    df_metadata=df_metadata.reset_index(drop=True).reset_index().rename(columns={'index':'id'})
    df_metadata.to_csv('insert_meta.csv', header=None, index=False)
    files_dict['Metadata']='insert_meta.csv'

    #si conteggia il numero di righe sulla base del quale suddividere i processi
    reads_count=samfile.count(until_eof=True)
    read_count_cpu=reads_count//(used_cpus)
    start=0
    end=read_count_cpu

    #in questa fase si lanciano i processi relativi al popolamento dei file Read e Tags
    workers=[]
    counter=0
    for file_name,sam_file in input_dict.items():
        if file_name != 'samfile_1':

            p=mp.Process(name=file_name,target=InsertReadTagsRows,args=(sam_file,start,end,counter)) 
            p.start()
            workers.append(p)
            start=start+read_count_cpu
            if str(files_number-1) in file_name:
                end=reads_count
            else:
                end=end+read_count_cpu

            counter+=1


    for worker in workers:
        worker.join() 

    interval_time=time.time()
    interval=interval_time-begin
    #print('\nFiles have been created in '+str(interval/60)+'\n')

    #in questa fase vengono lanciati in parallelo gli inserimenti nel database
    insert_workers=[]
    for table,csv_file in files_dict.items():
        p=mp.Process(name=table,target=processor,args=(csv_file,table))
        p.start()
        insert_workers.append(p)
        
    for worker in insert_workers:
        worker.join()

    stop=time.time()
    duration=(stop-begin)/60

    if performance==False:
        print('\nCongratulations! All the records of '+input_file+' file has been inserted in the database.')
        print('Overall time: '+str(duration))
    else:
        print(duration)

if __name__ == "__main__":
    main(sys.argv[1:])