#importazione pacchetto
import re
import sqlite3
from DBfunctions import *
import pysam
from tqdm import tqdm

def SamFileRowsList(filesList,IndexList): #aggiungere gestione index
    InsertList=[]
    AlignmentList=[]
    file_num=len(filesList)
    for file in range(0,file_num):
        samfile=pysam.AlignmentFile(filesList[file],"rb",filepath_index=IndexList[file])
        AlignmentList.append(samfile)
        description=samfile.description
        if description:
            InsertList.append((filesList[file],IndexList[file],description))
        else:
            InsertList.append((filesList[file],IndexList[file],None))
    return InsertList, AlignmentList


#Le funzioni seguenti dovranno essere inserite in un ciclo sui file multipli inseriti
#questa funzione deve essere eseguita per 3 scopi:
#   1.estrarre l'id massimo dalla tabella RS_Dictionary prima di creare la lista della tabella RS_Dictionary
#   2.estrarre l'id corrispondente al file in considerazione dopo aver popolato la tabella SamFile
#   3.estrarre l'id massimo dalla tabella Read prima di creare la lista della tabella Read 



#estrazione header
def RS_DictionaryRowsList(sam_file,samfile_id,max_id):

    try:    
        sq_header=sam_file.header['SQ']
    except KeyError:
        sq_header=None

    InsertList=[]
    SN_dict={}

    if sq_header is not None:
        for i in sq_header:
            try:
                name=str(i['SN'])
                maxid=max_id+1
                SN_dict[name]=maxid
            except:
                name=None

            try:
                length=int(i['LN'])
            except:
                length=None
            
            InsertList.append((name,length,samfile_id))

    return InsertList, SN_dict


def ReadGroupRowsList(sam_file,samfile_id):
    
    try:
        rg_header=sam_file.header['RG']
    except KeyError:
        rg_header=None

    InsertList=[]
    
    if rg_header is not None:
        for i in rg_header:
            try:
                gruoup_id=str(i['ID'])
            except:
                gruoup_id=None

            try:
                sample=int(i['SM'])
            except:
                sample=None

            try:
                platform=int(i['PL'])
            except:
                platform=None
            
            InsertList.append((gruoup_id,sample,platform,samfile_id))

    return InsertList


def MetadataRowsList(sam_file,samfile_id):

    main_tags=['HD','PG','CO']

    InsertList=[]

    for main_tag in main_tags:

        try:
            header=sam_file.header[main_tag]

            for line in header:
                tags=line.keys()

                for tag in tags:
                    value=line[tag]
                    InsertList.append((main_tag,tag,value,samfile_id))

        except:
            header=None

    return InsertList


#NB add to SN_dict '*':None

def ReadRowsList(sam_file,samfile_id,SN_dict,max_id,key):
    
    ReadInsertList=[]
    TagsInsertList=[]

    if key != '*':
        rname=str(SN_dict[key])
    else:
        rname= None

    read_id=max_id

    for read in sam_file.fetch(key): #remembre that * is not included
        read_dict=read.to_dict()

        qname=str(read_dict['name'])
        flag=int(read_dict['flag'])
        pos=int(read_dict['ref_pos'])
        mapq=int(read_dict['map_quality'])
        cigar=str(read_dict['cigar'])
        rnext=str(read_dict['next_ref_name'])
        pnext=int(read_dict['next_ref_pos'])
        tlen=int(read_dict['length'])
        seq=str(read_dict['seq'])
        qual=str(read_dict['qual'])
        group_id=None  #this is fo now....
        samfileid=samfile_id

        ReadInsertList.append((qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,group_id,samfileid))

        read_id += 1
        tags=read_dict['tags']
        tags_string=''

        if tags != []:
            for tag in tags:
                tags_string=tags_string+re.sub(':(.*?):',':',tag)+';'
            tags_string=tags_string[:-1]

            TagsInsertList.append((tags_string,read_id))

    return ReadInsertList, TagsInsertList


def insertRecordsIntoTable(TableName,rowsList,connection,record_len=0):

    if len(rowsList)==0 and TableName not in ['Read','Tags']:
        print('WARNING: The list of rows to populate '+TableName+' table is empty. Insertion will be skipped \n')
    
    else:

        cursor=connection.cursor()

        if TableName not in ['Read','Tags']:
            print('Ready to populate '+TableName+' table')
        else:
            record_len+=len(rowsList)
    
        str_info="pragma table_info('"+TableName+"')"
        cursor.execute(str_info)
        columns_tuple=tuple([i[1] for i in cursor.fetchall() if i[1]!='id'])
        cursor.close()
        
        len_tuple=len(columns_tuple)
        columns_string=str(columns_tuple).replace("u'","").replace("'","")
        values_string="("+','.join("?"*len_tuple)+")"
    
        insertQuery="INSERT INTO "+TableName+" "+columns_string+" VALUES "+values_string+";"
    
        try:
            cursor=connection.cursor()
            cursor.executemany(insertQuery,rowsList)
            connection.commit()

            if TableName not in ['Read','Tags']:
                print(str(cursor.rowcount)+" records successfully inserted into "+TableName+"\n")
    
            cursor.close()
        
        except sqlite3.Error as error:
            print("Insertion into "+TableName+" has failed: "+str(error))

    return record_len
    
#iterare sull'indice?

#inserire barra di progressione



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
    input_files=[input('Please enter the single .bam file name or the list of files: \n')]
    index_files=[input('Please enter the corrisponding .bam.bai index file or list of files: \n')]
    connection=DBConnection('PysamDB.db','OPEN')

    InsertList_sam, AlignmentList=SamFileRowsList(input_files,index_files)
    insertRecordsIntoTable('SamFile',InsertList_sam,connection)

    

    for num,file in enumerate(AlignmentList):

        file_id=GetMaxIdFromTable(InsertList_sam[num][0],'SamFile','SAMFILE_ID',connection)
        max_RSDictionary_id=GetMaxIdFromTable(InsertList_sam[num][0],'RS_Dictionary','MAX',connection)

        InsertList_RS, SN_dict=RS_DictionaryRowsList(file,file_id,max_RSDictionary_id)
        insertRecordsIntoTable('RS_Dictionary',InsertList_RS,connection)
        
        InsertList_RG=ReadGroupRowsList(file,file_id)
        insertRecordsIntoTable('ReadGroup',InsertList_RG,connection)
        
        InsertList_Meta=MetadataRowsList(file,file_id)
        insertRecordsIntoTable('Metadata',InsertList_Meta,connection)

        Index_stat=file.get_index_statistics()
        Index_stat_len=len(Index_stat)

        print("\nReady to populate Read table")

        pbar = tqdm(total=Index_stat_len)
        Read_len=0
        Tags_len=0

        for index in Index_stat:
            max_Read_id=GetMaxIdFromTable(file,'Read','MAX',connection)
            InsertList_Read,InsertList_Tags=ReadRowsList(file,file_id,SN_dict,max_Read_id,index[0])

            Read_len=insertRecordsIntoTable('Read',InsertList_Read,connection,Read_len)
            Tags_len=insertRecordsIntoTable('Tags',InsertList_Tags,connection,Tags_len)
            
            pbar.update(n=1)

        print("\n "+str(Read_len)+" records succesfully inserted into Reads")
        print("\n "+str(Tags_len)+" records succesfully inserted into Tags")

        print('\n Congratulations! All the records of '+InsertList_sam[num][0]+' file has been inserted in the database.')


if __name__ == "__main__":
    main()