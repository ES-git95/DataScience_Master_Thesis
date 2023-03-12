import re
import pandas as pd
import math

def main():
    file1 = open("mpileup_V5.txt",'r')
    file2 = open("mpileup_samtools.txt",'r')

    file1_lines = file1.readlines()
    file2_lines = file2.readlines()

    with open('V5_adj.txt', 'w') as f1,open('samtools_adj.txt', 'w') as f2:

        len_file1=math.ceil(len(file1_lines))
        len_file2=math.ceil(len(file2_lines))

        for i in range(len_file1):
            listed1=re.sub(' +', ' ',file1_lines[i])
            f1.write(listed1)

        for i in range(len_file2):
            listed2=re.sub(' +', ' ',file2_lines[i])
            f2.write(listed2)

    file1.close()
    file2.close()

    df_1=pd.read_csv('V5_adj.txt', sep=' ',header=None,names=['rname1','pos1','N1','num1','seq1','qual1'])
    df_2=pd.read_csv('samtools_adj.txt', sep=' ',header=None,names=['rname2','pos2','N2','num2','seq2','qual2'])
    df_result = pd.merge(df_1, df_2, how="outer",left_on='pos1', right_on='pos2')

    df_result_excluded=df_result[(df_result.pos1.isnull())|(df_result.pos2.isnull()) ]
    df_result_final=df_result.drop(df_result_excluded.index)

    correct=0
    with open('unmatched_lines.txt', 'w') as output:
        for index, row in df_result_final.iterrows():
            compare1=''.join(sorted(row[['num1','seq1','qual1']].apply(lambda x: ''.join(sorted(str(x))))))
            compare2=''.join(sorted(row[['num2','seq2','qual2']].apply(lambda x: ''.join(sorted(str(x))))))

            if compare1==compare2:
                correct+=1
            else:
                len1=len(compare1)
                len2=len(compare2)
                if len1 < len2:
                    len_diff=len2-len1
                    compare1=compare1 + ' '*len_diff
                if len1>len2:
                    len_diff=len1-len2
                    compare2=compare2 + ' '*len_diff
                else:
                    pass
                
                difference1 = ''.join(x for x, y in zip(list(compare1), list(compare2)) if x != y)
                difference2 = ''.join(y for x, y in zip(list(compare1), list(compare2)) if x != y)

                output.write("index: {} \tdifference: {} {} \tcompare: {} {}".format(index,difference1,difference2,compare1,compare2))

        print('\noverall correctness: {}'.format(correct/len(df_result_final)))

if __name__ == "__main__":
    main()