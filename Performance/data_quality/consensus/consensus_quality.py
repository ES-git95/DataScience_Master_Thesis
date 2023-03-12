import random
from Bio import Align
import os


def main():
    aligner=Align.PairwiseAligner()
    
    nb_interval=200
    length=20000
    limit_low=0
    limit_high=32078900
    
    
    
    for x in range(nb_interval):
        number = random.randint(limit_low,limit_high-length)
        min_number=number
        max_number= number + length
    
        os.system("python3 consensus_v5.py -r chr3R:"+str(min_number)+"-"+str(max_number)+" >> consensus_V5.txt")
        os.system("samtools consensus ENCFF469GFN_sorted.bam -r chr3R:"+str(min_number)+"-"+str(max_number)+" --mode simple --call-frac 0.55 >> samtools_consensus.txt")
    
        with open('samtools_consensus.txt','r') as f_read1,open('samtools_consensus_1.txt','w') as f_write1:
            f_write1.write("".join(line.strip() for line in f_read1)  )
    
        with open('consensus_V5.txt','r') as f_read2,open('consensus_V5_1.txt','w') as f_write2:
            f_write2.write("".join(line.strip() for line in f_read2)  )
    
        with open('consensus_V5_1.txt', 'r') as file_v5,  open('samtools_consensus_1.txt', 'r') as file_samtools,open('consensus_output_1.txt', 'a') as output:
            variableV5=file_v5.read().rstrip()
            variableSAM=file_samtools.read().rstrip()
            result=aligner.score(variableV5,variableSAM)
            output.write("min: {}\tmax: {}\tresult: {} \tperc: {}\n".format(min_number,max_number,result,result/length))
    
        os.system("rm consensus_V5.txt")
        os.system("rm samtools_consensus.txt")
        os.system("rm consensus_V5_1.txt")
        os.system("rm samtools_consensus_1.txt")


if __name__ == "__main__":
    main()

