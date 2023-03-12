import random
from Bio import Align
import os


def main():
    
    nb_interval=100
    length=20000
    limit_low=0
    limit_high=32078900
    
    
    
    for x in range(nb_interval):
        number = random.randint(limit_low,limit_high-length)
        min_number=number
        max_number= number + length
    
        os.system("python3 coverage_v5.py -r chr3R:"+str(min_number)+"-"+str(max_number)+" >> coverage_V5.txt")
        os.system("samtools coverage ENCFF469GFN_sorted.bam -r chr3R:"+str(min_number)+"-"+str(max_number)+" >> samtools_coverage.txt")
    
        with open('coverage_V5.txt', 'r') as file_v5,  open('samtools_coverage.txt', 'r') as file_samtools,open('coverage_output.txt', 'a') as output:
            variableV5=file_v5.read().rstrip()
            variableSAM=file_samtools.read().rstrip()
            output.write("V5\t"+variableV5+"\n")
            output.write("samtools\t"+variableSAM+"\n")
    
        os.system("rm coverage_V5.txt")
        os.system("rm samtools_coverage.txt")


if __name__ == "__main__":
    main()

