# 
# This require that the Bio package is downloaded on your computer
# Download at: https://biopython.org/wiki/Download
# 

from Bio import SeqIO
import numpy as np

#reading in file and splitting each row into list of lists
split_seq = []
for seq_record in SeqIO.parse('Total.0416.2020.clean.fa','fasta'):
    row = seq_record.seq
    split_seq.append([*row])

print(split_seq)
