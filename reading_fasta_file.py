# 
# This requires that the Bio package is downloaded on your computer
# Download at: https://biopython.org/wiki/Download
# 

from Bio import SeqIO
import numpy as np

records = list(SeqIO.parse('Total.0416.2020.clean.fa','fasta'))

n_taxa = len(records)
sequence_length = len(records[0].seq)
sequences = []
sequence_labels = []
msa = np.zeros((n_taxa, sequence_length), dtype = "uint8")

for i in range(n_taxa):
    label = records[i].id
    sequence = str(records[i].seq)
    sequence_labels.append(label)
    sequences.append(sequence)
    msa[i] = np.fromstring(sequence, dtype = "uint8")

print(sequences, sequence_length, sequence_labels, msa)

'''
phylip_file = open('real_test2.phy')
phylip_header = phylip_file.readline().strip().split()

n_taxa = int(phylip_header[0])
sequence_length = int(phylip_header[1])
sequences = []
sequence_labels = []
msa = np.zeros((n_taxa, sequence_length), dtype = "uint8")

for i in range(n_taxa):
    line = phylip_file.readline()
    label, sequence = line.strip().split()
    sequence_labels.append(label)
    msa[i] = np.fromstring(sequence, dtype = "uint8")
    print(msa[i], type(msa[i]))
    sequences.append(sequence)

print(sequences, sequence_length, sequence_labels, msa)
'''