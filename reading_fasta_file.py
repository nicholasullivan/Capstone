# 
# This requires that the Bio package is downloaded on your computer
# Download at: https://biopython.org/wiki/Download
# 

from Bio import SeqIO, AlignIO
import numpy as np
import os
from Bio.Align.Applications import MuscleCommandline
import subprocess






filepath = 'C:\\Users\\keerp\\Documents\\Data Science Practicum\\ECH2.align.capstone.muscle'

muscle_cline = subprocess.Popen(['muscle', '-in', filepath, '-quiet', '-stdout'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = muscle_cline.communicate()
records = list(SeqIO.parse(filepath,'fasta'))
alignment = AlignIO.read(filepath, "fasta")