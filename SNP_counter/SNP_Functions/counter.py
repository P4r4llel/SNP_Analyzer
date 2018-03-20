# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 14:53:57 2018

@author: V3rt1
COUNTER
"""

import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

input_file = "reference.fasta"

record_iterator = SeqIO.parse(input_file, "fasta")
output_sequences = []

while True:
    record = next(record_iterator, None)
    if record is None:
        break
    record.seq.alphabet = generic_dna
    A = str(record.seq.count("A"))
    T = str(record.seq.count("T"))
    G = str(record.seq.count("G"))
    C = str(record.seq.count("C"))
    print (record.seq[3:8])

print ("A = " + A + "; T = " + T + "; G = " + G + "; C = " + C)

# Test repo bullshit