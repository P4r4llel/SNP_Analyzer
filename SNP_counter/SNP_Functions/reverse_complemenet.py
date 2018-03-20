# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 14:07:07 2018

@author: V3rt1
"""

import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

input_file = "original_rasmus.fasta"

record_iterator = SeqIO.parse(input_file, "fasta")
output_sequences = []

while True:
    record = next(record_iterator, None)
    if record is None:
        break
    record.seq.alphabet = generic_dna
    new_record = record.reverse_complement()
    new_record.id = record.id + "-RC"
    output_sequences.append(new_record)

    
SeqIO.write(output_sequences, input_file + "RC.fasta", "fasta")
    
    