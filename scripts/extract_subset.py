#!/usr/bin/env python
#--------------------------------------------------------------------------------------------------------------------
# Script that is used to retrieve a subset of sequences from a multifasta. As input the function takes the 
# headerslines of the desired sequences, the complete set of headerlines from the original multifasta file, the 
# multifasta file and the desired name for the file to contain the subset of sequences.
#--------------------------------------------------------------------------------------------------------------------
# Import required packages
import re
import sys
from Bio import SeqIO

# Import files
file1 = open(sys.argv[1], 'r')
protein_fasta_headers = file1.readlines()
file1.close()

file2 = open(sys.argv[2], 'r')
nucleotide_fasta_headers = file2.readlines()
file2.close()

# Retrive the nucleotide fasta sequences coresponding to the protein fasta headers
with open(sys.argv[3]) as original_fasta, open(sys.argv[4], 'w') as corrected_fasta:
    records = SeqIO.parse(original_fasta, 'fasta')
    i = 0
    for record in records:
        if nucleotide_fasta_headers[i] not in protein_fasta_headers:
            i += 1
            continue
        else:
            SeqIO.write(record, corrected_fasta, 'fasta')
            i += 1
            continue
