# Importing necessary packages, including Streamlit, BiopPython, and Math
import streamlit as st
import Bio
import math
from math import log
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
from Bio import SeqIO
from Bio import Entrez

# Adds a title
st.write("""
# DNA Sequence Alignment Application

This app aligns two sequences
""")


# This cell has a Gap Function.
# The gap function returns a score of 0 if there is no gap,
# if there is 1 gap it returns -2,
# and otherwise it returns an a negative score based on how many gaps are present.


def gap_function(x, y):
    if y == 0:  # no gap
        return 0
    elif y == 1:  # gap open penalty
        return -2
    return - (2 + y / 4.0 + log(y) / 2.0)


###########################################
# This cell has a Global Alignment function.
# It uses pairwise2 algorithm to align two sequences.
# It takes in sequence A, sequence B, and a gap function and returns an alignment and a score.
# Match score represents score increase due to matching characters.
# Mismatch score represents score decrease due to mismatch characters.
# By default match score is set to 5, and mismatch score is set to -4.
# The 'format_alignment()' function formats the alignment prettily into a string.

# check if it's a valid sequence
DNA_Nucleotides = ['A', 'C', 'G', 'T']
def validate_seq(seq):
    """
    Check the sequence to make sure it is a valid DNA string
    """

    if seq == None:
        return False
    if not isinstance(seq, str):
        return False

    tmpseq = seq
    for nuc in tmpseq:
        if nuc not in DNA_Nucleotides:
            return False

    return tmpseq


def global_alignment_with_validation(SeqA, SeqB, gap_function, match_score=5, mismatch_score=-4):
  if validate_seq(SeqA) and validate_seq(SeqB):
    alignment = pairwise2.align.globalmc(SeqA, SeqB, match_score, mismatch_score,
                                       gap_function, gap_function)
    return(format_alignment(*alignment[0]))
  else:
    return "Please enter a valid DNA sequence."


def global_alignment(SeqA, SeqB, gap_function, match_score=5, mismatch_score=-4):
    alignment = pairwise2.align.globalmc(SeqA, SeqB, match_score, mismatch_score,
                                         gap_function, gap_function)

    return (format_alignment(*alignment[0]))

#############################################
#############################################
# Fetching Data from NCBI
#############################################

# This function downloads records from genbank
def retrieveRecord(filename):
    if not os.path.isfile(filename):

        # Downloading...
        filename_short = filename.replace(".fasta", "")
        Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
        net_handle = Entrez.efetch(
            db="nucleotide", id=filename_short, rettype="fasta", retmode="text"
        )
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        # print("Saved")
        # record = SeqIO.read(filename, "fasta")
        return None


# This function uses SeqIO.read to read the record into an object, record.
def returnRecord(filename):
    record = SeqIO.read(filename, "fasta")
    return record

def retrieve_return_record(filename):
    retrieveRecord(filename)
    record = SeqIO.read(filename, "fasta")
    return record
def retrieve_return_record_validation(filename):
    try:
        retrieveRecord(filename)
        record = SeqIO.read(filename, "fasta")
    except Exception:
        record = "Please enter a valid sequence"
    return record

def get_seq_with_validation(record):
    try:
        seq = record.seq
    except Exception:
        seq = ""
    return seq


    #Seq_A_record = retrieve_return_record_validation(Code_1)
    #Seq_B_record = retrieve_return_record_validation(Code_2)

#    Seq_A_record_description = Seq_A_record
 #   Seq_B_record_description = Seq_B_record

    # Seq_A_record_seq_full = get_seq_with_validation(Seq_A_record)
    # Seq_B_record_seq_full = get_seq_with_validation(Seq_B_record)
##########################################################

# Input DNA Code 1
# Input DNA code 1 index
# Input DNA Code 2
# Input DNA code 2 index
Code_1 = "MG762674"  # rosettus bat coronavirus
Code_1_index_1 = 0
Code_1_index_2 = 100
Code_2 = "NC_019843"  # Middle Eastern coronavirus strain
Code_2_index_1 = 0
Code_2_index_2 = 100

col1, col2 = st.columns(2)

# Input DNA code 1  and DNA code 2
with col1:
    st.sidebar.header('Input code 1 and code 2')
    Code_1 = st.text_input('Enter DNA code 1 as text (example: MG762674 Rousettus bat coronavirus HKU9 isolate)')
    st.write('Sequence A is:', Code_1)
with col2:
    Code_2 = st.text_input('Enter DNA code 2 as text (example: NC_019843 Middle East respiratory syndrome-related coronavirus isolate HCoV-EMC/2012)')
    st.write('Sequence B is :', Code_2)

# Input DNA code 1 index 1 and DNA code 2 index 1
with col1:
    st.sidebar.header('Input DNA sequence A Indices')
    Code_1_index_1 = st.text_input('Enter DNA code 1 index 1 (example: 0)' )
    st.write('DNA sequence A Index 1 is:', Code_1_index_1)
    Code_1_index_2 = st.text_input('Enter DNA code 1 index 2 (example: 100)')
    st.write('DNA sequence A Index 2 is:', Code_1_index_2)

with col2:
    st.sidebar.header('Input DNA sequence B Indices')
    Code_2_index_1 = st.text_input('Enter DNA code 2 index 1 (example: 0)')
    st.write('DNA sequence B Index 1 is:', Code_2_index_1)
    Code_2_index_2 = st.text_input('Enter DNA code 2 index 2 (example: 100)')
    st.write('DNA sequence B Index 2 is:', Code_2_index_2)
# Input DNA code 1 index 2 and DNA code 2 index 2


# Output Alignment with codes
# Output Alignment
Seq_A_record = retrieve_return_record_validation(Code_1)
Seq_B_record = retrieve_return_record_validation(Code_2)

Seq_A_record_description = Seq_A_record
Seq_B_record_description = Seq_B_record


#Seq_A_record_seq_full = Seq_A_record.seq
#Seq_B_record_seq_full = Seq_B_record.seq
#Seq_A_record_seq_full = retrieve_return_record_validation(Code_1).seq
#Seq_B_record_seq_full = retrieve_return_record_validation(Code_2).seq
Seq_A_record_seq_full = get_seq_with_validation(Seq_A_record)
Seq_B_record_seq_full = get_seq_with_validation(Seq_B_record)

Seq_A_record_seq_indexed = Seq_A_record_seq_full[int(Code_1_index_1): int(Code_1_index_2)]
Seq_B_record_seq_indexed = Seq_B_record_seq_full[int(Code_2_index_1): int(Code_2_index_2)]

#st.write('Sequence A description:', Seq_A_record_description)
st.write('Sequence A description:')
st.write(Seq_A_record_description)
st.write('Sequence B description:')
st.write(Seq_B_record_description)

# Output Alignment Using Codes
st.subheader('Alignment using Codes')
alignment_using_codes = global_alignment(Seq_A_record_seq_indexed, Seq_B_record_seq_indexed, gap_function)
# st.subheader(alignment)
#st.text('Aligned Sequence using codes', alignment_using_codes)
#st.text(alignment_using_codes)
st.text(alignment_using_codes)


