#Importing necessary packages, including Streamlit, BiopPython, and Math
import streamlit as st
import Bio
import math
from math import log
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
# Importing necessary packages, including Streamlit, BiopPython, and Math
import os
from Bio import SeqIO
from Bio import Entrez


#Adds a title
st.write("""
# DNA Sequence Alignment Application

This app aligns two sequences. It also computes GC content, and Amino acid frequency.
""")

#This cell has a Gap Function. 
#The gap function returns a score of 0 if there is no gap, 
#if there is 1 gap it returns -2, 
#and otherwise it returns an a negative score based on how many gaps are present.


def gap_function(x, y): 
  if y == 0: # no gap
    return 0
  elif y== 1: # gap open penalty
    return -2
  return - (2 + y/4.0 + log(y)/2.0)

###########################################
#This cell has a Global Alignment function. 
#It uses pairwise2 algorithm to align two sequences. 
#It takes in sequence A, sequence B, and a gap function and returns an alignment and a score. 
#Match score represents score increase due to matching characters. 
#Mismatch score represents score decrease due to mismatch characters. 
#By default match score is set to 5, and mismatch score is set to -4. 
#The 'format_alignment()' function formats the alignment prettily into a string.

def global_alignment(SeqA, SeqB, gap_function, match_score=5, mismatch_score=-4):
  
  alignment = pairwise2.align.globalmc(SeqA, SeqB, match_score, mismatch_score,
                                       gap_function, gap_function)
  
  return(format_alignment(*alignment[0]))


# source :https://python.plainenglish.io/bioinformatics-in-python-dna-toolkit-part-1-validating-and-counting-nucleotides-1cd1037ef06f
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
  SeqA = SeqA.upper()
  SeqB = SeqB.upper()
  if validate_seq(SeqA) and validate_seq(SeqB):
    alignment = pairwise2.align.globalmc(SeqA, SeqB, match_score, mismatch_score,
                                       gap_function, gap_function)
    return(format_alignment(*alignment[0]))
  else:
    return "Please enter a valid DNA sequence."

# This function uses SeqIO.read to read the record into an object, record.
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

def returnRecord(filename):
    record = SeqIO.read(filename, "fasta")
    return record

def retrieve_return_record(filename):
    retrieveRecord(filename)
    record = SeqIO.read(filename, "fasta")
    return record


def gc_content_percentage(seq):
  """GC Content in a DNA/RNA sequence"""
  GC_content = round((seq.count('C') + seq.count('G')) / len(seq) *100)
  return str(GC_content) + "%"

def at_content_percentage(seq):
  """AT Content in a DNA/RNA sequence"""
  AT_content = round((seq.count('A') + seq.count('T')) / len(seq) *100)
  return str(AT_content) + "%"

DNA_Codons_full_name = {
    # 'M' - START, '_' - STOP
    "GCT": "Alanine, Ala, A", "GCC": "Alanine, Ala, A", "GCA": "Alanine, Ala, A", "GCG": "Alanine, Ala, A",
    "TGT": "Cysteine, Cys, C", "TGC": "Cysteine, Cys, C",
    "GAT": "Aspartic Acid, Asp, D", "GAC": "Aspartic Acid, Asp, D",
    "GAA": "Glutamic Acid, Glut, E", "GAG": "Glutamic Acid, Glut, E",
    "TTT": "Phenylalanine, Phe, F", "TTC": "Phenylalanine, Phe, F",
    "GGT": "Glycine, Gly, G", "GGC": "Glycine, Gly, G", "GGA":"Glycine, Gly, G", "GGG":"Glycine, Gly, G",
    "CAT": "Histidine, His, H", "CAC": "Histidine, His, H",
    "ATA": "Isoleucine, Ile, I", "ATT": "Isoleucine, Ile, I", "ATC": "Isoleucine, Ile, I",
    "AAA": "Lysine, Lys, K", "AAG": "Lysine, Lys, K",
    "TTA": "Leucine, Leu, L", "TTG": "Leucine, Leu, L", "CTT": "Leucine, Leu, L", "CTC": "Leucine, Leu, L", "CTA": "Leucine, Leu, L", "CTG":"Leucine, Leu, L",
    "ATG": "Methionine, Met, M",
    "AAT": "Asparagine, Asn, N", "AAC": "Asparagine, Asn, N",
    "CCT": "Proline, Pro, P", "CCC": "Proline, Pro, P", "CCA":"Proline, Pro, P", "CCG":"Proline, Pro, P",
    "CAA": "Glutamine, Gln, Q", "CAG": "Glutamine, Gln, Q",
    "CGT": "Arginine, Arg, R", "CGC": "Arginine, Arg, R", "CGA": "Arginine, Arg, R", "CGG": "Arginine, Arg, R", "AGA": "Arginine, Arg, R", "AGG": "Arginine, Arg, R",
    "TCT": "Serine, Ser, S", "TCC": "Serine, Ser, S", "TCA": "Serine, Ser, S", "TCG": "Serine, Ser, S", "AGT": "Serine, Ser, S", "AGC": "Serine, Ser, S",
    "ACT": "Threonine, Thr, T", "ACC": "Threonine, Thr, T", "ACA": "Threonine, Thr, T", "ACG": "Threonine, Thr, T",
    "GTT": "Valine, Val, V", "GTC": "Valine, Val, V", "GTA": "Valine, Val, V", "GTG": "Valine, Val, V",
    "TGG": "Tryptophan, Trp, W",
    "TAT": "Tyrosine, Tyr, Y", "TAC": "Tyrosine, Tyr, Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}
def translate_seq_full_name(seq, init_pos=0):
  """Translates a DNA sequence into an aminoacid sequence"""
  return [
          DNA_Codons_full_name[seq[pos:pos + 3]]
          for pos in range(init_pos, len(seq) - 2, 3)
  ]
def peptide_frequency_full_name(seq):
  """
  Count the codons in a given sequence.
  Return a dictionary.
  """
  tmpFreqDict = {
    # 'M' - START, '_' - STOP
    "Alanine, Ala, A": 0,
    "Cysteine, Cys, C": 0,
    "Aspartic Acid, Asp, D": 0,
    "Glutamic Acid, Glut, E": 0,
    "Phenylalanine, Phe, F": 0,
    "Glycine, Gly, G": 0,
    "Histidine, His, H": 0,
    "Isoleucine, Ile, I": 0,
    "Lysine, Lys, K": 0,
    "Leucine, Leu, L": 0,
    "Methionine, Met, M": 0,
    "Asparagine, Asn, N": 0,
    "Proline, Pro, P": 0,
    "Glutamine, Gln, Q": 0,
    "Arginine, Arg, R": 0,
    "Serine, Ser, S": 0,
    "Threonine, Thr, T": 0,
    "Valine, Val, V": 0,
    "Tryptophan, Trp, W": 0,
    "Tyrosine, Tyr, Y": 0,
    "_": 0
  }
  tmpPeptide = translate_seq_full_name(seq)

  for pep in tmpPeptide:
    tmpFreqDict[pep] += 1
  return tmpFreqDict

def global_alignment_no_string(SeqA, SeqB, gap_function=gap_function, match_score=5, mismatch_score=-4):
    alignment = pairwise2.align.globalmc(SeqA, SeqB, match_score, mismatch_score,
                                         gap_function, gap_function)

    return (alignment[0])

def return_alignment_score_ratio(SeqA, SeqB):
  """Takes Sequence A and computes an alignment against itself.
  Next, it aligns Sequence A and Sequence B and computes the score.
  Finally, it computes a ratio of the ideal alignment to the actual alignment"""
  ideal_alignment = global_alignment_no_string(SeqA, SeqA)
  ideal_score = ideal_alignment.score
  actual_alignment = global_alignment_no_string(SeqA, SeqB)
  actual_score = actual_alignment.score
  actual_to_ideal_ratio = actual_score / ideal_score
  return actual_score, ideal_score, actual_to_ideal_ratio

import numpy as np
import matplotlib.pyplot as plt
#sizes = np.array([actual_score, ideal_score])

def plot_alignment_score(SeqA, SeqB):
  actual_score, ideal_score, actual_to_ideal_ratio = return_alignment_score_ratio(SeqA, SeqB)
  labels_alignment_score = 'Actual Alignment Score, Points', 'Potential Score, Points'
  potential_score = ideal_score - actual_score
  sizes_alignment_scores = [actual_score, potential_score]
  fig1, ax1 = plt.subplots()
  def absolute_value(val):
    a  = np.round(val/100.*sum(sizes_alignment_scores), 0)
    return a
  ax1.pie(sizes_alignment_scores, labels=labels_alignment_score, autopct = absolute_value,
          shadow=True, startangle=90)
  ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

  plt.show()

# Input DNA Sequence A and B
st.sidebar.header('Input sequence 1 and Sequence 2')
SeqA = 'AAATTT'
SeqB = 'ATAATT'
#st.set_option('deprecation.showPyplotGlobalUse', False)
col1, col2 = st.columns(2)
with col1:
  SeqA = st.text_input('Enter DNA sequence 1 as text')
  st.write('Sequence A is:', SeqA)
with col2:
  SeqB = st.text_input('Enter DNA sequence 2 as text')
  st.write('Sequence B is :', SeqB)

#Output AlignmentS
st.subheader('Alignment') 
#alignment = global_alignment(SeqA, SeqB, gap_function)
alignment = global_alignment_with_validation(SeqA, SeqB, gap_function)

#st.subheader(alignment)

if st.button(alignment):

    st.text(alignment)
else:
    st.write("Click Button above to align sequences")
#st.text(alignment)

st.subheader('Alignment quality:',)

if SeqA and SeqB:
    st.pyplot(plot_alignment_score(SeqA, SeqB))
else:
    st.pyplot(None)


Code = st.text_input('Enter DNA string to get GC content')
st.write('GC content in DNA sequence is:',)

if Code:
    GC_content = gc_content_percentage(Code)
else:
    GC_content = '0%'
st.sidebar.header('Input DNA sequence code to retreive sequence')
st.markdown(GC_content)

Code = st.text_input('Enter DNA string to get amino acid frequencies')
st.write('Amino acid frequency in this DNA sequence is:',)
peptide_frequency = peptide_frequency_full_name(Code)
st.sidebar.header('Input DNA sequence code to retreive sequence')
st.write(peptide_frequency)

Code = st.text_input('Enter DNA code to pull from NCBI (example: Homo sapiens insulin BC005255)')
st.write('DNA sequence A Index 1  (Example: BC005255):', Code)
Seq_C_record = retrieve_return_record(Code)
Seq_C = Seq_C_record.seq
st.sidebar.header('Input DNA sequence code to retreive sequence')

if Seq_C:
    st.markdown(Seq_C)
else:
    Code = BC005255
    Seq_C_record = retrieve_return_record(Code)
    Seq_C = Seq_C_record.seq
    st.markdown(Seq_C)