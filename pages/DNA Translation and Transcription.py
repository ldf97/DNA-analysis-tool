#Importing necessary packages, including Streamlit, BiopPython, and Math
import streamlit as st
import Bio
import math
from math import log
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment




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

DNA_Codons = {
    # 'M' - START, '_' - STOP = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
def translate_seq(seq, init_pos=0):
  """Translates a DNA sequence into an aminoacid sequence"""
  return [
          DNA_Codons[seq[pos:pos + 3]]
          for pos in range(init_pos, len(seq) - 2, 3)
  ]
DNA_Nucleotides = ['A', 'C', 'G', 'T']
DNA_ReverseComplement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}


def reverse_complement(seq):
  """
  Swapping adenine with thymine and guanine with cytosine.
  Reversing newly generated string
  """
  return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]

def transcription(seq):
  """
  DNA -> RNA Transcription.
  Replacing Thymine with Uracil
  """
  return seq.replace("T", "U")

amino_acid_list =['Alanine, Ala, A',
            'Cysteine, Cys, C',
            'Aspartic Acid, Asp, D',
            'Glutamic Acid, Glut, E',
            'Phenylalanine, Phe, F',
            'Glycine, Gly, G',
            'Histidine, His, H',
            'Isoleucine, Ile, I',
            'Lysine, Lys, K',
            'Leucine, Leu, L',
            'Methionine, Met, M',
            'Asparagine, Asn, N',
            'Proline, Pro, P',
            'Glutamine, Gln, Q',
            'Arginine, Arg, R',
            'Serine, Ser, S',
            'Threonine, Thr, T',
            'Valine, Val, V',
            'Tryptophan, Trp, W',
            'Tyrosine, Tyr, Y'
]
##Adds a title
#Adds a title

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


# amino acid groupings
amino_acids_charged_side_chains = ['Arginine, Arg, R', 'Histidine, His, H', 'Lysine, Lys, K', 'Aspartic Acid, Asp, D',
                                   'Glutamic Acid, Glut, E']
amino_acids_positively_charged = ['Arginine, Arg, R', 'Histidine, His, H', 'Lysine, Lys, K']
amino_acids_negatively_charged = ['Aspartic Acid, Asp, D', 'Glutamic Acid, Glut, E']
amino_acids_with_polar_uncharged_side_chains = ['Serine, Ser, S', 'Threonine, Thr, T', 'Asparagine, Asn, N',
                                                'Glutamine, Gln, Q']
amino_acids_special_cases = ['Cysteine, Cys, C', 'Glycine, Gly, G', 'Proline, Pro, P']
amino_acids_with_hydrophobic_side_chains = ['Alanine, Ala, A', 'Valine, Val, V', 'Isoleucine, Ile, I',
                                            'Methionine, Met, M', 'Phenylalanine, Phe, F', 'Tyrosine, Tyr, Y',
                                            'Tryptophan, Trp, W']


def amino_acid_type_frequency(amino_acid_list):
    """
    Count the amino acids in a given class.
    Return a counts for each type in this order: polar uncharged side chain, charged side chain,
    positively charged side chain, negatively charged side chain, special cases, hydrophobic side chains
    """
    tmp_amino_acids_charged_side_chains_dict = {'Arginine, Arg, R': 0, 'Histidine, His, H': 0, 'Lysine, Lys, K': 0,
                                                'Aspartic Acid, Asp, D': 0, 'Glutamic Acid, Glut, E': 0}
    tmp_amino_acids_positively_charged_dict = {'Arginine, Arg, R': 0, 'Histidine, His, H': 0, 'Lysine, Lys, K': 0}
    tmp_amino_acids_negatively_charged_dict = {'Aspartic Acid, Asp, D': 0, 'Glutamic Acid, Glut, E': 0}
    tmp_amino_acids_with_polar_uncharged_side_chains_dict = {'Serine, Ser, S': 0, 'Threonine, Thr, T': 0,
                                                             'Asparagine, Asn, N': 0, 'Glutamine, Gln, Q': 0}
    tmp_amino_acids_special_cases_dict = {'Cysteine, Cys, C': 0, 'Glycine, Gly, G': 0, 'Proline, Pro, P': 0}
    tmp_amino_acids_with_hydrophobic_side_chains_dict = {'Alanine, Ala, A': 0, 'Valine, Val, V': 0,
                                                         'Isoleucine, Ile, I': 0, 'Leucine, Leu, L': 0,
                                                         'Methionine, Met, M': 0, 'Phenylalanine, Phe, F': 0,
                                                         'Tyrosine, Tyr, Y': 0, 'Tryptophan, Trp, W': 0}
    for amino in amino_acid_list:
        if amino in tmp_amino_acids_with_polar_uncharged_side_chains_dict:
            tmp_amino_acids_with_polar_uncharged_side_chains_dict[amino] += 1
        if amino in tmp_amino_acids_charged_side_chains_dict:
            tmp_amino_acids_charged_side_chains_dict[amino] += 1
        if amino in tmp_amino_acids_positively_charged_dict:
            tmp_amino_acids_positively_charged_dict[amino] += 1
        if amino in tmp_amino_acids_negatively_charged_dict:
            tmp_amino_acids_negatively_charged_dict[amino] += 1
        if amino in tmp_amino_acids_special_cases_dict:
            tmp_amino_acids_special_cases_dict[amino] += 1
        if amino in tmp_amino_acids_with_hydrophobic_side_chains_dict:
            tmp_amino_acids_with_hydrophobic_side_chains_dict[amino] += 1

    return sum(tmp_amino_acids_with_polar_uncharged_side_chains_dict.values()), \
           sum(tmp_amino_acids_charged_side_chains_dict.values()), \
           sum(tmp_amino_acids_positively_charged_dict.values()), \
           sum(tmp_amino_acids_negatively_charged_dict.values()), \
           sum(tmp_amino_acids_special_cases_dict.values()), \
           sum(tmp_amino_acids_with_hydrophobic_side_chains_dict.values())


# https://stackoverflow.com/questions/44493417/pandas-dataframe-bar-plot-plot-bars-different-colors-from-specific-colormap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_amino_class_frequency(amino_acid_list):
    num_aminos = len(amino_acid_list)
    amino_acid_class_frequencies = amino_acid_type_frequency(amino_acid_list)
    amino_class_dataframe = pd.DataFrame(amino_acid_class_frequencies)
    fig, ax = plt.subplots()
    amino_class_dataframe.T.plot(ax=ax, kind='bar', xlabel='Amino amino acid class', ylabel='Amino acid count',
                                 color=['Red', 'Orange', 'Yellow', 'Green', 'Blue', 'Indigo'])
    # amino_class_dataframe.plot(ax=ax, kind='bar',xlabel='Amino amino acid class', ylabel='Amino acid count', color={'Polar Uncharged Side Chains':'Green', '1':'Blue', '2':'Red', '3':'Yellow', '4':'Purple', '5':'Orange'})
    ax.set_xticks(np.arange(len(amino_acid_class_frequencies)))
    # labels = ['Polar Uncharged Side Chains', 'Charged Side Chains','Positively Charged Sidechains', 'Negatively Charged Side Chains','Special Case Side Chains','Hydrophobic Side Chains']
    # ax.set_xticklabels(labels)
    ax = plt.gca()
    # ax.set_xlim([xmin, xmax])
    ax.set_xlim(-0.25, 0.25)
    ax.set_xticks([-0.25, -0.125, 0, 0.125, 0.25])
    ymin = 0
    ymax = num_aminos
    ax.set_ylim([ymin, ymax])
    ax.legend(["Polar Uncharged -R Amino Acids", "Total Charged -R Amino Acids", "Positively Charged -R Amino Acids",
               "Negatively Charged -R Amino Acids", "Special Case -R Amino Acids", "Hydrophobic -R Amino Acids"]);

    plt.show()


st.write("""
# DNA Transcription and Translation Tool

This tool Transcribes DNA to mRNA and Translates DNA strands into peptides.
""")
# Input DNA Sequence A and B
st.sidebar.header('Input sequence 1 and Sequence 2')
SeqA = 'AAATTT'
SeqB = 'ATAATT'

col1, col2 = st.columns(2)

SeqA = st.text_input('Enter DNA sequence to Transcribe')
st.write('DNA sequence transcribed to RNA is:')
st.markdown(transcription(SeqA))



SeqB = st.text_input('Enter DNA sequence to Translate')
st.write('DNA sequence translated to peptides is :',)
st.markdown(translate_seq(SeqA))

st.markdown(amino_acid_list)

Code = st.text_input('Enter DNA string to get amino acid frequencies')
st.write('Amino acid frequency in this DNA sequence is:',)
peptide_frequency = peptide_frequency_full_name(Code)
st.sidebar.header('Input DNA sequence code to retreive sequence')
st.write(peptide_frequency)

st.subheader('Amino Acid type frequencies:',)
amino_acid_list = translate_seq_full_name(SeqA)

st.pyplot(plot_amino_class_frequency(amino_acid_list))



st.write('Serine, threonine, glutamine, and asparagine are polar but neutral (uncharged) amino acids. \
These side chains can form multiple hydrogen bonds, so they prefer to project into the aqueous phase. \
If they are on the inside of the protein they are hydrogen-bonded to other buried polar groups.')