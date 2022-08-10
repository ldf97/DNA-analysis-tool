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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import py3Dmol
import urllib.request
import sys

from stmol import showmol

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



# rendering proteins
datadir = os.getcwd()
def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None

class Atom(dict):

    def __init__(self, line):
        self["type"] = line[0:6].strip()
        self["idx"] = line[6:11].strip()
        self["name"] = line[12:16].strip()
        self["resname"] = line[17:20].strip()
        self["resid"] = int(int(line[22:26]))
        self["x"] = float(line[30:38])
        self["y"] = float(line[38:46])
        self["z"] = float(line[46:54])
        self["sym"] = line[76:78].strip()

    def __str__(self):
        line = list(" " * 80)

        line[0:6] = self["type"].ljust(6)
        line[6:11] = self["idx"].ljust(5)
        line[12:16] = self["name"].ljust(4)
        line[17:20] = self["resname"].ljust(3)
        line[22:26] = str(self["resid"]).ljust(4)
        line[30:38] = str(self["x"]).rjust(8)
        line[38:46] = str(self["y"]).rjust(8)
        line[46:54] = str(self["z"]).rjust(8)
        line[76:78] = self["sym"].rjust(2)
        return "".join(line) + "\n"

class Molecule(list):
    def __init__(self, file):
        for line in file:
            if "ATOM " in line: #or "HETATM " in line: #this was only way I could get it to skip to where the ATOMS were located in pdb txt file, so commented out 'or 'hetatm' and added space next to 'atom'
                self.append(Atom(line))

    def __str__(self):
        outstr = ""
        for at in self:
            outstr += str(at)

        return outstr

datadir = os.getcwd()
def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None


def download_pdb_and_visualize_protein(pdbcode, datadir=datadir):
    download_pdb(pdbcode, datadir)

    with open(str(pdbcode + ".pdb")) as ifile:
        mol = Molecule(ifile)

    view = py3Dmol.view(width=400, height=300)
    view.addModelsAsFrames(str(mol))
    view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    view.zoomTo()
    #view.show()
    return view

#width=400
#height=300
#showmol(view, height=height, width=width)

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


def plot_nucleotide_frequency(seq):
    num_nucleotides = len(seq)
    nucleotide_frequencies = nucleotide_frequency(seq)
    nucleotides_dataframe = pd.DataFrame(nucleotide_frequencies, index=['0'])
    fig, ax = plt.subplots()
    nucleotides_dataframe.plot(ax=ax, kind='bar', xlabel='Nucleotide', ylabel='Sequence Length')
    ax = plt.gca()
    # ax.set_xlim([xmin, xmax])
    ymin = 0
    ymax = num_nucleotides
    ax.set_ylim([ymin, ymax])

    plt.show()

def nucleotide_frequency(seq):
  """
  Count the nucleotides in a given sequence.
  Return a dictionary.
  """
  tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
  for nuc in seq:
    tmpFreqDict[nuc] += 1
  return tmpFreqDict

def plot_aminoacid_frequency(seq):
    amino_frequencies = peptide_frequency_full_name(seq)
    num_aminoacids = len(amino_frequencies)
    amino_dataframe = pd.DataFrame(amino_frequencies, index=['0'])
    fig, ax = plt.subplots()
    amino_dataframe.plot(ax=ax, kind='bar', xlabel='Amino Acids', ylabel='Total Number of amino acids')
    ax = plt.gca()
    # ax.set_xlim([xmin, xmax])
    ymin = 0
    ymax = num_aminoacids
    ax.set_ylim([ymin, ymax])

    plt.show()
def gc_content(seq):
  """GC Content in a DNA/RNA sequence"""
  GC_content = round((seq.count('C') + seq.count('G')) / len(seq) *100)
  return (GC_content)

def at_content(seq):
  """AT Content in a DNA/RNA sequence"""
  AT_content = round((seq.count('A') + seq.count('T')) / len(seq) *100)
  return (AT_content)
def plot_gc_content(seq):
  GC_content = gc_content(seq)
  AT_content = at_content(seq)
  labels_GC_content = 'GC content', 'AT content'
  sizes_GC_content = [GC_content, AT_content]
  fig1, ax1 = plt.subplots()
  ax1.pie(sizes_GC_content, labels=labels_GC_content, autopct='%1.1f%%',
          shadow=True, startangle=90)
  ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

  plt.show()



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

#amino acid groupings
amino_acids_charged_side_chains = ['Arginine, Arg, R', 'Histidine, His, H', 'Lysine, Lys, K', 'Aspartic Acid, Asp, D', 'Glutamic Acid, Glut, E']
amino_acids_positively_charged = ['Arginine, Arg, R', 'Histidine, His, H', 'Lysine, Lys, K']
amino_acids_negatively_charged = ['Aspartic Acid, Asp, D', 'Glutamic Acid, Glut, E']
amino_acids_with_polar_uncharged_side_chains = ['Serine, Ser, S', 'Threonine, Thr, T', 'Asparagine, Asn, N', 'Glutamine, Gln, Q']
amino_acids_special_cases = ['Cysteine, Cys, C', 'Glycine, Gly, G', 'Proline, Pro, P']
amino_acids_with_hydrophobic_side_chains = ['Alanine, Ala, A', 'Valine, Val, V', 'Isoleucine, Ile, I','Methionine, Met, M','Phenylalanine, Phe, F','Tyrosine, Tyr, Y','Tryptophan, Trp, W']

def amino_acid_type_frequency(amino_acid_list):
  """
  Count the amino acids in a given class.
  Return a counts for each type in this order: polar uncharged side chain, charged side chain,
  positively charged side chain, negatively charged side chain, special cases, hydrophobic side chains
  """
  tmp_amino_acids_charged_side_chains_dict = {'Arginine, Arg, R': 0, 'Histidine, His, H': 0, 'Lysine, Lys, K': 0, 'Aspartic Acid, Asp, D': 0, 'Glutamic Acid, Glut, E':0}
  tmp_amino_acids_positively_charged_dict = {'Arginine, Arg, R':0, 'Histidine, His, H':0, 'Lysine, Lys, K':0}
  tmp_amino_acids_negatively_charged_dict = {'Aspartic Acid, Asp, D':0, 'Glutamic Acid, Glut, E':0}
  tmp_amino_acids_with_polar_uncharged_side_chains_dict = {'Serine, Ser, S':0, 'Threonine, Thr, T':0, 'Asparagine, Asn, N':0, 'Glutamine, Gln, Q':0}
  tmp_amino_acids_special_cases_dict = {'Cysteine, Cys, C':0, 'Glycine, Gly, G':0, 'Proline, Pro, P':0}
  tmp_amino_acids_with_hydrophobic_side_chains_dict = {'Alanine, Ala, A':0, 'Valine, Val, V':0, 'Isoleucine, Ile, I':0, 'Leucine, Leu, L':0, 'Methionine, Met, M':0,'Phenylalanine, Phe, F':0,'Tyrosine, Tyr, Y':0,'Tryptophan, Trp, W':0}
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


  return sum(tmp_amino_acids_with_polar_uncharged_side_chains_dict.values()),\
  sum(tmp_amino_acids_charged_side_chains_dict.values()),\
  sum(tmp_amino_acids_positively_charged_dict.values()),\
  sum(tmp_amino_acids_negatively_charged_dict.values()),\
  sum(tmp_amino_acids_special_cases_dict.values()),\
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


width=400
height=300


pdbcode = st.text_input('Enter PDB code as text to visualize a protein of interest. (e.g. Human insulin 3I40, or Water Structure of Hydrophobic Plant Protein 1CRN)')
if pdbcode:
    view = download_pdb_and_visualize_protein(pdbcode)
else:
    view = None
st.write('PDB is:', pdbcode)
if pdbcode:
    showmol(view, height=height, width=width)
else:
    st.write('Enter a code above to pull a structure file from Protein Database. (e.g. Human insulin 3I40, or Water Structure of Hydrophobic Plant Protein 1CRN)')












