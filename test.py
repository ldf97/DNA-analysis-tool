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
# DNA Sequence Analysis and Alignment Application

This app aligns two sequences
""")


