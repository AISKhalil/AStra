#!/usr/bin/env python3

#Aneuploidy Spectrum Analysis (Digital Karyotyping) for Rapid Authentication of Cell Lines.
__Name__	= "AStra"
__Author__      = "Ahmed Khalil"
__Email__       = "ahmed.khalil.bioinformatics@gmail.com"
__URL__		= "https://github.com/AISKhalil/AStra"
__Software__    = "Python 3"

"""
AStra (Aneuploid Spectrum (detection) through read depth analysis) is 
a Python-based software for capturing the aneuploidy profile spectrum of cancer genomes 
without explicit assumption of the ploidy level of the input sequencing data.
"""

import sys
import os
import argparse
import pysam    
import numpy as np
import xlsxwriter 
import matplotlib.pyplot as plt
from AStra import AStra

#Making sure you are running a version of python that works with this script.
if sys.version_info[0] != 3 or sys.version_info[1] < 6 or sys.version_info[2] < 5:
    print("This script requires Python version 3.6.5 or higher within major version 3")
    sys.exit(1)

# Input parsing
parser = argparse.ArgumentParser(description='AStra: Aneuploidy Spectrum Analysis as a Primer for Copy Number Studies of Cancer Cells', add_help=True)
parser.add_argument('-b','--bam', dest='bam', metavar='mappingGenome.bam', type=str, help='BAM file used to get allele frequencies', required=True)
parser.add_argument('-f','--fa' , dest='genomeFastaFile', metavar='hg19.ucsc.fa', type=str, help='Fasta file of the reference genome', required=True)
parser.add_argument('-o','--out', dest='outputDirectory', metavar='folder', type=str, help='Folder to keep all output files', required=True)

# Get information from the argparse (arguments)
args = parser.parse_args()
bam = args.bam
genomeFastaFile = args.genomeFastaFile
outputDirectory=args.outputDirectory

# Check if bam index is present; if not, create it
bamindexname = args.bam + ".bai"
if os.path.isfile(bamindexname):
	print("BAM index present...")
else:
	print("No BAM index is available, indexing it...")
	pysam.index(args.bam)

	
# Executing AStra on the input bam 
x = AStra(bam,  genomeFastaFile, outputDirectory)
x.runAStra()

# additional plots of the aneuploidy profile and spectrum
CNmax = 8 
delCN = 1.5
ampCN = 2.5
histogramBins = 200
x.saveGenomeAneuploidy(CNmax, delCN, ampCN)
x.saveGenomeAneuploidy2(CNmax, delCN, ampCN)	
x.saveGenomeAneuploidy3(CNmax, delCN, ampCN)
x.saveGenomeAneuploidy4()
x.saveGenome()
x.saveHistogram(histogramBins)