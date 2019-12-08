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
import pysam    
import numpy as np
import xlsxwriter 
import matplotlib.pyplot as plt
from AStra import AStra

# Inputs
genomeFastaFile = 'hg19.ucsc.fa'
outputDirectory = './AStraResults'
BamList = ['file1.bam','file2.bam','file3.bam','file4.bam']

# Output sheet
workbook = xlsxwriter.Workbook(outputDirectory + '/aneuploidySpectrum.xlsx') 
worksheet = workbook.add_worksheet() 

## Fields of the output excel file:
#1 ploidy number
#2 copy number reference
#3 CE
#4 CS 
#5 RD-median
#6 RD-median/CN ref
#7-17 aneuploidy spectrum
#18 ploidy state
#19 ploidy model
#20-25: CE for m1-m6

# Parameters of AStra's plots
CNmax = 8 
delCN = 1.5
ampCN = 2.5
histogramBins = 200

i = 1
for bam in BamList:
	print(bam)
	bamindexname = bam + ".bai"
	#
	if os.path.isfile(bamindexname):
		print("BAM index present...")
	else:
		print("No BAM index is available, indexing it...")
		pysam.index(args.bam)
	#
	x = AStra(mybam,  genomeFastaFile, outputDirectory)
	x.runAStra()	
	x.saveGenomeAneuploidy(CNmax,delCN,ampCN)
	x.saveGenomeAneuploidy2(CNmax,delCN,ampCN)
	x.saveGenomeAneuploidy3(CNmax,delCN,ampCN)
	x.saveGenomeAneuploidy4()
	x.saveGenome()
	x.saveHistogram(histogramBins)
	#
	spectrum = x.getAneuploidySpectrum()
	worksheet.write_row('A'+str(i),spectrum)
	i = i+1
	#
	del x
####
workbook.close() 





