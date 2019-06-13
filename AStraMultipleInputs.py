#!/usr/bin/env python3

#Aneuploidy Spectrum Analysis as a Primer for Copy Number Studies of Cancer Cells
__Name__		= "AStra"
__Author__      = "Ahmed Khalil"
__Email__       = "ahmed.khalil.bioinformatics@gmail.com"
__URL__			= "https://github.com/AISKhalil/AStra"
__Software__    = "Python 3"

"""AStra (Aneuploid Spectrum (detection) through read depth analysis) is 
	a Python-based software for capturing the aneuploidy spectrum of cancer genomes 
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

## Aneuploidy-Spectrum with fields:
# 1 nearest ploidy
# 2 number of reads
# 3-12: CE for model1-model10
# 13: HS
# 14: CN State0 percentageDS
# 15: Median error
# 16: Median correction factor
# 17-26: Ploidy-spectrum (percentages of genome-segments per CN-state)
# 27:36: number of segments per each CN-state
workbook = xlsxwriter.Workbook(outputDirectory + '/aneuploidySpectrum.xlsx') 
worksheet = workbook.add_worksheet() 

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
	x = AStra(bam,  genomeFastaFile, outputDirectory)
	x.ploidyEstimatorPipeline()	
	x.saveGenome()
	x.saveHistogram(200)
	#
	spectrum = x.ploidySpectrum()
	worksheet.write_row('A'+str(i),spectrum)
	i = i+1
	#
	del x
####
workbook.close() 





