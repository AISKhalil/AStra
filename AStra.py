#!/usr/bin/env python3

#Aneuploidy Spectrum Analysis as a Primer for Copy Number Studies of Cancer Cells
__Name__		= "AStra"
__Author__      = "Ahmed Khalil"
__Email__       = "ahmed.khalil.bioinformatics@gmail.com"
__URL__			= "https://github.com/AISKhalil/AStra"
__Software__    = "Python 3"

"""
AStra (Aneuploid Spectrum (detection) through read depth analysis) is 
a Python-based software for capturing the aneuploidy spectrum of cancer genomes 
without explicit assumption of the ploidy level of the input sequencing data.

Input data: BAM/SAM formatted-file.


"""
import warnings
import os
import six
import gc
#
import math
import numpy as np
import scipy as sp
import scipy.signal as signal
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import h5py
import xlsxwriter
#
import pysam
import pysamstats
#
import ruptures as rpt
from ruptures.utils import pairwise
from itertools import cycle
#
from textwrap import dedent
import pickle
import numexpr
import math
import logging
log = logging.getLogger(__name__)





class AStra(object):
	"""Base class to estimate the aneuploidy spectrum. """
	#---Attributes---#
	inputFile  = ''		      	  #Input bam-formatted file
	genomeFastaFile = ''	      #Genome sequence folder
	outputFolder = '.'	          #Output folder
	#
	binSize = 200000              #bin-size
	minMAPQ = 1                   #minimum MAPQ of a read to be kept for downstream analysis
	aneuploidySpectrumMethod = 1  #the method used for computing the aneuploidy spectrum- 1:segment-wise, 2:RD-wise
	#
	readCounts = dict()           #RD-signal
	gcPercent  = dict()           #GC-Percentage/bin [0-100%] used for filtering-out low-mappability, centromeres, and telomeres. 
	FIndices   = dict()	          #Indices of the filtered-indices in the chromosome coordinates.
	chrLengths = dict()	          #Dict of chromosome lengths (in bps).
	chrNames = np.array([])       #Array of chromosome names, from the bam file.
	genomeRD = np.array([])
	readDepthMedian = {}		  #Median of the RD signal.
	#
	chrSegments = dict()          #Segments per chromosome: start, end, width (in bins), and CN.
	#
	centralizationErrors = {}	  #the centralization-errors.	
	copyNumberReference = {}      #Reference copy-number.
	minimumCE = {}				  #Centralization-error corresponding to the copy number reference
	ploidyModel = {}			  #Model that achieve the minimum centralization error.
	ploidySegments = dict()       #The Anueploidy profile.
	#
	ploidyLevel = {}              #The ploidy of the cell.
	ploidyNumber = {} 			  #The ploidy number of the cell.

	
	

	
	def __init__(self, inputFile, fastaFile, outputFolder):
		"""
		Initialization
		"""
		self.inputFile = inputFile         #Input bam-formatted file
		self.genomeFastaFile =fastaFile    #Genome sequence folder
		self.outputFolder  = outputFolder  #Output folder
		#
		if not os.path.exists(self.outputFolder):
		    os.makedirs(self.outputFolder)

			
			


	def readsCounting(self):
		"""
		counting the number of reads per bin using pysam
		"""
		mybam = pysam.AlignmentFile(self.inputFile)
		mybamHeader = mybam.header
		#
		self.chrNames   = np.array([])
		self.chrLengths = dict()

		for i in range(0,len(mybamHeader['SQ'])):
			chrName = mybamHeader['SQ'][i]['SN']
			#print(chrName)
			chrNameList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
			if(chrName in chrNameList):
				self.chrNames = np.append(self.chrNames, chrName)
				self.chrLengths[chrName]	= mybam.lengths[i]

		#
		for chrom in self.chrNames:
			print('Reading ' + chrom + ' ...')
			coverageObject = pysamstats.stat_coverage_binned(self.inputFile, self.genomeFastaFile, chrom= chrom, window_size=self.binSize, window_offset=0, min_mapq=self.minMAPQ, no_dup = True)
			#
			coverageCount  = np.array([]) 
			gcRatio       = np.array([])
			for rec in coverageObject: 
				#rec format is "rec: {'pos': 100000, 'gc': 0, 'chrom': 'chr21', 'reads_pp': 0, 'reads_all': 0}"
				coverageCount = np.append(coverageCount, rec['reads_all'])
				gcRatio = np.append(gcRatio, rec['gc'])
			#
			self.readCounts[chrom] = coverageCount
			self.gcPercent[chrom]  = gcRatio			





	def RDfiltering(self, chrom):
		"""
		filtering the RD signal to remove centromeres,telomeres, and low-mappability bins. We filter them based on their GC-percentages.
		"""
		chrRDsignal = self.readCounts[chrom]
		chrGcCount  = self.gcPercent[chrom]
		#
		cond1 = (chrGcCount  > 30)
		cond2 = (chrGcCount  <= 70)
		cond3 = (chrRDsignal != 0)
		mask  = np.bitwise_and(cond1, cond2)
		mask  = np.bitwise_and(mask, cond3)
		#
		chrFIndex    = np.where(mask == True)
		self.FIndices[chrom] = chrFIndex[0]
		#
		chrFRDsignal = chrRDsignal[mask]
		return chrFRDsignal




	def RDsegmentation(self, chrom):
		"""
		segmentation of the RD signal using the change-points detection
		"""
		chrLength   = self.chrLengths[chrom]
		binSize    	= self.binSize
		chrLengthInBins = math.ceil(chrLength/binSize)
		#
		RDsignal = self.readCounts[chrom]

		# Ruptures segmentations
		model = "rbf" #Kernelized mean change
		algo = rpt.Pelt(model=model, min_size= 3, jump=1).fit(RDsignal)
		my_bkps = algo.predict(pen = 5)

		#Segments
		bkps = np.array(my_bkps)
		segmentsStart = np.concatenate(([0], bkps[0:-1])) #In filtered-coordinates
		segmentsEnd   = bkps - 1 #In filtered-coordinates
		segmentsWidth = segmentsEnd - segmentsStart + 1
		noSegments    = len(segmentsWidth)
		#
		segmentsRD    = np.zeros(noSegments)
		for i in range(0,noSegments):
			segmentData = RDsignal[segmentsStart[i]:segmentsEnd[i]+1]
			segmentsRD[i] = np.median(segmentData)

		#Concatenation
		segmentsStart = np.reshape(segmentsStart,(noSegments,1))
		segmentsEnd   = np.reshape(segmentsEnd,(noSegments,1))
		segmentsWidth = np.reshape(segmentsWidth,(noSegments,1))
		segmentsRD    = np.reshape(segmentsRD,(noSegments,1))
		segmentsData  = np.concatenate((segmentsStart, segmentsEnd, segmentsWidth, segmentsRD), axis=1)

		return segmentsData





	def computeMultimodalProbability(self, modelNumber, x):
		"""
		create the multimodal probability.
		"""
		sigma = 0.5/3 #standard deviation
		a1 = 1
		#
		if(modelNumber == 'model1'):
			c = 1.0/(a1)			
			RDDist = a1*stats.norm.pdf(x,2,sigma) 
		elif(modelNumber == 'model2'):
			c = 1.0/(a1)
			RDDist = a1*stats.norm.pdf(x,3,sigma)
		elif(modelNumber == 'model3'):
			c = 1.0/(a1)
			RDDist = a1*stats.norm.pdf(x,4,sigma) 
		elif(modelNumber == 'model4'):
			c = 1.0/(a1 + a1)
			RDDist = a1*stats.norm.pdf(x,2,sigma) + a1*stats.norm.pdf(x,3,sigma) 
		elif(modelNumber == 'model5'):
			c = 1.0/(a1 + a1)
			RDDist = a1*stats.norm.pdf(x,3,sigma) + a1*stats.norm.pdf(x,4,sigma) 
		elif(modelNumber == 'model6'):
			c = 1.0/(a1 + a1 + a1)
			RDDist = a1*stats.norm.pdf(x,2,sigma) + a1*stats.norm.pdf(x,3,sigma) + a1*stats.norm.pdf(x,4,sigma) 
		elif(modelNumber == 'model7'):
			c = 1.0/(a1 + a1 + a1 + a1 + a1 + a1 + a1 + a1 + a1)
			RDDist = a1*stats.norm.pdf(x,2,sigma) + a1*stats.norm.pdf(x,3,sigma) + a1*stats.norm.pdf(x,4,sigma) + a1*stats.norm.pdf(x,5,sigma) + a1*stats.norm.pdf(x,6,sigma) + a1*stats.norm.pdf(x,7,sigma) + a1*stats.norm.pdf(x,8,sigma) + a1*stats.norm.pdf(x,9,sigma) + a1*stats.norm.pdf(x,10,sigma) 
		#
		prob = c*RDDist
		#
		return prob





	def estimateCopyNumberReference(self, genomeFData, ploidyModel):
		"""
		estimate the copy number reference based on the assumed ploidy model
		"""
		##
		genomeFDataLen = len(genomeFData)
		genomeFBData = np.array([]) 
		#
		n = 1
		if(n==1):
			genomeFBData = genomeFData
		else:
			i = 0
			while ((i+n) < genomeFDataLen):
				genomeFBData = np.append(genomeFBData, np.mean(genomeFData[i:i+n]))
				i = i + n
		#
		minValue = min(genomeFBData) 
		maxValue = max(genomeFBData)
		#
		f,l = np.histogram(genomeFBData, 2000, (minValue, maxValue))
		l = l[:-1] + np.diff(l)/2.0
		#
		scanInterval = l
		noScanItems = len(scanInterval)
		inIntervalSum = np.zeros(noScanItems)
		inIntervalSum = inIntervalSum.reshape(-1,1)
		for j in range(noScanItems):
			CN = scanInterval[j] #Corresponding RD of CN = 2
			nScanInterval = scanInterval*2.0/CN
			correlationData = self.computeMultimodalProbability(ploidyModel, nScanInterval)*f
			inIntervalSum[j] = sum(correlationData)  
		#
		chrCNIndex = np.where(inIntervalSum == max(inIntervalSum))[0]
		arr = scanInterval[chrCNIndex]
		estimatedCNReference = np.mean(arr)
		#
		return estimatedCNReference





	def segmentsMerging(self, copyNumberReference):
		"""
		estimate the chromosomal aneuploidy of each chromosome
		"""
		# RDs & widths of genome segments
		mergedChrSegments = dict()
		for chrom in self.chrSegments.keys():

			##########################################################
			#-------------------- Segmentation ----------------------#
			# Chromosome-length
			chrLength       = self.chrLengths[chrom]
			binSize    	    = self.binSize
			chrLengthInBins = math.ceil(chrLength/binSize)

			# RD-signal
			chrRDSignal = self.readCounts[chrom]
			chrSegments = self.chrSegments[chrom]

			chrSegmentsStart = chrSegments[:,0].astype(int)
			chrSegmentsEnd   = chrSegments[:,1].astype(int)
			chrSegmentsWidth = chrSegments[:,2]
			chrSegmentsRD    = chrSegments[:,3]

			# CNs of chromosome segments
			chrSegmentsCN = chrSegmentsRD *2/copyNumberReference
			noSegments = len(chrSegmentsCN)
			#########################################################
			#########################################################
			##Merging subsequent long segments of same CN state
			chrSegmentsStart0 = chrSegmentsStart
			chrSegmentsEnd0   = chrSegmentsEnd
			chrSegmentsCN0    = chrSegmentsCN

			# states of gneome segments
			chrSegmentsStates = np.round(chrSegmentsCN0)

			# chromosome Aneuploidy
			changePointsIndices      = np.where(chrSegmentsStates[:-1] != chrSegmentsStates[1:])[0]
			chrNSegmentsStartIndices = np.append(0, changePointsIndices+1)			
			chrNSegmentsEndIndices   = np.append(changePointsIndices,len(chrSegmentsStates)-1)			

			# New-segments
			chrSegmentsStart1 = chrSegmentsStart0[chrNSegmentsStartIndices]
			chrSegmentsEnd1   = chrSegmentsEnd0[chrNSegmentsEndIndices]
			noSegments1       = len(chrSegmentsStart1)
			#
			chrSegmentsCN1    = []
			chrSegmentsRD1    = []
			for i in range(0, noSegments1):
				newSegmentRD      = np.median(chrRDSignal[chrSegmentsStart1[i]:chrSegmentsEnd1[i]+1]) 
				newSegmentCN      = newSegmentRD *2/copyNumberReference	
				chrSegmentsRD1    = np.append(chrSegmentsRD1, newSegmentRD)			
				chrSegmentsCN1    = np.append(chrSegmentsCN1, newSegmentCN)

			#######################################################
			##------------------------ Output  ------------------##
			chrSegmentsWidth1 = chrSegmentsEnd1 - chrSegmentsStart1 + 1
			segmentsStart     = np.reshape(chrSegmentsStart1,(noSegments1,1))
			segmentsEnd       = np.reshape(chrSegmentsEnd1,(noSegments1,1))
			segmentsWidth     = np.reshape(chrSegmentsWidth1,(noSegments1,1))
			segmentsRD        = np.reshape(chrSegmentsRD1,(noSegments1,1))
			segmentsCN        = np.reshape(chrSegmentsCN1,(noSegments1,1))

			segmentsData   = np.concatenate((segmentsStart, segmentsEnd, segmentsWidth, segmentsRD, segmentsCN), axis=1)
			mergedChrSegments[chrom] = segmentsData

		return mergedChrSegments




	def segmentsMerging2(self, copyNumberReference):
		"""
		estimate the chromosomal aneuploidy of each chromosome
		"""
		# RDs & widths of genome segments
		mergedChrSegments = dict()
		for chrom in self.chrSegments.keys():

			##########################################################
			#-------------------- Segmentation ----------------------#
			# Chromosome-length
			chrLength       = self.chrLengths[chrom]
			binSize    	    = self.binSize
			chrLengthInBins = math.ceil(chrLength/binSize)

			# RD-signal
			chrRDSignal = self.readCounts[chrom]
			chrSegments = self.chrSegments[chrom]

			chrSegmentsStart = chrSegments[:,0].astype(int)
			chrSegmentsEnd   = chrSegments[:,1].astype(int)
			chrSegmentsWidth = chrSegments[:,2]
			chrSegmentsRD    = chrSegments[:,3]

			# CNs of chromosome segments
			chrSegmentsCN = chrSegmentsRD *2/copyNumberReference
			noSegments = len(chrSegmentsCN)

			#######################################################
			##------------------------ Output  ------------------##
			segmentsStart     = np.reshape(chrSegmentsStart,(noSegments,1))
			segmentsEnd       = np.reshape(chrSegmentsEnd,(noSegments,1))
			segmentsWidth     = np.reshape(chrSegmentsWidth,(noSegments,1))
			segmentsRD        = np.reshape(chrSegmentsRD,(noSegments,1))
			segmentsCN        = np.reshape(chrSegmentsCN,(noSegments,1))

			segmentsData   = np.concatenate((segmentsStart, segmentsEnd, segmentsWidth, segmentsRD, segmentsCN), axis=1)
			mergedChrSegments[chrom] = segmentsData

		return mergedChrSegments




	def computeCentralizationError(self, estimatedCNReference, chrSegments):
		"""
		compute the centralization error using the estimated copy number reference from the genomic segments.
		"""
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		for chrom in chrSegments.keys():
				chrSegmentsData  = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])


		minWidthMask        = (genomeSegmentsWidth > 0)
		genomeSegmentsWidth = genomeSegmentsWidth[minWidthMask]		
		genomeSegmentsRD    = genomeSegmentsRD[minWidthMask]				

		# CNs of genome segments
		genomeSegmentsCN = genomeSegmentsRD *2/estimatedCNReference

		# states of gneome segments
		genomeSegmentsStates = np.round(genomeSegmentsCN)

		##############################################################################
		CEDict = dict()
		for i in range(0,20):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsCNi       = genomeSegmentsCN[CNiMask]
			genomeSegmentsStatesCNi = genomeSegmentsStates[CNiMask]
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			if(len(genomeSegmentsCNi) > 0):
				genomeSegmentsErrorCNi = np.abs(genomeSegmentsCNi - genomeSegmentsStatesCNi)
				centralizationErrorCNi = np.sum(genomeSegmentsErrorCNi*genomeSegmentsWidthCNi)
			else:
				centralizationErrorCNi = 0.0
			#print(centralizationErrorCNi)
			CEDict[i] = centralizationErrorCNi
		#-------------- Centralization-Error --------------#		
		CWDict = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1, 12:1, 13:1, 14:1, 15:1, 16:1, 17:1, 18:1, 19:1, 20:1}
		centralizationError = 0
		for i in CEDict.keys():
			centralizationError = centralizationError + CEDict[i]*CWDict[i]
		#
		return centralizationError




	def computeCentralizationError2(self, estimatedCNReference):
		"""
		compute the centralization error using the estimated copy number reference from the RD signal.
		"""
		print('CN = ' + str(estimatedCNReference))
		##
		genomeRD  = np.array([]) 
		for chrom in self.readCounts.keys():
			genomeRD = np.append(genomeRD, self.readCounts[chrom])

		# CNs of genome segments
		genomeCN = genomeRD*2/estimatedCNReference

		# states of gneome segments
		genomeStates = np.round(genomeCN)

		##############################################################################
		##############################################################################
		CEDict = dict()
		for i in range(0,20):
			CNiMask = (genomeStates == i)
			genomeSegmentsCNi       = genomeCN[CNiMask]
			genomeSegmentsStatesCNi = genomeStates[CNiMask]
			if(len(genomeSegmentsCNi) > 0):
				genomeSegmentsErrorCNi = np.abs(genomeSegmentsCNi - genomeSegmentsStatesCNi)
				centralizationErrorCNi = np.sum(genomeSegmentsErrorCNi)
			else:
				centralizationErrorCNi = 0.0
			#
			CEDict[i] = centralizationErrorCNi
		#-------------- Centralization-Error --------------#
		CWDict = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1, 12:1, 13:1, 14:1, 15:1, 16:1, 17:1, 18:1, 19:1, 20:1}		
		centralizationError = 0
		for i in CEDict.keys():
			centralizationError = centralizationError + CEDict[i]*CWDict[i]
		print('Error = ' + str(centralizationError))	

		return centralizationError



		

	def findNearestPloidy(self):
		"""
		find the nearest ploidy level based on the majority
		"""
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		chrSegments = self.ploidySegments
		for chrom in chrSegments.keys():
				chrSegmentsData     = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])

		# CNs of genome segments
		genomeSegmentsCN = genomeSegmentsRD *2/self.copyNumberReference
		genomeSegmentsStates = np.round(genomeSegmentsCN)
		#
		totalSegmentsWidth = np.sum(genomeSegmentsWidth)
		#
		spectrum = []
		for i in range(1,9):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
			iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
			spectrum.append(iSegmentsRatio)
		#
		ploidyIndex = spectrum.index(max(spectrum))
		if(ploidyIndex == 0):
			ploidyLevel = 'haploid'			
		elif(ploidyIndex == 1):			
			ploidyLevel = 'diploid'
		elif(ploidyIndex == 2):
			ploidyLevel = 'triploid'
		elif(ploidyIndex == 3):
			ploidyLevel = 'tetraploid'
		elif(ploidyIndex == 4):			
			ploidyLevel = 'pentaploid'
		elif(ploidyIndex == 5):
			ploidyLevel = 'hexaploid'
		elif(ploidyIndex == 6):
			ploidyLevel = 'heptaploid'
		elif(ploidyIndex == 7):
			ploidyLevel = 'octaploid'
		#
		return ploidyLevel
	



	def findPloidyNumber(self):
		"""
		find the ploidy number as the average of the observed copy number states across the entire genome.
		"""
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		chrSegments = self.ploidySegments
		for chrom in chrSegments.keys():
				chrSegmentsData     = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])

		# CNs of genome segments
		genomeSegmentsCN = genomeSegmentsRD *2/self.copyNumberReference
		genomeSegmentsStates = np.round(genomeSegmentsCN)
		#
		CNiMask = (genomeSegmentsStates != 0)
		genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
		genomeSegmentsStatesCNi = genomeSegmentsStates[CNiMask]
		#
		totalSegmentsWeightedWidth = np.sum(genomeSegmentsWidthCNi*genomeSegmentsStatesCNi)
		totalSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
		
		ploidyNumber = totalSegmentsWeightedWidth/totalSegmentsWidth 
		#
		return ploidyNumber



		
	def computeCS(self):
		"""
		Compute the centralization score.
		"""
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		chrSegments = self.ploidySegments
		for chrom in chrSegments.keys():
				chrSegmentsData     = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])
		
		# CNs of genome segments
		genomeSegmentsCN     = genomeSegmentsRD *2/self.copyNumberReference
		genomeSegmentsStates = np.round(genomeSegmentsCN)

		# Centralization score	
		nearStatesMask  = (abs(genomeSegmentsCN - genomeSegmentsStates) <= 0.25)
		nearStatesWidth = np.sum(genomeSegmentsWidth[nearStatesMask])
		totalSegmentsWidth   = np.sum(genomeSegmentsWidth)
		CS = nearStatesWidth*100/totalSegmentsWidth		
		#	
		return CS

		
		
		
	############################################################################################################################################	
	############################################################################################################################################
	############################################################################################################################################
	def writeAneuploidyResults(self):
		"""
		outputs the chromosomal aneuploidy profile.
		"""
		if not os.path.exists(self.outputFolder):
			os.makedirs(self.outputFolder)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		outputFilePath     = self.outputFolder + '/' + inputFileNameNoExt +  '_Chromosomal_Anueploidy.bed' 

		
		########### Chromosome Anueploidy ##########
		# bed-file header of Chromosomal Anueploidy
		fileID = open(outputFilePath, 'w')
		fileID.write('track type=narrowPeak visibility=3 description="Chromosomal Anueploidy (col5: copy-number state, col7: copy-number)" \n');

		
		# bed-file data
		for chrom in self.ploidySegments.keys():
			#
			chrName = chrom
			#
			chrPloidyData  = self.ploidySegments[chrom]
			chrPloidyStart = chrPloidyData[:,0]*self.binSize
			chrPloidyEnd   = ((chrPloidyData[:,1]+1)*self.binSize)-1
			#
			chrPloidyCN    = chrPloidyData[:,4]
			chrPloidyState = np.round(chrPloidyCN)
			#
			chrPloidyStart = chrPloidyStart.astype(int)
			chrPloidyEnd   = chrPloidyEnd.astype(int)
			chrPloidyState = chrPloidyState.astype(int)
			#
			noSegments = len(chrPloidyStart)
			for i in range(0,noSegments):
				lineData = chrName + '\t' + str(chrPloidyStart[i]) + '\t' + str(min(chrPloidyEnd[i], self.chrLengths[chrName]-1)) + '\t' + '.' + '\t' + str(chrPloidyState[i]) + '\t' + '.' + '\t' + str(chrPloidyCN[i]) + '\t' + '-1' +'\t' + '-1' +'\t' + '-1' + '\n' 
				fileID.write(lineData)
		#
		fileID.close()
		



		
	def writeAneuploidySpectrum(self, ploidyModelsError):
		"""
		outputs the aneuploidy spectrum 
		"""
		if not os.path.exists(self.outputFolder):
			os.makedirs(self.outputFolder)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		outputFilePath    = self.outputFolder + '/' + inputFileNameNoExt +  '_Aneuploidy_Spectrum.bed' 	

		######################################### Main features ######################################
		# text-file of the ploidy-result1
		fileID2 = open(outputFilePath, 'w')
		#
		fileID2.write('Ploidy number is ' + str(self.ploidyNumber) +'\n')
		fileID2.write('\n')
		#
		fileID2.write('Copy number reference = ' + str(self.copyNumberReference) +'\n')
		fileID2.write('Centralization error = ' + str(self.minimumCE) +'\n')
		fileID2.write('Centralization score = ' + str(self.CS) +' %' +'\n')
		fileID2.write('\n')
		#
		fileID2.write('RD signal median = ' + str(self.readDepthMedian) +'\n')
		fileID2.write('Median/CN reference = ' + str(self.readDepthMedian/self.copyNumberReference) +'\n')

		#################################### Aneuploidy spectrum #####################################		
		if(self.aneuploidySpectrumMethod == 1):
			# Segment-wise
			genomeSegmentsRD    = np.array([]) 
			genomeSegmentsWidth = np.array([])		
			chrSegments = self.ploidySegments
			for chrom in chrSegments.keys():
				chrSegmentsData     = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])
		elif(self.aneuploidySpectrumMethod == 2):
			# Segment-wise
			genomeSegmentsRD    = np.array([]) 			
			for chrom in self.readCounts.keys():
				genomeSegmentsRD    = np.append(genomeSegmentsRD, self.readCounts[chrom])
			genomeSegmentsWidth = np.ones(len(genomeSegmentsRD))

		# CNs of genome segments
		genomeSegmentsCN = genomeSegmentsRD *2/self.copyNumberReference
		# states of gneome segments
		genomeSegmentsStates = np.round(genomeSegmentsCN)
		#
		totalSegmentsWidth = np.sum(genomeSegmentsWidth)
		#
		for i in range(0,10):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
			iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
			fileID2.write('CN = ' + str(i) + ' : ' + str(iSegmentsRatio) +' %' + '\n')
		#
		CNiMask = (genomeSegmentsStates >= 10)
		genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
		iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
		iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
		fileID2.write('CN >= 10 : ' + str(iSegmentsRatio) +' %' + '\n')
		fileID2.write('\n')

		
		################################# Additional features #######################################
		fileID2.write('Ploidy state is ' + self.ploidyLevel +'\n')
		fileID2.write('Ploidy model is ' + self.ploidyModel +'\n')		
		#
		i = 1;
		for error in ploidyModelsError:
			fileID2.write('CE (model' + str(i) + ')    = ' + str(error) +'\n')
			i = i+1
		fileID2.write('\n')
		#
		fileID2.write('\n')
		fileID2.write('\n')
		fileID2.close()


		


	def getAneuploidySpectrum(self):
		"""
		return the information of the aneuploidy spectrum
		"""

		spectrum = []
		
		################################# Main features ###################################
		spectrum.append(self.ploidyNumber)
		spectrum.append(self.copyNumberReference)
		spectrum.append(self.minimumCE)
		spectrum.append(self.CS)
		spectrum.append(self.readDepthMedian)
		spectrum.append(self.readDepthMedian/self.copyNumberReference)
		
		############################# Aneuploidy spectrum ################################
		if(self.aneuploidySpectrumMethod == 1):
			# Segment-wise
			genomeSegmentsRD    = np.array([]) 
			genomeSegmentsWidth = np.array([])		
			chrSegments = self.ploidySegments
			for chrom in chrSegments.keys():
				chrSegmentsData     = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])
		elif(self.aneuploidySpectrumMethod == 2):
			# Segment-wise
			genomeSegmentsRD    = np.array([]) 			
			for chrom in self.readCounts.keys():
				genomeSegmentsRD    = np.append(genomeSegmentsRD, self.readCounts[chrom])
			genomeSegmentsWidth = np.ones(len(genomeSegmentsRD))

		#
		genomeSegmentsCN = genomeSegmentsRD *2/self.copyNumberReference
		genomeSegmentsStates = np.round(genomeSegmentsCN)
		totalSegmentsWidth = np.sum(genomeSegmentsWidth)
		#		
		for i in range(0,10):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
			iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
			spectrum.append(iSegmentsRatio)
		#
		CNiMask = (genomeSegmentsStates >= 10)
		genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
		iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
		iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
		spectrum.append(iSegmentsRatio)		
		
		################################# Additional features ############################
		if(self.ploidyLevel == 'haploid'):
			ploidyState = 1		
		elif(self.ploidyLevel == 'diploid'):
			ploidyState = 2		
		elif(self.ploidyLevel == 'triploid'):
			ploidyState = 3
		elif(self.ploidyLevel == 'tetraploid'):
			ploidyState = 4	
		elif(self.ploidyLevel == 'pentaploid'):
			ploidyState = 5		
		elif(self.ploidyLevel == 'hexaploid'):
			ploidyState = 6
		elif(self.ploidyLevel == 'heptaploid'):
			ploidyState = 7				
		elif(self.ploidyLevel == 'octaploid'):
			ploidyState = 8				
		else:
			raise ValueError('Error: genome-wide ploidy level is not from {haploid, diploid, triploid, tetraploid, pentaploid, hexaploid, heptaploid, octaploid}')	
		####	
		if(self.ploidyModel == 'model1'):
			ploidyModel = 1		
		elif(self.ploidyModel == 'model2'):
			ploidyModel = 2		
		elif(self.ploidyModel == 'model3'):
			ploidyModel = 3
		elif(self.ploidyModel == 'model4'):
			ploidyModel = 4	
		elif(self.ploidyModel == 'model5'):
			ploidyModel = 5		
		elif(self.ploidyModel == 'model6'):
			ploidyModel = 6
		elif(self.ploidyModel == 'model7'):
			ploidyModel = 7				
		else:
			raise ValueError('Error: genome-wide ploidy Model is not from model 1 to model 7')	
		####
		spectrum.append(ploidyState)
		spectrum.append(ploidyModel)
		# Models Centralization errors
		for error in self.centralizationErrors:
			spectrum.append(error)

		###############
		return spectrum








		
	############################################################################################################################################	
	############################################################################################################################################
	############################################################################################################################################
	def genomeDataForPlot(self):
		"""
		return the RD signal and segments of the genome.
		"""		
		genomeRDsignal     = np.array([])
		genomeBkps         = list()
		genomeChrStarts    = [0]
		startIndex         = 0
		genomeLengthInBins = 0
		
		###########################
		noChrs = len(self.chrNames)
		chrNames = []
		for i in range(1,noChrs):
			chrNames = chrNames + ['chr'+str(i)]
		chrNames = chrNames + ['chrX']

		chrIndex = 1
		for chrName in chrNames:	
			chrRDSignal    = self.readCounts[chrName]*2/self.copyNumberReference
			#
			chrPloidyData  = self.ploidySegments[chrName]
			chrPloidyStart = chrPloidyData[:,0].astype(int)
			chrPloidyStart = chrPloidyStart + startIndex
			#
			chrLength       = self.chrLengths[chrName]
			binSize    	    = self.binSize
			chrLengthInBins = math.ceil(chrLength/binSize)
			genomeLengthInBins = genomeLengthInBins + chrLengthInBins
			#
			if(chrIndex == 1):
				if(len(chrPloidyStart) > 1):
					bkps = list(chrPloidyStart[1:])	
				else:	
					bkps = []
			elif(chrIndex == noChrs):
					bkps = np.concatenate((chrPloidyStart[0:], [genomeLengthInBins]))
					bkps = list(bkps)	
			else:	
					bkps = list(chrPloidyStart)			
			#
			genomeRDsignal  = np.concatenate((genomeRDsignal,chrRDSignal))
			genomeBkps      = genomeBkps +  bkps
			startIndex      = startIndex + chrLengthInBins
			genomeChrStarts = genomeChrStarts + [startIndex]
			chrIndex        = chrIndex + 1
		####
		out = dict()
		out['a'] = genomeRDsignal
		out['b'] = genomeBkps
		out['c'] = genomeChrStarts
		return out	
		#---------------------------------------------------------------------------------------------------#




	def plotChrAneuploidy(self, chrName):
		"""
		plot the resulted segments of a chromosome.
		"""
		chrRDSignal    = self.readCounts[chrName]*2/self.copyNumberReference
		#
		chrPloidyData  = self.ploidySegments[chrName]
		chrPloidyStart = chrPloidyData[:,0].astype(int)
		#
		chrLength       = self.chrLengths[chrName]
		binSize    	    = self.binSize
		chrLengthInBins = math.ceil(chrLength/binSize)
		if(len(chrPloidyStart) > 1):
			bkps = np.concatenate((chrPloidyStart[1:], [chrLengthInBins]))
			bkps = list(bkps)	
		else:	
			bkps = [chrLengthInBins]
		
		######### Plotting ##########
		CNmax = 10
		signal = np.clip(chrRDSignal,0,CNmax)
		true_chg_pts = bkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#4286f4", "#f44174"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]
		for axe, sig in zip(axarr, signal.T):
			color_cycle = cycle(COLOR_CYCLE)
			# plot s
			axe.plot(range(n_samples), sig, linestyle = 'None',marker = '.', markersize = 2)

			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)

			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				axe.axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)
		#
		plt.title(chrName, fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		plt.show()
		#--------#




	def saveChrAneuploidy(self, chrName):
		"""
		save the resulted segments of a chromosome.
		"""
		chrRDSignal    = self.readCounts[chrName]*2/self.copyNumberReference
		#
		chrPloidyData  = self.ploidySegments[chrName]
		chrPloidyStart = chrPloidyData[:,0].astype(int)
		#
		chrLength       = self.chrLengths[chrName]
		binSize    	    = self.binSize
		chrLengthInBins = math.ceil(chrLength/binSize)
		if(len(chrPloidyStart) > 1):
			bkps = np.concatenate((chrPloidyStart[1:], [chrLengthInBins]))
			bkps = list(bkps)	
		else:	
			bkps = [chrLengthInBins]
		
		######### Plotting ##########
		CNmax = 10
		signal = np.clip(chrRDSignal,0,CNmax)
		true_chg_pts = bkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#4286f4", "#f44174"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]
		for axe, sig in zip(axarr, signal.T):
			color_cycle = cycle(COLOR_CYCLE)
			# plot s
			axe.plot(range(n_samples), sig, linestyle = 'None',marker = '.', markersize = 2)

			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)

			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				axe.axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)
		#
		plt.title(chrName, fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)

		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]		
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_' + chrName + '_AneuploidySegments.png')
		plt.close()
		#---------------------------------------------------------------------------------------------------#




	def plotGenomeAneuploidy(self, CNmax, delTh, ampTh):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']
		genomeStarts   = genomeData['c']
		
		############# Plotting #########
		signal = np.clip(genomeRDsignal,0,CNmax)
		true_chg_pts  = genomeBkps
		true_chg_pts2 = genomeStarts
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#afafaf", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			bkps = [0] + sorted(true_chg_pts)
			startIndices = bkps[0:-1]
			endIndices   = bkps[1:]
			noGenomeSegments = len(startIndices)
			for i in range(0,noGenomeSegments):
				start = startIndices[i]
				end   = endIndices[i]
				segmentCN = np.median(sig[start:end])
				if(segmentCN <= delTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#010101")
				elif(segmentCN >= ampTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#ff0000")					
				else:
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#0000ff")					

			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = sorted(true_chg_pts2)
			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				axe.axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)
		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		plt.show()
		#--------#




	def saveGenomeAneuploidy(self, CNmax, delTh, ampTh):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']
		genomeStarts   = genomeData['c']
		
		############# Plotting #########
		signal = np.clip(genomeRDsignal,0,CNmax)
		true_chg_pts  = genomeBkps
		true_chg_pts2 = genomeStarts
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#afafaf", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			bkps = [0] + sorted(true_chg_pts)
			startIndices = bkps[0:-1]
			endIndices   = bkps[1:]
			noGenomeSegments = len(startIndices)
			for i in range(0,noGenomeSegments):
				start = startIndices[i]
				end   = endIndices[i]
				segmentCN = np.median(sig[start:end])
				if(segmentCN <= delTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#010101")
				elif(segmentCN >= ampTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#ff0000")					
				else:
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#0000ff")					

			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = sorted(true_chg_pts2)
			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				axe.axvspan(max(0, start - 0.5), end - 0.5, facecolor=col, alpha=alpha)
		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_GenomeAneuploidyWithChrsMarkers.png')
		plt.close()
		#---------#	






	def plotGenomeAneuploidy2(self, CNmax, delTh, ampTh):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']
		
		############# Plotting #########
		signal = np.clip(genomeRDsignal,0,CNmax)
		true_chg_pts = genomeBkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#5f5f5f", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)
			startIndices = bkps[0:-1]
			endIndices   = bkps[1:]
			noGenomeSegments = len(startIndices)
			for i in range(0,noGenomeSegments):
				start = startIndices[i]
				end   = endIndices[i]
				segmentCN = np.median(sig[start:end])
				if(segmentCN <= delTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#010101")
				elif(segmentCN >= ampTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#ff0000")					
				else:
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#0000ff")					
		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])		
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		plt.show()
		#--------#



	def saveGenomeAneuploidy2(self, CNmax, delTh, ampTh):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']
		
		############# Plotting #########
		signal = np.clip(genomeRDsignal,0,CNmax)
		true_chg_pts = genomeBkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#5f5f5f", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)
			startIndices = bkps[0:-1]
			endIndices   = bkps[1:]
			noGenomeSegments = len(startIndices)
			for i in range(0,noGenomeSegments):
				start = startIndices[i]
				end   = endIndices[i]
				segmentCN = np.median(sig[start:end])
				if(segmentCN <= delTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#010101")
				elif(segmentCN >= ampTh):
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#ff0000")					
				else:
					axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = "#0000ff")					
		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_GenomeAneuploidy.png')
		plt.close()
		#---------#	



	def plotGenomeAneuploidy3(self, CNmax, delTh, ampTh):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']
		
		############# Plotting #########
		CMin  = 0.5
		signal = np.clip(genomeRDsignal,CMin,CNmax)
		true_chg_pts = genomeBkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#5f5f5f", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)
			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				segmentCN = np.median(sig[start:end])
				if(segmentCN <= delTh):
					axe.plot(list(range(start, end)), np.log2(sig[start:end]/2), linestyle = 'None',marker = '.', markersize = 2, color = "#010101")
				elif(segmentCN >= ampTh):
					axe.plot(list(range(start, end)), np.log2(sig[start:end]/2), linestyle = 'None',marker = '.', markersize = 2, color = "#ff0000")					
				else:
					axe.plot(list(range(start, end)), np.log2(sig[start:end]/2), linestyle = 'None',marker = '.', markersize = 2, color = "#0000ff")					
		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('log2(CN/2)', fontweight="bold", fontsize=12)
		plt.show()
		#--------#




	def saveGenomeAneuploidy3(self, CNmax, delTh, ampTh):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']
		
		############# Plotting #########
		CMin  = 0.5
		signal = np.clip(genomeRDsignal,CMin,CNmax)
		true_chg_pts = genomeBkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#5f5f5f", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)
			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				segmentCN = np.median(sig[start:end])
				if(segmentCN <= delTh):
					axe.plot(list(range(start, end)), np.log2(sig[start:end]/2), linestyle = 'None',marker = '.', markersize = 2, color = "#010101")
				elif(segmentCN >= ampTh):
					axe.plot(list(range(start, end)), np.log2(sig[start:end]/2), linestyle = 'None',marker = '.', markersize = 2, color = "#ff0000")					
				else:
					axe.plot(list(range(start, end)), np.log2(sig[start:end]/2), linestyle = 'None',marker = '.', markersize = 2, color = "#0000ff")					
		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('log2(CN/2)', fontweight="bold", fontsize=12)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_GenomeAneuploidy_logScale.png')
		plt.close()
		#---------#		




	def plotGenomeAneuploidy4(self):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']

		############# Plotting #########
		CNmax = 8
		signal = np.clip(genomeRDsignal,0,CNmax)
		true_chg_pts = genomeBkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#5f5f5f", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)
			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = col)


		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])		
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		plt.show()
		#--------#

		


	def saveGenomeAneuploidy4(self):
		"""
		plot the resulted segments of the genome.
		"""
		genomeData     = self.genomeDataForPlot()
		genomeRDsignal = genomeData['a']
		genomeBkps     = genomeData['b']

		############# Plotting #########
		CNmax = 8
		signal = np.clip(genomeRDsignal,0,CNmax)
		true_chg_pts = genomeBkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#5f5f5f", "#010101"]
		figsize = (10, 3 * n_features)  # figure size
		alpha = 0.2  # transparency of the colored background
		#
		fig, axarr = plt.subplots(n_features, figsize=figsize, sharex=True)
		if n_features == 1:
			axarr = [axarr]

		for axe, sig in zip(axarr, signal.T):
			#
			color_cycle = cycle(COLOR_CYCLE)
			# color each (true) regime
			bkps = [0] + sorted(true_chg_pts)
			for (start, end), col in zip(pairwise(bkps), color_cycle): 
				axe.plot(list(range(start, end)), sig[start:end], linestyle = 'None',marker = '.', markersize = 2, color = col)


		#
		plt.xlim([0, len(signal)])
		plt.ylim([0, CNmax])
		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_GenomeAneuploidy_Segments.png')
		plt.close()
		#---------#	
		

	def plotGenome(self):
		"""
		plot the RD-signal of the genome.
		"""
		genomeRD  = np.array([]) 
		for chrom in self.readCounts.keys():
			genomeRD = np.append(genomeRD, self.readCounts[chrom])
		#
		genomeCN = genomeRD*2/self.copyNumberReference
		######### Plotting ##########
		CNmax = 8
		signal = np.clip(genomeCN,0,CNmax)
		#
		plt.figure()
		plt.plot(signal,'r', linestyle = 'None',marker = '.', markersize = 1)
		plt.plot([0,len(signal)],[1, 1],'k--')				
		plt.plot([0,len(signal)],[2, 2],'k--')
		plt.plot([0,len(signal)],[3, 3],'k--')				
		plt.plot([0,len(signal)],[4, 4],'k--')
		plt.plot([0,len(signal)],[5, 5],'k--')

		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)
		plt.show()
		#--------#


	def saveGenome(self):
		"""
		plot the RD-signal of the genome.
		"""
		genomeRD  = np.array([]) 
		for chrom in self.readCounts.keys():
			genomeRD = np.append(genomeRD, self.readCounts[chrom])
		#
		genomeCN = genomeRD*2/self.copyNumberReference
		######### Plotting ##########
		CNmax = 8
		signal = np.clip(genomeCN,0,CNmax)
		#
		plt.figure()
		plt.plot(signal,'r', linestyle = 'None',marker = '.', markersize = 1)
		plt.plot([0,len(signal)],[1, 1],'k--')				
		plt.plot([0,len(signal)],[2, 2],'k--')
		plt.plot([0,len(signal)],[3, 3],'k--')				
		plt.plot([0,len(signal)],[4, 4],'k--')
		plt.plot([0,len(signal)],[5, 5],'k--')

		plt.title('Genome', fontweight="bold", fontsize=12)
		plt.xlabel('Bin number', fontweight="bold", fontsize=12)
		plt.ylabel('Copy number', fontweight="bold", fontsize=12)


		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_GenomeRD.png')
		plt.close()
		#---------#

		


	def plotHistogram(self, nbins):
		"""
		plot the histogram of the genome RD-signal.
		"""
		genomeRD  = self.genomeRD
		genomeRD = genomeRD[genomeRD < (self.copyNumberReference*4)]

		#################### Plotting ##################
		plt.figure()
		n = plt.hist(genomeRD,nbins)
		plt.title('Genome RD-frequency')
		plt.xlabel('Reads/bin')
		plt.ylabel('Frequency')
		####
		plt.plot([(self.readDepthMedian), (self.readDepthMedian)],[0,max(n[0])], 'r', lineWidth = 2)
		####
		plt.plot([self.copyNumberReference, self.copyNumberReference],[0,max(n[0])],'k--')
		####
		plt.plot([(0.5*self.copyNumberReference), (0.5*self.copyNumberReference)],[0,max(n[0])],'k--')				
		plt.plot([(1.5*self.copyNumberReference), (1.5*self.copyNumberReference)],[0,max(n[0])],'k--')
		plt.plot([(2*self.copyNumberReference), (2*self.copyNumberReference)],[0,max(n[0])],'k--')
		plt.plot([(2.5*self.copyNumberReference), (2.5*self.copyNumberReference)],[0,max(n[0])],'k--')
		plt.plot([(3*self.copyNumberReference), (3*self.copyNumberReference)],[0,max(n[0])],'k--')		
		####
		plt.show()
		#--------#



	def saveHistogram(self,nbins):
		"""
		plot the histogram of the genome RD-signal.
		"""
		genomeRD  = self.genomeRD
		genomeRD = genomeRD[genomeRD < (self.copyNumberReference*4)]

		#################### Plotting ##################
		plt.figure()
		n = plt.hist(genomeRD,nbins)
		plt.title('Genome RD-frequency')
		plt.xlabel('Reads/bin')
		plt.ylabel('Frequency')
		####
		plt.plot([(self.readDepthMedian), (self.readDepthMedian)],[0,max(n[0])], 'r', lineWidth = 2)
		####
		plt.plot([self.copyNumberReference, self.copyNumberReference],[0,max(n[0])],'k--')
		####
		plt.plot([(0.5*self.copyNumberReference), (0.5*self.copyNumberReference)],[0,max(n[0])],'k--')				
		plt.plot([(1.5*self.copyNumberReference), (1.5*self.copyNumberReference)],[0,max(n[0])],'k--')
		plt.plot([(2*self.copyNumberReference), (2*self.copyNumberReference)],[0,max(n[0])],'k--')
		plt.plot([(2.5*self.copyNumberReference), (2.5*self.copyNumberReference)],[0,max(n[0])],'k--')
		plt.plot([(3*self.copyNumberReference), (3*self.copyNumberReference)],[0,max(n[0])],'k--')		
		####
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_' + str(nbins) + 'bin_GenomeHistogram.png')
		plt.close()		
		#--------#






	###################################################################################################################################	
	###################################################################################################################################
	###################################################################################################################################
	def runAStra(self):
		"""
		the pipeline for estimating the aneuploidy spectrum of the WGS data 
		and generating the output files using our 6 models.	
		"""

		############################## RD-calculation ############################
		self.readsCounting()
		# Filtering & Segmentation
		genomeRD  = np.array([]) 
		for chrom in self.readCounts.keys():
			chrFRDsignal = self.RDfiltering(chrom)
			segmentsData = self.RDsegmentation(chrom) 
			self.chrSegments[chrom] = segmentsData
			genomeRD = np.append(genomeRD, chrFRDsignal)
		#
		self.genomeRD = genomeRD
		self.readDepthMedian = np.median(genomeRD)


		############################# Coarse scanning ############################
		ploidyModels = ['model1', 'model2', 'model3', 'model4', 'model5','model6']
		ploidySegments = dict()
		copyNumberReference = []
		ploidyModelsError   = []
		#
		for ploidy in ploidyModels:

			print('Model = ' + ploidy)

			##################### Level#1: Coarse CN reference ###################
			coarseCNReference = self.estimateCopyNumberReference(genomeRD, ploidy)
			print('Coarse CN reference = ' + str(coarseCNReference))		

			##################### Level#2: Fine CN reference ##################### 
			candidatesCNReference = np.linspace(1.9,2.1,100)*coarseCNReference/2
			#
			fineCopyNumberReference = []
			candidateCNRefError = []
			#
			for estimatedCNReference in candidatesCNReference:
				# copy-number estimation
				fineCopyNumberReference.append(estimatedCNReference)
				finePloidySegments = self.segmentsMerging(estimatedCNReference)

				# centralization-error
				centralizationError  = self.computeCentralizationError(estimatedCNReference, finePloidySegments)					
				candidateCNRefError.append(centralizationError)
			#
			CNRefIndex  = candidateCNRefError.index(min(candidateCNRefError))
			CNRefError  = candidateCNRefError[CNRefIndex]
			CNReference = fineCopyNumberReference[CNRefIndex]
			#
			ploidySegments[ploidy] = self.segmentsMerging(CNReference)
			copyNumberReference.append(CNReference)
			ploidyModelsError.append(CNRefError)
			#
			print(ploidy + ' CN reference = ' + str(CNReference))		
			print(ploidy + ' CE = ' + str(CNRefError))		

		# Final CN reference
		ploidyIndex = ploidyModelsError.index(min(ploidyModelsError))
		finalCE              = ploidyModelsError[ploidyIndex]
		finalCNReference     = copyNumberReference[ploidyIndex]
		finalPloidyModel     = ploidyModels[ploidyIndex]
		finalPloidySegments  = ploidySegments[finalPloidyModel]
		#
		print('\n')
		print('Final CN reference = ' + str(finalCNReference))
		print('Final CE           = ' + str(finalCE))			
		
		######################### Aneuploidy Spectrum ##########################
		self.copyNumberReference  = finalCNReference
		self.ploidySegments       = finalPloidySegments
		self.minimumCE            = finalCE
		self.ploidyModel          = finalPloidyModel
		self.centralizationErrors = ploidyModelsError
		#
		# Ploidy
		self.ploidyLevel    = self.findNearestPloidy()
		self.ploidyNumber   = self.findPloidyNumber()
		self.CS             = self.computeCS()
		print('\n')
		print('Ploidy number is ' + str(self.ploidyNumber))
		print('Centralization score is ' + str(self.CS) + ' %')
		print('Ploidy state is ' + self.ploidyLevel)
		############################ Output Files #############################
		self.writeAneuploidyResults()
		self.writeAneuploidySpectrum(ploidyModelsError)






	###################################################################################################################################	
	###################################################################################################################################
	###################################################################################################################################
	def runAStraLite1(self):
		"""
		the pipeline for estimating the aneuploidy spectrum of the WGS data 
		and generating the output files using only one distribution.	
		"""

		############################## RD-calculation ############################
		self.readsCounting()
		# Filtering & Segmentation
		genomeRD  = np.array([]) 
		for chrom in self.readCounts.keys():
			chrFRDsignal = self.RDfiltering(chrom)
			segmentsData = self.RDsegmentation(chrom) 
			self.chrSegments[chrom] = segmentsData
			genomeRD = np.append(genomeRD, chrFRDsignal)
		#
		self.genomeRD = genomeRD
		self.readDepthMedian = np.median(genomeRD)


		############################# Coarse scanning ############################
		ploidyModels = ['model7']
		ploidySegments = dict()
		copyNumberReference = []
		ploidyModelsError   = []
		#
		for ploidy in ploidyModels:

			print('Model = ' + ploidy)

			##################### Level#1: Coarse CN reference ######################
			coarseCNReference = self.estimateCopyNumberReference(genomeRD, ploidy)
			print('Coarse CN reference = ' + str(coarseCNReference))		

			##################### Level#2: Fine CN reference ##################### 
			fineCopyNumberReference = []
			candidateCNRefError = []
			#
			estimatedCNReference = coarseCNReference
			# copy-number estimation
			fineCopyNumberReference.append(estimatedCNReference)
			finePloidySegments = self.segmentsMerging(estimatedCNReference)
			# centralization-error
			centralizationError  = self.computeCentralizationError(estimatedCNReference, finePloidySegments)					
			candidateCNRefError.append(centralizationError)

			#
			CNRefIndex  = candidateCNRefError.index(min(candidateCNRefError))
			CNRefError  = candidateCNRefError[CNRefIndex]
			CNReference = fineCopyNumberReference[CNRefIndex]
			#
			ploidySegments[ploidy] = self.segmentsMerging(CNReference)
			copyNumberReference.append(CNReference)
			ploidyModelsError.append(CNRefError)
			#
			print(ploidy + ' CN reference = ' + str(CNReference))		
			print(ploidy + ' CE = ' + str(CNRefError))		

		# Final CN reference
		ploidyIndex = ploidyModelsError.index(min(ploidyModelsError))
		finalCE              = ploidyModelsError[ploidyIndex]
		finalCNReference     = copyNumberReference[ploidyIndex]
		finalPloidyModel     = ploidyModels[ploidyIndex]
		finalPloidySegments  = ploidySegments[finalPloidyModel]
		#
		print('\n')
		print('Final CN reference = ' + str(finalCNReference))
		print('Final CE           = ' + str(finalCE))			
		
		######################### Aneuploidy Spectrum ##########################
		self.copyNumberReference  = finalCNReference
		self.ploidySegments       = finalPloidySegments
		self.minimumCE            = finalCE
		self.ploidyModel          = finalPloidyModel
		self.centralizationErrors = ploidyModelsError
		#
		# Ploidy
		self.ploidyLevel    = self.findNearestPloidy()
		self.ploidyNumber   = self.findPloidyNumber()
		self.CS             = self.computeCS()
		print('\n')
		print('Ploidy number is ' + str(self.ploidyNumber))
		print('Centralization score is ' + str(self.CS) + ' %')
		print('Ploidy state is ' + self.ploidyLevel)
		
		############################ Output Files #############################
		self.writeAneuploidyResults()
		self.writeAneuploidySpectrum(ploidyModelsError)





	###################################################################################################################################	
	###################################################################################################################################
	###################################################################################################################################
	def runAStraLite2(self):
		"""
		the pipeline for estimating the aneuploidy spectrum of the WGS data 
		and generating the output files.	
		"""

		############################## RD-calculation ############################
		self.readsCounting()
		# Filtering & Segmentation
		genomeRD  = np.array([]) 
		for chrom in self.readCounts.keys():
			chrFRDsignal = self.RDfiltering(chrom)
			segmentsData = self.RDsegmentation(chrom) 
			self.chrSegments[chrom] = segmentsData
			genomeRD = np.append(genomeRD, chrFRDsignal)
		#
		self.genomeRD = genomeRD
		self.readDepthMedian = np.median(genomeRD)


		##################### Level#1: Coarse CN reference #################
		coarseCNReference = self.readDepthMedian
		print('Coarse CN reference = ' + str(coarseCNReference))		


		##################### Level#2: Fine CN reference ################### 
		candidatesCNReference = np.linspace(0.25,4,400)*coarseCNReference
		#
		fineCopyNumberReference = []
		candidateCNRefError = []
		#
		for estimatedCNReference in candidatesCNReference:
			# copy-number estimation
			fineCopyNumberReference.append(estimatedCNReference)
			finePloidySegments = self.segmentsMerging(estimatedCNReference)
			# centralization-error
			centralizationError  = self.computeCentralizationError(estimatedCNReference, finePloidySegments)					
			candidateCNRefError.append(centralizationError)
		#
		CNRefIndex       = candidateCNRefError.index(min(candidateCNRefError))
		finalCE          = candidateCNRefError[CNRefIndex]
		finalCNReference = fineCopyNumberReference[CNRefIndex]
		#
		finalPloidySegments   = self.segmentsMerging(finalCNReference)
		finalPloidyModel = 'model7'

		# Final CN reference
		print('\n')
		print('Final CN reference = ' + str(finalCNReference))
		print('Final CE           = ' + str(finalCE))			
		

		######################### Aneuploidy Spectrum ##########################
		self.copyNumberReference  = finalCNReference
		self.ploidySegments       = finalPloidySegments
		self.minimumCE            = finalCE
		self.ploidyModel          = finalPloidyModel
		self.centralizationErrors = [finalCE]
		#
		# Ploidy
		self.ploidyLevel    = self.findNearestPloidy()
		self.ploidyNumber   = self.findPloidyNumber()
		self.CS             = self.computeCS()
		print('\n')
		print('Ploidy number is ' + str(self.ploidyNumber))
		print('Centralization score is ' + str(self.CS) + ' %')
		print('Ploidy state is ' + self.ploidyLevel)
		

		############################ Output Files #############################
		ploidyModelsError = finalPloidyModel
		self.writeAneuploidyResults()
		self.writeAneuploidySpectrum(ploidyModelsError)
		#
		out = dict()
		out['CNR'] = fineCopyNumberReference
		out['CE']  = candidateCNRefError
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		outputFilePath     = self.outputFolder + '/' + inputFileNameNoExt +  '_CNRef_CE_Values.xlsx' 
		#
		workbook = xlsxwriter.Workbook(outputFilePath) 
		worksheet = workbook.add_worksheet() 
		worksheet.write_column('A1',fineCopyNumberReference)
		worksheet.write_column('B1',candidateCNRefError)
		workbook.close()

		##########
		return out