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
	binSize = 100000              #bin-size
	minMAPQ = 1                   #minimum MAPQ of a read to be kept for downstream analysis
	#
	readCounts = dict()           #RD-signal
	gcPercent  = dict()           #GC-Percentage/bin [0-100%] used for filtering-out low-mappability, centromeres, and telomeres. 
	FIndices   = dict()	          #Indices of the filtered-indices in the chromosome coordinates.
	chrLengths = dict()	          #Dict of chromosome lengths (in bps).
	chrNames = np.array([])       #Array of chromosome names, from the bam file.
	genomeRD = np.array([])
	readDepthMedian = {}		  #Median of the RD signal.
	#
	segmentsFrequency = 1000000   #the expected length of segments (in bps)
	chrSegments = dict()          #Segments per chromosome: start, end, width (in bins), and CN.
	#
	centralizationErrors = {}	  #the centralization-errors.	
	copyNumberReference = {}      #Reference copy-number.
	ploidyModel = {}			  #Model that achieve the minimum centralization error.
	ploidySegments = dict()       #the Anueploidy profile.
	#
	ploidyLevel = {}              #the Ploidy of the cell.


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
			if((chrName != 'chrY') and (chrName != 'chrM')):
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

		#Number of expected segments
		avgSize = round(self.segmentsFrequency/self.binSize)

		# Ruptures segmentations
		model = "rbf" #Kernelized mean change
		algo = rpt.Pelt(model=model, min_size= 5, jump=1).fit(RDsignal)
		my_bkps = algo.predict(pen = 3)


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
		a2 = 0.5
		a3 = 0.25
		a4 = 0.125
		a5 = 0.0625
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
			c = 1.0/(a1 + a2 + a3 + a4 + a5)			
			RDDist = a1*stats.norm.pdf(x,2,sigma) + a2*stats.norm.pdf(x,3,sigma) + a3*stats.norm.pdf(x,4,sigma) + a4*stats.norm.pdf(x,5,sigma) + a5*stats.norm.pdf(x,6,sigma)  
		elif(modelNumber == 'model8'):
			c = 1.0/(a2 + a1 + a2 + a3 + a4)
			RDDist = a2*stats.norm.pdf(x,2,sigma) + a1*stats.norm.pdf(x,3,sigma) + a2*stats.norm.pdf(x,4,sigma) + a3*stats.norm.pdf(x,5,sigma) + a4*stats.norm.pdf(x,6,sigma) 
		elif(modelNumber == 'model9'):
			c = 1.0/(a3 + a2 + a1 + a2 + a3)
			RDDist = a3*stats.norm.pdf(x,2,sigma) + a2*stats.norm.pdf(x,3,sigma) + a1*stats.norm.pdf(x,4,sigma) + a2*stats.norm.pdf(x,5,sigma) + a3*stats.norm.pdf(x,6,sigma) 
		elif(modelNumber == 'model10'):
			c = 1.0/(a1 + a1 + a1 + a1 + a1)
			RDDist = a1*stats.norm.pdf(x,2,sigma) + a1*stats.norm.pdf(x,3,sigma) + a1*stats.norm.pdf(x,4,sigma)	+ a1*stats.norm.pdf(x,5,sigma) + a1*stats.norm.pdf(x,6,sigma)			
		prob = c*RDDist
		#
		return prob




	def estimateCopyNumberReference(self, genomeFData, ploidyModel):
		"""
		estimate the copy number reference based on the assumed ploidy model
		"""
		cond1 = (genomeFData <= np.percentile(genomeFData,99))
		cond2 = (genomeFData > np.percentile(genomeFData,1))
		filteringMask  = np.bitwise_and(cond1, cond2)
		genomeFData = genomeFData[filteringMask]
		##
		genomeFDataLen = len(genomeFData)
		genomeFBData = np.array([]) 
		#
		n = round(self.segmentsFrequency/self.binSize)
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
			segmentsStart  = np.reshape(chrSegmentsStart1,(noSegments1,1))
			segmentsEnd    = np.reshape(chrSegmentsEnd1,(noSegments1,1))
			segmentsWidth  = np.reshape(chrSegmentsWidth1,(noSegments1,1))
			segmentsRD     = np.reshape(chrSegmentsRD1,(noSegments1,1))
			segmentsCN     = np.reshape(chrSegmentsCN1,(noSegments1,1))

			segmentsData   = np.concatenate((segmentsStart, segmentsEnd, segmentsWidth, segmentsRD, segmentsCN), axis=1)
			mergedChrSegments[chrom] = segmentsData



		return mergedChrSegments






	def computeCentralizationError(self, estimatedCNReference, chrSegments):
		"""
		compute the centralization error using the estimated copy number reference
		"""
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		for chrom in chrSegments.keys():
				chrSegmentsData  = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])


		#minSegmentWidth = round(self.segmentsFrequency/self.binSize)
		minWidthMask    = (genomeSegmentsWidth > 0)
		genomeSegmentsWidth = genomeSegmentsWidth[minWidthMask]		
		genomeSegmentsRD   = genomeSegmentsRD[minWidthMask]				

		# CNs of genome segments
		genomeSegmentsCN = genomeSegmentsRD *2/estimatedCNReference
		print('CN = ' + str(estimatedCNReference))

		# states of gneome segments
		genomeSegmentsStates = np.round(genomeSegmentsCN)

		##############################################################################
		##############################################################################
		CEDict = dict()
		for i in range(0,10):
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
		####
		CWDict = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1}
		#-------------- Centralization-Error --------------#		
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
		for i in range(2,5):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
			iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
			spectrum.append(iSegmentsRatio)
		#
		ploidyIndex = spectrum.index(max(spectrum))
		if(ploidyIndex == 0):
			ploidyLevel = 'diploid'
		elif(ploidyIndex == 1):
			ploidyLevel = 'triploid'
		elif(ploidyIndex == 2):
			ploidyLevel = 'tetraploid'

		######	
		return ploidyLevel
	







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
		outputFilePath2    = self.outputFolder + '/' + inputFileNameNoExt +  '_Ploidy_Results.bed' 	


		######################################################## Chromosome Anueploidy #########################################################
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
		#------------#


	def writeAneuploidySpectrum(self, ploidyModelsError):
		"""
		outputs the aneuploidy spectrum 
		"""
		if not os.path.exists(self.outputFolder):
			os.makedirs(self.outputFolder)
		#
		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		outputFilePath    = self.outputFolder + '/' + inputFileNameNoExt +  '_Ploidy_Results.bed' 	


		########################################## Ploidy Spectrum ######################################

		# text-file of the ploidy-result1
		fileID2 = open(outputFilePath, 'w')
		#
		fileID2.write('Ploidy-Level = ' + self.ploidyLevel +'\n')
		fileID2.write('Copy number reference = ' + str(self.copyNumberReference) +'\n')
		#
		fileID2.write('\n')
		#
		i = 1;
		for error in ploidyModelsError:
			fileID2.write('CE (model' + str(i) + ')    = ' + str(error) +'\n')
			i = i+1
		fileID2.write('\n')
		fileID2.write('\n')


		################################################################################################
		################################ Final-segments Statistics #####################################		
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
		# states of gneome segments
		genomeSegmentsStates = np.round(genomeSegmentsCN)
		##
		totalSegmentsWidth = np.sum(genomeSegmentsWidth)
		for i in range(0,10):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
			iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
			fileID2.write('CN = ' + str(i) + ' : ' + str(iSegmentsRatio) +' %' +' in '+ str(np.sum(CNiMask)) + ' segments' + '\n')
		##
		fileID2.write('\n')
		fileID2.write('\n')


		################################################################################################
		########################################### Errors #############################################
		##
		nearStatesMask  = (abs(genomeSegmentsCN - genomeSegmentsStates) <= 0.25)
		nearStatesWidth = np.sum(genomeSegmentsWidth[nearStatesMask])
		nearStatesRatio = nearStatesWidth*100/totalSegmentsWidth
		##
		nearBoundariesMask  = (abs(genomeSegmentsCN - genomeSegmentsStates) > 0.25)
		nearBoundariesWidth = np.sum(genomeSegmentsWidth[nearBoundariesMask])
		nearBoundariesRatio     = nearBoundariesWidth*100/totalSegmentsWidth
		##		
		fileID2.write('Homogenous score : ' + str(nearStatesRatio) +' %' +'\n')
		##
		########################
		########################
		genomeRD = self.genomeRD
		genomeCN = genomeRD *2/self.copyNumberReference
		genomeStates = np.round(genomeCN)
		genomeWidth = len(genomeRD)
		##
		state0Mask  = (genomeStates == 0)
		state0Width = np.sum(state0Mask)
		state0Ratio = state0Width*100/genomeWidth
		##
		state1Mask  = (genomeStates == 1)
		state1Width = np.sum(state1Mask)
		state1Ratio = state1Width*100/genomeWidth

		##		
		fileID2.write('Dispersion score: ' + str(state0Ratio) +' %' +'\n')
		##
		###############################################################################################
		####################################### Median/ploidy ratio ###################################
		if(self.ploidyLevel == 'diploid'):
			peakValue = self.copyNumberReference
		elif(self.ploidyLevel == 'triploid'):
			peakValue = self.copyNumberReference*1.5
		elif(self.ploidyLevel == 'tetraploid'):
			peakValue = self.copyNumberReference*2
		else:
			raise ValueError('Error: ploidyLevel should be {diploid, triploid, or tetraploid}')	
		ploidyMedianError = 100*abs(peakValue - self.readDepthMedian)/peakValue
		####
		fileID2.write('Median error:' + str(ploidyMedianError) +' %' +'\n')
		fileID2.write('Median-Peak ratio:' + str(self.readDepthMedian/peakValue) +' %' +'\n')
		fileID2.write('\n')
		fileID2.write('\n')



		##############
		fileID2.close()
		#------------#


	def ploidySpectrum(self):
		"""
		return the ploidy spectrum info
		"""

		spectrum = []
		####################################### Ploidy & CE ###################################
		if(self.ploidyLevel == 'diploid'):
			ploidyNumber = 2
		elif(self.ploidyLevel == 'triploid'):
			ploidyNumber = 3
		elif(self.ploidyLevel == 'tetraploid'):
			ploidyNumber = 4
		else:
			raise ValueError('Error: ploidyLevel should be {diploid, triploid, or tetraploid}')	
		spectrum.append(ploidyNumber)
		#
		genomeRD = []
		for chrom in self.readCounts.keys():
			chrRD = self.readCounts[chrom]
			genomeRD = np.append(genomeRD, chrRD)
		spectrum.append(np.sum(genomeRD))	
		#
		for error in self.centralizationErrors:
			spectrum.append(error)

		#######################################################################################
		########################################## Errors #####################################	
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		chrSegments = self.ploidySegments
		for chrom in chrSegments.keys():
				chrSegmentsData  = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])

		# CNs of genome segments
		genomeSegmentsCN = genomeSegmentsRD *2/self.copyNumberReference
		# states of gneome segments
		genomeSegmentsStates = np.round(genomeSegmentsCN)
		#
		totalSegmentsWidth = np.sum(genomeSegmentsWidth)
		for i in range(0,10):
			CNiMask = (genomeSegmentsStates == i)
			genomeSegmentsWidthCNi  = genomeSegmentsWidth[CNiMask]
			iSegmentsWidth = np.sum(genomeSegmentsWidthCNi)
			iSegmentsRatio = iSegmentsWidth*100/totalSegmentsWidth
		##
		nearStatesMask  = (abs(genomeSegmentsCN - genomeSegmentsStates) <= 0.25)
		nearStatesWidth = np.sum(genomeSegmentsWidth[nearStatesMask])
		nearStatesRatio = nearStatesWidth*100/totalSegmentsWidth
		##
		nearBoundariesMask  = (abs(genomeSegmentsCN - genomeSegmentsStates) > 0.25)
		nearBoundariesWidth = np.sum(genomeSegmentsWidth[nearBoundariesMask])
		nearBoundariesRatio     = nearBoundariesWidth*100/totalSegmentsWidth
		##		
		spectrum.append(nearStatesRatio)


		########################
		genomeRD = self.genomeRD
		genomeCN = genomeRD *2/self.copyNumberReference
		genomeStates = np.round(genomeCN)
		genomeWidth = len(genomeRD)
		##
		state0Mask  = (genomeStates == 0)
		state0Width = np.sum(state0Mask)
		state0Ratio = state0Width*100/genomeWidth
		##
		state1Mask  = (genomeStates == 1)
		state1Width = np.sum(state1Mask)
		state1Ratio = state1Width*100/genomeWidth
		##	
		spectrum.append(state0Ratio)	


		##################################
		if(self.ploidyLevel == 'diploid'):
			peakValue = self.copyNumberReference
			ploidyNumber = 2
		elif(self.ploidyLevel == 'triploid'):
			peakValue = self.copyNumberReference*1.5
			ploidyNumber = 3
		elif(self.ploidyLevel == 'tetraploid'):
			peakValue = self.copyNumberReference*2
			ploidyNumber = 4
		else:
			raise ValueError('Error: ploidyLevel should be {diploid, triploid, or tetraploid}')	
		ploidyMedianError = 100*abs(peakValue - self.readDepthMedian)/peakValue
		spectrum.append(ploidyMedianError)	
		spectrum.append(self.readDepthMedian/peakValue)


		#######################################################################################
		######################################## Spectrum #####################################	
		genomeSegmentsRD    = np.array([]) 
		genomeSegmentsWidth = np.array([]) 
		
		# RDs & widths of genome segments
		chrSegments = self.ploidySegments
		for chrom in chrSegments.keys():
				chrSegmentsData  = chrSegments[chrom]
				genomeSegmentsWidth = np.append(genomeSegmentsWidth, chrSegmentsData[:,2])
				genomeSegmentsRD    = np.append(genomeSegmentsRD, chrSegmentsData[:,3])

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
			spectrum.append(iSegmentsRatio)

		for i in range(0,10):
			CNiMask = (genomeSegmentsStates == i)
			spectrum.append(np.sum(CNiMask))

		######	
		return spectrum




	def plotAneuploidy(self, chrName):
		"""
		plot resulted-segments of a chromosome
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
		
		#############################
		######### Plotting ##########
		CNmax = 10
		signal = np.clip(chrRDSignal,0,CNmax)
		true_chg_pts = bkps
		#
		if signal.ndim == 1:
			signal = signal.reshape(-1, 1)
		n_samples, n_features = signal.shape 

		COLOR_CYCLE = ["#4286f4", "#f44174"]
		figsize = (10, 2 * n_features)  # figure size
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
		plt.title(chrName)
		plt.xlabel('Bin number')
		plt.ylabel('Copy number')
		plt.show()
		#--------#



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
		CNmax = 10
		signal = np.clip(genomeCN,0,CNmax)
		#
		plt.figure()
		plt.plot(signal,'r', linestyle = 'None',marker = '.', markersize = 1)
		plt.plot([0,len(signal)],[1, 1],'k--')				
		plt.plot([0,len(signal)],[2, 2],'k--')
		plt.plot([0,len(signal)],[3, 3],'k--')				
		plt.plot([0,len(signal)],[4, 4],'k--')
		plt.plot([0,len(signal)],[5, 5],'k--')

		plt.title('Genome')
		plt.xlabel('Bin number')
		plt.ylabel('Copy number')

		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
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
		CNmax = 10
		signal = np.clip(genomeCN,0,CNmax)
		#
		plt.figure()
		plt.plot(signal,'r', linestyle = 'None',marker = '.', markersize = 1)
		plt.plot([0,len(signal)],[1, 1],'k--')				
		plt.plot([0,len(signal)],[2, 2],'k--')
		plt.plot([0,len(signal)],[3, 3],'k--')				
		plt.plot([0,len(signal)],[4, 4],'k--')
		plt.plot([0,len(signal)],[5, 5],'k--')

		plt.title('Genome')
		plt.xlabel('Bin number')
		plt.ylabel('Copy number')

		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_GenomeRD.png')
		#--------#



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
		if(self.ploidyLevel == 'diploid'):
			plt.plot([self.copyNumberReference, self.copyNumberReference],[0,max(n[0])],'k', lineWidth = 2)
		elif(self.ploidyLevel == 'triploid'):
			plt.plot([1.5*self.copyNumberReference, 1.5*self.copyNumberReference],[0,max(n[0])],'k', lineWidth = 2)
		elif(self.ploidyLevel == 'tetraploid'):
			plt.plot([2*self.copyNumberReference, 2*self.copyNumberReference],[0,max(n[0])],'k', lineWidth = 2)
		else:
			raise ValueError('Error: ploidyLevel should be {diploid, triploid, or tetraploid}')	


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
		if(self.ploidyLevel == 'diploid'):
			plt.plot([self.copyNumberReference, self.copyNumberReference],[0,max(n[0])],'k', lineWidth = 2)
		elif(self.ploidyLevel == 'triploid'):
			plt.plot([1.5*self.copyNumberReference, 1.5*self.copyNumberReference],[0,max(n[0])],'k', lineWidth = 2)
		elif(self.ploidyLevel == 'tetraploid'):
			plt.plot([2*self.copyNumberReference, 2*self.copyNumberReference],[0,max(n[0])],'k', lineWidth = 2)
		else:
			raise ValueError('Error: ploidyLevel should be {diploid, triploid, or tetraploid}')	

		inputFileName      = os.path.basename(self.inputFile)
		inputFileNameNoExt = os.path.splitext(inputFileName)[0]
		plt.savefig(self.outputFolder + '/' + inputFileNameNoExt + '_' + str(nbins) + 'bin_GenomeHistogram.png')
		plt.close()		
		#--------#





	########################################################################################################################
	########################################################################################################################
	########################################################################################################################
	def ploidyEstimatorPipeline(self):
		"""
		the complete pipeline for estimating the ploidy level of the sequencing data 
		and write the output files.	
		"""
		# RD-calculation
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

		# Anueploidy Spectrum 
		ploidyModels = ['model1', 'model2', 'model3', 'model4', 'model5','model6', 'model7', 'model8', 'model9', 'model10']
		ploidySegments = dict()
		copyNumberReference = []
		ploidyModelsError   = []
		for ploidy in ploidyModels:
			print(ploidy)
			# copy-number estimation
			estimatedCNReference = self.estimateCopyNumberReference(genomeRD, ploidy)	
			copyNumberReference.append(estimatedCNReference)
			# segments merging based on the estimated copy-number reference
			ploidySegments[ploidy] = self.segmentsMerging(estimatedCNReference)
			# centralization-error
			centralizationError  = self.computeCentralizationError(estimatedCNReference, ploidySegments[ploidy])
			ploidyModelsError.append(centralizationError)
		#
		ploidyIndex = ploidyModelsError.index(min(ploidyModelsError))
		self.copyNumberReference = copyNumberReference[ploidyIndex]
		self.ploidyModel  = ploidyModels[ploidyIndex]
		self.ploidySegments = ploidySegments[self.ploidyModel]
		self.centralizationErrors = ploidyModelsError

		# Ploidy level
		self.ploidyLevel    = self.findNearestPloidy()
		print(self.ploidyLevel)
		print(self.ploidyModel)
		# Output
		self.writeAneuploidyResults()
		self.writeAneuploidySpectrum(ploidyModelsError)
		


