# AStra: Aneuploidy Spectrum Analysis as a Primer for Copy Number Profiling of Cancer Cells 

**Motivation:** Hyperploidy and segmental aneuploidy are hallmarks of cancer cells due to chromosome segregation errors and genomic instability. In such a situation, accurate aneuploidy profiling of cancer data is critical for calibration of copy-number detection tools. Additionally, cancer data may represent cell population which suffers from different levels of clonal heterogeneity. The degree of heterogeneity adversely affects the segregation of the genome into integral copy number states. This, in turn, strongly influences the reliability of this data for aneuploidy profiling and copy number analysis.

**Results:** We have developed `AStra` framework for aneuploidy profiling of cancer NGS (next-generation sequencing) data and assessing their suitability for copy number analysis without any prior knowledge of the input cell line. AStra estimates the best-fit aneuploidy profile as the spectrum with most genomic segments around integral copy number (CN) states. After defining the best-fit spectrum, we compute the homogeneity score (HS) that measures the localization degree of genomic regions around CN states.

**Availability and implementation:** `AStra` is an open source software implemented in Python and is available at https://github.com/AISKhalil/AStra.

**For more information please visit [AStra Manuscript](https://www.biorxiv.org/content/10.1101/639294v1?rss=1).**


## Requirements

- python >=3.6.5
- pysam  >=0.12

## Installation

### Download AStra from github

```bash
cd ~
git clone https://github.com/AISKhalil/AStra.git
```

### Python virtual environment and python dependencies

Make sure pip3 is installed. If not, then please check your distribution's documentation and install it. In Ubuntu 18.04 LTS you can do:

```bash
sudo apt update
sudo apt install python3-pip
pip3 install --upgrade pip
```

Go to the AStra folder, create and activate a new virtual environment

```bash
cd AStra
python3 -m venv AStraPythonEnv
source AStraPythonEnv/bin/activate
```

Install dependencies: please make sure you have install all the libraries/headers. In Ubuntu 18.04 LTS you can do:

```bash
sudo apt install python3-tk
sudo apt install zlib1g zlib1g-dev zlib1g-dbg
pip3 install wheel
pip3 install numpy
pip3 install scipy
pip3 install pandas
pip3 install cython
pip3 install h5py
pip3 install matplotlib
pip3 install ruptures
pip3 install scikit-learn
pip3 install hmmlearn
pip3 install pysam
pip3 install pysamstats
pip3 install numexpr
pip3 install xlsxwriter
```


## AStra execution 

`AStra` is developed as a Python-class. Therefore, we added two scripts as simpler interfaces for AStra: **AStraSingleInput** and **AStraMultipleInputs**. In order to use `AStra`, you need first to activate the python virtual environment that you created before:

```bash
cd AStra
source AStraPythonEnv/bin/activate
```

### 1) Running AStra on single input:

```bash
python ./AStraSingleInput.py -b input.bam -f hg19.ucsc.fa -o AStraResults
```

There are three required parameters: the input BAM file (`-b` or `--bam`), the reference genome (`-f` or `--fa`), and the output folder (`-o` or `--out`) that will be used to generate the output files.

```bash
######
usage: AStraSingleInput.py -b input.bam -f hg19.ucsc.fa -o AStraResults

arguments:
  -b, --bam             the input sorted BAM file. If not sorted, you can use samtools to sort it ("samtools sort input.bam > input.sorted.bam").
  -f, --fa              the fasta file of the reference genome. For human hg19, you can download from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz.
  -o, --out             the folder path to keep the results of AStra.
```

### 2) Running AStra on multiple bam files:

```bash
python ./AStraMultipleInputs.py
```

You can modify this script to add the path of many BAM files. Besides to the default output files, this script outputs an excel file contains the complete profile of all BAM files (profile of each BAM file in a row). Each row contains fields of 1: nearest ploidy, 2: number of reads, 3-12: Centralization error (CE) for model1-model10, 13: HS, 14: CN State 0 percentage, 15: Median error, 16: Median correction factor, 17-26: Anueploidy spectrum (percentages of genome-segments per CN-state), and 27-36: number of segments per each CN-state. This eases the analysis of large number of BAM files.


## AStra outputs 

`AStra` generates many output files providing the detailed characterization of the aneuploidy profile 
of the input cell line. Results of the analyzed cell lines were uploaded in **AStraResults** folder.

>    **a. Aneuploidy characterization:** a text file containing the important features of aneuploidy profile of the input BAM
>    such as nearest ploidy, copy number reference, CE for each model, anueploidy spectrum, homogeneity score, CN State 0 percentage, Median 
>    error, and Median correction factor.
   
>    **b. Aneuploidy profile:** a narrowPeak BED format file of the approximated segmental anueploidy of the complete genome
>    for UCSC Genome Browser.

>    **c. Coverage plot of the genome:** a figure of the RD signal of the genome after setting the CN reference. 
>    For example, coverage plot of VCaP cell line is:
![VCaP coverage plot](/AStraResults/VCap_ENCFF273KUS_ENCFF466WDC_merged_Input_hg19_CK_bowtie2_default_rmdup.removed.blackList.Region.bedtools_GenomeRD.png)
 
>    **d. Anueploidy spectrum:** a figure of the frequency distribution of the input sequencing data. In the frequency distribution,
>    the black line denotes the ploidy CN state (major peak) where the red line denotes the global median. The dotted black lines 
>    denote the other CN states. For example, anueploidy spectrum of VCaP cell line is:
![VCap anueploidy spectrum](/AStraResults/VCap_ENCFF273KUS_ENCFF466WDC_merged_Input_hg19_CK_bowtie2_default_rmdup.removed.blackList.Region.bedtools_200bin_GenomeHistogram.png)
 
 
