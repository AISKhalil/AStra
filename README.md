# AStra: Aneuploidy Spectrum (Digital Karyotyping) for Rapid Authentication of Cell Lines

**Motivation:** The widespread concern about genetic drift and cross-contamination of cell lines calls for a pressing need for their authentication as a basic quality control for scientific data. The current genetic techniques for authentication are time-consuming and require specific documentary standard and laboratory protocols. Given the fact that whole-genome sequencing (WGS) data are readily available, read depth (RD)-based computational analyses has allowed the estimation of genetic profiles of cell lines.

**AStra** is a Python-based software AStra for de novo estimation of the genome-wide aneuploidy profile from raw WGS reads. 
It is a simple, user-friendly, and free tool that provides the elementary information about the chromosomal aneuploidy for 
cell line authentication. **AStra** provides an analytical and visualization platform for rapid and easy comparison between different 
cell lines/strains.  We recommend **AStra** for rapid first-pass quality assessment of scientific data that employ cancer cell lines.

## Contents
- [Requirements](#requirement)
- [Download](#Download)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)

### <a name="requirement"></a>Requirements

- python >=3.6.5
- pysam  >=0.12


### <a name="Download"></a>Download

```bash
cd ~
git clone https://github.com/AISKhalil/AStra.git
```

### <a name="installation"></a>Installation

User can create "Python virtual environment" and install the required libraries using the following commands.

- Make sure pip3 is installed. If not, then please check your distribution's documentation and install it. In Ubuntu 18.04 LTS you can do:

```bash
sudo apt update
sudo apt install python3-pip
pip3 install --upgrade pip
```

- Go to the AStra folder, create and activate a new virtual environment

```bash
cd AStra
sudo apt-get install python3-venv
python3 -m venv AStraPythonEnv
source AStraPythonEnv/bin/activate
```

- Install dependencies: please make sure you have install all the libraries/headers. In Ubuntu 18.04 LTS you can do:

```bash
sudo apt install python3-tk
sudo apt install zlib1g zlib1g-dev zlib1g-dbg liblzma-dev libbz2-dev
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


### <a name="usage"></a>Usage

`AStra` is developed as a Python-class. Therefore, we added two scripts as simpler interfaces for AStra: **AStraSingleInput** and **AStraMultipleInputs**. In order to use `AStra`, you need first to activate the python virtual environment that you created before:

```bash
cd AStra
source AStraPythonEnv/bin/activate
```

#### 1) Running AStra on single input:

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

#### 2) Running AStra on multiple bam files:

```bash
python ./AStraMultipleInputs.py
```

You can modify this script to add the path of many BAM files. Besides to the default output files, this script outputs an excel file contains the complete profile of all BAM files (profile of each BAM file in a row). Each row contains fields of 1: whole-genome ploidy number, 2: copy number reference (CNR), 3: centralization error (CE), 4: centralization score, 5: read depeth(RD)-median, 6: RD-median/CNR, 7: total counts of NGS reads, 8-18: anueploidy spectrum (percentages of genome-segments per CN-state), 19: whole-genome ploidy state, 20: AStra model, and 21-26: CEs for all models (m1-m6). This eases the analysis of large number of BAM files.


### <a name="output"></a>Output

**AStra** generates many output files providing the detailed characterization of the aneuploidy profile 
of the input cell line. Results of the analyzed cell lines were uploaded in `AStraResults/` folder.

>    **a. Aneuploidy characterization:** a text file containing the important features of aneuploidy profile of the input
>    BAM such as whole-genome ploidy number, copy number reference (CNR), centralization error (CE), centralization score, 
>    read depth(RD)-median, RD-median/CNR, total counts of NGS reads, and anueploidy spectrum (percentages of genome->
>    segments per CN-state).
   
>    **b. Chromosomal aneuploidy:** a narrowPeak BED format file of the approximated segmental anueploidy of the 
>    complete genome for UCSC Genome Browser.

>    **c. Aneuploidy profile:** a figure of the RD signal of the genome after setting the CN reference. 
>    For example, aneuploidy profile of A427 cell line is:
![VCaP coverage plot](/AStraResults/A427_GenomeAneuploidyWithChrsMarkers.png)
 
>    **d. Anueploidy spectrum:** a figure of the frequency distribution of the input sequencing data. In the frequency distribution,
>    The red line denotes the RD median whereas the dotted black lines denote the CN states. 
>    For example, anueploidy spectrum of A427 cell line is:
![VCap anueploidy spectrum](/AStraResults/A427_200bin_GenomeHistogram.png)
 
 
