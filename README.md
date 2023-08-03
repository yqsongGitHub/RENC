
<!-- README.md is generated from README.Rmd. Please edit that file -->

</div>

<h2 align="center"> RENC: Ranking of target genes based on ENhancer-promoter Contacts </h2>

<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/Schematic.png">
</p>   

------------------------------------------------------------------------

<!-- badges: start -->
<!-- badges: end -->
## Introduction 
Our objective was to identify the target genes associated with putative enhancers located within duplication hotspots. These hotspots often contain multiple putative enhancer elements, forming complex enhancer-promoter pairs that pose challenges in pinpointing the actual targets of duplication events. To address this, we developed an innovative methodology called Ranking of target genes based on ENhancer-promoter Contacts (RENC). RENC enables the definition of enhancer-gene relationships for duplication hotspots, considering both the enhancer and the gene's perspectives.

To determine the enhancer-gene relationships, we employed HiChIP paired-end tags (PETs) specifically for enhancers within each duplication hotspot. This allowed us to measure the amount of enhancer activity delivered to each candidate gene, referred to as hotspot-delivered enhancer activity. Subsequently, we evaluated the relative contribution of hotspot-delivered enhancer activity for each candidate gene, as a percentage of the overall activity delivered by all enhancers interacting with that gene.To prioritize target genes for enhancers within the duplication hotspots, we multiplied these two factors—hotspot-delivered enhancer activity and the relative contribution—to assign a ranking. Additionally, we applied weights to the enhancers associated with each target gene based on their respective contributions.

Through the utilization of RENC, our methodology offers a comprehensive framework for identifying and prioritizing target genes in duplication hotspots, considering enhancer-gene relationships from multiple perspectives.<br />

RENC is designed with respect reference to [BEDTools](http://quinlanlab.org/tutorials/bedtools/bedtools.html) for command-line style programming, we provide a shell script and user guide via this github page.

------
## Requirements
For each cell-type, the inputs to the RENC methodology are:

* Required Inputs
 	* 1. BED files for the genomic coordinate of duplication or region 
 	* 2. BEDPE files for HiChIP
 	* 3. BED files for transcription start site（TSS）of protein coding genes
 	* 4. BED files for gene body of protein coding genes
* Optional Inputs
 	* 5. genes or the genomic coordinate in duplications to search; e.g.,"C4BPB","chr1:201970000-202085000"or "C4BPB;chr1:201970000-202085000"  

------

## Basic Usage and Quick Guide
### Example and reference data background introduction

Example data for testing is available at [RENC/example](https://github.com/yqsongGitHub/RENC/tree/main/example/input). 
 * 1. BED file for the genomic coordinate of duplication or region were from the duplication hotspots of Stomach adenocarcinoma (STAD) in the Pan-Cancer Atlas of Whole Genomes (PCAWG) project.
 * 2. BEDPE file for HiChIP were from HiChIP experiments mapped to hg19 for chromosome 1 in STAD (Stomach adenocarcinoma) cell line. Only intra-chromosomal PETs were kept.

Reference data is available at [RENC/reference](https://github.com/yqsongGitHub/RENC/tree/main/reference/hg19).
 * 3. BED file for transcription start site（TSS）of protein coding genes were download from  [UCSC](https://genome.ucsc.edu/)
 * 4. BED file for gene body of protein coding genes were download from  [UCSC](https://genome.ucsc.edu/)

### Installation
#### 1. Create Conda environment
By creating a Conda environment, you can have separate environments for different projects or purposes, ensuring that the packages and dependencies within each environment do not conflict with each other.
```
conda create -n RENC
```
#### 2. Activate your Conda environment
Once the environment is created, you can activate it and start installing packages or running programs within that environment, ensuring that the installed packages and dependencies are isolated from your system's global environment.
```
conda activate RENC
```
#### 3. Install bedtools using the Conda package manager
To install the bedtools as a dependency, once the installation is finished, you can start using bedtools within your Conda environment.
```
conda install -c bioconda bedtools
# If the installation is successful, the script below should run without any issues
bedtools -h
```
#### 4. Download RENC script via GitHub
It is recommended to download using ***git clone***. However, in case of any network issues, you can download directly by clicking the "Download ZIP" button on the GitHub website. By downloading the script, you acquire a local copy of the file that you can use and modify as needed on your own machine.
```
git clone https://github.com/yqsongGitHub/RENC.git
cd RENC
```
Then, you can use the following example command lines to test the functionality of RENC after installation.

### Routine analysis 1: Get all the enhancer-gene regulatory relationships from the input BEDPE file and BED file.
```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_STAD_chr1.bed \
-c ./example/input/AGS_STAD_chr1.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed 
```
Please note, in RENC, all files using the reference genome must be consistent, with the example using hg19. If there are files using the hg38 reference genome, the [LiftOver](https://liftover.broadinstitute.org/) tool can be used to convert the hg38 genome to hg19.

The informative output is a .RENC.txt file with annotation of information as follows.  

column | name | explanation
------ | ---- | ------------
1th | chr1 | chromosomal for the loop first anchor
2th | start1 | genomic coordinate of the start site for the first anchor
3th | end1 | genomic coordinate of the end site for the first anchor
4th | chr2 | chromosomal for the loop second anchor
5th | start2 | genomic coordinate of the start site for the second anchor
6th | end2 | genomic coordinate of the end site for the second anchor
7th | pets | observed paired-end tags number linking the two anchors
8th | chr | chromosomal for the duplication
9th | start | genomic coordinate of the start site for the duplication
10th | end | genomic coordinate of the end site for the duplication
11th | dup_sum_pets | the amount of duplicated enhancer activity delivered to the candidate gene
12th | gene_ID_sum_pets | the amount of all enhancers activity delivered to the candidate gene
13th | percentage | the relative contribution of the duplicated enhancer activity as percentage of activity delivered from all enhancers interacting with the candidate gene
14th | score | RENC score for the enhancer interacting with the gene is calculated by multiplying the amount of duplicated enhancer activity (dup_sum_pets) by the relative contribution of the duplicated enhancer activity (percentage)
15th | rank_gene | the rank of the target genes for the enhancers in the duplication hotspot
16th | gene | the specific symbol or identifier associated with the candidate gene
17th | rank_enhancer | the rank of the enhancers in the duplication hotspot delivered to the candidate gene
18th | enhancer | genomic coordinate of the enhancer


### Routine analysis 2: Get the target gene from input genomic coordinate of the duplicated enhancer.

If using multiple regions as input, please separate them with a semicolon (;). If no results are found, return a placeholder (.).

```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_STAD_chr1.bed \
-c ./example/input/AGS_STAD_chr1.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed \
-s  "chr1:201970000-202085000;chr1:33135000-33370000"
```

The informative output is a .RENC.search.txt file with annotation of information as follows.  
column | name | explanation
------ | ---- | ------------
1th | enhancer | the input genomic coordinate of the duplicated enhancer
2th | gene | the target gene. If no results are found, return a placeholder (.)

### Routine analysis 3: Get genomic coordinate of the duplicated enhancer from input target gene.

If using multiple regions as input, please separate them with a semicolon (;). If no results are found, return a placeholder (.).

```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_STAD_chr1.bed \
-c ./example/input/AGS_STAD_chr1.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed \
-s  "C4BPB;ELF3"
```
The informative output is a .RENC.search.txt file with annotation of information as follows.  
column | name | explanation
------ | ---- | ------------
1th | enhancer | genomic coordinate of the duplicated enhancer. If no results are found, return a placeholder (.)
2th | gene | the input target gene

------
## RENC Main Functions

Run ***bash RENC.sh*** or ***bash RENC.sh -h*** can show the main functions of RENC.sh with short descriptions and examples.     

```
Usage:   
  bash RENC.sh [-r <RSCRIPT>][-d <DUPLICATION_FILE>] [-c <HICHIP_FILE>][-t <TSS_FILE>][-g <GENEBODY_FILE>]
  -r Rscript	Rscript file [required].
  -d DUPLICATION_FILE	Duplication or region file [required].
  -c HICHIP_FILE	HiChIP file from Hichipper pipeline or other bedpe file [required].
  -t TSS_FILE	TSS file for protein coding genes or other TSS file [required].
  -g GENEBODY_FILE	Gene body file for protein coding genes [required].
  -s GENE_ENH	Gene or region in Duplication to search; ALL gene and region by default. e.g.,"C4BPB","chr1:201970000-202085000" or "C4BPB:chr1:201970000-202085000"

Example:     
  bash ./src/RENC.sh -r ./src/RENC.R  
  -d ./example/input/Duplication_STAD_chr1.bed  
  -c ./example/input/AGS_STAD_chr1.bedpe  
  -t ./reference/hg19/RefSeq_proteinCoding.tss.bed  
  -g ./reference/hg19/RefSeq_proteinCoding.body.bed 

```

--------
## Questions & Answers
Please submit a github issue with any questions or if you experience any issues/bugs.
- Create a new [issue](https://github.com/yqsongGitHub/RENC/issues).

-------
## Warning   
No tested bugs needed warning. 

--------
## Citation
Coming soon