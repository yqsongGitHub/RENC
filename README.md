
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

The formula is as follows：
<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/Formula.png" alt="Formula" width="50%" height="50%">
</p>   

$E_i$ is the number of PETs connecting the enhancer i to a candidate gene’s promoter. $E_{in}$ is all enhancers within the hotspot connecting to the candidate gene promoter. $E_{out}$ is all enhancers outside the hotspot connecting to the candidate gene promoter. $\sum_{i\in E_{in}} E_i$ is the total number of PETs connecting the enhancers within the hotspot to a candidate gene’s promoter. $\sum_{i\in E_{out}} E_i$ is the total number of PETs connecting the enhancers outside the hotspot to a candidate gene’s promoter.

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
17th | ID | the transcript ID or identifier associated with the candidate gene
18th | rank_enhancer | the rank of the enhancers in the duplication hotspot delivered to the candidate gene
19th | enhancer | genomic coordinate of the enhancer


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

------
## KLF5 example
### Step1:  Get all the enhancer-gene regulatory relationships from the input BEDPE file and BED file. 
The ***Duplication_Squamous.bed*** file contains the genomic coordinates of duplication hotspots or regions associated with squamous cell carcinoma. The ***SqCC_BICR31.bedpe*** file, derived from HiChIP experiments, is mapped to hg19 in the BICR31 (squamous cell carcinoma) cell line. The ***output*** provides all the enhancer-gene regulatory relationships and the RENC score for each candidate gene, including KLF5.
```
bash ./src/RENC.sh -r ./src/RENC.R  \
-d ./example/input/Duplication_Squamous.bed \
-c ./example/input/SqCC_BICR31.bedpe  \
-t ./reference/hg19/RefSeq_proteinCoding.tss.bed \
-g ./reference/hg19/RefSeq_proteinCoding.body.bed
```

### Step2: Visualize the region near KLF5.
#### Step2.1：Set global variables and load R packages.
We used R version 4.2.2 and the [gTrack](https://github.com/mskilab-org/gTrack) package for visualization.
We also need to load some functions.
```
library(gTrack)
library(gUtils)
library(rtracklayer)
source("./src/gTrack.R")
source("./src/function.R")
options(scipen = 100)
options(warn = -1)
```

#### Step2.2：Set global variables and load R packages.
Here we load some necessary files, including the ***SqCC_BICR31_RENC.txt*** and ***SqCC_BICR31_tss.bed*** files generated in step 1.
```
gt.ge <- readRDS("./data/gt.ge.rds")
genes <- read.table("./reference/hg19/RefSeq_proteinCoding.body.bed")
tss <- read.table("./reference/hg19/RefSeq_proteinCoding.tss.bed")
chip <- './data/GSM2356644_H3K27ac-BICR31.bw'
ins <- read.table("./data/BICR31_KLF5.bedpe")
dup <- read.table("./data/dup.Squamous.merged.1mb.sort.genomecov")
hs <- read.table("./example/input/Duplication_Squamous.bed")
pps <- read.table("./example/output/SqCC_BICR31_RENC.txt",header = T)
bp <- read.table("./example/output/SqCC_BICR31_tss.bed")
```

#### Step2.3：Visualize the region near KLF5
We selected the region from 73,205,000 to 74,345,000 on chromosome 13, which includes H3K27ac HiChIP signal, positions of duplication hotspots, duplication events observed in squamous cancer, H3K27ac ChIP-seq signal, RENC scores prioritizing genes that are more likely to be activated by enhancers within the first duplication hotspot presented in the window, and relative contributions of the enhancers to KLF5.

Here is the plot of KLF5 that we created. Further editing and enhancement using [Adobe Illustrator](https://www.adobe.com/products/illustrator/free-trial-download.html) is required.
<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/BICR31_KLF5.png" width="80%" height="80%">
</p>  

```
Gene <- "KLF5"
Chr <- "chr13"
Chr.start <-  73205000
Chr.end <- 74345000
Height <- 3
Ygap <- 0.8

hs <- hs[hs$V1=="chr13" & hs$V2>73210000 & hs$V3 < 74340000,]
hs1 <- hs[1,] 
pps1 <- pps[pps$chr == hs1$V1 & pps$start == hs1$V2 & pps$end == hs1$V3,]
bp <- find_enh_anchor(bp)
hara <- draw_gtrack(ins,chr=Chr,chr.start = Chr.start,chr.end=Chr.end)
dat.gt <- hara[[1]]
wins <- hara[[2]]
ind <- which(names(gt.ge@data[[1]]) %in% genes$V4)
gt.ge@data[[1]] <- gt.ge@data[[1]][ind] 
gt.ge@formatting$height <- Height+2
gt.ge@formatting$ygap <- Ygap

# H3K27ac ChIP-seq signal 
chip.gr.score <- gTrack(chip,col = '#858480', name = 'chip',bars = TRUE,height = Height+0.5, ygap = Ygap)

# Duplication hotspot
hs.gr <- GRanges(seqnames = Rle(hs$V1), 
                 ranges = IRanges(start = hs$V2,
                                  end = hs$V3))
hs.gr <- gTrack(hs.gr,col = '#198330', name = '',bars = TRUE,height = Height/2, ygap =Ygap/2)

# Duplication events 
dup.gr <- GRanges(seqnames = Rle(dup$V1), 
                  ranges = IRanges(start = dup$V2,end = dup$V3),
                  score =dup$V4)
dup.gr <- gTrack(dup.gr, y.field = 'score',
                 col = '#858480', name = 'dup',bars = TRUE,height = Height+0.5, ygap =Ygap, y0=1)

# RENC scores 
pps1 <- pps1[,c("ID","score")]
pps1 <- merge(tss,pps1,by.x="V5",by.y="ID",all.x=T)
pps1$score[is.na(pps1$score)] <- 0
pps1 <- pps1[!(duplicated(pps1$V4) & pps1$score==0),]
pps.gr1 <- GRanges(seqnames = Rle(pps1$V1), 
                   ranges = IRanges(start = pps1$V2,end = pps1$V3),
                   score =pps1$score)
pps.gr1 <- gTrack(pps.gr1, y.field = 'score',
                  col = '#858480', name = 'RENC',circles = TRUE,height = Height+0.5, ygap = Ygap)

# Enhancer contributions
bp <- bp[bp$gene==Gene,]
enr.gr <- GRanges(seqnames = Rle(bp$chr), 
                  ranges = IRanges(start = bp$start,end = bp$end),
                  pets =bp$pets)
enr.gr <- gTrack(enr.gr, y.field = 'pets',
                 col = '#858480', name = 'contrib.',circles = TRUE, height = Height+1,ygap = Ygap)

pdf("./Images/BICR31_KLF5.pdf",height = 12,width = 8)
plot(c(enr.gr,pps.gr1,chip.gr.score,dup.gr,hs.gr,gt.ge,dat.gt),wins)
dev.off()
```

### Step3 : Enhance the image using Adobe Illustrator.
In step 2, we obtained the plot of KLF5. We used [Adobe Illustrator](https://www.adobe.com/products/illustrator/free-trial-download.html) to enhance the representation of protein-coding genes and annotated the previously publishedfunctional enhancers for KLF5.

<p align="center">
<img align="center" src="https://github.com/yqsongGitHub/RENC/blob/main/Images/BICR31_KLF5_ai.png" width="80%" height="80%">
</p>  

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