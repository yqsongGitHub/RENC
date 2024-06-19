library(gTrack)
library(gUtils)
library(rtracklayer)
source("./src/gTrack.R")
source("./src/function.R")
options(scipen = 100)
options(warn = -1)
#gt.ge <- track.gencode()
#saveRDS(gt.ge,"./input/gt.ge.rds")
gt.ge <- readRDS("./data/gt.ge.rds")
genes <- read.table("./reference/hg19/RefSeq_proteinCoding.body.bed")
tss <- read.table("./reference/hg19/RefSeq_proteinCoding.tss.bed")
chip <- './data/GSM2356644_H3K27ac-BICR31.bw'
ins <- read.table("./data/BICR31_KLF5.bedpe")
dup <- read.table("./data/dup.Squamous.merged.1mb.sort.genomecov")
hs <- read.table("./example/input/Duplication_Squamous.bed")
pps <- read.table("./example/output/SqCC_BICR31_RENC.txt",header = T)
bp <- read.table("./example/output/SqCC_BICR31_tss.bed")

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

# add chip 
chip.gr.score <- gTrack(chip,col = '#858480', name = 'chip',bars = TRUE,height = Height+0.5, ygap = Ygap)

# add hotspot
hs.gr <- GRanges(seqnames = Rle(hs$V1), 
                 ranges = IRanges(start = hs$V2,
                                  end = hs$V3))
hs.gr <- gTrack(hs.gr,col = '#198330', name = '',bars = TRUE,height = Height/2, ygap =Ygap/2)

# add dup 
dup.gr <- GRanges(seqnames = Rle(dup$V1), 
                  ranges = IRanges(start = dup$V2,end = dup$V3),
                  score =dup$V4)
dup.gr <- gTrack(dup.gr, y.field = 'score',
                 col = '#858480', name = 'dup',bars = TRUE,height = Height+0.5, ygap =Ygap, y0=1)

# add score 
pps1 <- pps1[,c("ID","score")]
pps1 <- merge(tss,pps1,by.x="V5",by.y="ID",all.x=T)
pps1$score[is.na(pps1$score)] <- 0
pps1 <- pps1[!(duplicated(pps1$V4) & pps1$score==0),]
pps.gr1 <- GRanges(seqnames = Rle(pps1$V1), 
                   ranges = IRanges(start = pps1$V2,end = pps1$V3),
                   score =pps1$score)
pps.gr1 <- gTrack(pps.gr1, y.field = 'score',
                  col = '#858480', name = 'RENC',circles = TRUE,height = Height+0.5, ygap = Ygap)

# add enhancer
bp <- bp[bp$gene==Gene,]
enr.gr <- GRanges(seqnames = Rle(bp$chr), 
                  ranges = IRanges(start = bp$start,end = bp$end),
                  pets =bp$pets)
enr.gr <- gTrack(enr.gr, y.field = 'pets',
                 col = '#858480', name = 'contrib.',circles = TRUE, height = Height+1,ygap = Ygap)

pdf("./Images/BICR31_KLF5.pdf",height = 12,width = 8)
plot(c(enr.gr,pps.gr1,chip.gr.score,dup.gr,hs.gr,gt.ge,dat.gt),wins)
dev.off()








