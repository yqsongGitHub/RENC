#!/usr/bin/env Rscript

args <-  commandArgs(trailingOnly=TRUE)
tssfile <- args[1]
enhg <- args[2]

# function 
filter_func <- function(hichip2,hichip3){
  hichip2$fea <- paste(hichip2[,9],hichip2[,10],hichip2[,11],hichip2[,23],
                       hichip2[,24],hichip2[,7],hichip2[,1],hichip2[,2],
                       hichip2[,3],hichip2[,4],hichip2[,5],hichip2[,6],
                       sep="_")
  hichip3$fea <- paste(hichip3[,9],hichip3[,10],hichip3[,11],hichip3[,23],
                       hichip3[,24],hichip3[,7],hichip3[,1],hichip3[,2],
                       hichip3[,3],hichip3[,4],hichip3[,5],hichip3[,6],
                       sep="_")
  hichip3 <- hichip3[!hichip3$fea %in% hichip2$fea,]
  hichip3 <- hichip3[,!colnames(hichip3) %in% "fea"]
  return(hichip3)
}

cal_hichip_by_gene <- function(hichip.tss,hichip1,hichip2,hichip3,tss,fc.dat=NULL,selectmax=TRUE){
  hichip11 <- unique(hichip1[,c(9:11,16,17,7,1:6)])
  hichip21 <- unique(hichip2[,c(9:11,23,24,7,1:6)])
  hichip31 <- unique(hichip3[,c(9:11,23,24,7,1:6)])
  colnames(hichip11) <- c("chr","start","end","gene","ID","pets","chr1","start1","end1","chr2","start2","end2")
  colnames(hichip21) <- c("chr","start","end","gene","ID","pets","chr1","start1","end1","chr2","start2","end2")
  colnames(hichip31) <- c("chr","start","end","gene","ID","pets","chr1","start1","end1","chr2","start2","end2")
  hichip11$type  <- "non_gene"
  hichip21$type  <- "gene"
  if(nrow(hichip31)==0){
    hichip31$type  <- NULL
  }else{
    hichip31$type  <- "gene_tss_out_dup"
  }
  hichip1 <- rbind(hichip11,hichip21)
  hichip1 <- rbind(hichip1,hichip31)
  hichip1 <- unique(hichip1)
  hichip1 <- hichip1[!hichip1$pets==".",]
  hichip1$fea <- paste(hichip1[,1],hichip1[,2],hichip1[,3],sep="_")
  hichip1$gene_ID <- paste0(hichip1[,4],"_",hichip1[,5])
  hichip1$fea_gene_ID <- paste0(hichip1$fea,"_",hichip1$gene_ID)
  hichip1$pets <- as.numeric(hichip1$pets)
  hichip2 <- lapply(split(hichip1,hichip1$fea_gene_ID), function(x){
    rt <- unique(x)
    rt$sum <- sum(x[,6])
    return(rt)
  })
  hichip2 <-  do.call(rbind,hichip2)
  hichip.tss <- unique(hichip.tss[,c(1:8,12,13)])
  hichip.tss$gene_ID <- paste0(hichip.tss[,9],"_",hichip.tss[,10])
  pets.sum <- lapply(split(hichip.tss,hichip.tss$gene_ID), function(x){
    sum(x[,7])
  })
  pets.sum <- do.call(rbind,pets.sum)
  # merge
  hichip2 <- merge(hichip2,pets.sum,by.x= "gene_ID",by.y="row.names")
  hichip2$percentage <- hichip2[,17]/hichip2[,18]
  if(selectmax==T){
    hichip2 <- lapply(split(hichip2,hichip2$gene), function(x){
      x[x[,18]==max(x[,18]),]
    })
    hichip2 <- do.call(rbind,hichip2)
  }
  colnames(hichip2)[17:18] <- c("dup_sum_pets","gene_ID_sum_pets")
  hichip2$contribution <- hichip2$pets/hichip2$gene_ID_sum_pets
  # rank gene by dup
  hichip3 <- hichip2
  hichip3$score <- hichip3$dup_sum_pets * hichip3$percentage
  hichip3 <- lapply(split(hichip3,hichip3$fea), function(x){
    x$rank_gene <- as.numeric(factor(rank(-x$score)))
    return(x)
  })
  hichip3 <- do.call(rbind,hichip3)
  # rank enhancer by gene 
  hichip4 <- lapply(split(hichip3,list(hichip3$fea,hichip3$gene)), function(x){
    x$rank_enhancer <- as.numeric(factor(rank(-x$pets)))
    return(x)
  })
  hichip4 <- do.call(rbind,hichip4)
  hichip5 <- merge(hichip4,tss,by.x="ID",by.y="V5")
  inanchor1 <- hichip5$start1 < hichip5$V2 & hichip5$end1 > hichip5$V3
  for (aa in inanchor1) {
    if(aa){
      hichip5$enhancer <-  paste0(hichip5$chr1,":",hichip5$start1,"-",hichip5$end1)
    }else{
      hichip5$enhancer <-  paste0(hichip5$chr2,":",hichip5$start2,"-",hichip5$end2)
    }
  }
  if(!is.null(fc.dat)){
    hichip5 <- merge(hichip5,fc.dat,by.x="gene",by.y ="gene")
  }
  hichip5 <- hichip5[,c(8:13,7,3:5,17:19,21:22,6,23,29)]
  return(hichip5)
}

sel_gene_enh <- function(RENC,tss,fea){
  fea <- as.character(fea)
  re0 <- data.frame()
  for(i in fea){
    if(grepl("^chr",i)){
      fea1 <- gsub(":|-","_",i)
      RENC$fea <- paste(RENC$chr,RENC$start,RENC$end,sep="_")
      re1 <- RENC[RENC$fea %in% fea1 & RENC$rank_gene==1,]
      
      if(nrow(re1)==0){
        re2 <- data.frame("enhancer"=i,"gene"=".")
      }else{
        re2 <- data.frame("enhancer"=i,"gene"= paste(unique(re1$gene),collapse = ";"))
      }
      re0 <- rbind(re0,re2)
    }else{
      re1 <- RENC[RENC$gene %in% i & RENC$rank_enhancer==1,]
      if(nrow(re1)==0){
        re2 <- data.frame("enhancer"=".","gene"=i)
      }else{
        re2 <- data.frame("enhancer"= paste(unique(re1$enhancer),collapse = ";"),"gene"=i)
      }
      re0 <- rbind(re0,re2)
    }
  }
  return(re0)
}

# begin
# if(!dir.exists("./RENC")){
#   dir.create("./RENC")
# }
files <- gsub("_tss.*","",list.files("./output/bedpe_tss/"))
tss <- read.table(tssfile)
tss <- tss[tss$V1 %in% c(paste0("chr",c(1:22,"X","Y"))),]
for (i in files){
  cellname <-i
  file1 <- paste0("./output/no_gene/",cellname,"_no_gene.bedpe")
  file2 <- paste0("./output/gene_tss_in_dup/",cellname,"_gene_tss_in_dup.bedpe")
  file3 <- paste0("./output/gene_tss_out_dup/",cellname,"_gene_tss_out_dup.bedpe")
  hichip.tss.file <- paste0("./output/bedpe_tss/",cellname,"_tss.bed")
  genelist_output <- paste0(cellname,"_RENC.txt")
  readornot <- (file.size(file1)==0 & file.size(file2) ==0 & (file3)==0)||file.size(hichip.tss.file)==0 
  if(!readornot){
    print(cellname)
    if(file.size(file1)==0){
      hichip1 <- as.data.frame(matrix(nrow = 1,ncol = 25,"."))
    }else{
      hichip1 <- read.table(file1)
    }
    
    if(file.size(file2)==0){
      hichip2 <- as.data.frame(matrix(nrow = 1,ncol = 25,"."))
    }else{
      hichip2 <- read.table(file2)
    }
    
    if(file.size(file3)==0){
      hichip3 <- as.data.frame(matrix(nrow = 1,ncol = 25,"."))
    }else{
      hichip3 <- read.table(file3)
      hichip3 <- filter_func(hichip2,hichip3)
    }
    hichip.tss <- read.table(hichip.tss.file)
    RENC <- cal_hichip_by_gene(hichip.tss,hichip1,hichip2,hichip3,tss,selectmax = F)
    write.table(RENC,genelist_output,row.names = F,quote = F,col.names = T,sep = "\t")
  }
  
  if(!enhg=="ALL_RESULT"){
    enhgs_output <- paste0(cellname,"_RENC_search.txt")
    RENC <- read.table(genelist_output,header = T)
    enhgs <- unlist(strsplit(enhg,";"))
    enhgs_result <- sel_gene_enh(RENC,tss,enhgs)
    write.table(enhgs_result,enhgs_output,row.names = F,quote = F,col.names = T,sep = "\t")
  }
}

