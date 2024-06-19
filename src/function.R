draw_gtrack <- function(ins,chr,chr.start,chr.end){
  ins <- ins[ins$V2>chr.start &ins$V3<chr.end,]
  ins <-  ins[ins$V5>chr.start &ins$V6<chr.end,]
  chr.start <- min(ins$V2)
  chr.end <- max(ins$V6)
  ins$fea1 <- paste(ins$V1,ins$V2,ins$V3,sep="_")
  ins$fea2 <- paste(ins$V4,ins$V5,ins$V6,sep="_")
  if(ncol(ins)==8){
    ins1 <- lapply(split(ins,list(ins$fea1,ins$fea2)), nrow)
    ins1 <- do.call(rbind,ins1)
    ins1 <- data.frame(ins1)
    ins1$fea1 <- unlist(lapply(strsplit(rownames(ins1),"\\."),function(x){x[1]}))
    ins1$fea2 <- unlist(lapply(strsplit(rownames(ins1),"\\."),function(x){x[2]}))
    ins1 <- ins1[!ins1$ins1==0,]
  }else if(ncol(ins)==9){
    ins1 <- ins[,7:9]
  }
  fea.all<- unique(c(ins1$fea1,ins1$fea2))
  fea.all <- fea.all[order(fea.all,decreasing = T)]
  mt <- matrix(0, nrow = length(fea.all), ncol = length(fea.all))
  rownames(mt) <- fea.all
  colnames(mt) <- fea.all
  for (i in 1:nrow(ins1)){
    NN1 <- which(rownames(mt)==ins1$fea1[i])
    NN2 <- which(colnames(mt)==ins1$fea2[i])
    mt[NN1,NN2] <- ins1[,1][i]
  }
  ss.chr <- unname(sapply(fea.all,function(x){strsplit(x,"_")[[1]][1]}))
  ss.start <- as.numeric(unname(sapply(fea.all,function(x){strsplit(x,"_")[[1]][2]})))
  ss.end <- as.numeric(unname(sapply(fea.all,function(x){strsplit(x,"_")[[1]][3]})))
  gr <- GRanges(seqnames = Rle(ss.chr), 
                ranges = IRanges(start = ss.start,
                                 end = ss.end))
  chr.range <- paste0(chr,":",min(ss.start),"-",max(ss.end))
  wins <- GRanges(chr.range)
  #dat.gt <- gTrack(gr, mdata = mt,colormaps = c("white","red"),cmap.min=1,cmap.max = 8)
  dat.gt <- gTrack(gr, mdata = mt,colormaps = c("white","#CD0D0C"),height = 10,ygap = 0.5)
  return(list("dat.gt"=dat.gt,
              "wins"=wins))
  #plot(dat.gt,wins)
}
gtrack_merge_plot <- function(Gene="KLF5",Chr="chr13",Chr.start = 73200000,Chr.end=74350000,
                              Height=2,Ygap=0.5,filename="HARA_KLF5_addgene_chip_dup_pps_enh1"){
  # add gene
  hara <- draw_gtrack(ins,chr=Chr,chr.start = Chr.start,chr.end=Chr.end)
  dat.gt <- hara[[1]]
  wins <- hara[[2]]
  ind <- which(names(gt.ge@data[[1]]) %in% genes$V4)
  gt.ge@data[[1]] <- gt.ge@data[[1]][ind] 
  gt.ge@formatting$height <- Height
  gt.ge@formatting$ygap <- Ygap
  # add chip 
  chip.gr.score <- gTrack(chip,col = 'grey', name = 'chip',bars = TRUE,height = Height, ygap = Ygap)
  #plot(c(chip.gr.score,gt.ge,dat.gt),wins)
  # add hotspot
  hs.gr <- GRanges(seqnames = Rle(hs$V1), 
                   ranges = IRanges(start = hs$V2,
                                    end = hs$V3))
  hs.gr <- gTrack(hs.gr,col = 'grey', name = '',bars = TRUE,height = Height/4, ygap =Ygap/4)
  # add dup 
  dup.gr <- GRanges(seqnames = Rle(dup$V1), 
                    ranges = IRanges(start = dup$V2,
                                     end = dup$V3),
                    score =dup$V4)
  dup.gr <- gTrack(dup.gr, y.field = 'score',
                   col = 'grey', name = 'dup',bars = TRUE,height = Height, ygap =Ygap)
  #plot(c(dup.gr,chip.gr.score,gt.ge,dat.gt),wins)
  # add pps score 
  pps1 <- pps[,c("ID","score")]
  pps1 <- merge(tss,pps1,by.x="V5",by.y="ID",all.x=T)
  pps1$score <- log2(pps1$score+0.01)
  pps1$score[is.na(pps1$score)] <- 0
  pps1 <- pps1[!(duplicated(pps1$V4) & pps1$score==0),]
  pps.gr <- GRanges(seqnames = Rle(pps1$V1), 
                    ranges = IRanges(start = pps1$V2,
                                     end = pps1$V3),
                    score =pps1$score)
  pps.gr <- gTrack(pps.gr, y.field = 'score',
                   col = 'grey', name = 'PPS',circles = TRUE,height = Height, ygap = Ygap)
  #plot(c(pps.gr,dup.gr,chip.gr.score,gt.ge,dat.gt),wins)
  # add enhancer
  bp1$inanchor1 <- bp1$V10 >bp1$V2 & bp1$V11 < bp1$V3
  bp1$inanchor2 <- bp1$V10 >bp1$V5 & bp1$V11 < bp1$V6
  bp1.1 <- bp1[bp1$inanchor1==1,]
  bp1.1 <- bp1.1[,c(4:6,7,12)]
  colnames(bp1.1) <- c("chr","start","end","pets","gene")
  bp1.2 <- bp1[bp1$inanchor2==1,]
  bp1.2 <- bp1.2[,c(1:3,7,12)]
  colnames(bp1.2) <- c("chr","start","end","pets","gene")
  bp2 <- unique(rbind(bp1.1,bp1.2))
  bp2 <- bp2[bp2$gene==Gene,]
  enr.gr <- GRanges(seqnames = Rle(bp2$chr), 
                    ranges = IRanges(start = bp2$start,
                                     end = bp2$end),
                    pets =bp2$pets)
  enr.gr <- gTrack(enr.gr, y.field = 'pets',
                   col = 'grey', name = 'pets',circles = TRUE, height = Height,ygap = Ygap)
  #plot(c(enr.gr,pps.gr,dup.gr,chip.gr.score,gt.ge,dat.gt),wins)
  pdf(paste0("./fig/",filename,".pdf"),height = 12,width = 8)
  plot(c(enr.gr,pps.gr,dup.gr,hs.gr,chip.gr.score,gt.ge,dat.gt),wins)
  dev.off()
}
find_enh_anchor <- function(bp1){
  bp1$inanchor1 <- bp1$V10 >bp1$V2 & bp1$V11 < bp1$V3
  bp1$inanchor2 <- bp1$V10 >bp1$V5 & bp1$V11 < bp1$V6
  bp1.1 <- bp1[bp1$inanchor1==1,]
  bp1.1 <- bp1.1[,c(4:6,7,12)]
  colnames(bp1.1) <- c("chr","start","end","pets","gene")
  bp1.2 <- bp1[bp1$inanchor2==1,]
  bp1.2 <- bp1.2[,c(1:3,7,12)]
  colnames(bp1.2) <- c("chr","start","end","pets","gene")
  bp <- unique(rbind(bp1.1,bp1.2))
  return(bp)
}
