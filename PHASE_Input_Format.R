library(snpStats)
library(data.table)
plinky <- "/exports/igmm/eddie/haley-soay/Imputation/sheep_geno_imputed_Plates1to91_20220608/Input/Plates_1-2_HD_QC2"
plinkfile <- read.plink(fam=paste0(plinky,".fam"),bed=paste0(plinky,".bed"),bim=paste0(plinky,".bim"))

snps <- as.data.frame(plinkfile$genotypes)
for (i in 1:26){chrs<-plinkfile$map[which(plinkfile$map$chromosome == i),]
bottom_output <- snps[,which(colnames(snps) %in% chrs$snp.name)]
bottom_output <- data.frame(apply(bottom_output,2,function(x) as.numeric(as.character(x))),check.names=F,row.names=row.names(bottom_output))
bottom_output[bottom_output == 0] <- NA
bottom_output <- bottom_output -1
bottom_output <- bottom_output[rep(seq_len(nrow(bottom_output)),each=3),]
bottom_output[seq(3,nrow(bottom_output),3),] <- floor(bottom_output[seq(3,nrow(bottom_output),3),]/2)
bottom_output[bottom_output == 2] <- 1
bottom_output[seq(1,nrow(bottom_output),3),1] <- paste0("#",rownames(bottom_output[seq(1,nrow(bottom_output),3),]))
bottom_output[seq(1,nrow(bottom_output),3),2:ncol(bottom_output)] <- ""
bottom_output[is.na(bottom_output)] <- "?"

top_output <- as.data.frame(t(data.frame(A=c(nrow(plinkfile$fam),rep("",nrow(plinkfile$map[which(plinkfile$map$chromosome == i),]))),
                         B=c(nrow(plinkfile$map[which(plinkfile$map$chromosome == i),]),rep("",nrow(plinkfile$map[which(plinkfile$map$chromosome == i),]))),
                         C=c("P",plinkfile$map[which(plinkfile$map$chromosome == i),"position"]),
                         D=c(rep("S",nrow(plinkfile$map[which(plinkfile$map$chromosome == i),])),""),check.rows=FALSE)))
colnames(output) <- 1:ncol(output)
colnames(bottom_output) <- 1:ncol(bottom_output)

top_output <- rbindlist(list(top_output,bottom_output),fill=TRUE)
top_output[is.na(top_output)] <- ""

write.table(top_output,paste0("Chr",i,".inp"),col.names=FALSE,row.names=FALSE,quote=FALSE)

}
