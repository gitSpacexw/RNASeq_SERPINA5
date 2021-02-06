### perform cqn normalization
### generate glm offset matrix for performing DE analysis using edgeR

library(cqn)
options(stringsAsFactors=F)

## read raw count table and GC information
raw        <- read.delim("./data/RawCount_55sample.txt", header=T, sep="\t")
df.GC      <- read.delim("./data/Ensembl_hg19.Gene_featuers.GC_Content.txt", header=T, sep="\t")
df.GC      <- df.GC[match(raw$GeneID, df.GC$GeneID), ]

## cqn normalization
mm.raw           <- as.matrix(raw[, 7:ncol(raw)])
rownames(mm.raw) <- raw$GeneID
idx.0            <- which(rowSums(mm.raw) < 1)

set.seed(12357)
cqn_rr    <-  cqn(mm.raw[-idx.0, ], lengths=raw$CodingLength[-idx.0], x=df.GC$GC_Content[-idx.0])
CQN       <-  cqn_rr$y + cqn_rr$offset
glmOffset <-  cqn_rr$glm.offset

## write CQN values and GLM offset values to external files
#df.CQN          <- cbind(raw[-idx.0, 1:6], CQN)
#df.glmOffset    <- cbind(raw[-idx.0, 1:6], glmOffset)
write.table(CQN,       file="./CQN.txt",        sep="\t", quote=F, row.names=T, col.names=T)
write.table(glmOffset, file="./GLM_Offset.txt", sep="\t", quote=F, row.names=T, col.names=T)
