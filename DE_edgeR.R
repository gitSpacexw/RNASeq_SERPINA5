library(edgeR)
library(matrixStats)

options(stringsAsFactors=F)

## read in sample information, raw count and glm offset table
sampleInfor   <- read.delim("./data/sampleInfor.txt", header=T, sep="\t")
raw           <- read.delim("./data/RawCount_55sample.txt", header=T, sep="\t")
CQN           <- read.delim("./data/CQN.txt",  header=T, sep="\t")
glmOffset     <- read.delim("./data/GLM_Offset.txt", header=T, sep="\t")
rownames(raw) <- raw$GeneID

## match tables with the samples ordered in sampleInfor and genes
mm.glmOffset <- as.matrix(glmOffset[, match(sampleInfor$SampleID, colnames(glmOffset))])
mm.CQN       <- as.matrix(CQN[, match(sampleInfor$SampleID, colnames(CQN))])
mm.raw       <- as.matrix(raw[rownames(CQN), match(sampleInfor$SampleID, colnames(raw))])

## DE analysis using edgeR
group               <- factor( gsub("-",".", sampleInfor$AD_type))
RIN                 <- sampleInfor$RIN

dge                 <- DGEList(counts=mm.raw, group=group)
design              <- model.matrix(~ 0+group+RIN)
colnames(design)    <- c("Ctrl", "HpSp","Limbic", "Typical", "RIN")
dge$offset          <- mm.glmOffset


trend               <- estimateGLMTrendedDisp(dge, design )
tag                 <- estimateGLMTagwiseDisp(trend, design )
fit                 <- glmFit(tag, design, dispersion=tag$tagwise.dispersion)

my.contrast  <- makeContrasts(TypicalVsCtrl=Typical-Ctrl, LimbicVsHpSp=Limbic-HpSp, levels=design)
lrt.Typical_Ctrl    <- glmLRT(fit, contrast=my.contrast[, "TypicalVsCtrl"])
lrt.Limbic_HpSp     <- glmLRT(fit, contrast=my.contrast[, "LimbicVsHpSp"])

## prepare output files
x.Typical_Ctrl   <- topTags(lrt.Typical_Ctrl,  n=nrow(CQN))
x.Limbic_HpSp    <- topTags(lrt.Limbic_HpSp,   n=nrow(CQN))
df.Typical_Ctrl  <- data.frame(raw[match(rownames(x.Typical_Ctrl), raw$GeneID), 1:6], x.Typical_Ctrl)
df.Limbic_HpSp   <- data.frame(raw[match(rownames(x.Limbic_HpSp),  raw$GeneID), 1:6], x.Limbic_HpSp)

## annotate genes that reliably expressed (mean cqn >0) in either groups
df.Typical_Ctrl$mean.Typical <-rowMeans(mm.CQN[df.Typical_Ctrl$GeneID, subset(sampleInfor,AD_type=="Typical")$SampleID])
df.Typical_Ctrl$mean.Ctrl    <-rowMeans(mm.CQN[df.Typical_Ctrl$GeneID, subset(sampleInfor,AD_type=="Ctl")$SampleID])
df.Limbic_HpSp$mean.Limbic   <-rowMeans(mm.CQN[df.Limbic_HpSp$GeneID, subset(sampleInfor, AD_type=="Limbic")$SampleID])
df.Limbic_HpSp$mean.HpSp     <-rowMeans(mm.CQN[df.Limbic_HpSp$GeneID, subset(sampleInfor, AD_type=="HpSp")$SampleID])
df.Typical_Ctrl$expressed <- ifelse(df.Typical_Ctrl$mean.Typical>0 | df.Typical_Ctrl$mean.Ctrl>0, 1, 0)
df.Limbic_HpSp$expressed  <- ifelse(df.Limbic_HpSp$mean.Limbic>0   | df.Limbic_HpSp$mean.HpSp>0,  1, 0)


## calculate FDR for genes that are reliably expressed
df.Typical_Ctrl$FDR.expressed <- rep(NA, nrow(df.Typical_Ctrl))
df.Limbic_HpSp$FDR.expressed  <- rep(NA, nrow(df.Limbic_HpSp))
df.Typical_Ctrl$FDR.expressed[which(df.Typical_Ctrl$expressed==1)] <-
                     p.adjust(subset(df.Typical_Ctrl, expressed==1)$PValue, method = "BH")
df.Limbic_HpSp$FDR.expressed[which(df.Limbic_HpSp$expressed==1)]   <-
                     p.adjust(subset(df.Limbic_HpSp, expressed==1)$PValue,  method = "BH")

### prepare output data frame
write.table(df.Typical_Ctrl, file="./data/DE_Typical_Ctrl.txt", sep="\t", quote=F, row.names=F, col.names=T)
write.table(df.Limbic_HpSp,  file="./data/DE_Limbic_HpSp.txt",   sep="\t", quote=F, row.names=F, col.names=T)

