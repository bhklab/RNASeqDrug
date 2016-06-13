load("data/PSets/CCLE_isoforms.RData")
load("data/PSets/GDSC.RData")
load("data/PSets/GRAY_isoforms.RData")

source("code/foo.R")
mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")

library(VennDiagram)
pdf(file.path(file.path(path.diagrams, "gray_ccle_gdsc_celllines.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(CCLE@cell$cellid), 
                                          area2 = length(GDSC@cell$cellid),
                                          area3 = length(GRAY@cell$cellid),
                                          n12 = length(intersect(CCLE@cell$cellid, GDSC@cell$cellid)),
                                          n13 = length(intersect(CCLE@cell$cellid, GRAY@cell$cellid)),
                                          n23 = length(intersect(GDSC@cell$cellid, GRAY@cell$cellid)),
                                          n123 = length(PharmacoGx::intersectList(CCLE@cell$cellid, GDSC@cell$cellid, GRAY@cell$cellid)),
                                          category = c("CCLE", "GDSC", "GRAY"),
                                          col = mycol[1:3],
                                          fill = mycol[1:3],
                                          margin=0.10,
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()


pdf(file.path(file.path(path.diagrams, "gray_ccle_gdsc_drugs.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(rownames(CCLE@drug)), 
                                          area2 = length(rownames(GDSC@drug)),
                                          area3 = length(rownames(GRAY@drug)),
                                          n12 = length(intersect(rownames(CCLE@drug), rownames(GDSC@drug))),
                                          n13 = length(intersect(rownames(CCLE@drug), rownames(GRAY@drug))),
                                          n23 = length(intersect(rownames(GDSC@drug), rownames(GRAY@drug))),
                                          n123 = length(PharmacoGx::intersectList(rownames(CCLE@drug), rownames(GDSC@drug), rownames(GRAY@drug))),
                                          category = c("CCLE", "GDSC", "GRAY"),
                                          col = mycol[1:3],
                                          fill = mycol[1:3],
                                          margin=0.10,
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()

##Supplementary Figure 6
if(file.exists(file.path(path.data, "rnaseq_microarray.RData"))){
  load(file.path(path.data, "rnaseq_microarray.RData"))
}else{
  cells <- intersect(pData(CCLE@molecularProfiles$rnaseq)[,"cellid"], pData(CCLE@molecularProfiles$rna)[,"cellid"])
  rnaseq.samples <- rownames(pData(CCLE@molecularProfiles$rnaseq))[match(cells, pData(CCLE@molecularProfiles$rnaseq)[,"cellid"])]
  rna.samples <- rownames(pData(CCLE@molecularProfiles$rna))[match(cells, pData(CCLE@molecularProfiles$rna)[,"cellid"])]
  
  
  features <- intersect(fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"], fData(CCLE@molecularProfiles$rna)[,"EnsemblGeneId"])
  rnaseq.features <- rownames(fData(CCLE@molecularProfiles$rnaseq))[match(features, fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"])]
  rna.features <- rownames(fData(CCLE@molecularProfiles$rna))[match(features, fData(CCLE@molecularProfiles$rna)[,"EnsemblGeneId"])]
  ccle.rnaseq.microarray.cor <- cor(exprs(CCLE@molecularProfiles$rna)[rna.features, rna.samples], exprs(CCLE@molecularProfiles$rnaseq)[rnaseq.features, rnaseq.samples], use="pairwise.complete.obs", method="spearman")
  save(ccle.rnaseq.microarray.cor, file=file.path(path.data, "rnaseq_microarray.RData"))
}

pdf(file.path(path.diagrams, "rnaseq_microarray.pdf"), height=7, width=7)
par(mar=c(9,5,5,2))
boxplot(cbind("Identical"=diag(ccle.rnaseq.microarray.cor), "Different"=c(ccle.rnaseq.microarray.cor[upper.tri(ccle.rnaseq.microarray.cor)], ccle.rnaseq.microarray.cor[lower.tri(ccle.rnaseq.microarray.cor)])), las = 2, col = "gray", cex.lab=1, cex.axis=1, pch=19, ylab="spearman correlation", outpch=20, outcex=0.5)
dev.off()

##Supplementary Figure 7
gray.ccle.cells <- intersect(cellNames(GRAY), cellNames(CCLE))
gray.gdsc.cells <- intersect(cellNames(GRAY), cellNames(GDSC))
gdsc.ccle.cells <- intersect(cellNames(GDSC), cellNames(CCLE))

gray.sensitivity <- PharmacoGx::summarizeSensitivityProfiles(GRAY, sensitivity.measure = "auc_recomputed", drugs = c("Crizotinib", "Sorafenib"))
ccle.sensitivity <- PharmacoGx::summarizeSensitivityProfiles(CCLE, sensitivity.measure = "auc_recomputed", drugs = c("Crizotinib", "Sorafenib"))
gdsc.sensitivity <- PharmacoGx::summarizeSensitivityProfiles(GDSC, sensitivity.measure = "auc_recomputed", drugs = c("Crizotinib", "Sorafenib"))

myScatterPlot(Name=file.path(path.diagrams, "gray_ccle_sorafenib.pdf"), 
              x=ccle.sensitivity["Sorafenib", gray.ccle.cells], 
              y=gray.sensitivity["Sorafenib", gray.ccle.cells], 
              method="plain", 
              minp=10,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "gray_gdsc_sorafenib.pdf"), 
              x=gdsc.sensitivity["Sorafenib", gray.gdsc.cells], 
              y=gray.sensitivity["Sorafenib", gray.gdsc.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="GDSC",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "ccle_gdsc_sorafenib.pdf"), 
              x=ccle.sensitivity["Sorafenib", gdsc.ccle.cells], 
              y=gdsc.sensitivity["Sorafenib", gdsc.ccle.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GDSC")


myScatterPlot(Name=file.path(path.diagrams, "gray_ccle_crizotinib.pdf"), 
              x=ccle.sensitivity["Crizotinib", gray.ccle.cells], 
              y=gray.sensitivity["Crizotinib", gray.ccle.cells], 
              method="plain", 
              minp=10,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "gray_gdsc_crizotinib.pdf"), 
              x=gdsc.sensitivity["Crizotinib", gray.gdsc.cells], 
              y=gray.sensitivity["Crizotinib", gray.gdsc.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="GDSC",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "ccle_gdsc_crizotinib.pdf"), 
              x=ccle.sensitivity["Crizotinib", gdsc.ccle.cells], 
              y=gdsc.sensitivity["Crizotinib", gdsc.ccle.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GDSC")
myScatterPlotle.path(path.diagrams, "microarray_rnaseq_cor.pdf"), method="transparent")

## doubling time
xx <- read.csv(file.path(path.data, "GrowthCurve_dataCombined.csv"))
dec.batch <- xx[3:10,2:6]
colnames(dec.batch) <- xx[2,2:6]
rownames(dec.batch) <- xx[3:10,1]

nov.batch.1 <- xx[14:18,2:17]
colnames(nov.batch.1) <- xx[13,2:17]
rownames(nov.batch.1) <- xx[14:18,1]

nov.batch.2 <- xx[c(22:26, 28:29),2:15]
colnames(nov.batch.2) <- xx[21,2:15]
rownames(nov.batch.2) <- xx[c(22:26, 28:29),1]

nov.batch.3 <- xx[c(33:35, 38:41),2:16]
colnames(nov.batch.3) <- xx[32,2:16]
rownames(nov.batch.3) <- xx[c(33:35, 38:41),1]

nov.batch.4 <- xx[44:50,2:14]
colnames(nov.batch.4) <- xx[43,2:14]
rownames(nov.batch.4) <- xx[44:50,1]

doubling.time.dec.batch <- function(initial.conc, final.conc, culture.duration) {
  return((culture.duration*log10(2))/(log10(final.conc)- log10(initial.conc)))
}
dec.batch.dt <- sapply(1:ncol(dec.batch), function(x){doubling.time.dec.batch(as.numeric(dec.batch[1,x]), as.numeric(dec.batch[nrow(dec.batch), x]), 8)})
names(dec.batch.dt) <- colnames(dec.batch)
names(dec.batch.dt)[which(names(dec.batch.dt) == "HDQP1")] <- "HDQP-1"

nov.batch.1.dt <- sapply(1:ncol(nov.batch.1), function(x){doubling.time.dec.batch(as.numeric(nov.batch.1[1,x]), as.numeric(nov.batch.1[nrow(nov.batch.1), x]), 5)})
names(nov.batch.1.dt) <- colnames(nov.batch.1)

nov.batch.2.dt <- sapply(1:ncol(nov.batch.2), function(x){doubling.time.dec.batch(as.numeric(nov.batch.2[1,x]), as.numeric(nov.batch.2[nrow(nov.batch.2), x]), 8)})
names(nov.batch.2.dt) <- colnames(nov.batch.2)
names(nov.batch.2.dt)[which(names(nov.batch.2.dt) == "MFN223")] <- "MFM223"

nov.batch.3.dt <- sapply(1:ncol(nov.batch.3), function(x){doubling.time.dec.batch(as.numeric(nov.batch.3[1,x]), as.numeric(nov.batch.3[nrow(nov.batch.3), x]), 8)})
names(nov.batch.3.dt) <- colnames(nov.batch.3)
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "HS578T")] <- "Hs578T"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "JMT1")] <- "JIMT1"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "OCUB")] <- "OCUB1"

nov.batch.4.dt <- sapply(1:ncol(nov.batch.4), function(x){doubling.time.dec.batch(as.numeric(nov.batch.4[1,x]), as.numeric(nov.batch.4[nrow(nov.batch.4), x]), 7)})
names(nov.batch.4.dt) <- colnames(nov.batch.4)
xx <- sort(unionList(names(dec.batch.dt), names(nov.batch.1.dt), names(nov.batch.2.dt), names(nov.batch.3.dt), names(nov.batch.4.dt)))
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[\xa0]|[ ]"
yy <- tolower(gsub(badchars, "", xx))
yy <- gsub("mda", "mdamb", yy)
yy[which(yy == "mpe600")] <- "600mpe"
cells <- read.csv("~/Documents/DRUGNET/Curation/output/cell_annotation_all.csv", na.strings= c("", " "), stringsAsFactors=FALSE)
cellid <- cells[match(yy, cells$Ben_Neel.cellid), "unique.cellid"]
cellid[is.na(cellid)] <- xx[is.na(cellid)]

tt <- NULL
for(cell in xx){
  t <- NA
  if(cell %in% names(dec.batch.dt)) {
    x <- dec.batch.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.1.dt)) {
    x <- nov.batch.1.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.2.dt)) {
    x <- nov.batch.2.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.3.dt)) {
    x <- nov.batch.3.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.4.dt)) {
    x <- nov.batch.4.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  names(t) <- cell
  tt <- c(tt, round(t, digits=2))
}

dd <- cbind("Cell line"=cellid, "Doubling time"=tt)
xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "doubling_time.tex"), append=FALSE)

