load("data/PSets/CCLE_hs.RData")
load("data/PSets/GDSC.RData")
load("data/PSets/GRAY_hs.RData")
load("data/PSets/gCSI_hs.RData")
require(PharmacoGx)

source("code/foo.R")

##Supplementary Figure 4
##rnaseq/microarray concordance
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


###supplementary Figure 5
##intersection of cells/drugs between ccle, gdsc, gray
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

pdf(file.path(file.path(path.diagrams, "gray_ccle_gdsc_breast_celllines.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")]), 
                                          area2 = length(GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")]),
                                          area3 = length(GRAY@cell$cellid),
                                          n12 = length(intersect(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")], GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")])),
                                          n13 = length(intersect(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")], GRAY@cell$cellid)),
                                          n23 = length(intersect(GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")], GRAY@cell$cellid)),
                                          n123 = length(PharmacoGx::intersectList(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")], GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")], GRAY@cell$cellid)),
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
##inconsistency of pharmacological profiles for sorafenib and critozenib between 3 studies
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

#Supp figure 9
load(file.path(path.diagrams, "Biomarkers_uhn_status.RData"), verbose=TRUE)

results <- file.path(path.diagrams, "confidence_interval_difference.txt")
drug <- "lapatinib"
x <- which(biomarkers[[drug]][, "symbol"] == "ERBB2")
cutoff <- 0.1

gene <- biomarkers[[drug]][x, "symbol"] 
gene.id <- biomarkers[[drug]][x, "gene.id"]
best.isoform <- biomarkers[[drug]][x, "transcript.id"]
annot.isoforms <- annot.ensembl.all.transcripts[which(annot.ensembl.all.transcripts$gene_id == gene.id),]
annot.isoforms <- intersectList(rownames(annot.isoforms), colnames(ccle.isoforms.fpkm), colnames(gray.isoforms.fpkm), colnames(uhn.isoforms.fpkm))
annot.gene <- annot.ensembl.all.genes[which(annot.ensembl.all.genes$gene_id == gene.id), ]

cat(sprintf("%s, %s, %s, cutoff=%s\n", gene, best.isoform, drug, cutoff), file=results, append=TRUE)

## training 
ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.isoforms.fpkm))
ccle.breast.cells <- intersectList(rownames(CCLE@cell)[which(CCLE@cell$tissueid == "breast")], rownames(ccle.drug.sensitivity), rownames(ccle.isoforms.fpkm))
ccle.breast.cells <- ccle.breast.cells[which(!is.na(ccle.drug.sensitivity[ccle.breast.cells, drug]))]
ccle.auc <- ccle.drug.sensitivity[ccle.breast.cells, drug]
ccle.isoforms.models <- fnBuildLinearModel(cells=ccle.cells, sensitivity=ccle.drug.sensitivity[ccle.cells, drug], features=annot.isoforms, fpkm.matrix=ccle.isoforms.fpkm)

gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.isoforms.fpkm))
gdsc.breast.cells <- intersectList(rownames(CCLE@cell)[which(CCLE@cell$tissueid == "breast")],  rownames(gdsc.drug.sensitivity), rownames(ccle.isoforms.fpkm))
gdsc.breast.cells <- gdsc.breast.cells[which(!is.na(gdsc.drug.sensitivity[gdsc.breast.cells, drug]))]
gdsc.auc <- gdsc.drug.sensitivity[gdsc.breast.cells, drug]
gdsc.isoforms.models <- fnBuildLinearModel(cells=gdsc.cells, sensitivity=gdsc.drug.sensitivity[gdsc.cells, drug], features=annot.isoforms, fpkm.matrix=ccle.isoforms.fpkm)

training.isoforms.models <- fnNormalizeTraining(ccle.isoforms.models$simplified, gdsc.isoforms.models$simplified)
#training.isoforms.models <- ccle.isoforms.models$simplified

ccle.gene.model <- fnBuildLinearModel(cells=ccle.breast.cells, sensitivity=ccle.auc, features=annot.gene$gene_id, fpkm.matrix=ccle.genes.fpkm)
gdsc.gene.model <- fnBuildLinearModel(cells=gdsc.breast.cells, sensitivity=gdsc.auc, features=annot.gene$gene_id, fpkm.matrix=ccle.genes.fpkm)

training.gene.model <- fnNormalizeTraining(ccle.gene.model$simplified, gdsc.gene.model$simplified)
#training.gene.model <- ccle.gene.model$simplified
#training.models.isoforms <- training.models.isoforms[which(training.models.isoforms[,"pvalue"] < cutoff), , drop=FALSE]

expression <- cbind(ccle.isoforms.fpkm[ , annot.isoforms], "gene"=ccle.genes.fpkm[, annot.gene$gene_id])
colnames(expression)[ncol(expression)] <- annot.gene$gene_id

isoforms.ordered <- fnPlotHeatMap(sensitivity=ccle.auc, expression=expression, file.name=sprintf("%s_%s", drug, gene), cluster=TRUE, best.isoform=best.isoform)
isoforms.ordered <- fnPlotHeatMap.union.cells(ccle.sensitivity=ccle.auc, gdsc.sensitivity=gdsc.auc, expression=expression, file.name=sprintf("%s_%s", drug, gene), cluster=TRUE, best.isoform=best.isoform)
if(is.null(isoforms.ordered)) {
  isoforms.ordered <- annot.isoforms
}
fnPlotSensitivity(sensitivity=ccle.auc, file.name=sprintf("%s_%s", drug, gene), type="training")
fnPlotEstimates(isoform.models=training.isoforms.models, gene.model=training.gene.model, isoforms=isoforms.ordered, file.name=sprintf("%s_%s", drug, gene), cutoff=cutoff)
fnPlotPvalues(isoform.models=training.isoforms.models, gene.model=training.gene.model, isoforms=isoforms.ordered, file.name=sprintf("%s_%s", drug, gene), cutoff=cutoff)

cat(sprintf("training confidence interval difference pvalue: %s\n", fnConfidenceIntervalDifference(model1=training.isoforms.models[best.isoform,], 
                                                                                                   model2=training.gene.model, 
                                                                                                   x1=ccle.isoforms.models$complete[[best.isoform]]$model$`tt[, "exp"]`,
                                                                                                   x2=ccle.gene.model$complete[[1]]$model$`tt[, "exp"]`)), file=results, append=TRUE)

## in silico

gray.cells <- intersectList(rownames(gray.drug.sensitivity), rownames(gray.isoforms.fpkm))
gray.auc <- gray.drug.sensitivity[gray.cells, drug]

gray.isoforms.models <- fnBuildLinearModel(cells=rownames(gray.isoforms.fpkm), sensitivity=gray.auc, features=isoforms.ordered, fpkm.matrix=gray.isoforms.fpkm)
gray.gene.model <- fnBuildLinearModel(cells=rownames(gray.genes.fpkm), sensitivity=gray.auc, features=as.character(annot.gene$gene_id), fpkm.matrix=gray.genes.fpkm)

#in.silico.models.isoforms <- in.silico.models.isoforms[which(in.silico.models.isoforms[,"pvalue"] < cutoff & sign(in.silico.models.isoforms[,"estimate"]) == sign(training.models.isoforms[,"estimate"])), , drop=FALSE]
expression <- cbind(gray.isoforms.fpkm[ , isoforms.ordered], "gene"=gray.genes.fpkm[, as.character(annot.gene$gene_id)])
xx <- cor(expression[gray.cells,], expression[gray.cells,], use="pairwise", method="spearman")
fnCor(drug, gene, xx)
fnExp(drug, gene, gray.isoforms.fpkm[ , isoforms.ordered])
colnames(expression)[ncol(expression)] <- annot.gene$gene_id

xx <- fnPlotHeatMap(sensitivity=gray.auc, expression=expression, file.name=sprintf("%s_%s_in_silico", drug, gene), best.isoform=best.isoform)
fnPlotSensitivity(sensitivity=gray.auc, file.name=sprintf("%s_%s_in_silico", drug, gene))
fnPlotEstimates(isoform.models=gray.isoforms.models$simplified, gene.model=gray.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_silico", drug, gene), cutoff=cutoff)
fnPlotPvalues(isoform.models=gray.isoforms.models$simplified, gene.model=gray.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_silico", drug, gene), cutoff=cutoff)

cat(sprintf("in silico confidence interval difference pvalue: %s\n", fnConfidenceIntervalDifference(model1=gray.isoforms.models$simplified[best.isoform,], 
                                                                                                    model2=gray.gene.model$simplified, 
                                                                                                    x1=gray.isoforms.models$complete[[best.isoform]]$model$`tt[, "exp"]`,
                                                                                                    x2=gray.gene.model$complete[[1]]$model$`tt[, "exp"]`)), file=results, append=TRUE)

## in vitro

uhn.cells <- intersectList(rownames(uhn.drug.sensitivity), rownames(uhn.isoforms.fpkm))
uhn.auc <- uhn.drug.sensitivity[ , drug]

uhn.isoforms.models <- fnBuildLinearModel(cells=rownames(uhn.isoforms.fpkm), sensitivity=uhn.auc, features=isoforms.ordered, fpkm.matrix=uhn.isoforms.fpkm)
uhn.gene.model <- fnBuildLinearModel(cells=rownames(uhn.genes.fpkm), sensitivity=uhn.auc, features=as.character(annot.gene$gene_id), fpkm.matrix=uhn.genes.fpkm)

#in.vitro.models.isoforms <- in.vitro.models.isoforms[which(in.vitro.models.isoforms[,"pvalue"] < cutoff & sign(in.vitro.models.isoforms[,"estimate"]) == sign(in.silico.models.isoforms[,"estimate"])), , drop=FALSE]
expression <- cbind(uhn.isoforms.fpkm[ , isoforms.ordered], "gene"=uhn.genes.fpkm[, as.character(annot.gene$gene_id)])
xx <- expression[, c(best.isoform,"gene"), drop=FALSE]
fnPlotHeatMap(sensitivity=uhn.auc, file.name=sprintf("%s_%s_image.pdf", gene, drug), cluster=FALSE, expression=xx, best.isoform=best.isoform)
colnames(expression)[ncol(expression)] <- annot.gene$gene_id

xx <- fnPlotHeatMap(sensitivity=uhn.auc, expression=expression, file.name=sprintf("%s_%s_in_vitro", drug, gene), best.isoform=best.isoform)
fnPlotSensitivity(sensitivity=uhn.auc, file.name=sprintf("%s_%s_in_vitro", drug, gene))
fnPlotEstimates(isoform.models=uhn.isoforms.models$simplified, gene.model=uhn.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_vitro", drug, gene), cutoff=cutoff)
fnPlotPvalues(isoform.models=uhn.isoforms.models$simplified, gene.model=uhn.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_vitro", drug, gene), cutoff=cutoff)

cat(sprintf("in vitro confidence interval difference pvalue: %s\n", fnConfidenceIntervalDifference(model1=uhn.isoforms.models$simplified[best.isoform,], 
                                                                                                   model2=uhn.gene.model$simplified, 
                                                                                                   x1=uhn.isoforms.models$complete[[best.isoform]]$model$`tt[, "exp"]`,
                                                                                                   x2=uhn.gene.model$complete[[1]]$model$`tt[, "exp"]`)), file=results, append=TRUE)

#######
##Supplementary table 2
##validation rate of isformic biomarkers 
drugs <- unionList(intersectList(drugNames(GRAY),
                                 drugNames(CCLE),
                                 drugNames(GDSC)), 
                   intersectList(drugNames(gCSI),
                                 drugNames(CCLE),
                                 drugNames(GDSC)))
drugs <- sort(drugs)
xx <- matrix("-", ncol=6, nrow=length(drugs), dimnames=list(drugs, c("Compound", "Training", "Pre-validation Pan-Cancer","Training Breast Specific", "Pre-validation Breast", "Final validation")))
xx[drugs, "Compound"] <- drugs
load(file.path(path.diagrams, "all.biomarkers.RData"), verbose=TRUE)
xx[drugs, "Training"] <- sapply(all.biomarkers, function(x){length(which(x["type"]=="isoform"))})[drugs]
load(file.path(path.diagrams, "validated.biomarkers.gCSI.RData"), verbose=TRUE)
xx[intersect(drugs, drugNames(gCSI)), "Pre-validation Pan-Cancer"] <- sapply(validated.biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})[intersect(drugs, drugNames(gCSI))]
load(file.path(path.diagrams, "all.biomarkers.RData"), verbose=TRUE)
xx[drugs, "Training Breast Specific"] <- sapply(all.biomarkers, function(x){length(which(x["type"]=="isoform" & x["breast"] > 0.55))})[drugs]
load(file.path(path.diagrams, "Biomarkers_validated_breast_cindex_gray_pvalue.RData"), verbose=TRUE)
xx[intersect(drugs, drugNames(GRAY)), "Pre-validation Breast"] <- sapply(biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})[intersect(drugs, drugNames(GRAY))]
load(file.path(path.diagrams, "Biomarkers_uhn_status.RData"), verbose=TRUE)
xx[names(biomarkers), "Final validation"] <- sapply(biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})
xtable::print.xtable(xtable::xtable(xx, digits=0, align=c("l","l|","l","l||","l","l","l")), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "validation_rate.tex"), append=FALSE)
#######

#######
##Supplementary table 4
##known biomarkers
file.sensitivity.all.types <- "auc_recomputed_drug_association_mut_cnv.RData"
file.sensitivity <- "auc_recomputed_drug_association.RData"
load("data/training_ccle_gdsc_mut_cnv.RData")
annot.ensembl.all.genes.mut <- annot.ensembl.all.genes
annot.ensembl.all.isoforms.mut <- annot.ensembl.all.isoforms
load("data/training_ccle_gdsc.RData")
load(file.path(path.data, file.sensitivity.all.types), verbose=T)
drug.association.mut <- drug.association
drug.association.statistics.mut <- drug.association.statistics
drug.association.best.isoforms.mut <- drug.association.best.isoforms
load(file.path(path.data, file.sensitivity), verbose=T)
knownBiomarkersCheck <-
  function()
  {
    
    known.biomarkers <- read.csv(file="data/known.biomarkers.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, na.strings=c("", " ", "NA"))
    known.biomarkers <- known.biomarkers[which(!is.na(known.biomarkers$type) & known.biomarkers$type != "fusion"),1:3]
    #known.biomarkers <- known.biomarkers[which(known.biomarkers$type == "expression"),1:2]
    known.biomarkers <- cbind(known.biomarkers,  "gene cindex"=NA, "gene pvalue"=NA, "isoform"=NA, "isoform cindex"=NA, "isoform pvalue"=NA, "null model cindex"=NA)
    
    for(i in 1:nrow(known.biomarkers)) {
      feature <- rownames(annot.ensembl.all.genes)[which(annot.ensembl.all.genes$Symbol == known.biomarkers[i ,"gene"])]
      if(known.biomarkers[i,"type"]== "expression" & length(feature) > 0){
        known.biomarkers[i ,"gene cindex"] <- drug.association.statistics[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M2"]
        known.biomarkers[i ,"gene pvalue"] <- drug.association[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M2"]
        #known.biomarkers[i ,"gene R2"] <- drug.association.statistics[["r.squared"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M2"]
        
        known.biomarkers[i ,"isoform"] <- drug.association.best.isoforms[[feature]][[known.biomarkers[i,"drug"]]]
        
        known.biomarkers[i ,"isoform cindex"] <- drug.association.statistics[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M3B"]
        known.biomarkers[i ,"isoform pvalue"] <- drug.association[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M3B"]
        #known.biomarkers[i ,"isoform R2"] <- drug.association.statistics[["r.squared"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M3B"]
        
        known.biomarkers[i ,"null model cindex"] <- drug.association.statistics[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M0"]
        
        #xx <- rownames(ccle.drug.sensitivity)[which(!is.na(ccle.drug.sensitivity[,known.biomarkers[i,"drug"]]))]
        #xx <- intersect(xx, rownames(ccle.genes.fpkm))
        #known.biomarkers[i ,"mean expr"] <- mean(ccle.genes.fpkm[xx, feature], na.rm=TRUE)
        #breast.cells <- rownames(ccle.cell.profiles)[which(ccle.cell.profiles == "breast")]
        #breast.cells <- intersect(breast.cells, xx)
        #known.biomarkers[i ,"breast expr"] <- mean(ccle.genes.fpkm[breast.cells, feature], na.rm=TRUE)
        #known.biomarkers[i ,"mean isoform expr"] <- mean(ccle.isoforms.fpkm[xx, known.biomarkers[i ,"isoform"]], na.rm=TRUE)
        #known.biomarkers[i ,"breast isoform expr"] <- mean(ccle.isoforms.fpkm[breast.cells, known.biomarkers[i ,"isoform"]], na.rm=TRUE)
        
      }else{
        feature <- rownames(annot.ensembl.all.genes.mut)[which(annot.ensembl.all.genes.mut$Symbol == known.biomarkers[i ,"gene"])]
        
        known.biomarkers[i ,"null model cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M0"]
        if(known.biomarkers[i,"type"]== "expression"){
          known.biomarkers[i ,"gene cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M2"]
          known.biomarkers[i ,"gene pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M2"]
          
          known.biomarkers[i ,"isoform"] <- drug.association.best.isoforms.mut[[feature]][[known.biomarkers[i,"drug"]]]
          
          known.biomarkers[i ,"isoform cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M3B"]
          known.biomarkers[i ,"isoform pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M3B"]
        }else if(known.biomarkers[i,"type"]== "mutation"){
          known.biomarkers[i ,"gene cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "MM"]
          known.biomarkers[i ,"gene pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "MM"]
        }else if(known.biomarkers[i,"type"]== "amplification"){
          known.biomarkers[i ,"gene cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "MC"]
          known.biomarkers[i ,"gene pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "MC"]
        }
      }
    }
    known.biomarkers <- known.biomarkers[which(known.biomarkers$`gene cindex` != 0), ]
    #colnames(known.biomarkers)[1:3] <- capitalize(colnames(known.biomarkers)[1:3])
    xtable::print.xtable(xtable::xtable(known.biomarkers, digits=c(0, 0, 0, 0, 2, -1, 0, 2, -1, 2)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "known_biomarkers.tex"), append=FALSE)
    write.csv(known.biomarkers, file=file.path(path.diagrams, "known_biomarkers.csv"), row.names = F)
  }

knownBiomarkersCheck()
#######
#######
##Supplementary table 5
##Novel isoformic biomarkers
load(file.path(path.diagrams, "isoform_specific/Biomarkers_uhn_status.RData"), verbose=TRUE)
rr <- c("drug", "gene", "isoform", "Training cindex", "Training corrected pvalue", "Final validation cindex")
xx <- as.data.frame(matrix(NA, ncol=length(rr), nrow=sum(sapply(biomarkers,dim)[1,]), dimnames=(list(unlist(sapply(biomarkers,function(x){x["biomarker.id"]})), rr))))
i <- 0
drugs <- sort(names(biomarkers))
for(drug in drugs){
  rr <- order(biomarkers[[drug]]$UHN.cindex, decreasing=T)
  for(j in rr){
    i <- i + 1
    xx[i ,"drug"] <- drug
    isoform <- biomarkers[[drug]][j,"biomarker.id"]
    #feature <- biomarkers[[drug]][i, "gene.id"]
    feature <- annot.ensembl.all.isoforms[which(annot.ensembl.all.isoforms$EnsemblTranscriptId == isoform), "EnsemblGeneId"]
    
    xx[i ,"gene"] <- annot.ensembl.all.genes[feature, "Symbol"]
    #xx[i ,"gene cindex"] <- drug.association.statistics[["cindex"]][[feature]][[drug]]["mean", "M2"]
    #xx[i ,"gene pvalue"] <- drug.association[["cindex"]][[feature]][[drug]][1, "M2"]
    xx[i ,"isoform"] <- isoform
    xx[i ,"Training cindex"] <- biomarkers[[drug]][j,"cindex"]
    xx[i ,"Training corrected pvalue"] <- biomarkers[[drug]][j,"fdr"]
    xx[i ,"Final validation cindex"] <- round(as.numeric(biomarkers[[drug]][j,"UHN.cindex"]), digits=2)
    
  }
}
#colnames(known.biomarkers)[1:3] <- capitalize(colnames(known.biomarkers)[1:3])
xtable::print.xtable(xtable::xtable(xx, digits=c(0, 0, 0, 0, 2, -1, 2)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "novel_biomarkers.tex"), append=FALSE)

###
###Supplementary figure 13
###the heatmaps for the other drugs not included in uhn (not in Fig 4)
load("data/PSets/UHN_hs.RData")
source("code/foo.R")
source("code/foo_PreValidation.R")
source("code/foo_FinalValidation.R")
mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
red <- mycol[1]  
blue <- mycol[2]

load(file.path(path.diagrams, "Biomarkers_validated_breast_cindex_gray_pvalue.RData"), verbose=TRUE)
drugs <- setdiff(names(biomarkers), drugNames(UHN))
for(drug in drugs) {
  if(nrow(biomarkers[[drug]]) > 0){
    vtt <- biomarkers[[drug]]
    vtt <- vtt[which(vtt[,"type"] == "isoform"), , drop=F]
    vtt <- vtt[which(vtt[, "isoforms.no"] > 1), , drop=FALSE]
    vtt <- vtt[which(vtt[,"gray.specificity"] != "gene.specific"), , drop=F]
    gray.specificity <- vtt[ ,"gray.specificity"]
    xx <- fnFetchBiomarkers(top.significant.biomarkers=vtt, drug=drug, indices=1:nrow(vtt))
    xx <- do.call(rbind, xx)
    xx[,"short.label"] <- gsub(".ISO$","",xx[,"short.label"])
    xx <- apply(xx, 1, function(x){x})
    
    rr <- fnPlotAUCoverCellLinesGray(drug=drug, tissue.type="all", biomarkers=xx, suffix="isoform.specific", gray.specificity=gray.specificity)#, biomarkers.toPlot)
    if(all(!is.na(rr))){
      biomarkers.order <- rr$hv$rowInd
      fnPlotEffectSize(drug, biomarkers=xx, effect.size=effect.size, biomarkers.order) 
      xx <- do.call(rbind, xx)
      names(biomarkers.order) <- xx[biomarkers.order, "isoform.id"]
    }
    
  }
}
###supp File 3
##list of pre validated biomarkers in gray
load(file.path(path.diagrams, "Biomarkers_validated_breast_cindex_gray_pvalue.RData"), verbose=TRUE)
WriteXLS::WriteXLS("biomarkers", ExcelFileName=file.path(path.diagrams, "prevalidated_biomarkers_gray.xlsx"), row.names=TRUE)
load(file.path(path.diagrams, "all.biomarkers.RData"))
WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=file.path(path.diagrams, "predicted_biomarkers_training.xlsx"), row.names=TRUE)
load(file.path(path.diagrams, "isoform_specific/Biomarkers_uhn_status.RData"), verbose=TRUE)
WriteXLS::WriteXLS("biomarkers", ExcelFileName=file.path(path.diagrams, "validated_biomarkers_uhn.xlsx"), row.names=TRUE)
load(file.path(path.diagrams, "validated.biomarkers.gCSI.RData"), verbose=TRUE)
WriteXLS::WriteXLS("validated.biomarkers", ExcelFileName=file.path(path.diagrams, "validated_biomarkers_gCSI.xlsx"), row.names=TRUE)

###
###Check if the number of isoform specific pre validated biomarkers are greater than gene specific ones
gg <- NULL
ii <- NULL
for(drug in names(biomarkers)) {
  if(nrow(biomarkers[[drug]]) > 0){
    xx <- table(biomarkers[[drug]]$gray.specificity)
    gg <- c(gg, ifelse("gene.specific" %in% names(xx), xx["gene.specific"], 0))
    ii <- c(ii, sum(xx)-ifelse("gene.specific" %in% names(xx), xx["gene.specific"], 0))
  }
}
pp <- wilcox.test(ii, gg, paired=TRUE, alternative="greater")$p.value



####pre validation rates
load("result/auc_recomputed_ccle_gdsc/Biomarkers_uhn_status.RData")
