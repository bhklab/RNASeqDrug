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

