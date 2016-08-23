require(PharmacoGx) || stop("Library PharmacoGx is not available!")
require(Biobase) || stop("Library Biobase is not available!")
require(gdata) || stop("Library gdata is not available!")
require(genefu) || stop("Library genefu is not available!")
require(survcomp) || stop("Library survcomp is not available!")

source("code/foo_FinalValidation.R")
source("code/foo_PreValidation.R")
adjustment.method <- "fdr"

phenotype <- sensitivity.type <- "auc_recomputed"
path.data <- "data"
path.code <- file.path("code")
path.result <- file.path("result")
path.diagrams <- file.path(sprintf("result/%s", sensitivity.type))

validation.method <- "R2"

load(file.path("data/annotation.RData"), verbose = T)

load(file.path(path.data, "PSets/CCLE_isoforms.RData"))
load(file.path(path.data, "PSets/GDSC.RData"))
load(file.path(path.data, "PSets/GRAY_isoforms.RData"))
load(file.path(path.data, "PSets/UHN.RData"))

gdsc.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC, sensitivity.measure=sensitivity.type))
ccle.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure=sensitivity.type))
ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", fill.missing=FALSE)))
ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", fill.missing=FALSE)))


gray.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GRAY, sensitivity.measure=sensitivity.type))
gray.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="rnaseq", fill.missing=FALSE)))
gray.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="isoforms", fill.missing=FALSE)))

uhn.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=UHN, sensitivity.measure=sensitivity.type))
uhn.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=UHN, mDataType="rnaseq", fill.missing=FALSE)))
uhn.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=UHN, mDataType="isoforms", fill.missing=FALSE)))

mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
mycol3 <- RColorBrewer::brewer.pal(n=4, name="Set3")

## Heatmaps, effect sizes and sensitivity plots for final validation, just isoforms considered
cutoff <- 0.1
red <- mycol[1]  
blue <- mycol[2]
load(file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_2_%s.RData", validation.method)), verbose=TRUE)
drugs <- colnames(uhn.drug.sensitivity)
for(drug in drugs) {
  biomarkers[[drug]][, c("UHN.estimate", "UHN.pvalue", "UHN.R2", "gene.biotype")] <- NA
  biomarkers[[drug]][,"id"] <- biomarkers[[drug]]$gene.id
  biomarkers[[drug]][which(biomarkers[[drug]]$type == "isoform"),"id"] <- biomarkers[[drug]][which(biomarkers[[drug]]$type == "isoform"), "transcript.id"]
  tt <- biomarkers[[drug]]
  tt <- subset(tt,  tt$ccle > 0)
  vtt <- subset(tt, tt$validation.stat == "validated")
  vtt <- subset(vtt, vtt$gray.specificity == "isoform.specific")
  #vtt <- subset(vtt, vtt$type == "isoform")
  
  xx <- fnFetchBiomarkers(top.significant.biomarkers=vtt, drug=drug, indices=1:nrow(vtt))
  xx <- do.call(rbind, xx)
  xx[,"short.label"] <- gsub(".ISO$","",xx[,"short.label"])
  xx <- apply(xx, 1, function(x){x})
  rr <- fnPlotAUCoverCellLinesGray(drug=drug, tissue.type="all", biomarkers=xx, suffix="isoform.specific")#, biomarkers.toPlot)
  
  if(all(!is.na(rr))){
    biomarkers.order <- rr$hv$rowInd
    xx <- do.call(rbind, xx)
    names(biomarkers.order) <- xx[biomarkers.order, "isoform.id"]
  }
  
  exp.db <- NULL
  if(length(which(vtt$type == "isoform")) > 0) {
    exp.db <- uhn.isoforms.fpkm[ , vtt[which(vtt$type == "isoform"), "transcript.id"], drop=FALSE]
  }
  if(!is.null(exp.db)) {
    exp.db <- cbind(exp.db, uhn.genes.fpkm[ , vtt[which(vtt$type == "gene"), "gene.id"], drop=FALSE])
  } else {
    exp.db <- uhn.genes.fpkm[ , vtt[which(vtt$type == "gene"), "gene.id"], drop=FALSE]
  }
  coding.biotypes <- c("PRT", "AS", "prcTR")
  names(coding.biotypes) <- c("protein_coding", "antisense", "processed_transcript")
  bb <- vtt[match(colnames(exp.db), vtt$biomarker.id), "biotype"]
  bb <- sapply(bb, function(x){ifelse(x %in% names(coding.biotypes), coding.biotypes[x], x)})
  xx <- sprintf("%s (%s)", vtt[match(colnames(exp.db), vtt$biomarker.id), "symbol"], bb)
  names(xx) <- colnames(exp.db)
  sensitivity <- uhn.drug.sensitivity[!is.na(uhn.drug.sensitivity[ , drug]) , drug]
  exp.db <- exp.db[names(sensitivity), , drop=FALSE]
  uhn.models <- matrix(NA, ncol=4, nrow=ncol(exp.db), dimnames=list(colnames(exp.db), c("pvalue", "estimate", "R2", "gene_biotype")))
  for(marker in colnames(exp.db)) {
    uhn.model <- lm(sensitivity  ~ exp.db[, marker])
    uhn.pvalue <- 2 
    uhn.estimate <- uhn.R2 <- 0
    if(all(!is.na(uhn.model)) & !is.na(uhn.model$coefficients[2]))
    {
      uhn.models[marker,"pvalue"] <- summary(uhn.model)$coefficients[2,4]
      uhn.models[marker,"estimate"] <- summary(uhn.model)$coefficients[2,1]
      uhn.models[marker,"R2"] <- summary(uhn.model)$adj.r.squared
      uhn.models[marker,"gene_biotype"] <- annot.ensembl.all.transcripts[marker, "gene_biotype"]
    }
  }
  if(!is.null(rownames(uhn.models))) {
    biomarkers[[drug]][match(rownames(uhn.models), biomarkers[[drug]][,"id"]), c("UHN.estimate", "UHN.pvalue", "UHN.R2", "gene.biotype")] <- uhn.models[, c("estimate", "pvalue", "R2", "gene_biotype")]
    uhn.models <- uhn.models[which(as.numeric(uhn.models[,"pvalue"]) < cutoff & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
  } 
  message(drug)
  message(nrow(vtt))
  message(nrow(uhn.models))
  if(!is.null(rownames(uhn.models))) {
    exp.db <- exp.db[ ,rownames(uhn.models), drop=FALSE]
    oo <- order(sensitivity)
    exp.db <- exp.db[names(sensitivity)[oo], ,drop=FALSE]
    sensitivity <- sensitivity[oo]
    exp.col = NULL;
    
    for(i in 1:ncol(exp.db))
    {
      exp.db[,i] <- (exp.db[,i] - mean(exp.db[,i], na.rm = T))/sd(exp.db[,i], na.rm = T)
      exp.col <- union(exp.col, exp.db[,i])  
    }
    #exp.col <- as.vector(unique(exp.db))
    names(exp.col) = 1:length(exp.col)
    exp.col = data.frame("exp" = exp.col, "col" ="#000000" )
    exp.col = exp.col[order(exp.col[,"exp"]),]
    if(drug != "paclitaxel") {
      exp.col[,"col"] = colorRampPalette(c(blue, "white", red, red))(nrow(exp.col))
    }else{
      exp.col[,"col"] = colorRampPalette(c(blue, "white", red))(nrow(exp.col))
    }
    oo <- names(biomarkers.order)[sort(match(colnames(exp.db), names(biomarkers.order)))]
    exp.db <- exp.db[ , oo, drop=FALSE]
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_uhn.pdf", drug, sensitivity.type)), height=7, width=14)  
    #par(mfrow=c(2,1))
    #library("gplots")
    #gplots::heatmap.2(t(exp.db), Colv = NA, Rowv = T, col = exp.col[,"col"], scale = "row", trace = "none", dendrogram = "none", )
    colnames(exp.db) <- xx[colnames(exp.db)]
    if(ncol(exp.db) != 1){at.place=(0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1)}else{at.place=.5}
    hv <- NULL
    if(ncol(exp.db) == 1){
      par(mar=c(9, 2, 5, 15))
      par(oma=c(2,2,2,2))
      image(exp.db, col = exp.col[,"col"], axes = FALSE)
      grid(nx = nrow(exp.db), ny = ncol(exp.db), lty = 1)
      axis(4,at = 0, labels=colnames(exp.db), las=2, cex.axis = 2, tick = FALSE)
      #axis(4,at = (0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1), labels=sapply(biomarkers, function(x){ifelse(x[["gtex"]] == "tumor.specific", "*", NA)}, simplify = T), las=2, cex.axis = .8, tick = FALSE)
      
    }else{
      par(mar=c(2, 2, 2, 8))
      
      #par(mar=c(9, 2, 5, 15))
      #par(oma=c(2,2,2,2))
      image(exp.db, col = exp.col[,"col"], axes = FALSE)
      grid(nx = nrow(exp.db), ny = ncol(exp.db), lty = 1)
      axis(4, at=at.place, labels=colnames(exp.db), las=2, cex.axis = .8, tick = FALSE)
      # hv <- heatmap(t(exp.db), Colv = NA, Rowv = T, col = exp.col[,"col"], scale = "row", 
      #             labRow = colnames(exp.db),
      #             labCol = NA)
    }
    #image(exp.db, col = exp.col[,"col"], axes = FALSE)
    #grid(nx = nrow(exp.db), ny = ncol(exp.db), lty = 1)
    #axis(2,at = (0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1), labels=colnames(exp.db), las=2, cex.axis = .6, tick = FALSE)
    #axis(4,at = (0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1), labels=sapply(biomarkers, function(x){ifelse(x[["gtex"]] == "tumor.specific", "*", NA)}, simplify = T), las=2, cex.axis = .8, tick = FALSE)
    dev.off()
  
  
    my.xlim = c(1,length(sensitivity))
    my.ylim = range(sensitivity, na.rm=TRUE)
    
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_uhn_2.pdf", drug, sensitivity.type)), height=5, width=15)  
    par(mar=c(9,5,5,8))
    par(oma=c(2,2,2,2))
    plot(NA, xlim = my.xlim, ylim = my.ylim,ylab='',xlab='', axes = FALSE)
    axis(1,at = 1:length(sensitivity), labels=names(sensitivity), las = 2, cex.axis = 1.6, tck = -.05)
    axis(2,at = c(0,.1,.2,.3,.4,.5), labels=c(0,.1,.2,.3,.4,.5), cex.axis = 1, tck = -.02)
    
    box(lty = 1)
    points(1:length(sensitivity),sensitivity, pch = 20, cex = 1.5, col = "#663399")
    
    dev.off()
    
    pdf(file.path(path.diagrams,sprintf("%s_%s_uhn_estimates.pdf", drug, sensitivity.type)), height = 8, width = 2)
    xx <- as.numeric(uhn.models[oo, "estimate"])
    barplot(xx, horiz = T, col = sapply(xx , function(x){ifelse(x >= 0, mycol3[4], mycol3[3])}), yaxt='n')
    dev.off()
    
    pdf(file.path(path.diagrams,sprintf("%s_%s_uhn_pvalues.pdf", drug, sensitivity.type)), height = 8, width = 2)
    xx <- -log10(as.numeric(uhn.models[oo, "pvalue"])) * sign(as.numeric(uhn.models[oo, "estimate"]))
    
    barplot(xx, horiz = T, col = sapply(xx , function(x){ifelse(x >= 0, mycol3[4], mycol3[3])}), yaxt='n')
    dev.off()
  }else{
    message(sprintf("No final validated biomarker for drug: %s", drug))
  }
}
save(biomarkers, file=file.path(path.diagrams, "Biomarkers_uhn_status.RData"))
####################################################

## For each validated biomarkers this plot shows the exxpression of its gene and all the other isoforms in training, pre-validation and final-validation sets
results <- file.path(path.diagrams, "confidence_interval_difference.txt")
cutoff <- 0.1
for(drug in names(biomarkers)) {
  xx <- which(biomarkers[[drug]]$UHN.pvalue < cutoff)
  if (length(xx) > 0) {
    for(x in xx) {
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
    }
  }
}

#####################################################

## compare the expression of biomarkers in cell lines to tumour and normal samples
mycol <- RColorBrewer::brewer.pal(n=8, name="Set2")[c(6, 4, 1)]

xx <- intersect(rownames(ccle.isoforms.fpkm), rownames(CCLE@cell)[which(CCLE@cell$tissueid == "breast")])
ccle.isoforms.fpkm.breast <- ccle.isoforms.fpkm[xx, ]
load("data/TCGA_BRCA_isoforms.RData", verbose=T)
load("data/GTex_BR.RData", verbose=T)

normal.tcga <- rownames(pData(TCGA_BRCA_isoforms))[which(pData(TCGA_BRCA_isoforms)[, "sample_type"] == "NT")]
tumor.tcga <- rownames(pData(TCGA_BRCA_isoforms))[which(pData(TCGA_BRCA_isoforms)[, "sample_type"] == "TP")]

for(drug in names(biomarkers)) {
  xx <- which(biomarkers[[drug]]$UHN.pvalue < cutoff)
  if (length(xx) > 0) {
    for(x in xx) {
      ensembl.id <- biomarkers[[drug]][x, "transcript.id"]
      symbol <- biomarkers[[drug]][x, "symbol"]
      pdf(file = file.path(path.diagrams, sprintf("%s_%s_exp.pdf", drug, ensembl.id)), width=5.5, height=5.5)
      par(mar=c(7,6,2,2))
      #tt <- matrix(NA, ncol = 4, nrow = max((nrow(ccle.isoforms.fpkm.breast) + nrow(gray.isoforms.fpkm) + nrow(uhn.isoforms.fpkm)), ncol(exprs(GTex.BR$isoforms)), ncol(exprs(TCGA_BRCA_isoforms))))
      #colnames(tt) <- c("Cell lines","TCGA tumor", "TCGA normal", "GTex")
      tt <- matrix(NA, ncol = 3, nrow = max((nrow(ccle.isoforms.fpkm.breast) + nrow(gray.isoforms.fpkm) + nrow(uhn.isoforms.fpkm)), ncol(exprs(GTex.BR$isoforms)), ncol(exprs(TCGA_BRCA_isoforms))))
      colnames(tt) <- c("Cell lines","Tumour", "Healthy")
      tt[1:(nrow(ccle.isoforms.fpkm.breast) + nrow(gray.isoforms.fpkm) + nrow(uhn.isoforms.fpkm)),1] <- c(ccle.isoforms.fpkm.breast[ , ensembl.id], 
                                                                                                                       gray.isoforms.fpkm[ , ensembl.id],
                                                                                                                       uhn.isoforms.fpkm[ , ensembl.id])
      tt[1:length(tumor.tcga),2] <- exprs(TCGA_BRCA_isoforms)[ensembl.id, tumor.tcga]
      #tt[1:length(normal.tcga),3] <- exprs(TCGA_BRCA_isoforms)[ensembl.id, normal.tcga]
      
      #tt[1:ncol(exprs(GTex.BR$isoforms)),4] <- exprs(GTex.BR$isoforms)[ensembl.id, ]
      tt[1:ncol(exprs(GTex.BR$isoforms)),3] <- exprs(GTex.BR$isoforms)[ensembl.id, ]
      
      #boxplot(tt, col = mycol, main = sprintf("%s-%s (%s)", all.biomarkers[[i]][j,"symbol"], all.biomarkers[[i]][j,"type"], validation.stat), pch = 20)
      par(cex.lab=1.8)
      par(cex.axis=1.8)
      #boxplot(tt, pt.pch = 20, ylim=c(0,max(as.numeric(tt), na.rm=T)), cex.names=1.5, names=colnames(tt), las=2)
      invisible(genefu::boxplotplus2(t(tt[,2:3]), pt.pch = 20,
                                     .ylim=c(0,max(as.numeric(tt[,2:3]), na.rm=T)),
                                     pt.col=mycol[2:3],
                                     names=c("",""),#c(paste0("none (",length(res.nosignifseeds),")"), paste0("all (",length(res.allsignifseeds),")")),
                                     ylab="Expression",
                                     #xlab="signif seeds (# of targets)",
                                     #main=symbol,#sprintf("%s-%s (%s)", all.biomarkers[[i]][j,"symbol"], all.biomarkers[[i]][j,"type"], validation.stat),
                                     .las=2,
                                     pt.cex=1.5))
      text(1:2, par("usr")[3], labels=colnames(tt)[2:3], srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=1.5)
      
      #cells.tumor <- wilcox.test(tt[,1], tt[,2], alternative = "greater")$p.value
      if(!all(is.na(tt[,1])) & !all(is.na(tt[,3])))
      {
        #cells.normal <- t.test(tt[,1], tt[,3], alternative = "greater")$p.value
        cells.normal <- wilcox.test(tt[,1], tt[,3], alternative=ifelse(sign(biomarkers[[drug]][x, "estimate"]) > 0, "greater", "less"))$p.value
      }else{
        cells.normal <- 1
      }
#       if(!all(is.na(tt[,1])) & !all(is.na(tt[,4])))
#       {
#         #cells.gtex <- t.test(tt[,1], tt[,4], alternative = "greater")$p.value
#         cells.gtex <- wilcox.test(tt[,1], tt[,4], alternative = "greater")$p.value
#       }else{
#         cells.gtex <- 1
#       }
      
      
      if(!all(is.na(tt[,2])) & !all(is.na(tt[,3])))
      {
        #tumour.normal <- t.test(tt[,2], tt[,3], alternative = "greater")$p.value
        tumour.normal <- wilcox.test(tt[,2], tt[,3], alternative=ifelse(sign(biomarkers[[drug]][x, "estimate"]) > 0, "greater", "less"))$p.value
      }else{
        tumour.normal <- 1
      }
#       if(!all(is.na(tt[,2])) & !all(is.na(tt[,4])))
#       {
#         #tumour.gtex <- t.test(tt[,2], tt[,4], alternative = "greater")$p.value
#         tumour.gtex <- wilcox.test(tt[,2], tt[,4], alternative = "greater")$p.value
#       }else{
#         tumour.gtex <- 1
#       }
#       if(!all(is.na(tt[,1])) & !all(is.na(tt[,3])) & !all(is.na(tt[,4])))
#       {
#         cells <- survcomp::combine.test(c(cells.normal, cells.gtex))
#       }else{
#         cells <- 1
#       }
#       if(!all(is.na(tt[,2])) & !all(is.na(tt[,3])) & !all(is.na(tt[,4])))
#       {
#         tumours <- survcomp::combine.test(c(tumour.normal, tumour.gtex))
#       }else{
#         tumours <- 1
#       }
#       #test <- wilcox.test(tt[,2], tt[,3], alternative = "greater")$p.value
      cells <- cells.normal
      tumours <- tumour.normal
#      legend("topright", legend = c(sprintf("Cell Lines %s Healthy : %.1E", ifelse(sign(biomarkers[[drug]][x, "effect.size"]) > 0, ">", "<"), cells),
#                                    sprintf("Tumour %s Healthy : %.1E", ifelse(sign(biomarkers[[drug]][x, "effect.size"]) > 0, ">", "<"), tumours)), fill = mycol[1:2] , bty="n", cex=1.5)
      legend("topright", legend=sprintf("Tumour %s Healthy : %.1E", ifelse(sign(biomarkers[[drug]][x, "estimate"]) > 0, ">", "<"), tumours), bty="n", cex=1.5)
      #round(cells, digits = 25)
      dev.off()
    }
  }
}
#####################################################
