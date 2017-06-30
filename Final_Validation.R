require(PharmacoGx) || stop("Library PharmacoGx is not available!")
require(Biobase) || stop("Library Biobase is not available!")
require(gdata) || stop("Library gdata is not available!")
require(genefu) || stop("Library genefu is not available!")
require(survcomp) || stop("Library survcomp is not available!")

source("code/foo_FinalValidation.R")

load(file.path(path.data, "PSets/UHN_hs.RData"))
ccle.drug.sensitivity[which(is.nan(ccle.drug.sensitivity))] <- NA
gdsc.drug.sensitivity[which(is.nan(gdsc.drug.sensitivity))] <- NA
gray.drug.sensitivity[which(is.nan(gray.drug.sensitivity))] <- NA
drugs <- intersectList(colnames(ccle.drug.sensitivity), 
                       colnames(gdsc.drug.sensitivity),
                       colnames(gray.drug.sensitivity),
                       drugNames(UHN))
uhn.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=UHN, drugs=drugs, sensitivity.measure=sensitivity.type))
uhn.drug.sensitivity[which(is.nan(uhn.drug.sensitivity))] <- NA
genes <- colnames(ccle.genes.fpkm)
isoforms <- colnames(ccle.isoforms.fpkm)

uhn.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=UHN, mDataType="rnaseq", features=genes, fill.missing=FALSE)))
uhn.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=UHN, mDataType="isoforms", features=isoforms, fill.missing=FALSE)))
uhn.isoforms.fpkm[which(is.na(uhn.isoforms.fpkm))] <- 0
if(RNA_seq.normalize == TRUE)
{
  uhn.genes.fpkm <- log2(uhn.genes.fpkm + 1)
  uhn.isoforms.fpkm <- log2(uhn.isoforms.fpkm + 1)
}
uhn.cells <- intersect(rownames(uhn.genes.fpkm), rownames(uhn.drug.sensitivity))
uhn.drug.sensitivity <- uhn.drug.sensitivity[uhn.cells, ]
uhn.genes.fpkm <- uhn.genes.fpkm[uhn.cells, ]
uhn.isoforms.fpkm <- uhn.isoforms.fpkm[uhn.cells, ]
mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
mycol3 <- RColorBrewer::brewer.pal(n=4, name="Set3")

## Heatmaps, effect sizes and sensitivity plots for final validation, just isoforms considered
cutoff <- 0.1
red <- mycol[1]  
blue <- mycol[2]

load(file.path(path.diagrams, "validated.biomarkers.gray.RData"), verbose=TRUE)
for(drug in drugs) {
  validated.biomarkers[[drug]][, c("UHN.estimate", "UHN.pvalue", paste0("UHN.",effect.size), "gene.biotype")] <- NA
  validated.biomarkers[[drug]][,"id"] <- validated.biomarkers[[drug]][,"biomarker.id"]
  vtt <- validated.biomarkers[[drug]]
  vtt <- vtt[which(vtt[,"type"] == "isoform"), , drop=F]
  vtt <- vtt[which(vtt[, "isoforms.no"]>1), , drop=FALSE]
  gray.specificity <- vtt[ ,"gray.specificity"]
  xx <- fnFetchBiomarkers(top.significant.biomarkers=vtt, drug=drug, indices=1:nrow(vtt))
  xx <- do.call(rbind, xx)
  xx[,"short.label"] <- gsub(".ISO$","",xx[,"short.label"])
  xx <- apply(xx, 1, function(x){x})
  
  ###Figure 4
  ###heatmap of all pre validated biomarkers in GRAY for the drugs in common with UHN
  ###row labels are colored according to their specificity in GRAY
  rr <- fnPlotAUCoverCellLinesGray(drug=drug, tissue.type="all", biomarkers=xx, suffix="all.specificity", gray.specificity=gray.specificity)#, biomarkers.toPlot)
  ###
  vtt <- vtt[which(vtt[,"gray.specificity"] != "gene.specific"), , drop=F]
  gray.specificity <- vtt[ ,"gray.specificity"]
  rr <- fnPlotAUCoverCellLinesGray(drug=drug, tissue.type="all", biomarkers=xx, suffix="isorom.specific", gray.specificity=gray.specificity)#, biomarkers.toPlot)
  
  if(all(!is.na(rr))){
    biomarkers.order <- rr$hv$rowInd
    fnPlotEffectSize(drug, biomarkers=xx, effect.size=effect.size, biomarkers.order) 
    xx <- do.call(rbind, xx)
    names(biomarkers.order) <- xx[biomarkers.order, "isoform.id"]
  }
  
  exp.db <- uhn.isoforms.fpkm[ , vtt[, "transcript.id"], drop=FALSE]
  
  coding.biotypes <- c("PRT", "AS", "prcTR", "pseudogene")
  names(coding.biotypes) <- c("protein_coding", "antisense", "processed_transcript", "processed_pseudogene")
  bb <- vtt[match(colnames(exp.db), vtt$biomarker.id), "biotype"]
  bb <- sapply(bb, function(x){ifelse(x %in% names(coding.biotypes), coding.biotypes[x], x)})
  xx <- sprintf("%s (%s)", vtt[match(colnames(exp.db), vtt$biomarker.id), "symbol"], bb)
  names(xx) <- colnames(exp.db)
  sensitivity <- uhn.drug.sensitivity[!is.na(uhn.drug.sensitivity[ , drug]) , drug]
  exp.db <- exp.db[names(sensitivity), , drop=FALSE]
  uhn.models <- matrix(NA, ncol=4, nrow=ncol(exp.db), dimnames=list(colnames(exp.db), c("pvalue", "estimate", effect.size, "gene_biotype")))
  for(marker in colnames(exp.db)) {
    uhn.model <- lm(sensitivity  ~ exp.db[, marker])
    uhn.pvalue <- 2 
    uhn.estimate <- uhn.effect.size <- 0
    if(all(!is.na(uhn.model)) & !is.na(uhn.model$coefficients[2])) {
      uhn.models[marker,"pvalue"] <- summary(uhn.model)$coefficients[2,4]
      uhn.models[marker,"estimate"] <- summary(uhn.model)$coefficients[2,1]
      if(effect.size == "r.squared") {
        uhn.models[marker, effect.size] <- summary(uhn.model)$r.squared
      }
      if(effect.size == "cindex"){
        uhn.models[marker, effect.size] <- Hmisc::rcorr.cens(x=predict(uhn.model), S=sensitivity, outx=TRUE)[[1]]
      }
      uhn.models[marker,"gene_biotype"] <- annot.ensembl.all.isoforms[marker, "TranscriptBioType"]
    }
  }
  if(!is.null(rownames(uhn.models))) {
    validated.biomarkers[[drug]][match(rownames(uhn.models), validated.biomarkers[[drug]][,"id"]), c("UHN.estimate", "UHN.pvalue", paste0("UHN.",effect.size), "gene.biotype")] <- uhn.models[, c("estimate", "pvalue", effect.size, "gene_biotype")]
    uhn.models <- uhn.models[which(as.numeric(uhn.models[,"pvalue"]) < cutoff & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
    #uhn.models <- uhn.models[which(as.numeric(uhn.models[ ,effect.size]) > 0.55 & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
  } 
  message(drug)
  message(nrow(vtt))
  message(nrow(uhn.models))
}
#rr <- list();for(drug in drugs){rr[[drug]]<-biomarkers[[drug]][which(biomarkers[[drug]]$UHN.cindex > effect.size.cut.off),]}
rr <- list();
for(drug in drugs){
  rr[[drug]] <- validated.biomarkers[[drug]][which(validated.biomarkers[[drug]]$UHN.pvalue < cutoff & sign(as.numeric(validated.biomarkers[[drug]]$UHN.estimate)) == sign(as.numeric(as.numeric(validated.biomarkers[[drug]]$estimate)))),]
}
final.validated.biomarkers <- lapply(rr, function(x){if("UHN.cindex" %in% colnames(x)){x[order(x[, "UHN.cindex"], na.last=T, decreasing=T),]}else{x}})
save(final.validated.biomarkers, file=file.path(path.diagrams, "validated.biomarkers.uhn.RData"))

isoforms <- list()
for(drug in drugs) {
  #ii <- which.max(as.numeric(final.validated.biomarkers[[drug]]$UHN.cindex))
  ii <- 1
  symbol <- final.validated.biomarkers[[drug]][ii, "symbol"]
  gene <- final.validated.biomarkers[[drug]][ii, "gene.id"]
  best.isoform <- final.validated.biomarkers[[drug]][ii, "transcript.id"]
  message(sprintf("%s:  %s(%s)  %s  (%s) %s", drug, symbol, gene, final.validated.biomarkers[[drug]][ii, "transcript.id"], final.validated.biomarkers[[drug]][ii, "UHN.cindex"], final.validated.biomarkers[[drug]][ii, "biotype"]))
  #View(final.validated.biomarkers[[drug]][ii,])
  
  isoforms[[drug]] <- intersect(annot.ensembl.all.isoforms[which(annot.ensembl.all.isoforms[,"EnsemblGeneId"] == gene),"EnsemblTranscriptId"],
                                colnames(ccle.isoforms.fpkm))
  xx <- apply(gray.isoforms.fpkm[ , isoforms[[drug]], drop=FALSE], MARGIN=2, function(x){length(which(x!=0))} )
  isoforms[[drug]] <- names(xx)[which(xx!=0)]
  ##set igv orders
  if(drug=="AZD6244"){
    isoforms[[drug]] <- isoforms[[drug]][rev(c(3, 5, 2, 4, 7, 8, 9, 6, 1, 10))]
  }
  if(drug=="lapatinib"){
  #  isoforms[[drug]] <- isoforms[[drug]][rev(c(2, 4, 1, 5, 3))]
    isoforms[[drug]] <- isoforms[[drug]][rev(c(2, 1))]
  }
  if(drug=="Erlotinib"){
  #  isoforms[[drug]] <- isoforms[[drug]][rev(c(1, 5, 4, 6, 8, 2, 3, 7))]
    isoforms[[drug]] <- isoforms[[drug]][rev(c(2, 6, 1, 5, 3, 4, 7, 9, 10, 8))]
  }
  if(drug=="paclitaxel"){
    isoforms[[drug]] <- isoforms[[drug]][rev(c(1, 3, 7, 4, 5, 2, 6, 8))]
  }
  
  message(paste(isoforms[[drug]], collapse="  "))
  expression <- cbind(gray.isoforms.fpkm[ , isoforms[[drug]], drop=FALSE], "gene"=gray.genes.fpkm[, gene])
  xx <- cor(expression, expression, use="pairwise", method="spearman")
  ###Figure 5
  ###Correlation of the biomarker to the other alternatively spliced products of the corresponding gene
  fnCor(drug, gene=symbol, xx, isoforms=isoforms[[drug]], best.isoform=best.isoform)
  ###Expression of the biomarker along with all the other alternatively spliced products of the corresponding gene
  fnExp(drug, gene=symbol, exp=gray.isoforms.fpkm[ , isoforms[[drug]], drop=FALSE], best.isoform=best.isoform)
  exprs <- cbind(best.isoform=uhn.isoforms.fpkm[ , best.isoform, drop=FALSE], "gene"=uhn.genes.fpkm[, gene])
  ###Heatmap of the expression of the biomarker in UHN cell lines
  fnPlotHeatMap(sensitivity=uhn.drug.sensitivity[ , drug], file.name=sprintf("%s_%s", drug, symbol), cluster=FALSE, expression=exprs, best.isoform=best.isoform)
  ###Sensitivity of UHN cell lines to drug
  fnPlotSensitivity(sensitivity=uhn.drug.sensitivity[ , drug], file.name=sprintf("%s_%s", drug, symbol))
  ####
  
  sensitivity <- uhn.drug.sensitivity[!is.na(uhn.drug.sensitivity[ , drug]) , drug]
  uhn.model <- lm(sensitivity  ~ exprs[names(sensitivity), best.isoform])
  ###Supplementary 14
  ###Predicted AUC values aginst the actual AUC values of the drug
  myScatterPlot(Name=file.path("result/scatter", sprintf("%s_%s_aac.pdf", drug, symbol)), 
                x=predict(uhn.model), 
                y=sensitivity, 
                legend.label=sprintf("cindex=%s",round(as.numeric(final.validated.biomarkers[[drug]][ii, "UHN.cindex"]), digits=2)), 
                method="plain", 
                minp=10,
                xlim=c(0, 1), 
                ylim=c(0, 1),
                xlab="predicted",
                ylab="AAC",
                main=best.isoform)
  ###
  expression <- cbind(uhn.isoforms.fpkm[ , isoforms[[drug]], drop=FALSE], "gene"=uhn.genes.fpkm[, gene])
  ###Supplementary 15
  ###Expression of the biomarker along with all the other alternatively spliced products of the corresponding gene in UHN cells
  fnPlotHeatMap(sensitivity=sensitivity, file.name=sprintf("%s_%s", drug, symbol), cluster=FALSE, expression=expression, best.isoform=best.isoform)
  ###
  colnames(expression)[ncol(expression)] <- annot.gene$EnsemblGeneId

  
  #colnames(expression)[ncol(expression)] <- gene
}
