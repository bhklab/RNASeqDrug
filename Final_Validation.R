require(PharmacoGx) || stop("Library PharmacoGx is not available!")
require(Biobase) || stop("Library Biobase is not available!")
require(gdata) || stop("Library gdata is not available!")
require(genefu) || stop("Library genefu is not available!")
require(survcomp) || stop("Library survcomp is not available!")

source("code/foo_FinalValidation.R")
validation.method <- effect.size

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
load(file.path(path.diagrams, "Biomarkers_validated_breast_cindex_gray_pvalue.RData"), verbose=TRUE)
for(drug in drugs) {
  biomarkers[[drug]][, c("UHN.estimate", "UHN.pvalue", paste0("UHN.",effect.size), "gene.biotype")] <- NA
  biomarkers[[drug]][,"id"] <- biomarkers[[drug]][,"biomarker.id"]
  vtt <- biomarkers[[drug]]
  vtt <- vtt[which(vtt[,"type"] == "isoform"), , drop=F]
  vtt <- vtt[which(vtt[, "isoforms.no"]>1), , drop=FALSE]
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
  
  exp.db <- NULL
  if(length(which(vtt$type == "isoform")) > 0) {
    exp.db <- uhn.isoforms.fpkm[ , vtt[which(vtt$type == "isoform"), "transcript.id"], drop=FALSE]
  }
  if(!is.null(exp.db)) {
    exp.db <- cbind(exp.db, uhn.genes.fpkm[ , vtt[which(vtt$type == "gene"), "gene.id"], drop=FALSE])
  } else {
    exp.db <- uhn.genes.fpkm[ , vtt[which(vtt$type == "gene"), "gene.id"], drop=FALSE]
  }
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
    biomarkers[[drug]][match(rownames(uhn.models), biomarkers[[drug]][,"id"]), c("UHN.estimate", "UHN.pvalue", paste0("UHN.",effect.size), "gene.biotype")] <- uhn.models[, c("estimate", "pvalue", effect.size, "gene_biotype")]
    uhn.models <- uhn.models[which(as.numeric(uhn.models[,"pvalue"]) < cutoff & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
    #uhn.models <- uhn.models[which(as.numeric(uhn.models[ , effect.size]) > effect.size.cut.off & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
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
rr <- list();for(drug in drugs){rr[[drug]]<-biomarkers[[drug]][which(biomarkers[[drug]]$UHN.pvalue < cutoff & sign(as.numeric(biomarkers[[drug]]$UHN.estimate)) == sign(as.numeric(as.numeric(biomarkers[[drug]]$estimate)))),]}
biomarkers <- rr
save(biomarkers, file=file.path(path.diagrams, "Biomarkers_uhn_status.RData"))

load(file.path(path.diagrams, "Biomarkers_validated_breast_cindex_gray_pvalue.RData"), verbose=TRUE)
path.diagrams2 <- file.path(path.diagrams, "isoform_specific")
if(!file.exists(path.diagrams2)){dir.create(path.diagrams2)}
for(drug in drugs) {
  biomarkers[[drug]][, c("UHN.estimate", "UHN.pvalue", paste0("UHN.",effect.size), "gene.biotype")] <- NA
  biomarkers[[drug]][,"id"] <- biomarkers[[drug]][,"biomarker.id"]
  vtt <- biomarkers[[drug]]
  vtt <- vtt[which(vtt[,"type"] == "isoform"), , drop=F]
  vtt <- vtt[which(vtt[, "isoforms.no"]>1), , drop=FALSE]
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
  
  exp.db <- NULL
  if(length(which(vtt$type == "isoform")) > 0) {
    exp.db <- uhn.isoforms.fpkm[ , vtt[which(vtt$type == "isoform"), "transcript.id"], drop=FALSE]
  }
  if(!is.null(exp.db)) {
    exp.db <- cbind(exp.db, uhn.genes.fpkm[ , vtt[which(vtt$type == "gene"), "gene.id"], drop=FALSE])
  } else {
    exp.db <- uhn.genes.fpkm[ , vtt[which(vtt$type == "gene"), "gene.id"], drop=FALSE]
  }
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
    biomarkers[[drug]][match(rownames(uhn.models), biomarkers[[drug]][,"id"]), c("UHN.estimate", "UHN.pvalue", paste0("UHN.",effect.size), "gene.biotype")] <- uhn.models[, c("estimate", "pvalue", effect.size, "gene_biotype")]
    uhn.models <- uhn.models[which(as.numeric(uhn.models[,"pvalue"]) < cutoff & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
    #uhn.models <- uhn.models[which(as.numeric(uhn.models[ ,effect.size]) > 0.55 & sign(as.numeric(uhn.models[,"estimate"])) == sign(as.numeric(vtt[match(rownames(uhn.models), vtt$id), "estimate"]))), , drop=FALSE]
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
    pdf(file = file.path(path.diagrams2, sprintf("%s_%s_uhn.pdf", drug, sensitivity.type)), height=7, width=14)  
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
    
    pdf(file = file.path(path.diagrams2, sprintf("%s_%s_uhn_2.pdf", drug, sensitivity.type)), height=5, width=15)  
    par(mar=c(9,5,5,8))
    par(oma=c(2,2,2,2))
    plot(NA, xlim = my.xlim, ylim = my.ylim,ylab='',xlab='', axes = FALSE)
    axis(1,at = 1:length(sensitivity), labels=names(sensitivity), las = 2, cex.axis = 1.6, tck = -.05)
    axis(2,at = c(0,.1,.2,.3,.4,.5), labels=c(0,.1,.2,.3,.4,.5), cex.axis = 1, tck = -.02)
    
    box(lty = 1)
    points(1:length(sensitivity),sensitivity, pch = 20, cex = 1.5, col = "#663399")
    
    dev.off()
    
    pdf(file.path(path.diagrams2,sprintf("%s_%s_uhn_estimates.pdf", drug, sensitivity.type)), height = 8, width = 2)
    xx <- as.numeric(uhn.models[oo, "estimate"])
    barplot(xx, horiz = T, col = sapply(xx , function(x){ifelse(x >= 0, mycol3[4], mycol3[3])}), yaxt='n')
    dev.off()
    
    pdf(file.path(path.diagrams2,sprintf("%s_%s_uhn_pvalues.pdf", drug, sensitivity.type)), height = 8, width = 2)
    xx <- -log10(as.numeric(uhn.models[oo, "pvalue"])) * sign(as.numeric(uhn.models[oo, "estimate"]))
    
    barplot(xx, horiz = T, col = sapply(xx , function(x){ifelse(x >= 0, mycol3[4], mycol3[3])}), yaxt='n')
    dev.off()
  }else{
    message(sprintf("No final validated biomarker for drug: %s", drug))
  }
}
#rr <- list();for(drug in drugs){rr[[drug]]<-biomarkers[[drug]][which(biomarkers[[drug]]$UHN.cindex > effect.size.cut.off),]}
rr <- list();for(drug in drugs){rr[[drug]]<-biomarkers[[drug]][which(biomarkers[[drug]]$UHN.pvalue < cutoff & sign(as.numeric(biomarkers[[drug]]$UHN.estimate)) == sign(as.numeric(as.numeric(biomarkers[[drug]]$estimate)))),]}
biomarkers <- rr
save(biomarkers, file=file.path(path.diagrams2, "Biomarkers_uhn_status.RData"))

####################################################
isoforms <- list()
for(drug in drugs) {
  #xx <- which.min(as.numeric(biomarkers[[drug]]$UHN.pvalue))
  ii <- which.max(as.numeric(biomarkers[[drug]]$UHN.cindex))
  #if(drug=="paclitaxel"){
  #  ii <- order(as.numeric(biomarkers[[drug]]$UHN.cindex), decreasing=T)[2]
  #}
  gene <- biomarkers[[drug]][ii,"gene.id"]
  best.isoform <- biomarkers[[drug]][ii, "transcript.id"]
  message(sprintf("%s:  %s(%s)  %s  (%s) %s", drug, biomarkers[[drug]][ii, "symbol"], gene, biomarkers[[drug]][ii, "transcript.id"], biomarkers[[drug]][ii, "UHN.cindex"], biomarkers[[drug]][ii, "biotype"]))
  #View(biomarkers[[drug]][ii,])
  
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
  fnCor(drug, gene, xx, isoforms=isoforms[[drug]], best.isoform=best.isoform)
  fnExp(drug, gene, exp=gray.isoforms.fpkm[ , isoforms[[drug]], drop=FALSE], best.isoform=best.isoform)
  exprs <- cbind(best.isoform=uhn.isoforms.fpkm[ , best.isoform], "gene"=uhn.genes.fpkm[, gene])
  fnPlotHeatMap(sensitivity=uhn.drug.sensitivity[ , drug], file.name=sprintf("%s_%s", drug, gene), cluster=FALSE, expression=exprs, best.isoform=best.isoform)
  fnPlotSensitivity(sensitivity=uhn.drug.sensitivity[ , drug], file.name=sprintf("%s_%s_sensitivity", drug, gene))
  
  
  #colnames(expression)[ncol(expression)] <- gene
}
#########
## For each validated biomarkers this plot shows the expression of its gene and all the other isoforms in training, pre-validation and final-validation sets
results <- file.path(path.diagrams, "confidence_interval_difference.txt")
cutoff <- 0.1
for(drug in names(biomarkers)) {
  xx <- which(biomarkers[[drug]]$UHN.pvalue < cutoff)
  if (length(xx) > 0) {
    for(x in xx) {
      gene <- biomarkers[[drug]][x, "symbol"] 
      gene.id <- biomarkers[[drug]][x, "gene.id"]
      best.isoform <- biomarkers[[drug]][x, "transcript.id"]
      annot.isoforms <- annot.ensembl.all.isoforms[which(annot.ensembl.all.isoforms[,"EnsemblGeneId"] == gene.id),]
      annot.isoforms <- intersectList(rownames(annot.isoforms), colnames(ccle.isoforms.fpkm), colnames(gray.isoforms.fpkm), colnames(uhn.isoforms.fpkm))
      annot.gene <- annot.ensembl.all.genes[which(annot.ensembl.all.genes[,"EnsemblGeneId"] == gene.id), ]
      
      cat(sprintf("%s, %s, %s, cutoff=%s\n", gene, best.isoform, drug, cutoff), file=results, append=TRUE)

      ## training 
      ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.isoforms.fpkm))
      ccle.breast.cells <- intersectList(rownames(ccle.cell.profiles)[which(ccle.cell.profiles$tissueid == "breast")], rownames(ccle.drug.sensitivity), rownames(ccle.isoforms.fpkm))
      ccle.breast.cells <- ccle.breast.cells[which(!is.na(ccle.drug.sensitivity[ccle.breast.cells, drug]))]
      ccle.auc <- ccle.drug.sensitivity[ccle.breast.cells, drug]
      ccle.isoforms.models <- fnBuildLinearModel(cells=ccle.cells, sensitivity=ccle.drug.sensitivity[ccle.cells, drug], features=annot.isoforms, fpkm.matrix=ccle.isoforms.fpkm)
      
      gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.isoforms.fpkm))
      gdsc.breast.cells <- intersectList(rownames(ccle.cell.profiles)[which(ccle.cell.profiles$tissueid == "breast")],  rownames(gdsc.drug.sensitivity), rownames(ccle.isoforms.fpkm))
      gdsc.breast.cells <- gdsc.breast.cells[which(!is.na(gdsc.drug.sensitivity[gdsc.breast.cells, drug]))]
      gdsc.auc <- gdsc.drug.sensitivity[gdsc.breast.cells, drug]
      gdsc.isoforms.models <- fnBuildLinearModel(cells=gdsc.cells, sensitivity=gdsc.drug.sensitivity[gdsc.cells, drug], features=annot.isoforms, fpkm.matrix=ccle.isoforms.fpkm)
      
      training.isoforms.models <- fnNormalizeTraining(ccle.isoforms.models$simplified, gdsc.isoforms.models$simplified)
      #training.isoforms.models <- ccle.isoforms.models$simplified
      
      ccle.gene.model <- fnBuildLinearModel(cells=ccle.breast.cells, sensitivity=ccle.auc, features=annot.gene$EnsemblGeneId, fpkm.matrix=ccle.genes.fpkm)
      gdsc.gene.model <- fnBuildLinearModel(cells=gdsc.breast.cells, sensitivity=gdsc.auc, features=annot.gene$EnsemblGeneId, fpkm.matrix=ccle.genes.fpkm)
      
      training.gene.model <- fnNormalizeTraining(ccle.gene.model$simplified, gdsc.gene.model$simplified)
      #training.gene.model <- ccle.gene.model$simplified
      #training.models.isoforms <- training.models.isoforms[which(training.models.isoforms[,"pvalue"] < cutoff), , drop=FALSE]
      
      expression <- cbind(ccle.isoforms.fpkm[ , annot.isoforms], "gene"=ccle.genes.fpkm[, annot.gene$EnsemblGeneId])
      colnames(expression)[ncol(expression)] <- annot.gene$EnsemblGeneId
      
      #isoforms.ordered <- fnPlotHeatMap(sensitivity=ccle.auc, expression=expression, file.name=sprintf("%s_%s", drug, gene), cluster=TRUE, best.isoform=best.isoform)
      #isoforms.ordered <- fnPlotHeatMap.union.cells(ccle.sensitivity=ccle.auc, gdsc.sensitivity=gdsc.auc, expression=expression, file.name=sprintf("%s_%s", drug, gene), cluster=TRUE, best.isoform=best.isoform)
      isoforms.ordered <- NULL
      if(is.null(isoforms.ordered)) {
        isoforms.ordered <- annot.isoforms
      }
      #fnPlotSensitivity(sensitivity=ccle.auc, file.name=sprintf("%s_%s", drug, gene), type="training")
      #fnPlotEstimates(isoform.models=training.isoforms.models, gene.model=training.gene.model, isoforms=isoforms.ordered, file.name=sprintf("%s_%s", drug, gene), cutoff=cutoff)
      #fnPlotPvalues(isoform.models=training.isoforms.models, gene.model=training.gene.model, isoforms=isoforms.ordered, file.name=sprintf("%s_%s", drug, gene), cutoff=cutoff)
      
      cat(sprintf("training confidence interval difference pvalue: %s\n", fnConfidenceIntervalDifference(model1=training.isoforms.models[best.isoform,], 
                                                                                                         model2=training.gene.model, 
                                                                                                         x1=ccle.isoforms.models$complete[[best.isoform]]$model$`tt[, "exp"]`,
                                                                                                         x2=ccle.gene.model$complete[[1]]$model$`tt[, "exp"]`)), file=results, append=TRUE)
      
      ## in silico
      
      gray.cells <- intersectList(rownames(gray.drug.sensitivity), rownames(gray.isoforms.fpkm))
      gray.auc <- gray.drug.sensitivity[gray.cells, drug]

      gray.isoforms.models <- fnBuildLinearModel(cells=rownames(gray.isoforms.fpkm), sensitivity=gray.auc, features=isoforms.ordered, fpkm.matrix=gray.isoforms.fpkm)
      gray.gene.model <- fnBuildLinearModel(cells=rownames(gray.genes.fpkm), sensitivity=gray.auc, features=as.character(annot.gene$EnsemblGeneId), fpkm.matrix=gray.genes.fpkm)
      
      #in.silico.models.isoforms <- in.silico.models.isoforms[which(in.silico.models.isoforms[,"pvalue"] < cutoff & sign(in.silico.models.isoforms[,"estimate"]) == sign(training.models.isoforms[,"estimate"])), , drop=FALSE]
      #expression <- cbind(gray.isoforms.fpkm[ , isoforms.ordered], "gene"=gray.genes.fpkm[, as.character(annot.gene$EnsemblGeneId)])
      #xx <- cor(expression[gray.cells,], expression[gray.cells,], use="pairwise", method="spearman")
      #fnCor(drug, gene, xx)
      #fnExp(drug, gene, gray.isoforms.fpkm[ , isoforms.ordered, drop=FALSE])
      #colnames(expression)[ncol(expression)] <- annot.gene$EnsemblGeneId

      #xx <- fnPlotHeatMap(sensitivity=gray.auc, expression=expression, file.name=sprintf("%s_%s_in_silico", drug, gene), best.isoform=best.isoform)
      #fnPlotSensitivity(sensitivity=gray.auc, file.name=sprintf("%s_%s_in_silico", drug, gene))
      #fnPlotEstimates(isoform.models=gray.isoforms.models$simplified, gene.model=gray.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_silico", drug, gene), cutoff=cutoff)
      #fnPlotPvalues(isoform.models=gray.isoforms.models$simplified, gene.model=gray.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_silico", drug, gene), cutoff=cutoff)
      
      cat(sprintf("in silico confidence interval difference pvalue: %s\n", fnConfidenceIntervalDifference(model1=gray.isoforms.models$simplified[best.isoform,], 
                                                                                                          model2=gray.gene.model$simplified, 
                                                                                                          x1=gray.isoforms.models$complete[[best.isoform]]$model$`tt[, "exp"]`,
                                                                                                          x2=gray.gene.model$complete[[1]]$model$`tt[, "exp"]`)), file=results, append=TRUE)
      
      ## in vitro
      
      uhn.cells <- intersectList(rownames(uhn.drug.sensitivity), rownames(uhn.isoforms.fpkm))
      uhn.auc <- uhn.drug.sensitivity[ , drug]
      
      uhn.isoforms.models <- fnBuildLinearModel(cells=rownames(uhn.isoforms.fpkm), sensitivity=uhn.auc, features=isoforms.ordered, fpkm.matrix=uhn.isoforms.fpkm)
      uhn.gene.model <- fnBuildLinearModel(cells=rownames(uhn.genes.fpkm), sensitivity=uhn.auc, features=as.character(annot.gene$EnsemblGeneId), fpkm.matrix=uhn.genes.fpkm)
      
      #in.vitro.models.isoforms <- in.vitro.models.isoforms[which(in.vitro.models.isoforms[,"pvalue"] < cutoff & sign(in.vitro.models.isoforms[,"estimate"]) == sign(in.silico.models.isoforms[,"estimate"])), , drop=FALSE]
      #expression <- cbind(uhn.isoforms.fpkm[ , isoforms.ordered, drop=FALSE], "gene"=uhn.genes.fpkm[, as.character(annot.gene$EnsemblGeneId)])
      #xx <- expression[, c(best.isoform,"gene"), drop=FALSE]
      #fnPlotHeatMap(sensitivity=uhn.auc, file.name=sprintf("%s_%s_image.pdf", gene, drug), cluster=FALSE, expression=xx, best.isoform=best.isoform)
      #colnames(expression)[ncol(expression)] <- annot.gene$EnsemblGeneId
      
      #xx <- fnPlotHeatMap(sensitivity=uhn.auc, expression=expression, file.name=sprintf("%s_%s_in_vitro", drug, gene), best.isoform=best.isoform)
      #fnPlotSensitivity(sensitivity=uhn.auc, file.name=sprintf("%s_%s_in_vitro", drug, gene))
      #fnPlotEstimates(isoform.models=uhn.isoforms.models$simplified, gene.model=uhn.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_vitro", drug, gene), cutoff=cutoff)
      #fnPlotPvalues(isoform.models=uhn.isoforms.models$simplified, gene.model=uhn.gene.model$simplified, isoforms=isoforms.ordered, file.name=sprintf("%s_%s_in_vitro", drug, gene), cutoff=cutoff)
      
      cat(sprintf("in vitro confidence interval difference pvalue: %s\n", fnConfidenceIntervalDifference(model1=uhn.isoforms.models$simplified[best.isoform,], 
                                                                                                         model2=uhn.gene.model$simplified, 
                                                                                                         x1=uhn.isoforms.models$complete[[best.isoform]]$model$`tt[, "exp"]`,
                                                                                                         x2=uhn.gene.model$complete[[1]]$model$`tt[, "exp"]`)), file=results, append=TRUE)
    }
  }
}

