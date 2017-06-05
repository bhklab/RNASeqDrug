require(gplots) || stop("Library gplots is not available!")
fnValidation <- function(top.significant.biomarkers, validation.cut.off, validation.method=c(effect.size, "pvalue"), effect.size.cut.off=0.55, pvalue.cut.off=0.05, plot.create=FALSE, common.cells=FALSE) {
  rr <- list()
  for(drug in drugs) {
    if(!all(is.na(top.significant.biomarkers[[drug]]))) {
      top.significant.biomarkers.drug <- top.significant.biomarkers[[drug]]
      if(validation.method==effect.size) {
        #  top.significant.biomarkers.drug <- top.significant.biomarkers.drug[which(top.significant.biomarkers.drug$effect.size <= top.significant.biomarkers.drug$breast), , drop=FALSE]
        top.significant.biomarkers.drug <- top.significant.biomarkers.drug[which(top.significant.biomarkers.drug$breast >= effect.size.cut.off), , drop=FALSE]
      } else {
        top.significant.biomarkers.drug <- top.significant.biomarkers.drug[which(top.significant.biomarkers.drug[ ,paste0(tissue, "_pvalue")] < pvalue.cut.off), , drop=FALSE]
      }
      ##The validation can be done on the most significant biomarkers which the number of significant biomarkers is determined by validation.cut.off
      top.significant.biomarkers.drug.5 <- NULL
      specificities <- table(top.significant.biomarkers.drug$specificity)
      for(spec in names(specificities)) {
        top.significant.biomarkers.drug.5 <- rbind(top.significant.biomarkers.drug.5, subset(top.significant.biomarkers.drug, specificity == spec)[1:min(validation.cut.off, specificities[spec], na.rm=TRUE),])
      }
      top.significant.biomarkers.drug <- top.significant.biomarkers.drug.5[order(top.significant.biomarkers.drug.5$rank),]
      biomarkers.no <- nrow(top.significant.biomarkers.drug)
      rr[[drug]] <- top.significant.biomarkers.drug
      
      if(plot.create) {
        pdf(file=file.path(path.diagrams, sprintf("%s.pdf", drug)), height=4 * (plots.no + 1), width=16)
        par(mfrow=c((plots.no + 1), 4))
        
        if(training.type == "CCLE_GDSC") {
          fnPlotAUC(first="ccle", second="gdsc", drug, tissue.type="all", color="#666666")
          fnPlotAUC(first="ccle", second="gdsc", drug, tissue.type=tissue, color="#996666")
          fnPlotAUC(first="ccle", second="gray", drug, tissue.type="", color="#CC3399")
          fnPlotAUC(first="gray", second="gdsc", drug, tissue.type="", color="#663399")
        }else {
          fnPlotAUC(first="ccle", second="gray", drug, tissue.type="", color="#CC3399")
          plot.new()
          plot.new()
          plot.new()
        }
      }
      for(i in 1:biomarkers.no) {
        biomarker <- fnDefineBiomarker (gene.id=as.character(top.significant.biomarkers.drug[i, "gene.id"]), 
                                        isoform.id=as.character(top.significant.biomarkers.drug[i, "transcript.id"]), 
                                        estimate=as.numeric(top.significant.biomarkers.drug[i, "estimate"]),
                                        pvalue=as.numeric(top.significant.biomarkers.drug[i, adjustment.method]), 
                                        effect.size=as.numeric(top.significant.biomarkers.drug[i, effect.size]), 
                                        type=as.character(top.significant.biomarkers.drug[i, "type"]))
        f <- FALSE
        if(biomarker$type == "isoform") {f <- biomarker$id %in% colnames(gray.isoforms.fpkm) & !all(is.na(gray.isoforms.fpkm[,biomarker$id]))}
        if(biomarker$type == "gene") {f <- biomarker$id %in% colnames(gray.genes.fpkm) & !all(is.na(gray.genes.fpkm[,biomarker$id]))}
        if (f) {
          if(plot.create) {
            if(training.type == "CCLE_GDSC") {
              fnPlotEXP(first=paste0("ccle.", biomarker$type), second=paste0("gray.", biomarker$type), biomarker)    
              fnPlotCCLEGDSC(biomarker, tissue.type=tissue, drug=drug)
            }else {
              fnPlotCCLE(biomarker, tissue.type=tissue, drug=drug)
            }
          }
          if(common.cells) {
            gray <- fnValidateWithGray(biomarker, cells="common", color="#336600", drug=drug, plot.create, validation.method="pvalue", validation.cut.off=pvalue.cut.off)
          }else {
            gray <- fnValidateWithGray(biomarker, cells="all", color="#66CC00", drug=drug, plot.create, validation.method="pvalue", validation.cut.off=pvalue.cut.off)
          }
          rr[[drug]][i, "gray.n"] <- gray$n
          rr[[drug]][i, "gray.pvalue"] <- gray$pvalue
          rr[[drug]][i, "gray.estimate"] <- gray$estimate
          rr[[drug]][i, "gray.effect.size"] <- gray$effect.size
          rr[[drug]][i, "validation.stat"] <- gray$validation.stat
          rr[[drug]][i, "biotype"] <- annot.ensembl.all.genes[biomarker$gene.id, "GeneBioType"]
        }
      }
      if(plot.create) {
        dev.off()
      }
    }
  }
  return(rr)
}
fnValidateWithGray <- function(biomarker, cells=c("common","all"), my.xlim, my.ylim, color, drug, plot.create=FALSE, validation.method=c(effect.size, "pvalue"), validation.cut.off, expression.cut.off=0.1) {
  gray.sensitivity <- subset(gray.drug.sensitivity, !is.na(gray.drug.sensitivity[, drug]) , select=drug)
  intersected.celllines.gray <- intersect(rownames(gray.sensitivity), rownames(gray.genes.fpkm))
  if(cells=="common") {
    intersected.celllines.gray <-  intersect(intersected.celllines.gray, rownames(ccle.genes.fpkm))
  }
  if(biomarker$type=="gene") {
    sensitivity <- data.frame(cbind(gray.sensitivity[intersected.celllines.gray, drug], gray.genes.fpkm[intersected.celllines.gray, biomarker$gene.id]))
  }else {
    sensitivity <- data.frame(cbind(gray.sensitivity[intersected.celllines.gray, drug], gray.isoforms.fpkm[intersected.celllines.gray, biomarker$isoform.id]))
  }
  sensitivity <- sensitivity[complete.cases(sensitivity), ,drop=FALSE]
  colnames(sensitivity) <- c(phenotype, "exprs")
  
  gray.n <- nrow(sensitivity)
  gray.model <- lm(sensitivity[, phenotype]  ~ sensitivity[ , "exprs"])
  gray.pvalue <- 2 
  gray.estimate <- gray.effect.size <- 0
  if(all(!is.na(gray.model)) && !is.na(gray.model$coefficients[2]) && length(which(sensitivity[ , "exprs"]!=0)) > (expression.cut.off*nrow(sensitivity)))
  {
    gray.pvalue <- summary(gray.model)$coefficients[2, 4]
    gray.estimate <- summary(gray.model)$coefficients[2, 1]
    if(effect.size=="r.squared"){
      gray.effect.size <- summary(gray.model)$r.squared
    }else if(effect.size=="cindex"){
      #gray.effect.size <- survcomp::concordance.index(x=-predict(gray.model), surv.time=sensitivity[, phenotype], surv.event=rep(1, nrow(sensitivity)), na.rm=TRUE, outx=TRUE)[[1]]
      gray.effect.size <- Hmisc::rcorr.cens(x=predict(gray.model), S=sensitivity[, phenotype], outx=TRUE)[[1]]
    }
  }
  #lm compue two sided pvalue
  gray.pvalue <- gray.pvalue / 2
  if(plot.create) {
    my.xlim <- range(as.numeric(sensitivity[, "exprs"]), na.rm=TRUE)
    my.ylim <- range(as.numeric(sensitivity[, phenotype]), na.rm=TRUE)
    plot(NA, xlim=my.xlim, ylim=my.ylim, xlab=biomarker$label, ylab=phenotype, cex.lab=.7, cex.axis=.7)
    points(sensitivity[, "exprs"], sensitivity[, phenotype], pch=17, col=color)
    title(sprintf("GRAY-%s\n%s and %s", toupper(method), biomarker$type, gsub("drugid_", "", drug)), cex.main=.8)
    textxy(my.xlim[1], my.ylim[2] - .08, sprintf("In GRAY model (Linear Regression)\n pvalue=%s\n%s=%s\nestimate=%s", round(gray.pvalue, digits=10), effect.size, round(gray.effect.size, digits=2), round(gray.estimate, digits=2)), cex=0.6)
  }
  validation.stat <- "unvalidated"
  if(validation.method=="pvalue") {
    if(!is.na(gray.pvalue)) {
      if(as.numeric(gray.pvalue) < validation.cut.off & 
         sign(as.numeric(gray.estimate)) == sign(as.numeric(biomarker$estimate))) {validation.stat <- "validated"}
    }
  } else {
    if(!is.na(gray.effect.size)) {
      if(as.numeric(gray.effect.size) >= validation.cut.off & 
         sign(as.numeric(gray.estimate)) == sign(as.numeric(biomarker$estimate))) {validation.stat <- "validated"}
    }
  }
  return(list(pvalue=gray.pvalue, estimate=gray.estimate, n=gray.n, effect.size=gray.effect.size, validation.stat=validation.stat))
}
fnDefineBiomarker <- function(gene.id, isoform.id, estimate, pvalue, effect.size, type) {
  biomarker <- list()
  biomarker[["gene.id"]] <- gene.id
  biomarker[["isoform.id"]] <- isoform.id
  
  biomarker[["pvalue"]] <- pvalue
  biomarker[["effect.size"]] <- effect.size
  biomarker[["estimate"]] <- estimate
  biomarker[["type"]] <- type
  biomarker[["label"]] <- annot.ensembl.all.genes[which(annot.ensembl.all.genes[ , "EnsemblGeneId"] == gene.id), "Symbol"]
  biomarker[["short.label"]] <- biomarker[["label"]]
  
  if(biomarker[["type"]] == "gene") {
    biomarker[["id"]] <- biomarker[["gene.id"]]
  }else {
    biomarker[["label"]] <- sprintf("%s(%s)", biomarker[["isoform.id"]], biomarker[["label"]])
    biomarker[["short.label"]] <- sprintf("%s.ISO", biomarker[["short.label"]])
    biomarker[["id"]] <- biomarker[["isoform.id"]]
  }
  return(biomarker)
}
fnGtex <- function(all.biomarkers, validation.method=c(stat, "pvalue")) {
  ccle.breast.cells <- rownames(ccle.tissuetype)[which(ccle.tissuetype[ , "tissue.type"]=="breast")]
  ccle.breast.cells <- intersectList(ccle.breast.cells, rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  ccle.breast.genes <- ccle.genes.fpkm[ccle.breast.cells, , drop=FALSE]
  ccle.breast.isoforms <- ccle.isoforms.fpkm[ccle.breast.cells, , drop=FALSE]
  
  gray.breast.genes <- gray.genes.fpkm 
  gray.breast.isoforms <- gray.isoforms.fpkm 
  
  gtex.breast.genes <- t(exprs(GTex.BR$genes))
  gtex.breast.isoforms <-  t(exprs(GTex.BR$isoforms))
  
  breast.no <- max(max(nrow(gtex.breast.genes), nrow(ccle.breast.genes)), nrow(gray.breast.genes))
  mycol <- RColorBrewer::brewer.pal(n=8, name="Set2")[c(3, 4, 6)]
  
  for(i in 1:length(all.biomarkers)) {
    no <- nrow(all.biomarkers[[i]])
    if(no > 0) {
      nn <- ifelse(no >= 4, no%/%4, 1)
      pdf(file=file.path(path.diagrams, sprintf("%s_exp_%s.pdf", names(all.biomarkers)[i], validation.method)), width=14, height=3.5 * nn)
      par(mfrow=c(nn,4))
      for( j in 1:nrow(all.biomarkers[[i]])){
        validation.stat <- all.biomarkers[[i]][j, "validation.stat"]
        if(!is.na(all.biomarkers[[i]][j, "gray.pvalue"]))
        {
          if(all.biomarkers[[i]][j, "type"] == "gene") {
            ensembl.id <- rownames(annot.ensembl.all.genes)[which(annot.ensembl.all.genes[ , "EnsemblGeneId"] == all.biomarkers[[i]][j, "gene.id"])]
            all.biomarkers[[i]][j, "ccle"] <- median(ccle.breast.genes[, as.character(all.biomarkers[[i]][j, "gene.id"])])
            all.biomarkers[[i]][j, "gray"] <- median(gray.breast.genes[, as.character(all.biomarkers[[i]][j, "gene.id"])])
            if(ensembl.id %in% colnames(gtex.breast.genes)) {
              all.biomarkers[[i]][j, "gtex"] <- median(gtex.breast.genes[, ensembl.id])
            }else {
              all.biomarkers[[i]][j, "gtex"] <- NA
            }
            tt <- matrix(NA, ncol=3, nrow=breast.no)
            colnames(tt) <- c("CCLE", "GRAY", "GTEX")
            tt[1:nrow(ccle.breast.genes),1] <- ccle.breast.genes[, as.character(all.biomarkers[[i]][j, "gene.id"])]
            tt[1:nrow(gray.breast.genes),2] <- gray.breast.genes[, as.character(all.biomarkers[[i]][j, "gene.id"])]
            if(ensembl.id %in% colnames(gtex.breast.genes)) {
              tt[1:nrow(gtex.breast.genes),3] <- gtex.breast.genes[, ensembl.id]
            }
          }else{ ##isoform
            all.biomarkers[[i]][j, "ccle"] <- median(ccle.breast.isoforms[, as.character(all.biomarkers[[i]][j, "transcript.id"])])
            all.biomarkers[[i]][j, "gray"] <- median(gray.breast.isoforms[, as.character(all.biomarkers[[i]][j, "transcript.id"])])
            if(as.character(all.biomarkers[[i]][j, "transcript.id"]) %in% colnames(gtex.breast.isoforms)){
              all.biomarkers[[i]][j, "gtex"] <- median(gtex.breast.isoforms[, as.character(all.biomarkers[[i]][j, "transcript.id"])])
            }else{
              all.biomarkers[[i]][j, "gtex"] <- NA
            }
            tt <- matrix(NA, ncol=3, nrow=breast.no)
            colnames(tt) <- c("CCLE", "GRAY", "GTEX")
            tt[1:nrow(ccle.breast.isoforms),1] <- ccle.breast.isoforms[, as.character(all.biomarkers[[i]][j, "transcript.id"])]
            tt[1:nrow(gray.breast.isoforms),2] <- gray.breast.isoforms[, as.character(all.biomarkers[[i]][j, "transcript.id"])]
            if(as.character(all.biomarkers[[i]][j, "transcript.id"]) %in% colnames(gtex.breast.isoforms)){
              tt[1:nrow(gtex.breast.isoforms),3] <- gtex.breast.isoforms[, as.character(all.biomarkers[[i]][j, "transcript.id"])]
            }
          }
          #boxplot(tt, col=mycol, main=sprintf("%s-%s (%s)", all.biomarkers[[i]][j, "symbol"], all.biomarkers[[i]][j, "type"], validation.stat), pch=20)
          invisible(genefu::boxplotplus2(t(tt),
                                         .ylim=c(0,9.5),
                                         pt.col=mycol,
                                         #names=c(paste0("none (", length(res.nosignifseeds), ")"), paste0("all (", length(res.allsignifseeds), ")")),
                                         ylab="expression", 
                                         #xlab="signif seeds (# of targets)",
                                         main=sprintf("%s-%s (%s)", all.biomarkers[[i]][j, "symbol"], all.biomarkers[[i]][j, "type"], validation.stat),
                                         #.las=2,
                                         pt.cex=1))
          train <- test <- NA
          if(!all(is.na(tt[,3]))) {
            train <- wilcox.test(tt[,1], tt[,3], alternative="greater")$p.value
            test <- wilcox.test(tt[,2], tt[,3], alternative="greater")$p.value
          }
          legend("topright", legend=c(sprintf("CCLE > GTEX : %s", round(train, digits=25)), sprintf("GRAY > GTEX : %s", round(test, digits=25))), fill=mycol[1:2], bty="n")
          all.biomarkers[[i]][j, "Gtex.train.sit"] <-  train
          all.biomarkers[[i]][j, "Gtex.test.sit"] <- test
        }
      }
      dev.off()
    }
  }
  
  tt <- lapply(all.biomarkers, function(x){table(x$validation.stat, ifelse(x$Gtex.train.sit< 0.05, "Lower", "Greater"))})
  validated.gtex.lower <- sapply(tt, function(x){if(("validated" %in% rownames(x)) & "Lower" %in% colnames(x)) {x["validated", "Lower"]/sum(x[, "Lower"])}else{0}})
  validated.gtex.greater <- sapply(tt, function(x){if(("validated" %in% rownames(x)) & "Greater" %in% colnames(x)) {x["validated", "Greater"]/sum(x[, "Greater"])} else {0}})
  validated.all <- sapply(tt, function(x){if("validated" %in% rownames(x)) {sum(x["validated",])/sum(x)} else {0}})  
  biomarkers.no <- sapply(all.biomarkers, function(x){nrow(x)})
  isoformstt <- lapply(all.biomarkers, function(x){table(x$validation.stat, x$type)})
  validated.isoforms <- sapply(isoformstt, function(x){if(("validated" %in% rownames(x)) & ("isoform" %in% colnames(x))){x["validated", "isoform"]/sum(x[, "isoform"])} else {0}})
  validated.genes <- sapply(isoformstt, function(x){if(("validated" %in% rownames(x)) & ("gene" %in% colnames(x))) {x["validated", "gene"]/sum(x[, "gene"])} else{0}})
  
  WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=file.path(path.diagrams, sprintf("biomarkers.gtex.%s.xlsx", validation.method)))
  save(all.biomarkers, file=file.path(path.diagrams, sprintf("biomarkers.gtex.%s.xlsx", validation.method)))
  
  validation.res <- cbind("gtex.lower"=validated.gtex.lower, 
                          "gtex.greater"=validated.gtex.greater, 
                          "all"=validated.all, 
                          "isoform"=validated.isoforms,
                          "gene"=validated.genes)
  write.csv(cbind(validation.res,
                  "biomarkers.no"=biomarkers.no), file.path(path.diagrams, sprintf("validated.regarding.gtex.%s.csv", validation.method)))
  pdf(file.path(path.diagrams, sprintf("validated.biomarkers.proportion.%s.pdf", validation.method)), height=7 , width=7)
  mycol <- RColorBrewer::brewer.pal(n=8, name="Set2")[1:5]
  par(mar=c(7,5,5,6))
  par(oma=c(2,2,2,2))
  barplot(t(validation.res), beside=T, las =2, col=mycol, main="Validated breast biomarkers proportion in several categories", ylim=c(0,1))
  legend("topright", legend=colnames(validation.res), fill=mycol, bty="n")
  dev.off()
}
fnPlotCCLEGDSC <- function(biomarker, tissue.type, my.xlim, my.ylim, drug) {
  
  exp.type <- paste0(biomarker$type, ".exp")
  ccle.sensitivity <- subset(ccle.drug.sensitivity, !is.na(ccle.drug.sensitivity[, drug]), select=drug)
  gdsc.sensitivity <- subset(gdsc.drug.sensitivity, !is.na(gdsc.drug.sensitivity[, drug]), select=drug)
  
  unique.ccle.celllines <- setdiff(rownames(ccle.sensitivity), rownames(gdsc.sensitivity))
  unique.ccle.celllines <- intersectList(unique.ccle.celllines, rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  intersected.celllines <- intersectList(rownames(ccle.sensitivity), rownames(gdsc.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  if(biomarker$type == "gene"){
      unique <- cbind(ccle.sensitivity[unique.ccle.celllines, drug], ccle.genes.fpkm[unique.ccle.celllines, biomarker$gene.id], ccle.tissuetype[unique.ccle.celllines, "tissue.type"])
      sensitivity <- cbind(ccle.sensitivity[intersected.celllines,], gdsc.sensitivity[intersected.celllines,], ccle.genes.fpkm[intersected.celllines, biomarker$gene.id], ccle.tissuetype[intersected.celllines, "tissue.type"])
  }else{
      unique <- cbind(ccle.sensitivity[unique.ccle.celllines, drug], ccle.isoforms.fpkm[unique.ccle.celllines, biomarker$isoform.id], ccle.tissuetype[unique.ccle.celllines, "tissue.type"])
      sensitivity <- cbind(ccle.sensitivity[intersected.celllines,], gdsc.sensitivity[intersected.celllines,], ccle.isoforms.fpkm[intersected.celllines, biomarker$isoform.id], ccle.tissuetype[intersected.celllines, "tissue.type"])
  }
  colnames(unique) <- c(paste0("ccle.", phenotype), exp.type, "tissue.type")
  if(tissue.type != "all")
  {
    unique <- subset(unique, unique[, "tissue.type"] == tissue.type)
  }
  colnames(sensitivity) <- c(paste0("ccle.", phenotype), paste0("gdsc.", phenotype), exp.type, "tissue.type")
  if(tissue.type != "all")
  {
    sensitivity <- subset(sensitivity, sensitivity[, "tissue.type"] == tissue.type)
  }
  
  my.xlim <- c(min(min(as.numeric(sensitivity[, exp.type]), na.rm=TRUE), min(as.numeric(unique[, exp.type]), na.rm=TRUE), na.rm=TRUE), max(max(as.numeric(sensitivity[, exp.type])), max(as.numeric(unique[, exp.type]), na.rm=TRUE), na.rm=TRUE))
  my.ylim <- c(min(min(as.numeric(sensitivity[, paste0("ccle.", phenotype)]), na.rm=TRUE), min(as.numeric(unique[, paste0("ccle.", phenotype)]), na.rm=TRUE), na.rm=TRUE), max(max(as.numeric(sensitivity[, paste0("ccle.", phenotype)]), na.rm=TRUE), max(as.numeric(unique[, paste0("ccle.", phenotype)]), na.rm=TRUE), na.rm=TRUE))
  
  sensitivity <- cbind(sensitivity, "grey")
  colnames(sensitivity)[ncol(sensitivity)] <- "col"
  sensitivity <- sensitivity[order(sensitivity[, paste0("gdsc.", phenotype)]),]
  sensitivity[!is.na(sensitivity[, paste0("gdsc.", phenotype)]), "col"] <- colorRampPalette(c("blue", "light blue", "red"))(nrow(sensitivity[!is.na(sensitivity[, paste0("gdsc.", phenotype)]),]))
  
  
  plot(NA, xlim=my.xlim, ylim=my.ylim, xlab=biomarker$label, ylab=paste0("ccle.", phenotype), cex.lab=.7, cex.axis=.7)
  points(unique[, exp.type], unique[, paste0("ccle.", phenotype)], pch=17, col="grey")
  points(sensitivity[, exp.type], sensitivity[, paste0("ccle.", phenotype)], pch=19, col=sensitivity[, "col"])
  #lines(my.xlim, my.ylim, lty=2, col="grey")
  
  x <- c(round(my.xlim[1], digits=1) + 0.1 , round(my.xlim[1], digits=1) + 0.3, round(my.xlim[1], digits=1) + 0.3, round(my.xlim[1], digits=1) + 0.1)
  y <- c(round(my.ylim[2], digits=1) - 0.1, round(my.ylim[2], digits=1) - 0.1, round(my.ylim[2], digits=1), round(my.ylim[2], digits=1))
  #legend.gradient(cbind(x, y), cols=colorRampPalette(c("blue" , "light blue", "red"))( nrow(sensitivity) ), title="gdsc.AUC", limits=c("resistent", "sensitive"), cex=.7)
  
  si <- (nrow(sensitivity) - 2):nrow(sensitivity)
  #textxy(as.numeric(sensitivity[si, "isoform.exp"]), as.numeric(sensitivity[si, "ccle.AUC"]), rownames(sensitivity)[si], cex=0.7)
  ri <- 1:3
  #textxy(as.numeric(sensitivity[ri, "isoform.exp"]), as.numeric(sensitivity[ri, "ccle.AUC"]), rownames(sensitivity)[ri], cex=0.7)
  
  textxy(my.xlim[1] , my.ylim[2] - .08, sprintf("In CCLE/GDSC model \n adj pvalue=%s\n    %s=%s\n estimate=%s", round(biomarker$pvalue, digits=10), effect.size, round(biomarker$effect.size, digits=2), round(biomarker$estimate, digits=2)), cex=0.6)
  
  title(sprintf("CCLE/GDSC\n%s and %s in %s cell lines", biomarker$label, gsub("drugid_", "", drug), tissue.type), cex.main=.8)
  
}
fnPlotCCLE <- function(biomarker, tissue.type, my.xlim, my.ylim, drug) {
  
  exp.type <- paste0(biomarker$type, ".exp")
  ccle.sensitivity <- subset(ccle.drug.sensitivity, !is.na(ccle.drug.sensitivity[, drug]), select=drug)
  
  ccle.celllines <- rownames(ccle.sensitivity)
  ccle <- cbind(ccle.sensitivity[ccle.celllines, drug], ccle.genes.fpkm[ccle.celllines, biomarker$gene.id], ccle.isoforms.fpkm[ccle.celllines, biomarker$isoform.id], ccle.tissuetype[ccle.celllines, "tissue.type"])
  colnames(ccle) <- c(phenotype, "gene.exp", "isoform.exp", "tissue.type")
  if(tissue.type != "all")
  {
    ccle <- subset(ccle, ccle[, "tissue.type"] == tissue.type)
  }
  
  my.xlim <- range(as.numeric(ccle[, exp.type]), na.rm=T)
  my.ylim <- range(as.numeric(ccle[, phenotype]), na.rm=T)
  
  plot(NA, xlim=my.xlim, ylim=my.ylim, xlab=biomarker$label, ylab=phenotype, cex.lab=.7, cex.axis=.7)
  points(ccle[, exp.type], ccle[, phenotype], pch=19, col="blue")
  #lines(my.xlim, my.ylim, lty=2, col="grey")
  
  x <- c(round(my.xlim[1], digits=1) + 0.1 , round(my.xlim[1], digits=1) + 0.3, round(my.xlim[1], digits=1) + 0.3, round(my.xlim[1], digits=1) + 0.1)
  y <- c(round(my.ylim[2], digits=1) - 0.1, round(my.ylim[2], digits=1) - 0.1, round(my.ylim[2], digits=1), round(my.ylim[2], digits=1))
  #legend.gradient(cbind(x, y), cols=colorRampPalette(c("blue" , "light blue", "red"))( nrow(sensitivity) ), title="gdsc.AUC", limits=c("resistent", "sensitive"), cex=.7)
  
  textxy(my.xlim[1] , my.ylim[2] - .08, sprintf("In CCLE model \n adj pvalue=%s\n    %s=%s\n estimate=%s", round(biomarker$pvalue, digits=10), effect.size, round(biomarker$effect.size, digits=2), round(biomarker$estimate, digits=2)), cex=0.6)
  
  title(sprintf("CCLE\n%s and %s in %s cell lines", biomarker$label, gsub("drugid_", "", drug), tissue.type), cex.main=.8)
  
}
fnPlotAUC <- function(first, second, drug, tissue.type, color) {
  switch(first, "gray"={first.db <- gray.drug.sensitivity}, "ccle"={first.db <- ccle.drug.sensitivity}, "gdsc"={first.db <- gdsc.drug.sensitivity})
  switch(second, "gray"={second.db <- gray.drug.sensitivity}, "ccle"={second.db <- ccle.drug.sensitivity}, "gdsc"={second.db <- gdsc.drug.sensitivity})
  
  first.db <- first.db[complete.cases(first.db[, drug]),]
  second.db <- second.db[complete.cases(second.db[, drug]),]
  
  celllines <- intersect(rownames(first.db), rownames(second.db))
  sensitivity <- data.frame(cbind(first.db[celllines, drug], second.db[celllines, drug]), ccle.tissuetype[celllines, "tissue.type"])
  colnames(sensitivity) <- c(first, second, "tissue.type")
  rownames(sensitivity) <- celllines
  tissue.type.str <- ""
  if((tissue.type != "") & (tissue.type != "all"))
  {
    sensitivity <- sensitivity[which(sensitivity[ , "tissue.type"] == tissue.type),] 
    tissue.type.str <- sprintf("in %s cell lines", tissue.type)
  }
  
  my.xlim <- c(min(as.numeric(sensitivity[,1]), na.rm=TRUE), max(as.numeric(sensitivity[,1]), na.rm=TRUE))
  my.ylim <- c(min(as.numeric(sensitivity[,2]), na.rm=TRUE), max(as.numeric(sensitivity[,2]), na.rm=TRUE))
  
  plot(NA, xlim=my.xlim, 
       ylim=my.ylim, 
       xlab=colnames(sensitivity)[1], 
       ylab=colnames(sensitivity)[2],
       main=sprintf("Correlation between %s and %s for %s %s is %s", 
                    toupper(first), 
                    toupper(second) , 
                    gsub("drugid_", "", drug),
                    tissue.type.str, 
                    round(cor(sensitivity[,1], sensitivity[,2], method="pearson"), digits=2)), cex.main=.9)
  #grid(nx=sensitivity[,1], ny=sensitivity[,2], lty=1)
  points(sensitivity[,1], sensitivity[,2], pch=19, col=color)
}
fnPlotEXP <- function(first, second, biomarker) {
  switch(first, "ccle.gene"={first.db <- ccle.genes.fpkm}, "ccle.isoform"={first.db <- ccle.isoforms.fpkm})
  switch(second, "gray.gene"={second.db <- gray.genes.fpkm}, "gray.isoform"={second.db <- gray.isoforms.fpkm})
  
  celllines <- intersect(rownames(first.db), rownames(second.db))
  expression <- data.frame(cbind(first.db[celllines, biomarker$id], second.db[celllines, biomarker$id]))
  colnames(expression) <- c(paste0(biomarker$label, " in CCLE"), paste0(biomarker$label, " in Gray"))
  rownames(expression) <- celllines
  
  my.xlim <- c(min(as.numeric(expression[,1]), na.rm=TRUE), max(as.numeric(expression[,1]), na.rm=TRUE))
  my.ylim <- c(min(as.numeric(expression[,2]), na.rm=TRUE), max(as.numeric(expression[,2]), na.rm=TRUE))
  
  plot(NA, xlim=my.xlim, 
       ylim=my.ylim, 
       xlab=colnames(expression)[1], 
       ylab=colnames(expression)[2],
       main=sprintf("Correlation between expression of %s \nin CCLE vs GRAY dataset is %s", 
                    biomarker$label,
                    round(cor(expression[,1], expression[,2], method="pearson", use="pairwise.complete.obs"), digits=2)), 
       cex.main=.9)
  points(expression[,1], expression[,2], pch=19, col="#3333CC")
}
#tissue <- c("all", "tissue types")
fnFetchBiomarkers <- function(top.significant.biomarkers, drug, indices) {
  biomarkers <- list()
  j <- 1;
  for(i in indices)
  {
    common <- which(top.significant.biomarkers$symbol == as.character(top.significant.biomarkers[i,"symbol"]))
    if(length(common) == 1){delta.rank <- "-"}else{delta.rank <- as.character(i - common[which(common != i)])}
    #ifelse(top.significant.biomarkers[common[which(common != i)], "rank"] < top.significant.biomarkers[i, "rank"], as.character(top.significant.biomarkers[i,"delta.rank"]),as.character(-1*top.significant.biomarkers[i,"delta.rank"]))}
    biomarker <- fnDefineBiomarker2 (gene.id=as.character(top.significant.biomarkers[i, "gene.id"]), 
                                    isoform.id=as.character(top.significant.biomarkers[i, "transcript.id"]), 
                                    estimate=as.numeric(top.significant.biomarkers[i, "estimate"]),
                                    pvalue=as.numeric(top.significant.biomarkers[i, adjustment.method]), 
                                    gray.estimate=as.numeric(top.significant.biomarkers[i, "all.gray.estimate"]),
                                    gray.pvalue=as.numeric(top.significant.biomarkers[i, "all.gray.pvalue"]), 
                                    effect.size.value=as.numeric(top.significant.biomarkers[i, effect.size]), 
                                    type=as.character(top.significant.biomarkers[i, "type"]),
                                    biotype=as.character(top.significant.biomarkers[i, "biotype"]),
                                    gtex=ifelse(as.numeric(top.significant.biomarkers[i, "Gtex.train.sit"]) < 0.05, "tumor.specific","no.diff"),
                                    delta.rank=delta.rank,
                                    isoform.no=as.character(top.significant.biomarkers[i,"isoforms.no"]))
    biomarkers[[j]] <- biomarker
    j <- j + 1;
  }
  return(biomarkers)
}
fnPlotAUCoverCellLinesCCLE.GDSC <- function(drug, tissue.type, biomarkers, biomarkers.order) {
  ccle.sensitivity <- subset(ccle.drug.sensitivity, !is.na(ccle.drug.sensitivity[,drug]), select=drug)
  cells <- intersectList(rownames(ccle.sensitivity), rownames(gdsc.drug.sensitivity), rownames(ccle.tissuetype), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  sensitivity <- cbind(ccle.sensitivity[cells, ], gdsc.drug.sensitivity[cells, drug], ccle.tissuetype[cells, "tissue.type"])
  colnames(sensitivity) <- c("ccle.AUC","gdsc.AUC","tissue.type")
  rownames(sensitivity) <- cells
  
  tissue.type.str <- ""
  if(tissue.type != "all")
  {
    sensitivity <- sensitivity[which(sensitivity[ , "tissue.type"] == tissue.type),] 
    tissue.type.str <- sprintf("in %s cell lines", tissue.type)
  }
  
  sensitivity <- cbind(sensitivity, "col"="#000000")
  sensitivity <- sensitivity[order(sensitivity[,"gdsc.AUC"]),]
  sensitivity[!is.na(sensitivity[,"gdsc.AUC"]),"col"] <- colorRampPalette(c("blue","light blue","red"))(nrow(sensitivity[!is.na(sensitivity[,"gdsc.AUC"]),]))
  
  sensitivity <- sensitivity[order(sensitivity[,"ccle.AUC"]),]
  
  exp.db <- NULL;
  exp.col <- NULL;
  for (i in 1:length(biomarkers))
  {
    if (biomarkers[[i]]$type == "isoform")
    {
      exp.db <- cbind(exp.db, ccle.isoforms.fpkm[rownames(sensitivity),biomarkers[[i]]$isoform.id])
    }else{
      exp.db <- cbind(exp.db, ccle.genes.fpkm[rownames(sensitivity),biomarkers[[i]]$gene.id])      
    }
  }
  for(i in 1:ncol(exp.db))
  {
    exp.db[,i] <- (exp.db[,i] - mean(exp.db[,i], na.rm=T))/sd(exp.db[,i], na.rm=T)
    exp.col <- union(exp.col, exp.db[,i])  
  }
  exp.col <- sort(exp.col)
  exp.col <- data.frame("exp"=exp.col, "col" ="#000000" )
  exp.col[,"col"] <- colorRampPalette(c(blue,"white",red))(nrow(exp.col))
  
  colnames(exp.db) <- sapply(biomarkers, function(x){x[["short.label"]]})
  if(!is.null(biomarkers.order)) {exp.db <- exp.db[ , biomarkers.order, drop=FALSE]}
  #exp.db <- exp.db[-nrow(exp.db),]
  #sensitivity <- sensitivity[-nrow(sensitivity),]
  
  
  my.xlim <- c(1, nrow(sensitivity))
  my.ylim <- as.numeric(range(sensitivity[,1]))
  
  # color.scheme <- rgb( seq(0,1,length=256),  # Red
  #                     seq(0,1,length=256),  # Green
  #                     seq(1,0,length=256))  # Blue
  
  #exp.db <- exp.db[, biomarkers.toPlot]
  pdf(file=file.path(path.diagrams, sprintf("%s_%s_CCLE_GDSC.pdf", gsub("drugid_","",drug), phenotype)), height=7, width=14)  
  par(mar=c(5, 8, 2, 2))
  #par(mfrow=c(2,1))
  image(exp.db, col=exp.col[,"col"], axes=FALSE)
  grid(nx=nrow(exp.db), ny=ncol(exp.db), lty=1)
  if(ncol(exp.db) != 1){at.place=(0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1)}else{at.place=.5}
  axis(2,at=at.place, labels=colnames(exp.db), las=2, cex.axis=.6, tick=FALSE)
  #plot(1:nrow(sensitivity),sensitivity[,1], type="l", xlim=my.xlim, ylim=my.ylim, xlab="cell lines", ylab=colnames(sensitivity)[1],main=sprintf("%s AUC over %s cell lines in %s", toupper(dataset.name), tissue.type, gsub("drugid_","",drug)), cex.main=.9)
  dev.off()
  pdf(file=file.path(path.diagrams, sprintf("%s_%s_CCLE_GDSC_2.pdf", gsub("drugid_","",drug), phenotype)), height=5, width=12)    
  par(mar=c(12, 5, 2, 8))
  plot(NA, xlim=my.xlim, ylim=my.ylim, ylab='', xlab='', axes=FALSE)
  axis(1, at=1:nrow(sensitivity), labels=rownames(sensitivity), las=2, cex.axis=1.5, tck=-.05)
  axis(2, at=c(0, .1, .2, .3, .4, .5), labels=c(0, .1, .2, .3, .4, .5), cex.axis=1.5, tck=-.02)
  
  box(lty=1)
  points(1:nrow(sensitivity),sensitivity[,"ccle.AUC"], pch=20,col=sensitivity[,"col"])
  
  dev.off()
}
fnPlotAUCoverCellLinesCCLE.GDSC.union.cells <- function(drug, tissue.type, biomarkers, biomarkers.order) {
  ccle.sensitivity <- ccle.drug.sensitivity[which(!is.na(ccle.drug.sensitivity[ , drug])) , drug, drop=FALSE]
  gdsc.sensitivity <- gdsc.drug.sensitivity[which(!is.na(gdsc.drug.sensitivity[ , drug])) , drug, drop=FALSE]

  ccle.cells <- intersectList(rownames(ccle.sensitivity), rownames(ccle.tissuetype), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  gdsc.cells <- intersectList(rownames(gdsc.sensitivity),rownames(ccle.tissuetype), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  cells <- union(ccle.cells, gdsc.cells)
  
  sensitivity <- matrix(NA, ncol=5, nrow=length(cells), dimnames=list(cells, c("ccle.AUC", "gdsc.AUC", "avg.AUC", "union.AUC","tissue.type")))
  sensitivity[ccle.cells, "ccle.AUC"] <- ccle.sensitivity[ccle.cells, ]
  sensitivity[gdsc.cells, "gdsc.AUC"] <- gdsc.sensitivity[gdsc.cells, ]
  sensitivity[cells, "tissue.type"] <- ccle.tissuetype[cells, "tissue.type"]
  sensitivity[ , "avg.AUC"] <- rowMeans(cbind(as.numeric(sensitivity[, "ccle.AUC"]), as.numeric(sensitivity[, "gdsc.AUC"])), na.rm=TRUE)
  sensitivity[ , "union.AUC"] <- sensitivity[ , "ccle.AUC"]
  sensitivity[which(is.na(sensitivity[ , "union.AUC"])), "union.AUC"] <- sensitivity[which(is.na(sensitivity[ , "union.AUC"])) , "gdsc.AUC"]

  tissue.type.str <- ""
  if(tissue.type != "all")
  {
    sensitivity <- sensitivity[which(sensitivity[ , "tissue.type"] == tissue.type),] 
    tissue.type.str <- sprintf("in %s cell lines", tissue.type)
  }
  
  sensitivity <- cbind(sensitivity, "col"="#000000")
  sensitivity <- cbind(sensitivity, "pch"=20)
  sensitivity[which(is.na(sensitivity[ , "ccle.AUC"])), "pch"] <- 17
  sensitivity <- sensitivity[order(sensitivity[ , "gdsc.AUC"]), ]
  sensitivity[which(!is.na(sensitivity[ , "ccle.AUC"]) & !is.na(sensitivity[ , "gdsc.AUC"])), "col"] <- colorRampPalette(c("blue", "light blue", "red"))(nrow(sensitivity[which(!is.na(sensitivity[ , "ccle.AUC"]) & !is.na(sensitivity[ , "gdsc.AUC"])), ]))
  
  sensitivity <- sensitivity[order(sensitivity[ , "avg.AUC"]),]
  
  exp.db <- NULL;
  exp.col <- NULL;
  for (i in 1:length(biomarkers))
  {
    if (biomarkers[[i]]$type == "isoform")
    {
      exp.db <- cbind(exp.db, ccle.isoforms.fpkm[rownames(sensitivity), biomarkers[[i]]$isoform.id])
    }else{
      exp.db <- cbind(exp.db, ccle.genes.fpkm[rownames(sensitivity), biomarkers[[i]]$gene.id])      
    }
  }
  for(i in 1:ncol(exp.db))
  {
    exp.db[,i] <- (exp.db[,i] - mean(exp.db[,i], na.rm=T))/sd(exp.db[,i], na.rm=T)
    exp.col <- union(exp.col, exp.db[,i])  
  }
  exp.col <- sort(exp.col)
  exp.col <- data.frame("exp"=exp.col, "col" ="#000000" )
  exp.col[,"col"] <- colorRampPalette(c(blue, "white", red))(nrow(exp.col))
  
  colnames(exp.db) <- sapply(biomarkers, function(x){x[["short.label"]]})
  if(!is.null(biomarkers.order)) {exp.db <- exp.db[ , biomarkers.order, drop=FALSE]}
  #exp.db <- exp.db[-nrow(exp.db),]
  #sensitivity <- sensitivity[-nrow(sensitivity),]
  
  
  my.xlim <- c(1, nrow(sensitivity))
  my.ylim <- as.numeric(range(sensitivity[ ,"union.AUC"]))
  
  # color.scheme <- rgb( seq(0,1,length=256),  # Red
  #                     seq(0,1,length=256),  # Green
  #                     seq(1,0,length=256))  # Blue
  
  #exp.db <- exp.db[, biomarkers.toPlot]
  pdf(file=file.path(path.diagrams, sprintf("%s_%s_CCLE_GDSC_union.pdf", gsub("drugid_","",drug), phenotype)), height=7, width=14)  
  par(mar=c(5, 8, 2, 2))
  #par(mfrow=c(2,1))
  image(exp.db, col=exp.col[,"col"], axes=FALSE)
  grid(nx=nrow(exp.db), ny=ncol(exp.db), lty=1)
  if(ncol(exp.db) != 1){at.place=(0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1)}else{at.place=.5}
  axis(2,at=at.place, labels=colnames(exp.db), las=2, cex.axis=.6, tick=FALSE)
  #plot(1:nrow(sensitivity),sensitivity[,1], type="l", xlim=my.xlim, ylim=my.ylim, xlab="cell lines", ylab=colnames(sensitivity)[1],main=sprintf("%s AUC over %s cell lines in %s", toupper(dataset.name), tissue.type, gsub("drugid_","",drug)), cex.main=.9)
  dev.off()
  pdf(file=file.path(path.diagrams, sprintf("%s_%s_CCLE_GDSC_union_2.pdf", gsub("drugid_","",drug), phenotype)), height=5, width=12)    
  par(mar=c(12, 5, 2, 8))
  plot(NA, xlim=my.xlim, ylim=my.ylim, ylab='', xlab='', axes=FALSE)
  axis(1, at=1:nrow(sensitivity), labels=rownames(sensitivity), las=2, cex.axis=1.5, tck=-.05)
  axis(2, at=c(0, .1, .2, .3, .4, .5), labels=c(0, .1, .2, .3, .4, .5), cex.axis=1.5, tck=-.02)
  
  box(lty=1)
  points(1:nrow(sensitivity),sensitivity[ ,"union.AUC"], pch=as.numeric(sensitivity[ , "pch"]), col=sensitivity[ , "col"])
  
  dev.off()
}
fnPlotAUCoverCellLinesGray <- function(drug, tissue.type, biomarkers, gray.specificity=NULL, suffix)  {
  sensitivity <- subset(gray.drug.sensitivity, !is.na(gray.drug.sensitivity[,drug]), select=drug)
  sensitivity <- cbind(sensitivity , "col"="#000000")
  sensitivity <- sensitivity[order(as.numeric(sensitivity[,drug])),]
  sensitivity <- subset(sensitivity, rownames(sensitivity) %in% rownames(gray.genes.fpkm))
  exp.db <- NULL;
  exp.col <- NULL;
  GeneBioType <- NULL;
  for (i in 1:length(biomarkers))
  {
    if (biomarkers[[i]]$type == "isoform")
    {
      exp.db <- cbind(exp.db, gray.isoforms.fpkm[rownames(sensitivity),biomarkers[[i]]$isoform.id])  
      #      exp.db <- cbind(exp.db, (gray.isoforms.fpkm[rownames(sensitivity),biomarkers[[i]]$isoform.id] - 
      #                                apply(gray.isoforms.fpkm[rownames(sensitivity),], 1, mean, na.rm=T))/
      #                                apply(gray.isoforms.fpkm[rownames(sensitivity),], 1, sd, na.rm=T))
    }else{
      exp.db <- cbind(exp.db, gray.genes.fpkm[rownames(sensitivity),biomarkers[[i]]$gene.id])  
      #      exp.db <- cbind(exp.db, (gray.genes.fpkm[rownames(sensitivity),biomarkers[[i]]$gene.id] - 
      #                                apply(gray.genes.fpkm[rownames(sensitivity),], 1, mean, na.rm=T))/
      #                                apply(gray.genes.fpkm[rownames(sensitivity),], 1, sd, na.rm=T))
    }
  }

  for(i in 1:ncol(exp.db))
  {
    exp.db[,i] <- (exp.db[,i] - mean(exp.db[,i], na.rm=T))/sd(exp.db[,i], na.rm=T)
    exp.col <- union(exp.col, exp.db[,i])  
  }
  #exp.col <- as.vector(unique(exp.db))
  names(exp.col) <- 1:length(exp.col)
  exp.col <- data.frame("exp"=exp.col, "col" ="#000000" )
  exp.col <- exp.col[order(exp.col[,"exp"]),]
  exp.col[,"col"] <- colorRampPalette(c(blue, "white", red))(nrow(exp.col))
  rownames(exp.db) <- rownames(sensitivity)
  colnames(exp.db) <- sapply(biomarkers, function(x){x[["short.label"]]})
  
  quantile.range <- quantile(exp.db, probs = seq(0, 1, 0.01), na.rm=T)
  palette.breaks <- seq(quantile.range["0%"], quantile.range["100%"], 0.05)
  color.palette  <- colorRampPalette(c(blue, "white", red))(length(palette.breaks) - 1)
  if(drug=="lapatinib"){color.palette  <- colorRampPalette(c(blue, "white", red, red))(length(palette.breaks) - 1)}
  #if(drug=="Erlotinib"){color.palette  <- colorRampPalette(c(blue, "white", red, red))(length(palette.breaks) - 1)}
  if(drug=="paclitaxel"){color.palette  <- colorRampPalette(c(blue, blue, "white", red, red))(length(palette.breaks) - 1)}
  exp.db <- exp.db[-nrow(exp.db), , drop=FALSE]
  sensitivity <- sensitivity[-nrow(sensitivity), , drop=FALSE]
  
  
  my.xlim <- c(1,nrow(sensitivity))
  my.ylim <- as.numeric(range(sensitivity[,1]))
  
  #exp.db <- exp.db[, biomarkers.toPlot]
  #par(mfrow=c(2,1))
  #library("gplots")
  #gplots::heatmap.2(t(exp.db), Colv=NA, Rowv=T, col=exp.col[,"col"], scale="row", trace="none", dendrogram="none", )
  coding.biotypes <- c("PRT", "AS", "prcTR", "pseudogene", "ncRNA")
  names(coding.biotypes) <- c("protein_coding", "antisense", "processed_transcript", "processed_pseudogene", "3prime_overlapping_ncRNA")
  #xx <- sapply(biomarkers, function(x){
  #  ifelse(x[["gtex"]] == "tumor.specific",
  #         sprintf("%s (%s) *", x[["short.label"]], ifelse(x[["biotype"]] %in% names(coding.biotypes), coding.biotypes[x[["biotype"]]], x[["biotype"]])), 
  #         sprintf("%s (%s)", x[["short.label"]], ifelse(x[["biotype"]] %in% names(coding.biotypes), coding.biotypes[x[["biotype"]]], x[["biotype"]])))})
  xx <- sapply(biomarkers, function(x){
    sprintf("%s (%s)", x[["short.label"]], ifelse(x[["biotype"]] %in% names(coding.biotypes), coding.biotypes[x[["biotype"]]], x[["biotype"]]))})

  label.col <- rep("black", length(xx))
  if(!is.null(gray.specificity)) {
    for( i in 1:length(gray.specificity)){
      if(gray.specificity[i]=="isoform.specific") {
        label.col[i]<- "green4"
        }else if(gray.specificity[i]=="gene.specific") {
          label.col[i]<- "red"
          }
      }
    }

  if(ncol(exp.db) == 1){
    pdf(file=file.path(path.diagrams, sprintf("%s_%s_Gray_%s.pdf", gsub("drugid_","",drug), phenotype, suffix)), height=2, width=22)  
    par(mar=c(1, 1, 1, 17))
    par(oma=c(2,2,2,2))
    image(exp.db, col=exp.col[,"col"], axes=FALSE)
    grid(nx=nrow(exp.db), ny=ncol(exp.db), lty=1)
    axis(4, at=0, labels=xx, las=2, cex.axis=1.5, tick=FALSE, col=label.col)
    #axis(4,at=(0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1), labels=sapply(biomarkers, function(x){ifelse(x[["gtex"]] == "tumor.specific", "*", NA)}, simplify=T), las=2, cex.axis=.8, tick=FALSE)
    
  }else {
    ff <- FALSE
    if(length(xx) < 4){ff <- TRUE}
    pdf(file=file.path(path.diagrams, sprintf("%s_%s_Gray_%s.pdf", gsub("drugid_","",drug), phenotype, suffix)), height=ifelse(ff, 3, 9) , width=17)
    #par(mar=c(2, 2, 2, 8))
   
    if(drug=="paclitaxel"){par(oma=c(0,0,0,13))} else{ par(oma=c(0,0,0,10))}
    if(drug=="Nutlin-3"){cr <- 1.3} else{cr <- ifelse(ff, 2, 0.2 + 1.3/log10(ncol(exp.db)))}
    # hv <- gplots::heatmap.2(t(exp.db), Colv=NA, Rowv=T, dendrogram="none", col=exp.col[,"col"], scale="row", trace="none", key=FALSE,
    #                         labRow=xx, colRow=label.col, cexRow = 0.2 + 1.3/log10(ncol(exp.db)),
    #                         #labRow=biomarkers.toPlot,
    #                         labCol=NA)
    hv <- gplots::heatmap.2(t(exp.db), Colv=NA, Rowv=T, dendrogram="none", col=color.palette, breaks=palette.breaks, trace="none", key=FALSE,
                            labRow=xx, colRow=label.col, cexRow=cr,
                            #labRow=biomarkers.toPlot,
                            labCol=NA)
    
  }
  #image(exp.db, col=exp.col[,"col"], axes=FALSE)
  #grid(nx=nrow(exp.db), ny=ncol(exp.db), lty=1)
  #axis(2,at=(0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1), labels=colnames(exp.db), las=2, cex.axis=.6, tick=FALSE)
  #axis(4,at=(0:(ncol(exp.db) - 1))/(ncol(exp.db) - 1), labels=sapply(biomarkers, function(x){ifelse(x[["gtex"]] == "tumor.specific", "*", NA)}, simplify=T), las=2, cex.axis=.8, tick=FALSE)
  #mtext(side=4, line=-3, outer=F, text=xx[hv$rowInd], at=1:length(xx), las=2, col=label.col)
  dev.off()
  pdf(file=file.path(path.diagrams, sprintf("%s_%s_AUC_Gray_2_%s.pdf", gsub("drugid_","",drug), phenotype, suffix)), height=5, width=15)  
  par(mar=c(12,5,2,8))
  plot(NA, xlim=my.xlim, ylim=my.ylim,ylab='',xlab='', axes=FALSE)
  axis(1,at=1:nrow(sensitivity), labels=rownames(sensitivity), las=2, cex.axis=1.5, tck=-.05)
  axis(2,at=c(0,.1,.2,.3,.4,.5), labels=c(0,.1,.2,.3,.4,.5), cex.axis=1.5, tck=-.02)
  
  box(lty=1)
  points(1:nrow(sensitivity),sensitivity[,drug], pch=20, cex=1.5, col="#663399")
  
  dev.off()
  if(ncol(exp.db) > 1)
  {
    return (list("hv"=hv, "label.col"=label.col[which(label.col != "black")]))
  }else{
    return(NA)
  }
}
fnPlotEffectSize <- function (drug, biomarkers, effect.size, biomarkers.order)  {
  biomarkers.effect.size <-  sapply(biomarkers, function(x){x[effect.size]})
  biomarkers.effect.size <- do.call("c", biomarkers.effect.size)
  biomarkers.effect.size <- biomarkers.effect.size - .5
  biomarkers.estimate <-  sapply(biomarkers, function(x){x["estimate"]})
  biomarkers.estimate <- do.call("c", biomarkers.estimate)
  
  if(!is.null(biomarkers.order)) {biomarkers.effect.size <- as.numeric(biomarkers.effect.size[biomarkers.order]) * sign(biomarkers.estimate[biomarkers.order])}
  names(biomarkers.effect.size) <- NULL
  mycol3 <- RColorBrewer::brewer.pal(n=4, name="Set3")
  
  pdf(file.path(path.diagrams,sprintf("%s_%s_GRAY.pdf", drug, effect.size)), height=8, width=2)
  barplot(biomarkers.effect.size , horiz=T, col=sapply(biomarkers.effect.size , function(x){ifelse(x>=0,mycol3[4], mycol3[3])}),axes=FALSE)
  axis(1, at=seq(-.2, .2, .1), cex=.7)
  dev.off()
  
}
fnPlotLogPvalue <- function (drug, biomarkers, pvalue.col.name, effect.size.col, biomarkers.order)  {
  biomakers.label <- sapply(biomarkers, function(x){x["short.label"]})
  biomarkers.pvalue <-  sapply(biomarkers, function(x){x[pvalue.col.name]})
  names(biomarkers.pvalue) <- biomakers.label
  biomarkers.pvalue <- do.call("c", biomarkers.pvalue)
  
  biomarkers.estimate <-  sapply(biomarkers, function(x){x[effect.size.col]})
  names(biomarkers.estimate) <- biomakers.label
  biomarkers.estimate <- do.call("c", biomarkers.estimate)
  
  if(!is.null(biomarkers.order)) {biomarkers.pvalue <- -log10(as.numeric(biomarkers.pvalue[biomarkers.order])) * sign(biomarkers.estimate[biomarkers.order])}
  names(biomarkers.pvalue) <- NULL
  mycol3 <- RColorBrewer::brewer.pal(n=4, name="Set3")
  
  pdf(file.path(path.diagrams,sprintf("%s_%s.pdf", drug, pvalue.col.name)), height=8, width=2)
  barplot(biomarkers.pvalue , horiz=T, col=sapply(biomarkers.pvalue , function(x){ifelse(x>=0, mycol3[4], mycol3[3])}))
  dev.off()
  
}
fnDefineBiomarker2 <- function(gene.id, isoform.id, estimate, pvalue, gray.estimate, gray.pvalue, effect.size.value, type, biotype, gtex, delta.rank, isoform.no) {
  biomarker <- list()
  biomarker[["gene.id"]] <- gene.id
  biomarker[["isoform.id"]] <- isoform.id
  biomarker[["gtex"]] <- gtex
  
  biomarker[[effect.size]] <- effect.size.value
  biomarker[["estimate"]] <- estimate
  biomarker[["pvalue"]] <- pvalue
  biomarker[["gray.estimate"]] <- gray.estimate
  biomarker[["gray.pvalue"]] <- gray.pvalue
  biomarker[["type"]] <- type
  biomarker[["biotype"]] <- biotype
  biomarker[["label"]] <- annot.ensembl.all.genes[which(annot.ensembl.all.genes[,"EnsemblGeneId"] == gene.id),"Symbol"]
  biomarker[["short.label"]] <- biomarker[["label"]]
  biomarker[["delta.rank"]] <- delta.rank
  biomarker[["isoform.no"]] <- isoform.no
  
  if(biomarker[["type"]] == "gene"){
    biomarker[["id"]] <- biomarker[["gene.id"]]
  }else{
    biomarker[["label"]] <- sprintf("%s(%s)",biomarker[["isoform.id"]], biomarker[["label"]])
    biomarker[["short.label"]] <- sprintf("%s.ISO", biomarker[["short.label"]])
    biomarker[["id"]] <- biomarker[["isoform.id"]]
  }
  return(biomarker)
}
