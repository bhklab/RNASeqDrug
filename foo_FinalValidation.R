fnBuildLinearModel <- function(cells, sensitivity, features, fpkm.matrix) {
  complete.models <- list()
  models <- data.frame(matrix(NA, ncol=5, nrow=length(features), dimnames=list(features, c("pvalue", "estimate", "R2", "se", "n"))), stringsAsFactors=FALSE)
  for(feature in features) {
    if(!all(is.na(fpkm.matrix[cells, feature]))) {
      tt <- cbind("auc"=sensitivity[cells], "exp"=fpkm.matrix[cells, feature])
      tt <- tt[complete.cases(tt), , drop=FALSE]
      n <- nrow(tt)
      model <- lm(tt[ ,"auc"]  ~ tt[ ,"exp"])
      models[feature,"n"] <- n
      
      if(all(!is.na(model)) & !is.na(model$coefficients[2]))
      {
        models[feature,"pvalue"] <- summary(model)$coefficients[2,4]
        models[feature,"estimate"] <- summary(model)$coefficients[2,1]
        models[feature,"R2"] <- summary(model)$adj.r.squared
        models[feature,"se"] <- summary(model)$coefficients[2,2]
      }
      complete.models[[feature]] <- model
    }
  }
  return(list("simplified"=models, "complete"=complete.models))
}
fnNormalizeTraining <- function(ccle.model, gdsc.model) {
  xx <- data.frame(cbind("estimate"=ccle.model$estimate * ccle.model$n/(ccle.model$n + gdsc.model$n) + gdsc.model$estimate * gdsc.model$n/(ccle.model$n + gdsc.model$n),
                         "pvalue"=ccle.model$pvalue * ccle.model$n/(ccle.model$n + gdsc.model$n) + gdsc.model$pvalue * gdsc.model$n/(ccle.model$n + gdsc.model$n),
                         "R2"=ccle.model$R2 * ccle.model$n/(ccle.model$n + gdsc.model$n) + gdsc.model$R2 * gdsc.model$n/(ccle.model$n + gdsc.model$n),
                         "se"=ccle.model$se * ccle.model$n/(ccle.model$n + gdsc.model$n) + gdsc.model$se * gdsc.model$n/(ccle.model$n + gdsc.model$n)), stringsAsFactors=FALSE)
  rownames(xx) <- rownames(ccle.model)
  return(xx)
}
fnPlotHeatMap <- function(sensitivity, expression, file.name, cluster=FALSE, best.isoform) {
  mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
  red <- mycol[1]  
  blue <- mycol[2]
  
  sensitivity <- sensitivity[which(!is.na(sensitivity))]
  sensitivity <- sensitivity[order(sensitivity)]
  exp.db <- expression[names(sensitivity), ,drop=FALSE]
  xx <- which(apply(exp.db, MARGIN = 2, function(x){all(x == 0) | all(is.na(x))}))
  if(length(xx) > 0){
    exp.db <- exp.db[ ,-xx]
  }
  exp.col = NULL;
  
  for(i in 1:ncol(exp.db))
  {
    exp.db[,i] <- (exp.db[,i]-min(exp.db[,i], na.rm=T))/(max(exp.db[,i], na.rm=T)-min(exp.db[,i], na.rm=T))
    exp.col <- union(exp.col, exp.db[,i])  
  }
  #exp.col <- as.vector(unique(exp.db))
  names(exp.col) = 1:length(exp.col)
  exp.col = data.frame("exp" = exp.col, "col" ="#000000" )
  exp.col = exp.col[order(exp.col[,"exp"]),]
  exp.col[,"col"] = colorRampPalette(c(blue, "white",red))(nrow(exp.col))
  quantile.range <- quantile(exp.db, probs = seq(0, 1, 0.01), na.rm=T)
  palette.breaks <- seq(quantile.range["0%"], quantile.range["100%"], 0.1)
  color.palette  <- colorRampPalette(c(blue, "white", red))(length(palette.breaks) - 1)

  pdf(file = file.path(path.diagrams, sprintf("%s_isofoms.pdf", file.name)), height=9, width=17) 
  exp.db.iso <- exp.db[,-ncol(exp.db), drop=FALSE]
  colnames(exp.db.iso)[which(colnames(exp.db.iso) == best.isoform)] <- sprintf("%s***", best.isoform)
  xx <- which(apply(exp.db.iso, MARGIN=2, function(x){which(all(is.na(x)))}) == 1)
  if(length(xx) > 0) {exp.db.iso <- exp.db.iso[, -xx, drop=FALSE]}
  if(ncol(exp.db.iso) == 1) {
    par(mar=c(9, 2, 5, 15))
    par(oma=c(2,2,2,2))
#    image(exp.db.iso, col = exp.col[,"col"], axes = FALSE)
    image(exp.db.iso, col=color.palette, axes = FALSE)
    grid(nx = nrow(exp.db.iso), ny = ncol(exp.db.iso), lty = 1)
    axis(4,at = 0, labels=colnames(exp.db.iso), las=2, cex.axis = 2, tick = FALSE)
  }else {
    par(oma=c(0,0,0,13))
    if(cluster){
        # hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=T, dendrogram="none", col=exp.col[,"col"], scale="row", trace="none", key=FALSE,
        #                         labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
        #                         #labRow=biomarkers.toPlot,
        #                         labCol=NA)
        hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=T, dendrogram="none", col=color.palette, breaks=palette.breaks, trace="none", key=FALSE,
                                labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
                                #labRow=biomarkers.toPlot,
                                labCol=NA)
    }else{
        # hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=NA, dendrogram="none", col=exp.col[,"col"], scale="row", trace="none", key=FALSE,
        #                         labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
        #                         #labRow=biomarkers.toPlot,
        #                         labCol=NA)
      hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=NA, dendrogram="none", col=color.palette, breaks=palette.breaks, trace="none", key=FALSE,
                              labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
                              #labRow=biomarkers.toPlot,
                              labCol=NA)
    }
  }
  dev.off()
  
  pdf(file = file.path(path.diagrams, sprintf("%s_gene.pdf", file.name)), height=5, width=12) 
  exp.db.gene <- exp.db[,ncol(exp.db), drop=FALSE]
  par(mar=c(7,1,7,18))
  par(oma=c(2,2,2,2))
#  image(exp.db.gene, col = exp.col[,"col"], axes = FALSE)
  image(exp.db.gene, col=color.palette, axes = FALSE)
  grid(nx = nrow(exp.db.gene), ny = ncol(exp.db.gene), lty = 1)
  axis(4,at = 0, labels=colnames(exp.db.gene), las=2, cex.axis = 2, tick = FALSE)
  dev.off()
  if(cluster & ncol(exp.db.iso) > 1){
    return(gsub("[***]", "", colnames(exp.db.iso)[rev(hv$rowInd)]))
  }
}
fnPlotHeatMap.union.cells <- function(ccle.sensitivity, gdsc.sensitivity, expression, file.name, cluster=FALSE, best.isoform) {
  mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
  red <- mycol[1]  
  blue <- mycol[2]
  
  ccle.cells <- names(ccle.sensitivity)
  gdsc.cells <- names(gdsc.sensitivity)
  cells <- union(names(ccle.sensitivity), names(gdsc.sensitivity))
  
  sensitivity <- matrix(NA, ncol=4, nrow=length(cells), dimnames=list(cells, c("ccle", "gdsc", "avg", "union")))
  sensitivity[ccle.cells, "ccle"] <- ccle.sensitivity[ccle.cells]
  sensitivity[gdsc.cells, "gdsc"] <- gdsc.sensitivity[gdsc.cells]
  sensitivity[ , "avg"] <- rowMeans(cbind(as.numeric(sensitivity[, "ccle"]), as.numeric(sensitivity[, "gdsc"])), na.rm=TRUE)
  sensitivity[ , "union"] <- sensitivity[ , "ccle"]
  sensitivity[which(is.na(sensitivity[ , "union"])), "union"] <- sensitivity[which(is.na(sensitivity[ , "union"])) , "gdsc"]
  
  sensitivity <- cbind(sensitivity, "col"="#000000")
  sensitivity <- cbind(sensitivity, "pch"=20)
  sensitivity[which(is.na(sensitivity[ , "ccle"])), "pch"] <- 17
  sensitivity <- sensitivity[order(sensitivity[ , "gdsc"]), ]
  sensitivity[which(!is.na(sensitivity[ , "ccle"]) & !is.na(sensitivity[ , "gdsc"])), "col"] <- colorRampPalette(c("blue", "light blue", "red"))(nrow(sensitivity[which(!is.na(sensitivity[ , "ccle"]) & !is.na(sensitivity[ , "gdsc"])), ]))
  
  sensitivity <- sensitivity[order(sensitivity[ , "avg"]),]
  
  exp.db <- expression[rownames(sensitivity), ,drop=FALSE]
  xx <- which(apply(exp.db, MARGIN = 2, function(x){all(x == 0) | all(is.na(x))}))
  if(length(xx) > 0){
    exp.db <- exp.db[ ,-xx]
  }
  exp.col = NULL;
  
  for(i in 1:ncol(exp.db))
  {
    exp.db[,i] <- (exp.db[,i]-min(exp.db[,i], na.rm=T))/(max(exp.db[,i], na.rm=T)-min(exp.db[,i], na.rm=T))
    exp.col <- union(exp.col, exp.db[,i])  
  }
  #exp.col <- as.vector(unique(exp.db))
  names(exp.col) = 1:length(exp.col)
  exp.col = data.frame("exp" = exp.col, "col" ="#000000" )
  exp.col = exp.col[order(exp.col[,"exp"]),]
  exp.col[,"col"] = colorRampPalette(c(blue, "white",red))(nrow(exp.col))
  quantile.range <- quantile(exp.db, probs = seq(0, 1, 0.01), na.rm=T)
  palette.breaks <- seq(quantile.range["0%"], quantile.range["100%"], 0.1)
  color.palette  <- colorRampPalette(c(blue, "white", red))(length(palette.breaks) - 1)
  
  pdf(file = file.path(path.diagrams, sprintf("%s_isofoms_union.pdf", file.name)), height=9, width=17) 
  exp.db.iso <- exp.db[,-ncol(exp.db), drop=FALSE]
  colnames(exp.db.iso)[which(colnames(exp.db.iso) == best.isoform)] <- sprintf("%s***", best.isoform)
  xx <- which(apply(exp.db.iso, MARGIN=2, function(x){which(all(is.na(x)))}) == 1)
  if(length(xx) > 0) {exp.db.iso <- exp.db.iso[, -xx]}
  if(ncol(exp.db.iso) == 1) {
    par(mar=c(9, 2, 5, 15))
    par(oma=c(2,2,2,2))
#    image(exp.db.iso, col=exp.col[,"col"], axes = FALSE)
    image(exp.db.iso, col=color.palette, axes = FALSE)
    grid(nx = nrow(exp.db.iso), ny = ncol(exp.db.iso), lty = 1)
    axis(4,at = 0, labels=colnames(exp.db.iso), las=2, cex.axis = 2, tick = FALSE)
  }else {
    par(oma=c(0,0,0,10))
    if(cluster){
      # hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=T, dendrogram="none", col=exp.col[,"col"], scale="row", trace="none", key=FALSE,
      #                         labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
      #                         #labRow=biomarkers.toPlot,
      #                         labCol=NA)
      hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=T, dendrogram="none", col=color.palette, breaks=palette.breaks, trace="none", key=FALSE,
                              labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
                              #labRow=biomarkers.toPlot,
                              labCol=NA)
    }else{
      # hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=NA, dendrogram="none", col=exp.col[,"col"], scale="row", trace="none", key=FALSE,
      #                         labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
      #                         #labRow=biomarkers.toPlot,
      #                         labCol=NA)
      hv <- gplots::heatmap.2(t(exp.db.iso), Colv=NA, Rowv=NA, dendrogram="none", col=color.palette, breaks=palette.breaks, trace="none", key=FALSE,
                              labRow=colnames(exp.db.iso), cexRow = 0.2 + 1.3/log10(ncol(exp.db.iso)),
                              #labRow=biomarkers.toPlot,
                              labCol=NA)
    }
  }
  dev.off()
  
  pdf(file = file.path(path.diagrams, sprintf("%s_gene_union.pdf", file.name)), height=5, width=12) 
  exp.db.gene <- exp.db[,ncol(exp.db), drop=FALSE]
  par(mar=c(7,1,7,18))
  par(oma=c(2,2,2,2))
#  image(exp.db.gene, col = exp.col[,"col"], axes = FALSE)
  image(exp.db.gene, col=color.palette, axes = FALSE)
  grid(nx = nrow(exp.db.gene), ny = ncol(exp.db.gene), lty = 1)
  axis(4,at = 0, labels=colnames(exp.db.gene), las=2, cex.axis = 2, tick = FALSE)
  dev.off()
  
  my.xlim = c(1,nrow(sensitivity))
  my.ylim = range(as.numeric(sensitivity[,"union"]))
  
  pdf(file=file.path(path.diagrams, sprintf("%s_sensitivity_union.pdf", file.name)), height=5, width=15)    
  par(mar=c(12, 5, 2, 8))
  plot(NA, xlim=my.xlim, ylim=my.ylim, ylab='', xlab='', axes=FALSE)
  axis(1, at=1:nrow(sensitivity), labels=rownames(sensitivity), las=2, cex.axis=1.5, tck=-.05)
  axis(2, at=c(0, .1, .2, .3, .4, .5), labels=c(0, .1, .2, .3, .4, .5), cex.axis=1.5, tck=-.02)
  
  box(lty=1)
  points(1:nrow(sensitivity), as.numeric(sensitivity[ ,"union"]), pch=as.numeric(sensitivity[ , "pch"]), col=sensitivity[ , "col"])
  
  dev.off()
  
  if(cluster & ncol(exp.db.iso) > 1){
    return(gsub("[***]", "", colnames(exp.db.iso)[rev(hv$rowInd)]))
  }
}
fnPlotSensitivity <- function(sensitivity, file.name, type="validation") {
  sensitivity <- sensitivity[which(!is.na(sensitivity))]
  sensitivity <- sensitivity[order(sensitivity)]
  
  my.xlim = c(1,length(sensitivity))
  my.ylim = range(sensitivity)
  
  pdf(file = file.path(path.diagrams, sprintf("%s_sensitivity.pdf", file.name)), height=5, width=ifelse(length(sensitivity) > 20, 13, 13))  
  par(mar=c(12,8,1,1))
#  par(oma=c(2,2,2,2))
  plot(NA, xlim = my.xlim, ylim = my.ylim,ylab='',xlab='', axes = FALSE)
  axis(1,at = 1:length(sensitivity), labels=rep("", length(sensitivity)), tck = -.02)
  text(1:length(sensitivity), par("usr")[3], labels=names(sensitivity), srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=ifelse(length(sensitivity) > 20, 1.2, 1.8))
  axis(2,at = c(0,.1,.2,.3,.4,.5), labels=c(0,.1,.2,.3,.4,.5), cex.axis = 1, tck = -.02)
  
  box(lty = 1)
  if(type == "training"){
    xx <- rep(NA, length(sensitivity))
    names(xx) <- names(sensitivity)
    xx[intersect(names(sensitivity), rownames(gdsc.drug.sensitivity))] <- gdsc.drug.sensitivity[intersect(names(sensitivity), rownames(gdsc.drug.sensitivity)), drug]
    xx[which(!is.na(xx))] <- colorRampPalette(c("blue","light blue","red"))(length(which(!is.na(xx))))
    xx[is.na(xx)] <- "#000000"
    points(1:length(sensitivity), sensitivity, pch=19, cex=1.5, col=xx)
  }else{
    points(1:length(sensitivity), sensitivity, pch=19, cex=1.5, col="#663399")
  }
  
  dev.off()
}
fnPlotEstimates <- function(isoform.models, gene.model, isoforms, file.name, cutoff) {
  isoforms=rev(isoforms)
  mycol3 <- RColorBrewer::brewer.pal(n=4, name="Set3")
  pdf(file.path(path.diagrams,sprintf("%s_estimates.pdf", file.name)), height = 8, width = 2)
  #par(oma=c(0, 0, 0, 0))
  par(mar=c(3, 1, 1, 1))
  tt <- sapply(isoforms , function(x){ifelse(isoform.models[x, "pvalue"] < cutoff, isoform.models[x, "estimate"], 0)})
  if(!is.null(gene.model)){
    tt <- c("gene"=ifelse(gene.model$pvalue < cutoff, gene.model$estimate , 0), tt)
  }
  barplot(tt, horiz = T, col = sapply(tt , function(x){ifelse(x>=0,mycol3[4], mycol3[3])}), yaxt='n', cex.axis=1.2)
  dev.off()
}
fnPlotPvalues <- function(isoform.models, gene.model, isoforms, file.name, cutoff) {
  isoforms=rev(isoforms)
  mycol3 <- RColorBrewer::brewer.pal(n=4, name="Set3")
  pdf(file.path(path.diagrams,sprintf("%s_pvalues.pdf", file.name)), height = 8, width = 2)
  #par(oma=c(0, 0, 0, 0))
  par(mar=c(3, 1, 1, 1))
  #tt <- sapply(isoforms , function(x){ifelse(as.numeric(isoform.models[x, "pvalue"]) < cutoff, -log10(as.numeric(isoform.models[x, "pvalue"])) * sign(as.numeric(isoform.models[x, "estimate"])), 0)})
  tt <- sapply(isoforms , function(x){-log10(as.numeric(isoform.models[x, "pvalue"])) * sign(as.numeric(isoform.models[x, "estimate"]))})
  if(!is.null(gene.model)){
    #tt <- c("gene"= ifelse(as.numeric(gene.model$pvalue) < cutoff, -log10(as.numeric(gene.model$pvalue)) * sign(gene.model$estimate) , 0), tt)
    tt <- c("gene"= -log10(as.numeric(gene.model$pvalue)) * sign(gene.model$estimate), tt)
  }
  barplot(tt, horiz = T, col = sapply(tt , function(x){ifelse(x>=0, mycol3[4], mycol3[3])}), yaxt='n', cex.axis=1.2)
  abline(v=cutoff, lty=2)
  dev.off()
}
fnConfidenceIntervalDifference <- function(model1, model2, x1, x2){
  n <- length(x1)
  r <- cor(x1, x2[names(x1)], use="complete.obs", method="spearman")  
  t.stat <- (model1$estimate - model2$estimate) / sqrt(model1$se^2 + model2$se^2 - 2 * r * model1$se * model2$se)
  diff.ci.p <- pt(q=as.numeric(t.stat), df=n - 1, lower.tail=FALSE)
}
fnCor <- function(drug, gene, xx) {
  tt <- xx
  xx[lower.tri(xx)] <- NA
  diag(xx) <- NA
  quantile.range <- quantile(xx[xx > 0], probs = seq(0, 1, 0.01), na.rm=T)
  palette.breaks <- seq(quantile.range["0%"], quantile.range["100%"], 0.1)
  color.palette  <- colorRampPalette(c("white", red))(length(palette.breaks))
  
  if(any(xx < 0, na.rm=TRUE)){
    neg.quantile.range <- seq(min(xx[xx <= 0], na.rm=T), 0, 0.01)
    neg.palette.breaks <- seq(neg.quantile.range[1], neg.quantile.range[length(neg.quantile.range)], 0.1)
    neg.color.palette  <- colorRampPalette(c(blue, "white"))(length(neg.palette.breaks))
  }else {
    neg.palette.breaks <- NULL
    neg.color.palette <- NULL
  }
  neg.palette.breaks <- c(-1, neg.palette.breaks)
  
  cc <- colorRampPalette(c(blue, "white", red))(length(unique(as.vector(xx[!is.na(xx)]))))
  pdf(file.path(path.diagrams, sprintf("%s_%s_isoforms_gene_corr.pdf", drug, gene)), height=6, width=6)
  par(oma=c(0,0,6,6))
  hv <- gplots::heatmap.2(xx, Colv=NA, Rowv=NA, dendrogram="none", col=c(neg.color.palette, color.palette), breaks = c(neg.palette.breaks, palette.breaks), trace="none", key=FALSE,
                          labRow=colnames(xx)[1:ncol(xx)-1], cexRow = 0.1 + 1/log10(ncol(xx)-1),
                          #labRow=biomarkers.toPlot,
                          labCol=NA)
  if(gene == "TNKS1BP1")
  {
    mtext(colnames(xx)[2:ncol(xx)], side=3, las=2, at=seq(0.29, 0.29 + (ncol(xx) - 2) * 0.1, 0.1))
  }
  if(gene == "TGFA")
  {
    mtext(colnames(xx)[2:ncol(xx)], side=3, las=2, at=seq(0.28, 0.28 + (ncol(xx) - 2) * 0.08, 0.08))
  }
  if(gene == "DUOX1")
  {
    mtext(colnames(xx)[2:ncol(xx)], side=3, las=2, at=seq(0.24, 0.24 + (ncol(xx) - 2) * 0.057, 0.057))
  }
  if(gene == "HNRPDL")
  {
    mtext(colnames(xx)[2:ncol(xx)], side=3, las=2, at=seq(0.29, 0.29 + (ncol(xx) - 2) * 0.1, 0.1))
  }
  if(gene == "CPEB4")
  {
    mtext(colnames(xx)[2:ncol(xx)], side=3, las=2, at=seq(0.24, 0.24 + (ncol(xx) - 2) * 0.057, 0.057))
  }
  dev.off()
  if(gene == "TNKS1BP1")
  {
    tt <- tt[best.isoform, rev(c("ENST00000528882", "ENST00000358252", "ENST00000527207", "ENST00000532437", "ENST00000427750", "ENST00000530920", "ENST00000532273"))]
  }
  if(gene == "TGFA")
  {
    tt <- tt[best.isoform, rev(c("ENST00000295400", "ENST00000445399", "ENST00000474101", "ENST00000418333", "ENST00000450929", "ENST00000444975", "ENST00000394241", "ENST00000460808", "ENST00000419940"))]
  }
  if(gene == "DUOX1")
  {
    tt <- tt[best.isoform, rev(c("ENST00000389037", "ENST00000321429", "ENST00000558322", "ENST00000561220", "ENST00000561166", "ENST00000557893", "ENST00000559716", "ENST00000558991", "ENST00000559219", "ENST00000558446", "ENST00000559221", "ENST00000558744"))]
  }
  if(gene == "HNRPDL")
  {
    tt <- tt[best.isoform, rev(c("ENST00000295470", "ENST00000507721", "ENST00000502762", "ENST00000349655", "ENST00000514511"))]
  }
  if(gene == "CPEB4")
  {
    tt <- tt[best.isoform, rev(c("ENST00000265085", "ENST00000334035", "ENST00000520867", "ENST00000519835", "ENST00000522336", "ENST00000519467", "ENST00000519152", "ENST00000517880", "ENST00000518141", "ENST00000522344"))]
  }
  cc <- rep("gray", length(tt))
  cc[which(names(tt) == best.isoform)] <- "red"
  pdf(file.path(path.diagrams, sprintf("%s_%s_isoforms_gene_corr_bar_plot.pdf", drug, gene)), height=6, width=2)
  barplot(tt, horiz = T, col=cc , yaxt='n', cex.axis=1.2)
  abline(v=0.8, lty=2)
  dev.off()
}
fnExp <- function(drug, gene, exp) {
  tt <- apply(exp, MARGIN=2, function(x){mean(x)})
  tt <- tt/sum(tt, na.rm=T)
  tt[which(is.na(tt))] <- 0
  if(gene == "TNKS1BP1")
  {
    tt <- tt[rev(c("ENST00000528882", "ENST00000358252", "ENST00000527207", "ENST00000532437", "ENST00000427750", "ENST00000530920", "ENST00000532273"))]
  }
  if(gene == "TGFA")
  {
    tt <- tt[rev(c("ENST00000295400", "ENST00000445399", "ENST00000474101", "ENST00000418333", "ENST00000450929", "ENST00000444975", "ENST00000394241", "ENST00000460808", "ENST00000419940"))]
  }
  if(gene == "DUOX1")
  {
    tt <- tt[rev(c("ENST00000389037", "ENST00000321429", "ENST00000558322", "ENST00000561220", "ENST00000561166", "ENST00000557893", "ENST00000559716", "ENST00000558991", "ENST00000559219", "ENST00000558446", "ENST00000559221", "ENST00000558744"))]
  }
  if(gene == "HNRPDL")
  {
    tt <- tt[rev(c("ENST00000295470", "ENST00000507721", "ENST00000502762", "ENST00000349655", "ENST00000514511"))]
  }
  if(gene == "CPEB4")
  {
    tt <- tt[rev(c("ENST00000265085", "ENST00000334035", "ENST00000520867", "ENST00000519835", "ENST00000522336", "ENST00000519467", "ENST00000519152", "ENST00000517880", "ENST00000518141", "ENST00000522344"))]
  }
  cc <- rep("gray", length(tt))
  cc[which(names(tt) == best.isoform)] <- "red"
  pdf(file.path(path.diagrams, sprintf("%s_%s_isoforms_relative_exp.pdf", drug, gene)), height=6, width=2)
  barplot(tt, horiz = T, col=cc , yaxt='n', cex.axis=1.2)
  dev.off()
}

