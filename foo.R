compute.effect.size <- function(effect.size, yi, y.hat){
  if(effect.size == "r.squared"){
    RSS <- sum((yi - y.hat)^2)   
    TSS <- sum((yi - mean(yi))^2)
    return(1-(RSS/TSS))
  }else if(effect.size == "cindex"){
    #return(survcomp::concordance.index(x=-y.hat, surv.time=yi, surv.event=rep(1, length(yi)), na.rm=TRUE, outx=TRUE)[[1]])
    return(Hmisc::rcorr.cens(x=y.hat, S = yi, outx=TRUE)[[1]])  
    #survcomp::concordance.index(x=-y.hat, surv.time=yi, surv.event=rep(1, length(yi)), na.rm=TRUE, outx=TRUE)[[1]]
    #Hmisc::rcorr.cens(x=y.hat, S = yi, outx=TRUE)[[1]]
  }
}
fnOrderDrugs<- function(data,filename, ylab, main){
  mean.data = double()
  for(i in 1:ncol(data)){ mean.data[i] = mean(data[complete.cases(data[,i]),i])}
  data.ordered <- data[,order(mean.data)]
  data.order <- order(mean.data)
  
  pdf(file = filename, height=8, width=16)
  #A <- gsub("drugid_","",colnames(data.ordered))
  mp = boxplot(data.ordered, main = main, xaxt="n", ylab = ylab,cex.axis= 0.7, las=2,col = rainbow(ncol(data)))
  text(1:ncol(data), par("usr")[3], labels = gsub("drugid_","",colnames(data.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7)
  dev.off()
  return(list(order = data.order,ordered = data.ordered))
}
mystacked.barplot.simple <- function(Filename, data, main.label, cex=1.3){
  mycol <- c(RColorBrewer::brewer.pal(n=9, name="Set1"), RColorBrewer::brewer.pal(n=5, name="Set3"), "black")[c(1, 4, 2, 3, 5:15)]
  
  C = mycol[1:nrow(data)]#c(mycol[1],mycol[4],mycol[2])
  #C = c("purple","blue","red")
  
  pdf(file = Filename, height=9, width=14)
  #par(mar=c(8,5,2,2))
  par(mar=c(6.1, 4.1, 4.1, 12.1), xpd=TRUE)
  yr = c(0, 100)
  # mp <- barplot(data, col = C, ylim = yr, ylab = "Percentage", main = main.label,border=NA,axisnames = FALSE, space = 0.5)
  mp <- barplot(data, col = C, ylim = yr, ylab = "Percentage",border=NA,axisnames = FALSE, space = 0.5)
  
  text(mp, par("usr")[3], labels = colnames(data), srt = 50, adj = c(1.1,1.1), xpd = TRUE, cex=cex)
  #axis(1, at=seq(1,length(xx), by=1),labels=FALSE)
  #axis(2)
  legend("topright", inset=c(-0.2,0), legend = rownames(data), fill = C, bty="n")
  dev.off()
  #pdf(file = paste0(Filename,"_legend.pdf"), height=2, width=2)
  #barplot(NA, col = C, ylim = c(0,1) , yaxt = "n", border=NA,axisnames = FALSE, space = 0.5)
  #legend("topright", legend = rownames(data), fill = C, bty="n") 
  #dev.off()
}
mybarplot <- function(Filename, data, barsNo, groupNo, group.labels, ylab.label, legend.lables, main.label, yaxis = c("Regular", "Log"), cex=1.1){
  #compare all
  barplot.colors= RColorBrewer::brewer.pal(n=7, name="Set1")[c(2 , 1, 3:7)]#c("steelblue3","palegreen3", "mediumpurple3","darkorange3", "indianred3")
  Sp <- c(.75,.25,.25,.25,.25)
  yr <- c(0, max(data,na.rm=TRUE)*1.02)
  idx <- NULL
  k <- 0
  for(i in 1:round(ncol(data)/barsNo))
  {
    idx <- c(idx,k + round(barsNo/2))
    k <- k + barsNo
  }
  
  Groups.labels <- matrix(NA, nrow = barsNo, ncol = groupNo)
  Sp <- NULL
  legend.colors <- NULL
  for(i in 1:barsNo)
  {
    Groups.labels[i,] <- ""
    Sp <- c(Sp,.85)
  }
  Groups.labels[round(barsNo/2+1),] <- group.labels
  Sp[1] <- .15
  
  pdf(file = sprintf("%s.pdf", Filename) , height=7, width=14)
  par(mar=c(9,5,5,8), xpd=T)
  par(oma=c(2,2,2,2))
  temp <- NULL
  for(i in 1:barsNo) {
    temp <- rbind(temp, data[seq(i, length(data), by=barsNo)])
  }
  if(yaxis == "Log") {
    if(0 %in% temp){
      temp <- temp + 1
    }
    #mp <- barplot(log10(temp), beside=TRUE, col=barplot.colors[1:barsNo], ylim=log10(c(0, max(temp)) + 1), ylab=ylab.label, space=Sp, border=NA, axes=FALSE)
    mp <- barplot(log10(temp), beside=TRUE, col=barplot.colors[1:barsNo], ylim=log10(c(0, max(temp)) + 1), ylab=ylab.label, space=c(.15,.85), border=NA, axes=FALSE)
    text(mp, par("usr")[3], labels=Groups.labels, srt=50, adj=c(1.1,1.1), xpd=TRUE, cex=cex)
    magicaxis::magaxis(unlog = 'y', side=2, tcl=-.3)
    legend("topright", inset=c(-.05,0), legend=legend.lables, fill=barplot.colors[1:barsNo], bty="n")
  }else {
    
    #    mp <- barplot(temp, beside=TRUE, col=barplot.colors[1:barsNo], ylim=c(0, max(temp)*1.02), ylab=ylab.label, space=Sp, border=NA, axes=FALSE)
    mp <- barplot(temp, beside=TRUE, col=barplot.colors[1:barsNo], ylim=c(0, max(temp)*1.02), ylab=ylab.label, space=c(.15,.85), border=NA, axes=FALSE)
    
    text(mp, par("usr")[3], labels=Groups.labels, srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=1.1)
    magicaxis::magaxis(side=2, tcl=-.3)
    legend("topright", inset=c(-.05,0), legend=legend.lables, fill=barplot.colors[1:barsNo], bty="n")
  }
  
  dev.off()
}
mybarplot.gap <- function(Filename, data, barsNo, group.labels, ylab.label, legend.lables, main.label, breakpoint){
  #compare all
  barplot.colors= c("steelblue3","palegreen3", "mediumpurple3","darkorange3", "indianred3")
  Sp = c(.75,.25,.25,.25,.25)
  Groups.labels <- matrix("", nrow = 1, ncol = ncol(data))
  Groups.labels[1, seq(from=1,to=ncol(Groups.labels),by=barsNo)] = group.labels
  
  Max = max(data,na.rm=TRUE)
  yr = c(0, Max*1.02)
  if(breakpoint == "No")
  {
    lower = c(0,round(Max/10))
  }
  else
  {
    lower = c(0,breakpoint)
  }
  
  upper = c(round(Max/2),(Max*1.02))
  #gr = c(100, round(max(Sum,na.rm=TRUE)/2))
  y_outer=21
  lowspan=c(0,11)
  topspan=c(lowspan[2]+1,21)
  
  pdf(file = Filename, height=7, width=14)
  
  plot(c(0,1),c(0,y_outer),type='n',axes=FALSE,ylab=ylab.label,xlab='', cex=0.7)
  subplot(barplot(t(data), beside = TRUE, col = barplot.colors[1:barsNo], ylim=lower,xpd=FALSE, space=Sp[1:barsNo], names.arg = Groups.labels,las =3, cex.names=0.7,cex.axis=0.7),x=c(0,1),y=lowspan)
  subplot(barplot(t(data), beside = TRUE, col = barplot.colors[1:barsNo], ylim=upper,xpd=FALSE, space=Sp[1:barsNo], main = main.label,cex.axis=0.7), x=c(0,1),y=topspan)
  
  legend("topright", legend = legend.lables, fill = barplot.colors[1:barsNo], bty="n")
  
  
  dev.off()
}
fnComputeAssociateGenes.effect.size <- function(FDR_CutOff=0.01, effect.size_CutOff=0.55, signed) {
  N.Isoforms <- rownames(subset(isoforms_No_List, isoforms.NO > 1))
  One.Isoforms <- rownames(subset(isoforms_No_List, isoforms.NO == 1))
  M <- list()
  index <- NULL
  for(i in 1:length(Models)){
    index <- c(index, apply(combn(Models, i), MARGIN=2, paste, collapse="_"))
  }
  for(i in 1:length(drugs)){
    M[[i]] <- list()
    names(M)[i] <- drugs[i]
    for(j in 1:length(index)){
      M[[i]][[index[j]]] <- matrix(0 , ncol=2, nrow=2, dimnames=list(c("1.isoform","n.isoforms"), c("positive","negative")))
    }
  }
  for (drug in drugs){
    for(ss in names(M[[drug]])){
      xx <- paste(drug, ss, sep ="_")
      if(length(unlist(strsplit(ss, "_"))) == 1){
        M[[drug]][[ss]]["1.isoform", "positive"] <- length(which(FDR_List[One.Isoforms, xx] < FDR_CutOff & 
                                                                   statistics.matrix[One.Isoforms, xx] > effect.size_CutOff &
                                                                   estimate_List[One.Isoforms, xx] > 0))
        M[[drug]][[ss]]["1.isoform", "negative"] <- length(which(FDR_List[One.Isoforms, xx] < FDR_CutOff & 
                                                                   statistics.matrix[One.Isoforms, xx] > effect.size_CutOff &
                                                                   estimate_List[One.Isoforms, xx] < 0))
        M[[drug]][[ss]]["n.isoforms", "positive"] <- length(which(FDR_List[N.Isoforms, xx] < FDR_CutOff & 
                                                                    statistics.matrix[N.Isoforms, xx] > effect.size_CutOff &
                                                                    estimate_List[N.Isoforms, xx] > 0))
        M[[drug]][[ss]]["n.isoforms", "negative"] <- length(which(FDR_List[N.Isoforms, xx] < FDR_CutOff & 
                                                                    statistics.matrix[N.Isoforms, xx] > effect.size_CutOff &
                                                                    estimate_List[N.Isoforms, xx] < 0))
      } else{
        mm <- paste(drug, unlist(strsplit(ss, "_")), sep="_")
        M[[drug]][[ss]]["1.isoform", "positive"] <- length(intersectList(sapply(mm, function(x){which(FDR_List[One.Isoforms, x] < FDR_CutOff &
                                                                                                        statistics.matrix[One.Isoforms, x] > effect.size_CutOff &
                                                                                                        estimate_List[One.Isoforms, x] > 0)})))
        M[[drug]][[ss]]["1.isoform", "negative"] <- length(intersectList(sapply(mm, function(x){which(FDR_List[One.Isoforms, x] < FDR_CutOff &
                                                                                                        statistics.matrix[One.Isoforms, x] > effect.size_CutOff &
                                                                                                        estimate_List[One.Isoforms, x] < 0)})))
        
        M[[drug]][[ss]]["n.isoforms", "positive"] <- length(intersectList(sapply(mm, function(x){which(FDR_List[N.Isoforms, x] < FDR_CutOff &
                                                                                                         statistics.matrix[N.Isoforms, x] > effect.size_CutOff &
                                                                                                         estimate_List[N.Isoforms, x] > 0)})))
        M[[drug]][[ss]]["n.isoforms", "negative"] <- length(intersectList(sapply(mm, function(x){which(FDR_List[N.Isoforms, x] < FDR_CutOff &
                                                                                                         statistics.matrix[N.Isoforms, x] > effect.size_CutOff &
                                                                                                         estimate_List[N.Isoforms, x] < 0)})))
      }
    }
  }
  return(M)
}
associations.all.drugs <-function(model.rank, annot.ensembl.all.genes) {
  rr <- NULL
  col.names <- c(".estimate",".pvalue",paste0(".",adjustment.method))
  
  exp <- expand.grid(".pvalue", Models.names)
  ranks <- unlist(mapply(paste0, exp[,2], exp[,1], SIMPLIFY = F))
  
  col.names <- c(".estimate",".pvalue",paste0(".",adjustment.method))
  exp <- expand.grid(col.names, Models.names)
  models.col.names <- unlist(mapply(paste0, exp[,2], exp[,1], SIMPLIFY = F))
  
  rank <- ranks[which(Models == model.rank)]
  #map.Model <- c(NA,2,3,4,4)
  for(i in 1:length(drugs))
  {
    table_write <- NULL
    drug <- colnames(ccle.drug.sensitivity.ordered)[i]
    for(j in 1:length(Models))
    {
      index.fdr <- paste(gsub("drugid_","",drug),Models[j],sep ="_")
      table_write <- cbind(table_write, estimate_List[,index.fdr], Pvalues.Numinals[,index.fdr], FDR_List[,index.fdr])
    }
    colnames(table_write) <- models.col.names
    table_write_dt <- data.frame(table_write, stringsAsFactors = FALSE)
    table_write_dt[,"symbol"] <- annot.ensembl.all.genes[rownames(table_write_dt), "Symbol"]
    Genedf0 <-data.frame(ccle.primary.drug.sensitivity[,i])
    Genedf0[,"Tissue"] <- factor(ccle.primary.tissuetype[,"tissue.type"], ordered =FALSE)
    table_write_dt[,"n"] <- nrow(Genedf0[complete.cases(Genedf0),])
    table_write_dt[,"isoform"] <- best.isoforms.matrix[,gsub("drugid_","",drug)]
    table_write_dt[,"isoform.exp"] <- mean.ccle.isoforms.fpkm[best.isoforms.matrix[,gsub("drugid_","",drug)]]
    table_write_dt[,"isoform.no"] <-  isoforms_No_List[rownames(table_write_dt), 1]
    table_write_dt[,sprintf("gene.%s",effect.size)] <- statistics.matrix[,paste(gsub("drugid_","",drug),"M2",sep ="_")]
    table_write_dt[,sprintf("isoform.%s",effect.size)] <- statistics.matrix[,paste(gsub("drugid_","",drug),"M3B",sep ="_")]
    
    table_write_dt <- cbind(table_write_dt$symbol, table_write_dt$n, table_write_dt[,models.col.names], table_write_dt$isoform, table_write_dt$isoform.exp, table_write_dt[, sprintf("gene.%s",effect.size)], table_write_dt[,sprintf("isoform.%s",effect.size)])
    colnames(table_write_dt) <- c("symbol","n",models.col.names,"isoform","isoform.exp",sprintf("gene.%s",effect.size),sprintf("isoform.%s",effect.size))
    ### Order based on pvalues of best isoform method
    table_write_dt <- table_write_dt[order(table_write_dt[,rank]),]
    
    rr <- c(rr, list(table_write_dt))
  }
  names(rr) <- gsub("drugid_", "", colnames(ccle.drug.sensitivity.ordered))
  return(rr)
  
}
fnTop.significant.biomarkers <- function(associations, cut_off=.01, BioNo=50, rank.type=c("pvalue", "pvalue.adj")){
  rr <- list()
  
  for(i in 1:length(associations))
  {
    significant_table <- data.frame(matrix(NA, ncol=10, nrow=nrow(associations[[i]])))
    colnames(significant_table) <- c("symbol", "biomarker.id", "type", "specificity", "estimate", effect.size, adjustment.method, "pvalue", "rank", "delta.rank")
    rownames(significant_table) <- rownames(associations[[i]])
    significant_table$symbol <- associations[[i]]$symbol
    significant_table$biomarker.id <- associations[[i]]$symbol
    significant_table$type <- "gene"
    if(rank.type == "pvalue.adj")
    {
      significant_table[rownames(subset(associations[[i]], associations[[i]][, paste0("Genes.", adjustment.method)] < cut_off & associations[[i]][, paste0("Isoforms.", adjustment.method)] >= cut_off)), "specificity"] <- "gene.specific"
      significant_table[rownames(subset(associations[[i]], associations[[i]][, paste0("Genes.", adjustment.method)] < cut_off & associations[[i]][, paste0("Isoforms.", adjustment.method)] < cut_off)), "specificity"] <- "common"
    }else{
      significant_table[rownames(subset(associations[[i]], associations[[i]]$Genes.pvalue < cut_off & associations[[i]]$Isoforms.pvalue >= cut_off)), "specificity"] <- "gene.specific"
      significant_table[rownames(subset(associations[[i]], associations[[i]]$Genes.pvalue < cut_off & associations[[i]]$Isoforms.pvalue < cut_off)), "specificity"] <- "common"
    }
    
    
    significant_table$estimate <- associations[[i]]$Genes.estimate
    significant_table[,effect.size] <- associations[[i]][,sprintf("gene.%s",effect.size)]
    significant_table[, adjustment.method] <- associations[[i]][, paste0("Genes.", adjustment.method)]
    significant_table$pvalue <- associations[[i]]$Genes.pvalue
    
    significant_table$gene.id <- rownames(associations[[i]])
    significant_table$transcript.id <- associations[[i]]$isoform
    
    significant_table$isoforms.no <- isoforms_No_List[rownames(significant_table),1]
    
    temp <- data.frame(matrix(NA, ncol=10, nrow=nrow(associations[[i]])))
    colnames(temp) <- c("symbol","biomarker.id","type","specificity","estimate",effect.size,adjustment.method, "pvalue","rank","delta.rank")
    rownames(temp) <- rownames(associations[[i]])
    temp$symbol <- associations[[i]]$symbol
    temp$biomarker.id <- associations[[i]]$isoform
    temp$type <- "isoform"
    if(rank.type == "pvalue.adj")
    {
      temp[rownames(subset(associations[[i]], associations[[i]][, paste0("Isoforms.", adjustment.method)] < cut_off & associations[[i]][, paste0("Genes.", adjustment.method)] >= cut_off)), "specificity"] <- "isoform.specific"
      temp[rownames(subset(associations[[i]], associations[[i]][, paste0("Isoforms.", adjustment.method)] < cut_off & associations[[i]][, paste0("Genes.", adjustment.method)] < cut_off)), "specificity"] <- "common"
    }else{
      temp[rownames(subset(associations[[i]], associations[[i]]$Isoforms.pvalue < cut_off & associations[[i]]$Genes.pvalue >= cut_off)), "specificity"] <- "isoform.specific"
      temp[rownames(subset(associations[[i]], associations[[i]]$Isoforms.pvalue < cut_off & associations[[i]]$Genes.pvalue < cut_off)), "specificity"] <- "common"
    }
    
    temp$estimate <- associations[[i]]$Isoforms.estimate
    temp[,effect.size] <- associations[[i]][,sprintf("isoform.%s",effect.size)]
    temp[, adjustment.method] <- associations[[i]][, paste0("Isoforms.", adjustment.method)]
    temp$pvalue <- associations[[i]]$Isoforms.pvalue
    
    temp$isoforms.no <- isoforms_No_List[rownames(temp),1]
    temp$gene.id <- rownames(associations[[i]])
    temp$transcript.id <- associations[[i]]$isoform
    
    
    significant_table <- rbind(significant_table,temp)
    rownames(significant_table) <- 1:nrow(significant_table)
    #significant_table <- significant_table[order(significant_table[,adjustment.method], ifelse(abs(significant_table[,"estimate"])!=0, 1/abs(significant_table[,"estimate"]),100000)),]
    if(rank.type == "pvalue.adj")
    { 
      significant_table <- significant_table[order(significant_table[,adjustment.method]),]
    }else{
      significant_table <- significant_table[order(significant_table[,"pvalue"]),]
    }
    significant_table$rank <- 1:nrow(significant_table)    
    significant_table$delta.rank <- 0
    if(BioNo != "all")
    {
      for ( j in 1: BioNo)
      {
        temp <- subset(significant_table, significant_table$symbol == significant_table[j,"symbol"])
        significant_table[j,"delta.rank"] <- abs(temp[1,"rank"] - temp[2,"rank"])
      }
      if(rank.type == "pvalue.adj")
      { 
        significant_table <- significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE]
      }
      last <- ifelse(nrow(significant_table) < BioNo,  nrow(significant_table), BioNo)
    }else{
      significant.NO <- nrow(significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE])
      
      for ( j in 1: significant.NO)
      {
        temp <- subset(significant_table, significant_table$symbol == significant_table[j,"symbol"])
        significant_table[j,"delta.rank"] <- abs(temp[1,"rank"] - temp[2,"rank"])
      }
      if(rank.type == "pvalue.adj")
      { 
        significant_table <-  significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE]
      }
      
      last <- nrow(significant_table)
    }
    
    if(rank.type == "pvalue.adj")
    { 
      significant_table <-  significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE]
    }
    last <- ifelse(nrow(significant_table) < BioNo,  nrow(significant_table), BioNo)
    rr[[names(associations)[i]]] <- significant_table[1:last,]
  }
  return(rr)
}
fnidentify.tissue.specific.biomarkers <- function(biomarkers, boot=FALSE) {
  
  rr <- list()
  for(i in 1:length(biomarkers))
  {
    significant_table <- biomarkers[[i]]
    drug <- names(biomarkers)[i]
    Drugs_ToCheck <- drug
    if(!all(is.na(significant_table)))
    {
      for(j in 1:nrow(significant_table))
      {
        tissue.pvalue <- 1
        tissue.effect.size <- -100
        tissue.effect.size.boot <- -100
        gene.id <- as.character(significant_table[j, "gene.id"])
        transcript.id <- as.character(significant_table[j, "transcript.id"])
        switch(training.type, "CCLE_GDSC"= {
          if (significant_table[j, "type"] == "gene"){
            M0 <- list(ccle = fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle"), 
                       gdsc = fnCreateNullModel(drug = Drugs_ToCheck, assay = "gdsc"))
            weight <- list(ccle = M0$ccle$n / (M0$ccle$n + M0$gdsc$n), gdsc = M0$gdsc$n / (M0$ccle$n + M0$gdsc$n))
            
            M2 <- list(ccle = fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0$ccle, data = ccle.genes.fpkm[,gene.id]), 
                       gdsc = fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0$gdsc, data = ccle.genes.fpkm[, gene.id]))
            
            if(!is.null(M2$ccle) && !is.null(M2$gdsc)){if(!is.na(M2$ccle$coefficient) && !is.na(M2$gdsc$coefficient)){if((sign(M2$ccle$coefficient) != sign(M2$gdsc$coefficient)) || M2$ccle$coefficient == 0 || M2$gdsc$coefficient == 0){M2$ccle <- NULL; M2$gdsc <- NULL}}}
            if(!is.null(M2$ccle) && !is.null(M2$gdsc))
            {
              ccle.effect.size <- compute.effect.size(effect.size, yi=M2$ccle$dataset$drug, y.hat=fitted(M2$ccle$model))
              gdsc.effect.size <- compute.effect.size(effect.size, yi=M2$gdsc$dataset$drug, y.hat=fitted(M2$gdsc$model))
              
              if(nrow(summary(M2$ccle$model)$coefficients) == 2){ ccle.pvalue <- summary(M2$ccle$model)$coefficients[2, 4]} else{ ccle.pvalue <- 1}
              if(nrow(summary(M2$gdsc$model)$coefficients) == 2){ gdsc.pvalue <- summary(M2$gdsc$model)$coefficients[2, 4]} else{ gdsc.pvalue <- 1}
              
              tissue.effect.size <- ifelse(!is.na(ccle.effect.size) & !is.na(gdsc.effect.size), ccle.effect.size * weight$ccle + gdsc.effect.size * weight$gdsc, -100)
              tissue.pvalue <- ccle.pvalue * weight$ccle + gdsc.pvalue * weight$gdsc
              
              if(boot)
              {
                boot.models <- list("M2" = list(data = M2$ccle$dataset, formula = M2$ccle$formula, object = M2$ccle$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
                
                boot.models <- list("M2" = list(data = M2$gdsc$dataset, formula = M2$gdsc$formula, object = M2$gdsc$model))
                M_res.gdsc <- fnboot(models = boot.models, R = 100)
                
                tissue.effect.size.boot <- median(M_res.ccle[["M2"]]) * weight$ccle + median(M_res.gdsc[["M2"]]) * weight$gdsc
              }
            }
          }
          if (significant_table[j, "type"] == "isoform"){
            M0 <- list(ccle = fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle"), 
                       gdsc = fnCreateNullModel(drug = Drugs_ToCheck, assay = "gdsc"))
            weight <- list(ccle = M0$ccle$n / (M0$ccle$n + M0$gdsc$n), gdsc = M0$gdsc$n / (M0$ccle$n + M0$gdsc$n))
            
            M3B <- list(ccle = fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0$ccle, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id)), 
                        gdsc = fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0$gdsc, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id)))
            if(!is.null(M3B$ccle) && !is.null(M3B$gdsc)){if(!is.na(M3B$ccle$coefficient) && !is.na(M3B$gdsc$coefficient)){if((sign(M3B$ccle$coefficient) != sign(M3B$gdsc$coefficient)) || M3B$ccle$coefficient == 0 || M3B$gdsc$coefficient == 0){M3B$ccle <- NULL; M3B$gdsc <- NULL}}}
            if(!is.null(M3B$ccle) && !is.null(M3B$gdsc))
            {
              ccle.effect.size <- compute.effect.size(effect.size, yi=M3B$ccle$dataset$drug, y.hat=fitted(M3B$ccle$model))
              gdsc.effect.size <- compute.effect.size(effect.size, yi=M3B$gdsc$dataset$drug, y.hat=fitted(M3B$gdsc$model))
              
              if(nrow(summary(M3B$ccle$model)$coefficients) == 2){ ccle.pvalue <- summary(M3B$ccle$model)$coefficients[2, 4]} else{ ccle.pvalue <- 1}
              if(nrow(summary(M3B$gdsc$model)$coefficients) == 2){ gdsc.pvalue <- summary(M3B$gdsc$model)$coefficients[2, 4]} else{ gdsc.pvalue <- 1}
              
              tissue.effect.size <- ifelse(!is.na(ccle.effect.size) & !is.na(gdsc.effect.size), ccle.effect.size * weight$ccle + gdsc.effect.size * weight$gdsc, -100)
              tissue.pvalue <- ccle.pvalue * weight$ccle + gdsc.pvalue * weight$gdsc
              
              if(boot)
              {
                
                boot.models <- list("M3B" = list(data = M3B$ccle$dataset, formula = M3B$ccle$formula, object = M3B$ccle$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
                
                boot.models <- list("M3B" = list(data = M3B$gdsc$dataset, formula = M3B$gdsc$formula, object = M3B$gdsc$model))
                M_res.gdsc <- fnboot(models = boot.models, R = 100)
                
                tissue.effect.size.boot <- median(M_res.ccle[["M3B"]]) * weight$ccle + median(M_res.gdsc[["M3B"]]) * weight$gdsc
              }
            }
          }
        }, "CCLE" = {
          
          if (significant_table[j, "type"] == "gene"){
            M0 <- fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle")
            M2 <- fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0, data = ccle.genes.fpkm[,gene.id])
            if(!is.null(M2))
            {
              tissue.effect.size <- compute.effect.size(effect.size, M2$dataset$drug, fitted(M2$model))
              tissue.pvalue <- summary(M2$model)$coefficients[2, 4]
              if(nrow(summary(M2$model)$coefficients) == 2){ tissue.pvalue <- summary(M2$model)$coefficients[2, 4]} else{ tissue.pvalue <- 1}
              
              if(boot)
              {
                boot.models <- list("M2" = list(data = M2$dataset, formula = M2$formula, object = M2$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
                tissue.effect.size.boot <- median(M_res.ccle[["M2"]])
              }
              
            }
          }
          if (significant_table[j, "type"] == "isoform"){
            M0 <- fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle")
            M3B <- fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id))
            if(!is.null(M3B))
            {
              tissue.effect.size <- compute.effect.size(effect.size, M3B$dataset$drug, fitted(M3B$model))
              if(nrow(summary(M3B$model)$coefficients) == 2){ tissue.pvalue <- summary(M3B$model)$coefficients[2, 4]} else{ tissue.pvalue <- 1}
              
              if(boot)
              {
                boot.models <- list("M3B" = list(data = M3B$dataset, formula = M3B$formula, object = M3B$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)  
                
                tissue.effect.size.boot <- median(M_res.ccle[["M3B"]])
              }
            }
          }
        })
        significant_table[j, tissue] <- tissue.effect.size
        significant_table[j, paste0(tissue,"_pvalue")] <- tissue.pvalue
        if(boot){significant_table[j, paste0(tissue,"_boot")] <- tissue.effect.size.boot}
      }
    }
    rr[[drug]] <- significant_table
  }
  return(rr)
}
fnPercentageBiomarkersType <- function(associations) {
  percentage.biomarkers.type <- matrix(0, ncol=length(gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered))), nrow = 3)
  colnames(percentage.biomarkers.type) <- gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered))
  rownames(percentage.biomarkers.type) <- c("isoform.specific","common","gene.specific")
  for( i in gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)))#names(associations))
  {
    type.table <- table(associations[[i]][,"specificity"])
    if(!is.na(type.table["common"]))
    {
      type.table["common"] <- type.table["common"]/2
    }
    percentage.biomarkers.type["gene.specific",i] <- ifelse(is.na(round((type.table["gene.specific"]/sum(type.table))*100,digit=0)), 0, round((type.table["gene.specific"]/sum(type.table))*100,digit=0))
    percentage.biomarkers.type["isoform.specific",i] <- ifelse(is.na(round((type.table["isoform.specific"]/sum(type.table))*100,digit=0)), 0, round((type.table["isoform.specific"]/sum(type.table))*100,digit=0))
    percentage.biomarkers.type["common",i] <- ifelse(is.na(round((type.table["common"]/sum(type.table))*100,digit=0)), 0, round((type.table["common"]/sum(type.table))*100,digit=0))
    if(sum(percentage.biomarkers.type[,i]) != 0) {
      if(sum(percentage.biomarkers.type[,i]) < 100){percentage.biomarkers.type["common",i] <- percentage.biomarkers.type["common",i] + 100 - sum(percentage.biomarkers.type[,i])}
      if(sum(percentage.biomarkers.type[,i]) > 100){percentage.biomarkers.type["common",i] <- 100 - sum(percentage.biomarkers.type[c(1,3),i])}
    }
  }
  return(percentage.biomarkers.type)
}
Check.KnownAssociations <- function(associations) {
  #FirtsLetter <-  tolower(unlist(strsplit(Models.names[which(Models==rank.model)],split=""))[1])
  known.associations <- read.csv(file = file.path(path.data,"KnownAssociations.csv"),stringsAsFactors = FALSE)
  for (i in 1:nrow(known.associations))
  {
    if(known.associations[i,"drug"] %in% colnames(ccle.drug.sensitivity))
    {
      index <- which(associations[[known.associations[i,"drug"]]]$symbol == known.associations[i,"gene"])[1]
      known.associations[i,"rank"] <- index
      known.associations[i,adjustment.method] <- associations[[known.associations[i,"drug"]]][index, paste0("Genes.", adjustment.method)]
    }
  }
  write.csv(known.associations, file = file.path(path.diagrams,"KnownAssociations.csv"))
}
barplot.models <- function(model, isoforms_No=c("all", "1.isoform", "n.isoforms"), signed=c("all", "positive", "negative"), prototype, main.title, breakpoint,  yaxis=c("Regular", "Log"), cex=1.1) {
  isoforms_No <- match.arg(isoforms_No)
  signed <- match.arg(signed)
  barplot.matrix <- matrix(NA, nrow=1 , ncol=(length(drugs) * length(prototype)))
  models.drugs.names <- expand.grid(prototype, colnames(ccle.drug.sensitivity.ordered))
  colnames(barplot.matrix) <- paste(models.drugs.names[,2], models.drugs.names[,1], sep ="_")
  for(i in 1:length(drugs))
  {
    for(j in 1:length(prototype))
    {
      matrix.index <- paste(colnames(ccle.drug.sensitivity.ordered)[i], prototype[j], sep ="_")
      if(isoforms_No=="all" & signed=="all"){
        barplot.matrix[1, matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]])
      }else if(isoforms_No=="all"){
        barplot.matrix[1, matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]][, signed])
      }else if(signed=="all"){
        barplot.matrix[1, matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]][isoforms_No, ])
      }else{
        barplot.matrix[1, matrix.index] <- model[[DOI[i]]][[prototype[j]]][isoforms_No, signed]
      }
    }
  }
  
  Glabels <- colnames(ccle.drug.sensitivity.ordered)
  if(breakpoint == "Regular")
  {
    File <- file.path(path.diagrams, sprintf("BarPlot_%s-isoforms-%s-sign_%s", isoforms_No, signed, str_replace_all(main.title, "[^[:alnum:]]","")))    
    mybarplot(Filename=File, data=barplot.matrix, barsNo=length(prototype), groupNo=length(drugs), group.labels=Glabels, ylab.label="Number of Associated Genes", legend.lables=Models.names[prototype], main.label=main.title, yaxis=yaxis, cex=cex) 
  }else{
    File <- file.path(path.diagrams, sprintf("Gapped_BarPlot_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mybarplot.gap(Filename=File, data=barplot.matrix, barsNo=length(prototype), group.labels=Glabels, ylab.label="Number of Associated Genes", legend.lables=Models.names[prototype], main.label=main.title, breakpoint) 
  }
}
fnWilcox <- function(model, signed) {
  W <- matrix(NA,nrow=length(Models), ncol=length(Models), byrow=FALSE)
  rownames(W) <- Models
  colnames(W) <- Models
  
  if(signed)
  {
    test.matrix <- matrix(NA, ncol=length(Models) ,nrow=length(drugs) * 2)
    colnames(test.matrix) <- Models
    for(i in 1:length(Models))
    {
      k <- 1
      for (j in 1:length(drugs))
      {
        test.matrix[k,Models[i]] <- sum(model[[j]][[Models[i]]][,"positive"])
        test.matrix[k+1,Models[i]] <- sum(model[[j]][[Models[i]]][,"negative"])
        k <- k + 2
      }
    }
  }else{
    test.matrix <- matrix(NA, ncol=length(Models) ,nrow=length(drugs))
    colnames(test.matrix) <- Models
    for(i in 1:length(Models))
    {
      for (j in 1:length(drugs))
      {
        test.matrix[j,Models[i]] <- sum(model[[j]][[Models[i]]][,"positive"])
      }
    }
    
  }
  
  ##Wilcoxon test
  for(j in 1:length(Models))
  {
    for(k in 1:length(Models))
    {
      if(j != k)
      {
        #W[j,k] <- t.test((test.matrix[,j])^2 , (test.matrix[,k])^2, paired = TRUE, alternative = "less")$p.value
        W[j,k] <- wilcox.test((test.matrix[,j])^2 , (test.matrix[,k])^2, paired = TRUE, alternative = "less")$p.value
      }
      else
      {
        W[j,k] <- 1
      }
    }
  }
  return(list(comparison = W))
  
}
myScatterPlot <- function(Name, x, y, method=c("plain", "transparent", "smooth"), legend.label, transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...){
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  
  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             pdf(file = Name, height=7, width=7)
             plot(x=x, y=y, col=col, pch=pch, ...)
             abline(0, 1, lty=2, col="gray")
             if(missing(legend.label)){
               legend("topright", legend=sprintf("r=%.1e", cor(as.numeric(x),as.numeric(y), method="spearman")), cex=1.5, bty="n")
             }else{
               legend("topright", legend=legend.label, cex=1.5, bty="n")
             }
             dev.off()
           },
           "transparent"={
             pdf(file = Name, height=7, width=7)
             myrgb <- grDevices::col2rgb(col, alpha=FALSE) / 255
             plot(x=x, y=y, pch=pch, col=rgb(red=myrgb[1], green=myrgb[2], blue=myrgb[3], alpha=transparency, maxColorValue=1), ...)
             dev.off()
           },
           "smooth"={
             pdf(file = Name, height=7, width=7)
             smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
             #abline(0,max(y)/max(x),col="red")
             abline(lm(y~x), col = "red")
             #lines(x=c(0,0), y=c(max(x),max(y)), lty="solid", col="red")
             dev.off()
           }
    )
  }
}

