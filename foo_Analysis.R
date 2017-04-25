## return all the isoforms of GeneId according to ensembl.map.genes.isoforms annotation provided to the function
fnIsoforms_of_Gene <- function(Gene_Map = ensembl.map.genes.isoforms, GeneId) {
  return(data.frame(strsplit(Gene_Map[,as.character(GeneId)],",")))
}
## plot a bar plot of all modeles in which each model's bar height is showing the number of 
## significant biomarkers caught by that model
#breakpoint c("Regular" ,"No", any number) 
barplot.models <- function(model, isoforms_No = c("1.isoform","n.isoforms","all"), prototype, main.title, breakpoint,  yaxis = c("Regular", "Log"), cex=1.1) {
  models.NO <- length(prototype)
  barplot.matrix <- matrix(NA, nrow=1 , ncol=(drugs_No * models.NO))
  models.drugs.names <- expand.grid(prototype,gsub("drugid_", "", colnames(ccle.drug.sensitivity.ordered)))
  colnames(barplot.matrix) <- paste(models.drugs.names[,2], models.drugs.names[,1], sep ="_")
  for(i in 1:drugs_No)
  {
    for(j in 1:models.NO)
    {
      matrix.index <- paste(gsub("drugid_", "", colnames(ccle.drug.sensitivity.ordered)[i]), prototype[j], sep ="_")
      if(isoforms_No != "all")
      {
        barplot.matrix[1,matrix.index] <- model[[DOI[i]]][[prototype[j]]][isoforms_No, "positive"] + model[[DOI[i]]][[prototype[j]]][isoforms_No, "negative"]
      }else{
        barplot.matrix[1,matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]][, "positive"]) + sum(model[[DOI[i]]][[prototype[j]]][, "negative"])
      }
    }
  }
  
  Glabels <- gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered))
  if(breakpoint == "Regular")
  {
    File <- file.path(path.diagrams, sprintf("BarPlot_%s_%s.pdf", isoforms_No, str_replace_all(main.title, "[^[:alnum:]]","")))    
    mybarplot(Filename=File, data=barplot.matrix, barsNo=models.NO, groupNo=drugs_No, group.labels=Glabels, ylab.label="Number of Associated Genes", legend.lables=Models.names[prototype], main.label=main.title, yaxis=yaxis, cex=cex) 
  }else{
    File <- file.path(path.diagrams, sprintf("Gapped_BarPlot_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mybarplot.gap(Filename=File, data=barplot.matrix, barsNo=models.NO, group.labels=Glabels,ylab.label="Number of Associated Genes", legend.lables=Models.names[prototype], main.label=main.title, breakpoint) 
  }
}
## plot a bar plot of all modeles in which each model's bar is divided to two parts for 
## positive and negative significant biomarkers caught by that model
stacked.barplot.models <- function(model, isoforms_No = c("1.isoform","n.isoforms","all"), prototype, main.title, breakpoint) {
  models.NO <- length(prototype)
  stacked.barplot.matrix <- matrix(NA, nrow=2 , ncol=(models.NO * drugs_No))
  models.drugs.names <- expand.grid(gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)),prototype)
  colnames(stacked.barplot.matrix) <- paste(models.drugs.names[,1],models.drugs.names[,2],sep ="_")
  for(i in 1:drugs_No)
  {
    for(j in 1:models.NO)
    {
      matrix.index <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)[i]),prototype[j],sep ="_")
      if(isoforms_No != "all")
      {
        stacked.barplot.matrix[1,matrix.index] <- model[[DOI[i]]][[prototype[j]]][isoforms_No,"positive"]
        stacked.barplot.matrix[2,matrix.index] <- model[[DOI[i]]][[prototype[j]]][isoforms_No,"negative"]
      }else{
        stacked.barplot.matrix[1,matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]][,"positive"])
        stacked.barplot.matrix[2,matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]][,"negative"])            
      }
    }
  }
  
  Glabels <- gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered))
  Sign <- c(" Pos."," Neg.")
  exp <- expand.grid(Sign,Models.names[prototype])
  Legends <- unlist(mapply(paste0, exp[,2], exp[,1], SIMPLIFY = F))
  if(breakpoint == "Regular")
  {
    File <- file.path(path.diagrams, sprintf("Stacked_BarPlot_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mystacked.barplot.destiny(Filename = File, data = stacked.barplot.matrix, barsNo = models.NO, stackedNo = 2, groupNo = drugs_No, group.labels = Glabels, ylab.label = "Number of Associated Genes", legend.lables= Legends, main.label = main.title)
  }else{
    File <- file.path(path.diagrams, sprintf("Stacked_Gapped_BarPlot_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mystacked.gap.barplot(Filename = File, data = stacked.barplot.matrix, barsNo = models.NO, stackedNo = 2, groupNo = drugs_No, group.labels = Glabels, ylab.label = "Number of Associated Genes", legend.lables= Legends, main.label = main.title, breakpoint)
  }
  
}
## plot a bar plot of all the models in which each density of color of each bar is changing according to the 
## strenght of the stat (cinde or r squared) of the significant biomarkers in that model
#breakpoit Regular or the point in y axis where should be gapped
stacked.barplot.models.ranges <- function(model, isoforms_No = c("1.isoform","n.isoforms","all"), prototype, main.title, breakpoint, yaxis = c("Regular", "Log"), cex=1.1) {
  models.NO <- length(prototype)
  stacked.barplot.matrix <- matrix(NA, nrow=length(ranges.names) , ncol=(models.NO * drugs_No))
  models.drugs.names <- expand.grid(gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)),prototype)
  colnames(stacked.barplot.matrix) <- paste(models.drugs.names[,1],models.drugs.names[,2],sep ="_")
  for(i in 1:drugs_No)
  {
    for(j in 1:models.NO)
    {
      matrix.index <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)[i]),prototype[j],sep ="_")
      if(isoforms_No != "all")
      {
        for(l in 1:length(ranges.names))
        {
          stacked.barplot.matrix[length(ranges.names)-l+1,matrix.index] <- model[[DOI[i]]][[prototype[j]]][isoforms_No,ranges.names[l]]
        }
      }else{
        for(l in 1:length(ranges.names))
        {
          stacked.barplot.matrix[length(ranges.names)-l+1,matrix.index] <- sum(model[[DOI[i]]][[prototype[j]]][,ranges.names[l]])
        }
      }
    }
  }
  
  Glabels <- gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered))
  Legends <- Models.names[prototype]
  if(breakpoint == "Regular")
  {
    File <- file.path(path.diagrams, sprintf("Stacked_BarPlot_ranges_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    ylabels <- mystacked.barplot(Filename = File, data = stacked.barplot.matrix, barsNo = models.NO, stackedNo = length(ranges.names), groupNo = drugs_No, group.labels = Glabels, ylab.label = "Number of Associated Genes", legend.lables= Legends, main.label = main.title, yaxis, cex)
  }else{
    File <- file.path(path.diagrams, sprintf("Stacked_Gapped_BarPlot_ranges_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mystacked.gap.barplot(Filename = File, data = stacked.barplot.matrix, barsNo = models.NO, stackedNo = length(ranges.names), groupNo = drugs_No, group.labels = Glabels, ylab.label = "Number of Associated Genes", legend.lables= Legends, main.label = main.title, breakpoint)
  }
  return (ylabels)
}
fnVennDiagramTriple <- function(model, m3=c("M3","M3B","M3P")) {
  pdf(file = file.path(path.diagrams, paste0(Models.names[which(Models == m3)],"Microarray_RNASeq_Venn.pdf")),height=48, width=32)
  gl <- grid.layout(nrow=6, ncol=4)
  pushViewport(viewport(layout=gl))
  par(new=TRUE, fig=gridFIG())
  
  i <- 1
  k <- 1
  for(j in 1:drugs_No)
  {
    vp.1 <- viewport(layout.pos.col=i, layout.pos.row=k) 
    pushViewport(vp.1)
    if((i ==1) && (k == 1)){  par(new=TRUE, fig=gridFIG())}
    venn.plot <- draw.triple.venn(sum(model[[DOI[j]]][["M1"]]), sum(model[[DOI[j]]][["M2"]]), sum(model[[DOI[j]]][[m3]]), sum(model[[DOI[j]]][["M1_M2"]]), sum(model[[DOI[j]]][[paste0("M2_",m3)]]),sum(model[[DOI[j]]][[paste0("M1_",m3)]]), sum(model[[DOI[j]]][[paste0("M1_M2_",m3)]]), c(Models.names[which(Models=="M1")], Models.names[which(Models=="M2")], paste(Models.names[which(Models==m3)], gsub("drugid_","", colnames(ccle.drug.sensitivity)[DOI[j]]),sep="\n")),fill = c("blue", "red", "green"),lty = "blank",cex = 2,cat.cex = 1,cat.col = c("black", "black", "black"))
    grid.draw(venn.plot)
    popViewport()
    if(i < 4) { i <- i + 1} else{i <- 1; k <- k+1}
    #grid.newpage()
  }
  dev.off()
}
fnVennDiagramPair <- function(model, m3=c("M3","M3B","M3P")) {
  pdf(file = file.path(path.diagrams, paste0(Models.names[which(Models == m3)],"_RNASeq_Pair_Venn.pdf")),height=48, width=32)
  gl <- grid.layout(nrow=6, ncol=4)
  pushViewport(viewport(layout=gl))
  par(new=TRUE, fig=gridFIG())
  
  i <- 1
  k <- 1
  for(j in 1:drugs_No)
  {
    vp.1 <- viewport(layout.pos.col=i, layout.pos.row=k) 
    pushViewport(vp.1)
    if((i ==1) && (k == 1)){  par(new=TRUE, fig=gridFIG())}
    venn.plot <- draw.pairwise.venn(sum(model[[DOI[j]]][["M2"]]), sum(model[[DOI[j]]][[m3]]), sum(model[[DOI[j]]][[paste0("M2_",m3)]]), c(Models.names[which(Models=="M2")], paste(gsub("drugid_","", colnames(ccle.drug.sensitivity)[DOI[j]]),Models.names[which(Models==m3)], sep="\n")),fill = c("blue", "red"),lty = "blank",cex = 2,cat.cex = 1,cat.col = c("black", "black"))
    grid.draw(venn.plot)
    popViewport()
    if(i < 4) { i <- i + 1} else{i <- 1; k <- k+1}
    #grid.newpage()
  }
  dev.off()
}
fnFancystat<- function(FDR, combinations, signed) {
  ## Fancy Statistics
  ##Comparison of RNA-seq and Microarray regarding to the number of genes with FDR less than FDR cutoff
  statistics <- list()
  for (i in 1:length(combinations))
  {
    statistics[[combinations[i]]] <- matrix(NA, ncol=drugs_No, nrow=length(FDR), signed)
  }
  barplot.colors <- c("steelblue3","palegreen3", "mediumpurple3","darkorange3", "indianred3")
  Filename <- file.path(path.diagrams, paste0("Models_Compare.pdf"))
  pdf(file = Filename, height=7, width=14)
  
  for(k in 1:length(FDR))
  {
    fdr_cuttoff <- FDR[k]
    
    barplot.matrix <- matrix(0,ncol=drugs_No,nrow=length(combinations))
    colnames(barplot.matrix) <- gsub("drugid_","",colnames(ccle.drug.sensitivity))
    rownames(barplot.matrix) <- 1:length(combinations)
    for(i in 1:length(combinations))
    {
      m1 <- unlist(strsplit(combinations[i],"_"))[1]
      m2 <- unlist(strsplit(combinations[i],"_"))[2]
      rownames(barplot.matrix)[i] <- paste(Models.names[which(Models==m2)] , Models.names[which(Models==m1)] , sep= " vs ")
      for(j in 1: drugs_No)
      {
        if (signed)
        {
          index.fdr.m1 <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[j]), m1, sep ="_")
          estimate.m1 <- estimate_List[, index.fdr.m1]
          
          fdr.m1.pos <- FDR_List[names(estimate.m1[estimate.m1 >= 0]), index.fdr.m1]
          fdr.m1.pos <- fdr.m1.pos[fdr.m1.pos < fdr_cuttoff]
          
          fdr.m1.neg <- FDR_List[names(estimate.m1[estimate.m1 < 0]), index.fdr.m1]
          fdr.m1.neg <- fdr.m1.neg[fdr.m1.neg < fdr_cuttoff]
          
          index.fdr.m2 <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[j]), m2, sep ="_")
          estimate.m2 <- estimate_List[, index.fdr.m2]
          
          fdr.m2.pos <- FDR_List[names(estimate.m2[estimate.m2 >= 0]), index.fdr.m2]
          fdr.m2.pos <- fdr.m2.pos[fdr.m2.pos < fdr_cuttoff]
          
          fdr.m2.neg <- FDR_List[names(estimate.m1[estimate.m2 < 0]), index.fdr.m2]
          fdr.m2.neg <- fdr.m2.neg[fdr.m2.neg < fdr_cuttoff]
          
          #unique genes in the first model
          MM1 <- length(setdiff(names(fdr.m1.pos),names(fdr.m2.pos))) +  length(setdiff(names(fdr.m1.neg),names(fdr.m2.neg)))
          #unique genes in the second model
          MM2 <- length(setdiff(names(fdr.m2.pos),names(fdr.m1.pos))) +  length(setdiff(names(fdr.m2.neg),names(fdr.m1.neg)))
        }
        else
        {
          index.fdr.m1 <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[j]), m1, sep ="_")
          
          fdr.m1 <- FDR_List[, index.fdr.m1]
          fdr.m1 <- fdr.m1[fdr.m1 < fdr_cuttoff]
          
          index.fdr.m2 <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[j]), m2, sep ="_")
          
          fdr.m2 <- FDR_List[, index.fdr.m2]
          fdr.m2 <- fdr.m2[fdr.m2 < fdr_cuttoff]
          
          #unique genes in the first model
          MM1 <- length(setdiff(names(fdr.m1),names(fdr.m2))) 
          #unique genes in the second model
          MM2 <- length(setdiff(names(fdr.m2),names(fdr.m1)))
        }
        if ((MM1 != 0) || (MM2 != 0))
        {barplot.matrix[i,j] <- ((MM2 - MM1) / (MM1+ MM2))}                                                                              
      }
      statistics[[combinations[i]]][k,] <- barplot.matrix[i,]
    }
    barplot(barplot.matrix, beside=TRUE, col = barplot.colors[1:length(combinations)] , main = sprintf("%s cutoff = %s", adjustment.method, fdr_cuttoff), cex.names = 0.7, las=2)
    legend("topright", legend = rownames(barplot.matrix), fill = barplot.colors[1:length(combinations)], bty="n")
  }
  dev.off()
  return(statistics)  
}
fnCuttOffs<-function(cutoff_statistics, combination) {
  m1 <- unlist(strsplit(combination,"_"))[1]
  m2 <- unlist(strsplit(combination,"_"))[2]
  pdf(file = file.path(path.diagrams, paste0(paste(Models.names[which(Models==m2)] , Models.names[which(Models==m1)] , sep= "_vs_"),"_Cutoff_Statistics.pdf")), height=7, width=14)
  L <- list()
  for(i in 1:drugs_No)
  {
    L[[i]]<- cutoff_statistics[[combination]][,i]
  }
  plot.multi(FDR,L)
  title(paste(Models.names[which(Models==m2)] , Models.names[which(Models==m1)] , sep= " vs "))
  dev.off()
}
fnWilcox <- function(model, signed) {
  W <- matrix(NA,nrow=length(Models), ncol=length(Models), byrow=FALSE)
  rownames(W) <- Models
  colnames(W) <- Models
  
  if(signed)
  {
    test.matrix <- matrix(NA, ncol=length(Models) ,nrow=drugs_No * 2)
    colnames(test.matrix) <- Models
    for(i in 1:length(Models))
    {
      k <- 1
      for (j in 1:drugs_No)
      {
        test.matrix[k,Models[i]] <- sum(model[[j]][[Models[i]]][,"positive"])
        test.matrix[k+1,Models[i]] <- sum(model[[j]][[Models[i]]][,"negative"])
        k <- k + 2
      }
    }
  }else{
    test.matrix <- matrix(NA, ncol=length(Models) ,nrow=drugs_No)
    colnames(test.matrix) <- Models
    for(i in 1:length(Models))
    {
      for (j in 1:drugs_No)
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
##report all the associations in a list of objets where each object is for one drug
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
  for(i in 1:drugs_No)
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
    table_write_dt[,sprintf("gene.%s",stat)] <- statistics.matrix[,paste(gsub("drugid_","",drug),"M2",sep ="_")]
    table_write_dt[,sprintf("isoform.%s",stat)] <- statistics.matrix[,paste(gsub("drugid_","",drug),"M3B",sep ="_")]
    
    table_write_dt <- cbind(table_write_dt$symbol, table_write_dt$n, table_write_dt[,models.col.names], table_write_dt$isoform, table_write_dt$isoform.exp, table_write_dt[, sprintf("gene.%s",stat)], table_write_dt[,sprintf("isoform.%s",stat)])
    colnames(table_write_dt) <- c("symbol","n",models.col.names,"isoform","isoform.exp",sprintf("gene.%s",stat),sprintf("isoform.%s",stat))
    ### Order based on pvalues of best isoform method
    table_write_dt <- table_write_dt[order(table_write_dt[,rank]),]
    
    rr <- c(rr, list(table_write_dt))
  }
  names(rr) <- gsub("drugid_", "", colnames(ccle.drug.sensitivity.ordered))
  return(rr)
  
}
associations.all.drugs.significant <-function(rank, method=c("isoform","gene"), cutoff=.01, R2_cut=.6, exp_cut=1, annot.ensembl.all.genes) {
  rr <- NULL
  col.names = c(".estimate",".pvalue",paste0(".",adjustment.method))
  
  exp <- expand.grid(".pvalue", Models.names)
  ranks <- unlist(mapply(paste0, exp[,2], exp[,1], SIMPLIFY = F))
  
  col.names <- c(".estimate",".pvalue",paste0(".",adjustment.method))
  exp <- expand.grid(col.names, Models.names)
  models.col.names <- unlist(mapply(paste0, exp[,2], exp[,1], SIMPLIFY = F))
  
  #rank <- ranks[which(Models == model.rank)]
  #map.Model = c(NA,2,3,4,4)
  for(temp in 1:drugs_No)
  {
    i <- which(colnames(ccle.drug.sensitivity)==colnames(ccle.drug.sensitivity.ordered)[temp])
    
    table_write <- NULL
    for(j in 1:length(Models))
    {
      index.fdr <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),Models[j],sep ="_")
      table_write <- cbind(table_write, estimate_List[,index.fdr], Pvalues.Numinals[,index.fdr], FDR_List[,index.fdr])
    }
    colnames(table_write) <- models.col.names
    table_write_dt <- data.frame(table_write, stringsAsFactors = FALSE)
    table_write_dt[,"symbol"] <- annot.ensembl.all.genes[rownames(table_write_dt), "Symbol"]
    Genedf0 <-data.frame(ccle.primary.drug.sensitivity[,i])
    Genedf0[,"Tissue"] <- factor(ccle.primary.tissuetype[,"tissue.type"], ordered =FALSE)
    table_write_dt[,"n"] <- nrow(Genedf0[complete.cases(Genedf0),])
    table_write_dt[,"isoform"] <- best.isoforms.matrix[,gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]
    table_write_dt[,"isoform.exp"] <- mean.ccle.isoforms.fpkm[best.isoforms.matrix[,gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]]
    table_write_dt[,sprintf("gene.%s",stat)] <- statistics.matrix[,paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),"M2",sep ="_")]
    table_write_dt[,sprintf("isoform.%s",stat)] <- statistics.matrix[,paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),"M3B",sep ="_")]
    table_write_dt <- cbind(table_write_dt$symbol, table_write_dt$n, table_write_dt[,models.col.names], table_write_dt$isoform, table_write_dt$isoform.exp, table_write_dt[,sprintf("gene.%s",stat)], table_write_dt[,sprintf("isoform.%s",stat)])
    colnames(table_write_dt) <- c("symbol","n",models.col.names,"isoform","isoform.exp",sprintf("gene.%s",stat),sprintf("isoform.%s",stat))
    
    if(method == "isoform")
    {
      table_write_dt <- subset(table_write_dt,  
                                table_write_dt[,sprintf("%s.%s",Models.names[2],adjustment.method)]==1 & 
                                table_write_dt[,sprintf("%s.%s",Models.names[3],adjustment.method)] < cutoff & 
                                table_write_dt[,sprintf("isoform.%s",stat)] > R2_cut & 
                                #&  table_write_dt[,sprintf("%s.estimate",Models.names[3])] > 0 &  
                                table_write_dt[,"isoform.exp"] > exp_cut)
    }else if(method == "gene")
    {
      table_write_dt <- subset(table_write_dt, table_write_dt[,sprintf("%s.%s",Models.names[3],adjustment.method)]== 1  & table_write_dt[,sprintf("%s.%s",Models.names[2],adjustment.method)] < cutoff )
    }
    ### Order based on pvalues of best isoform method
    table_write_dt <- table_write_dt[order(table_write_dt[,rank], decreasing = TRUE),]
    rr <- c(rr, list(table_write_dt))
    #write.table(table_write_dt, file = file.path(path.result, paste0(colnames(ccle.drug.sensitivity)[i],".csv")), sep = "\t", quote = FALSE)
  }
  names(rr) <- gsub("drugid_", "", colnames(ccle.drug.sensitivity.ordered))
  return(rr)
  
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
Check.KnownAssociations.PLOS1 <- function(associations) {
  #FirtsLetter <-  tolower(unlist(strsplit(Models.names[which(Models==rank.model)],split=""))[1])
  known.associations <- read.csv(file = file.path(path.data,"Table_S4.csv"),stringsAsFactors = FALSE)
  for (i in 1:nrow(known.associations))
  {
    if(known.associations[i,"drug"] %in% colnames(ccle.drug.sensitivity))
    {
      if(known.associations[i,"drug"] %in% colnames(ccle.drug.sensitivity))
      {
        index <- which(associations[[known.associations[i,"drug"]]]$symbol == known.associations[i,"gene"])[1]
        known.associations[i,"rank"] <- index
        known.associations[i,adjustment.method] <- associations[[known.associations[i,"drug"]]][index, paste0("Genes.", adjustment.method)]
      }
    }
  }
  write.csv(known.associations, file = file.path(path.diagrams,"KnownAssociations_PLOS1.csv"))
}
fnComputeAssociateGenes <- function(FDR_CutOff = 0.01, signed) {
  N.Isoforms <- rownames(subset(isoforms_No_List,isoforms.NO > 1))
  One.Isoforms <- rownames(subset(isoforms_No_List,isoforms.NO == 1))
  M <- list()
  for(i in 1:drugs_No)
  {
    M[[i]] <- list()
    names(M)[i] <- gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])
    for(j in 1:length(Models))
    {
      M[[i]][[Models[j]]] <- matrix(0 , ncol = 2, nrow = 2)
      colnames(M[[i]][[Models[j]]]) <- c("positive","negative")
      rownames(M[[i]][[Models[j]]]) <- c("1.isoform","n.isoforms") 
      for(k in 1:j)
      {
        if(Models[j] != Models[k])
        {
          index.comb.2 <- paste(Models[k],Models[j],sep="_")
          M[[i]][[index.comb.2]] <- matrix(0 , ncol = 2, nrow = 2)
          colnames(M[[i]][[index.comb.2]]) <- c("positive","negative")
          rownames(M[[i]][[index.comb.2]]) <- c("1.isoform","n.isoforms")            
        }
      }
    }
    for(j in 3:length(Models))
    {
      index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[j],sep="_")
      M[[i]][[index.comb.3]] <- matrix(0 , ncol = 2, nrow = 2)
      colnames(M[[i]][[index.comb.3]]) <- c("positive","negative")
      rownames(M[[i]][[index.comb.3]]) <- c("1.isoform","n.isoforms")            
    }
  }
  index.fdr <- NULL
  if(signed)
  {
    for (i in 1:drugs_No)
    {
      One.Isoforms.k <- list()
      N.Isoforms.k <- list()
      
      for(k in 1:length(Models))
      {
        index.fdr[k] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),Models[k],sep ="_")
        N.Isoforms.common <- N.Isoforms
        
        One.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[One.Isoforms, index.fdr[k]], "sign" = estimate_List[One.Isoforms, index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["1.isoform","positive"] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0))
        M[[i]][[Models[k]]]["1.isoform","negative"] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0))
        
        
        N.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[N.Isoforms.common, index.fdr[k]], "sign" = estimate_List[N.Isoforms.common, index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["n.isoforms","positive"] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0))
        M[[i]][[Models[k]]]["n.isoforms","negative"] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0))
      }
      for(k in 1:length(Models))
      {
        for(k2 in 1:k)
        {
          index.comb.2 <- paste(Models[k2],Models[k],sep="_")
          if(Models[k] != Models[k2])
          {
            M[[i]][[index.comb.2]]["1.isoform","positive"] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0)),
                                                                              rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign > 0))))
            M[[i]][[index.comb.2]]["n.isoforms","positive"] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0)),
                                                                               rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign > 0))))
            
            M[[i]][[index.comb.2]]["1.isoform","negative"] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0)),
                                                                              rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign < 0))))
            M[[i]][[index.comb.2]]["n.isoforms","negative"] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0)),
                                                                               rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign < 0))))
            
          }
        }
      }
      for(k in 3:length(Models))
      {
        index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[k],sep="_")
        
        M[[i]][[index.comb.3]]["1.isoform","positive"] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign > 0)),
                                                                                    rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign > 0))),
                                                                          rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0))))
        M[[i]][[index.comb.3]]["n.isoforms","positive"] <-  length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign > 0)),
                                                                                      rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign > 0))),
                                                                            rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0))))
        
        M[[i]][[index.comb.3]]["1.isoform","negative"] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign < 0)),
                                                                                    rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign < 0))),
                                                                          rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0))))
        M[[i]][[index.comb.3]]["n.isoforms","negative"] <-  length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign < 0)),
                                                                                      rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign < 0))),
                                                                            rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0))))
        
      }
    }#drugs
  }
  else
  {
    for (i in 1:drugs_No)
    {
      One.Isoforms.k <- list()
      N.Isoforms.k <- list()
      for(k in 1:length(Models))
      {
        index.fdr[k] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),Models[k],sep ="_")
        N.Isoforms.common <- N.Isoforms
        One.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[One.Isoforms, index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["1.isoform","positive"] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff))
        
        N.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[N.Isoforms.common, index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["n.isoforms","positive"] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff))
        
      }
      for(k in 1:length(Models))
      {
        for(k2 in 1:k)
        {
          index.comb.2 <- paste(Models[k2],Models[k],sep="_")
          if(Models[k] != Models[k2])
          {
            M[[i]][[index.comb.2]]["1.isoform","positive"] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff)),
                                                                              rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff))))
            M[[i]][[index.comb.2]]["n.isoforms","positive"] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff)),
                                                                               rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff))))
            
          }
        }
      }
      for(k in 3:length(Models))
      {
        index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[k],sep="_")
        
        M[[i]][[index.comb.3]]["1.isoform","positive"] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff )),
                                                                                    rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff ))),
                                                                          rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff ))))
        M[[i]][[index.comb.3]]["n.isoforms","positive"] <-  length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff )),
                                                                                      rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff ))),
                                                                            rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff ))))
        
      }
    }#drugs
  }
  return(M)
}
fnComputeAssociateGenes.stat <- function(FDR_CutOff = 0.01, stat_CutOff = 0.6, signed) {
  N.Isoforms <- rownames(subset(isoforms_No_List,isoforms.NO > 1))
  One.Isoforms <- rownames(subset(isoforms_No_List,isoforms.NO == 1))
  M <- list()
  for(i in 1:drugs_No)
  {
    M[[i]] <- list()
    names(M)[i] <- gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])
    for(j in 1:length(Models))
    {
      M[[i]][[Models[j]]] <- matrix(0 , ncol = 2, nrow = 2)
      colnames(M[[i]][[Models[j]]]) <- c("positive","negative")
      rownames(M[[i]][[Models[j]]]) <- c("1.isoform","n.isoforms") 
      for(k in 1:j)
      {
        if(Models[j] != Models[k])
        {
          index.comb.2 <- paste(Models[k],Models[j],sep="_")
          M[[i]][[index.comb.2]] <- matrix(0 , ncol = 2, nrow = 2)
          colnames(M[[i]][[index.comb.2]]) <- c("positive","negative")
          rownames(M[[i]][[index.comb.2]]) <- c("1.isoform","n.isoforms")            
        }
      }
    }
    for(j in 3:length(Models))
    {
      index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[j],sep="_")
      M[[i]][[index.comb.3]] <- matrix(0 , ncol = 2, nrow = 2)
      colnames(M[[i]][[index.comb.3]]) <- c("positive","negative")
      rownames(M[[i]][[index.comb.3]]) <- c("1.isoform","n.isoforms")            
    }
  }
  index.fdr <- NULL
  if(signed)
  {
    for (i in 1:drugs_No)
    {
      One.Isoforms.k <- list()
      N.Isoforms.k <- list()
      for(k in 1:length(Models))
      {
        index.fdr[k] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),Models[k],sep ="_")
        N.Isoforms.common <- N.Isoforms
        
        One.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[One.Isoforms, index.fdr[k]], "stat" = statistics.matrix[One.Isoforms,index.fdr[k]], "sign" = estimate_List[One.Isoforms, index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["1.isoform","positive"] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff & sign > 0 ))
        M[[i]][[Models[k]]]["1.isoform","negative"] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff & sign < 0 ))
        
        
        N.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[N.Isoforms.common, index.fdr[k]], "stat" = statistics.matrix[N.Isoforms.common,index.fdr[k]],"sign" = estimate_List[N.Isoforms.common, index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["n.isoforms","positive"] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff & sign > 0 ))
        M[[i]][[Models[k]]]["n.isoforms","negative"] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff & sign < 0 ))
      }
      for(k in 1:length(Models))
      { 
        for(k2 in 1:k)
        {
          index.comb.2 <- paste(Models[k2],Models[k],sep="_")
          if(Models[k] != Models[k2])
          {          
            M[[i]][[index.comb.2]]["1.isoform","positive"] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff )),
                                                                              rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign > 0  & stat > stat_CutOff ))))
            M[[i]][[index.comb.2]]["n.isoforms","positive"] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0  & stat > stat_CutOff )),
                                                                               rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign > 0  & stat > stat_CutOff ))))
            
            M[[i]][[index.comb.2]]["1.isoform","negative"] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff )),
                                                                              rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign < 0  & stat > stat_CutOff ))))
            M[[i]][[index.comb.2]]["n.isoforms","negative"] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0  & stat > stat_CutOff )),
                                                                               rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & sign < 0  & stat > stat_CutOff ))))
            
          }
        }
      }
      for(k in 3:length(Models))
      {
        index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[k],sep="_")
        
        M[[i]][[index.comb.3]]["1.isoform","positive"] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff )),
                                                                                    rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff ))),
                                                                          rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff ))))
        M[[i]][[index.comb.3]]["n.isoforms","positive"] <-  length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff )),
                                                                                      rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff ))),
                                                                            rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign > 0 & stat > stat_CutOff ))))
        
        M[[i]][[index.comb.3]]["1.isoform","negative"] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff )),
                                                                                    rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff ))),
                                                                          rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff ))))
        M[[i]][[index.comb.3]]["n.isoforms","negative"] <-  length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff )),
                                                                                      rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff ))),
                                                                            rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & sign < 0 & stat > stat_CutOff ))))
        
      }
    }#drugs
  }
  else
  {
    for (i in 1:drugs_No)
    {
      One.Isoforms.k <- list()
      N.Isoforms.k <- list()
      
      
      for(k in 1:length(Models))
      {
        index.fdr[k] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),Models[k],sep ="_")
        N.Isoforms.common <- N.Isoforms
        
        One.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[One.Isoforms, index.fdr[k]], "stat" = statistics.matrix[One.Isoforms,index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["1.isoform","positive"] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff))
        
        N.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[N.Isoforms.common, index.fdr[k]], "stat" = statistics.matrix[N.Isoforms.common,index.fdr[k]]), stringsAsFactors = FALSE)
        M[[i]][[Models[k]]]["n.isoforms","positive"] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff))
      }
      for(k in 1:length(Models))
      {
        for(k2 in 1:k)
        {
          index.comb.2 <- paste(Models[k2],Models[k],sep="_")
          if(Models[k] != Models[k2])
          {
            M[[i]][[index.comb.2]]["1.isoform","positive"] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff)),
                                                                              rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & stat > stat_CutOff))))
            M[[i]][[index.comb.2]]["n.isoforms","positive"] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff)),
                                                                               rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & stat > stat_CutOff))))
          }
        }
      }
      for(k in 3:length(Models))
      {
        index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[k],sep="_")
        
        M[[i]][[index.comb.3]]["1.isoform","positive"] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & stat > stat_CutOff)),
                                                                                    rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & stat > stat_CutOff))),
                                                                          rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff))))
        M[[i]][[index.comb.3]]["n.isoforms","positive"] <-  length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & stat > stat_CutOff)),
                                                                                      rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & stat > stat_CutOff))),
                                                                            rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > stat_CutOff))))
      }
    }#drugs
  }
  return(M)
}
fnComputeAssociateGenes.stat.range <- function(FDR_CutOff = 0.01, biomarkers.sign = c("positive","negative","all")) {
  switch(biomarkers.sign, "positive"= {sign.condition = "sign > 0"}, "negative"= {sign.condition = "sign < 0"}, "all"= {sign.condition = "sign > -2"})
  
  N.Isoforms <- rownames(subset(isoforms_No_List,isoforms.NO > 1))
  One.Isoforms <- rownames(subset(isoforms_No_List,isoforms.NO == 1))
  M = list()
  for(i in 1:drugs_No)
  {
    M[[i]] <- list()
    names(M)[i] <- gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])
    for(j in 1:length(Models))
    {
      M[[i]][[Models[j]]] <- matrix(0 , ncol = length(ranges.names), nrow = 2)
      colnames(M[[i]][[Models[j]]]) <- ranges.names
      rownames(M[[i]][[Models[j]]]) <- c("1.isoform","n.isoforms") 
      for(k in 1:j)
      {
        if(Models[j] != Models[k])
        {
          index.comb.2 <- paste(Models[k],Models[j],sep="_")
          M[[i]][[index.comb.2]] <- matrix(0 , ncol = length(ranges.names), nrow = 2)
          colnames(M[[i]][[index.comb.2]]) <- ranges.names
          rownames(M[[i]][[index.comb.2]]) <- c("1.isoform","n.isoforms")            
        }
      }
    }
    for(j in 3:length(Models))
    {
      index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[j],sep="_")
      M[[i]][[index.comb.3]] <- matrix(0 , ncol = length(ranges.names), nrow = 2)
      colnames(M[[i]][[index.comb.3]]) <- ranges.names
      rownames(M[[i]][[index.comb.3]]) <- c("1.isoform","n.isoforms")            
    }
  }
  index.fdr <- NULL
  
  for (i in 1:drugs_No)
  {
    One.Isoforms.k <- list();
    N.Isoforms.k <- list();
    for(k in 1:length(Models))
    {
      index.fdr[k] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),Models[k],sep ="_")
      N.Isoforms.common <- N.Isoforms
      
      One.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[One.Isoforms, index.fdr[k]], "stat" = statistics.matrix[One.Isoforms,index.fdr[k]],"sign" = estimate_List[One.Isoforms, index.fdr[k]]), stringsAsFactors = FALSE)
      N.Isoforms.k[[Models[k]]] <- data.frame(cbind("fdr" = FDR_List[N.Isoforms.common, index.fdr[k]], "stat" = statistics.matrix[N.Isoforms.common,index.fdr[k]],"sign" = estimate_List[N.Isoforms.common, index.fdr[k]]), stringsAsFactors = FALSE)
      
      
      for(l in 1:length(ranges.names))
      {  
        M[[i]][[Models[k]]]["1.isoform",ranges.names[l]] <- nrow(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))
        M[[i]][[Models[k]]]["n.isoforms",ranges.names[l]] <- nrow(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))
      }
    }  
    for(k in 1:length(Models))
    {
      for(k2 in 1:k)
      {
        index.comb.2 <- paste(Models[k2],Models[k],sep="_")
        if(Models[k] != Models[k2])
        {
          for(l in 1:length(ranges.names))
          { 
            M[[i]][[index.comb.2]]["1.isoform",ranges.names[l]] <- length(intersect(rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition)))),
                                                                                   rownames(subset (One.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))))
            
            M[[i]][[index.comb.2]]["n.isoforms",ranges.names[l]] <- length(intersect(rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition)))),
                                                                                    rownames(subset (N.Isoforms.k[[Models[k2]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))))
          }
        }
      }
    }
    for(k in 3:length(Models))
    {
      index.comb.3 <- paste(paste(Models[1],Models[2],sep="_"),Models[k],sep="_")
      
      for(l in 1:length(ranges.names))
      { 
        M[[i]][[index.comb.3]]["1.isoform",ranges.names[l]] <- length(intersect(intersect(rownames(subset (One.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition)))),
                                                                                         rownames(subset (One.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))),
                                                                               rownames(subset (One.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))))
        
        M[[i]][[index.comb.3]]["n.isoforms",ranges.names[l]] <- length(intersect(intersect(rownames(subset (N.Isoforms.k[[Models[1]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition)))),
                                                                                          rownames(subset (N.Isoforms.k[[Models[2]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))),
                                                                                rownames(subset (N.Isoforms.k[[Models[k]]] , fdr < FDR_CutOff & stat > ranges[l+1] & stat <= ranges[l] & eval(parse(text=sign.condition))))))
      }
    }
  }#drugs
  
  return(M)
}
fnMarkersPercent <- function(model, m1, m2) {
  barplot.colors <- c("indianred3","steelblue3","palegreen3", "mediumpurple3","darkorange3")
  Filename <- file.path(path.diagrams, paste0("MarkersPercent.pdf"))
  pdf(file = Filename, height=7, width=14)
  barplot.matrix <- NULL
  for(j in 1:drugs_No)
  {
    barplot.matrix[j] <- (sum(model[[DOI[j]]][[m2]]))/(sum(model[[DOI[j]]][[m1]])+sum(model[[DOI[j]]][[m2]])-sum(model[[DOI[j]]][[paste0(paste0(m1,"_"),m2)]]))
  }
  mp <- barplot(barplot.matrix, beside=TRUE,space = .6, ylim = c(0,1), ylab = "percentage covered by isoform model", main = sprintf("Percent of markers covered by %s vs %s",Models.names[which(Models==m2)],Models.names[which(Models==m1)]), 
                col = rainbow(ncol(ccle.drug.sensitivity)),border=NA,xaxt='n',cex.main = .7, cex.lab = .7, cex.axis = .7)
  text(mp, par("usr")[3], labels = gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7)
  
  dev.off()
  
}
fnIsoformVSGeneR2 <- function(FDR_CutOff = 0.01, m1 ,m2) {
  index.fdr <- NULL
  exp <- expand.grid(c(m1,m2,paste(m1,m2,sep="_")),gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)))
  M1 <- matrix(NA, nrow(FDR_List), ncol = 3 * drugs_No)
  rownames(M1) <- rownames(FDR_List)
  colnames(M1) <- paste(exp[,1],exp[,2],sep="_")
  
  M2 <- matrix(NA, nrow(FDR_List), ncol = drugs_No)
  rownames(M2) <- rownames(FDR_List)
  colnames(M2) <- colnames(ccle.drug.sensitivity.ordered)
  
  M12 <- matrix(NA, nrow(FDR_List), ncol = drugs_No)  
  rownames(M12) <- rownames(FDR_List)
  colnames(M12) <- colnames(ccle.drug.sensitivity.ordered)
  
  k <- 1
  for (j in 1:drugs_No)
  {
    i <- DOI[j]
    index.fdr[1] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),m1,sep ="_")
    index.fdr[2] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),m2,sep ="_")
    
    
    T1 <- data.frame(cbind("fdr" = FDR_List[, index.fdr[1]], "stat" = statistics.matrix[,index.fdr[1]]), stringsAsFactors = FALSE)
    T1 <- rownames(subset (T1 , fdr < FDR_CutOff ))
    
    T2 <- data.frame(cbind("fdr" = FDR_List[, index.fdr[2]], "stat" = statistics.matrix[,index.fdr[2]]), stringsAsFactors = FALSE)
    T2 <- rownames(subset (T2 , fdr < FDR_CutOff))
    
    
    T12 <- intersect(T1, T2)
    
    M1[setdiff(T1,T12),k] <- statistics.matrix[setdiff(T1,T12), index.fdr[1]]    
    M1[setdiff(T2,T12),k + 1] <- statistics.matrix[setdiff(T2,T12), index.fdr[2]]    
    M1[T12,k + 2] <- statistics.matrix[T12, index.fdr[2]]    
    k <- k + 3
    #M2[setdiff(names(T2),T12),i] <- statistics.matrix[setdiff(names(T2),T12), index.fdr[2]]
    #M12[T12,i] <- statistics.matrix[T12, index.fdr[2]]
  }#drugs
  
  Filename <- file.path(path.diagrams, paste0("MarkersPercentBoxPlot.pdf"))
  pdf(file = Filename, height=7, width=14)
  boxplot(M1, main="Isoforms vs Genes", xaxt="n", ylab=stat,cex.axis= 0.7, las=2,col = c("blue","red","purple"))
  
  text(seq(2,ncol(ccle.drug.sensitivity.ordered)*3,by=3), par("usr")[3], labels = gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7) 
  legend("topleft", legend = c("Gene specific","Isoform specific","Common genes"), fill = c("blue","red","purple"), bty="n")
  
  dev.off()
  
  
}
fnIsoformVSGeneExp <- function(FDR_CutOff = 0.01, m1 ,m2) {
  index.fdr <- NULL
  exp <- expand.grid(c(m1,m2,paste(m1,m2,sep="_")),gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)))
  M1 <- matrix(NA, nrow(FDR_List), ncol = 3 * drugs_No)
  rownames(M1) <- rownames(FDR_List)
  colnames(M1) <- paste(exp[,1],exp[,2],sep="_")
  
  M2 <- matrix(NA, nrow(FDR_List), ncol = drugs_No)
  rownames(M2) <- rownames(FDR_List)
  colnames(M2) <- colnames(ccle.drug.sensitivity.ordered)
  
  M12 <- matrix(NA, nrow(FDR_List), ncol = drugs_No)  
  rownames(M12) <- rownames(FDR_List)
  colnames(M12) <- colnames(ccle.drug.sensitivity.ordered)
  
  k <- 1
  for (j in 1:drugs_No)
  {
    i <- DOI[j]
    
    index.fdr[1] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),m1,sep ="_")
    index.fdr[2] <- paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),m2,sep ="_")
    
    
    T1 <- data.frame(cbind("fdr" = FDR_List[, index.fdr[1]], "stat" = statistics.matrix[,index.fdr[1]], "isoform" =  best.isoforms.matrix[,gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]), stringsAsFactors = FALSE)
    T1 <- rownames(subset (T1 , fdr < FDR_CutOff & isoform != ""))
    
    T2 <- data.frame(cbind("fdr" = FDR_List[, index.fdr[2]], "stat" = statistics.matrix[,index.fdr[2]], "isoform" =  best.isoforms.matrix[,gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]), stringsAsFactors = FALSE)
    T2 <- rownames(subset (T2 , fdr < FDR_CutOff & isoform != ""))
    
    
    T12 <- intersect(T1, T2)
    
    if (length(setdiff(T1,T12))<= 1){M1[setdiff(T1,T12),k] <- median(ccle.isoforms.fpkm[,best.isoforms.matrix[setdiff(T1,T12),gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]])}else{M1[setdiff(T1,T12),k] = apply(ccle.isoforms.fpkm[  ,best.isoforms.matrix[setdiff(T1,T12),gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]], 2, FUN = median)}
    if (length(setdiff(T2,T12)) <= 1){M1[setdiff(T2,T12),k+1] <- median(ccle.isoforms.fpkm[,best.isoforms.matrix[setdiff(T2,T12),gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]])}else{ M1[setdiff(T2,T12),k+1] = apply(ccle.isoforms.fpkm[  ,best.isoforms.matrix[setdiff(T2,T12),gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]], 2, FUN = median)}
    if(length(T12) <=1){M1[T12,k+2] <- median(ccle.isoforms.fpkm[  ,best.isoforms.matrix[T12,gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]])}else{M1[T12,k+2] = apply(ccle.isoforms.fpkm[  ,best.isoforms.matrix[T12,gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])]], 2, FUN = median)}
    k <- k + 3
    #M2[setdiff(names(T2),T12),i] <- statistics.matrix[setdiff(names(T2),T12), index.fdr[2]]
    #M12[T12,i] <- statistics.matrix[T12, index.fdr[2]]
  }#drugs
  
  Filename <- file.path(path.diagrams, paste0("MarkersPercentBoxPlot_Expression.pdf"))
  pdf(file = Filename, height=7, width=14)
  boxplot(M1, main="Isoforms vs Genes", xaxt="n", ylab="Expression",cex.axis= 0.7, las=2,col = c("blue","red","purple"))
  
  text(seq(2,ncol(ccle.drug.sensitivity.ordered)*3,by=3), par("usr")[3], labels = gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7) 
  legend("topleft", legend = c("Gene Exp","Isoform Exp","Common (Isoform Exp)"), fill = c("blue","red","purple"), bty="n")
  
  dev.off()
  
  
}
fnSensitivity.tissue <- function (gene, isoform, drug) {
  sensitivity.ccle <- ccle.drug.sensitivity[,drug]
  sensitivity.ccle <- sensitivity[complete.cases(sensitivity.ccle)]
  tissue.types.ccle <- table(ccle.drug.tissuetype[names(sensitivity.ccle),])
  
  sensitivity.gdsc <- gdsc.drug.sensitivity[,drug]
  sensitivity.gdsc <- sensitivity[complete.cases(sensitivity.gdsc)]
  tissue.types.gdsc <- table(ccle.drug.tissuetype[names(sensitivity.gdsc),])
  
  tissue.types <- intersect(names(tissue.types.ccle[tissue.types.ccle >= 5]) & names(tissue.types.gdsc[tissue.types.gdsc >= 5]))
  
  models.tissues <- c("M0","M1","M2","M3B")
  M_Comp_cases <- list()
  M_res <- list()
  for(i in 1:length(models.tissues))
  {
    M_res[[i]] <- numeric()
    names(M_res)[i] <- models.tissues[i]
  }
  boot.models <- list()
  Pvalues.tissue <- list()
  statistics.tissue <- list()
  #layout(matrix(1:24,6,4))
  
  for(i in 1: length(tissue.types))
  {
    statistics.tissue[[i]] <- matrix(0,nrow=5, ncol=length(models.tissues), byrow=FALSE)
    rownames(statistics.tissue[[i]]) <- c("mean","median","min","max","var")
    colnames(statistics.tissue[[i]]) <- models.tissues
    names(statistics.tissue)[i] <- names(tissue.types)[i]
    
    Pvalues.tissue[[i]] <- matrix(NA,nrow=length(models.tissues), ncol=length(models.tissues), byrow=FALSE)
    rownames(Pvalues.tissue[[i]]) <- models.tissues
    colnames(Pvalues.tissue[[i]]) <- models.tissues
    names(Pvalues.tissue)[i] <- names(tissue.types)[i]
    
    tissue <- subset(ccle.drug.tissuetype, ccle.drug.tissuetype$tissue.type == names(tissue.types)[i])
    Genedf0 <-data.frame(sensitivity[intersect(rownames(tissue),names(sensitivity))])
    colnames(Genedf0)[1] <- drug
    
    linearFormula0 <- paste(drug,"~ 1")
    M0 <- lm (as.formula(linearFormula0) , Genedf0)
    Pvalues.tissue[[i]]["M0","M0"] <- (summary(M0))$coefficients[1,1]
    
    linearFormula <- paste(drug,"~ gene")
    GenedfM <- Genedf0
    GenedfM[,"gene"] <- ccle.drug.microarray.exp[rownames(Genedf0),gene]
    M1 <- lm (as.formula(linearFormula) , GenedfM)
    Pvalues.tissue[[i]]["M1","M0"] <- ifelse(nrow((summary(M1))$coefficients)==2,(summary(M1))$coefficients[2,1],0)
    
    Genedf <- Genedf0
    Genedf[,"gene"] <- ccle.genes.fpkm[rownames(Genedf0),gene]
    M2 <- lm (as.formula(linearFormula), Genedf)
    Pvalues.tissue[[i]]["M2","M0"] <- ifelse(nrow((summary(M2))$coefficients)==2,(summary(M2))$coefficients[2,1],0)
    
    linearFormula_Iso <- paste(drug,"~ isoform")
    Isoformsdf <- Genedf0
    Isoformsdf[,"isoform"] <- ccle.isoforms.fpkm[rownames(Genedf0),isoform]
    M3B <- lm (linearFormula_Iso, Isoformsdf)
    Pvalues.tissue[[i]]["M3B","M0"] <- ifelse(nrow((summary(M3B))$coefficients)==2,(summary(M3B))$coefficients[2,1],0)
    
    boot.models[["M0"]] <- list(data = Genedf0, formula = linearFormula0, object = M0)
    boot.models[["M1"]] <- list(data = GenedfM, formula = linearFormula, object = M1)
    boot.models[["M2"]] <- list(data = Genedf, formula = linearFormula, object = M2)
    boot.models[["M3B"]] <- list(data = Isoformsdf, formula = linearFormula_Iso, object = M3B)
    
    M_res <- fnboot(boot.models, isoforms.no = isoforms_No_List[gene,"isoforms.NO"], "adj.r.squared", R = 100)
    for(j in 1:length(models.tissues))
    {
      statistics.tissue[[i]]["mean",models.tissues[j]] <- mean(M_res[[models.tissues[j]]])
      statistics.tissue[[i]]["median",models.tissues[j]] <- median(M_res[[models.tissues[j]]])
      statistics.tissue[[i]]["min",models.tissues[j]] <- min(M_res[[models.tissues[j]]])
      statistics.tissue[[i]]["max",models.tissues[j]] <- max(M_res[[models.tissues[j]]])
      statistics.tissue[[i]]["var",models.tissues[j]] <- var(M_res[[models.tissues[j]]])
    }
    for(j in 1:(length(models.tissues) - 1))
    {
      for(k in (j + 1):length(models.tissues))
      {
        Pvalues.tissue[[i]][models.tissues[j],models.tissues[k]] <- wilcox.test(M_res[[models.tissues[k]]] , M_res[[models.tissues[j]]], paired = TRUE, alternative = "greater")$p.value
        if(is.na(Pvalues.tissue[[i]][models.tissues[j],models.tissues[k]])){Pvalues.tissue[[i]][models.tissues[j],models.tissues[k]] <- 1}
      }
      ### For Select Based
      Pvalues.tissue[[i]][models.tissues[j],"M3B"] <- min(isoforms_No_List[gene,"isoforms.NO"] * Pvalues.tissue[[i]][models.tissues[j],"M3B"], 1)
    }
  }
  return (list(p.values = Pvalues.tissue, statistics = statistics.tissue, tissues = tissue.types))
}
associations.all.drugs.tissues.significant <-function(associations,cut_off) {
  rrr <- NULL
  for(i in 1:length(associations))
  {
    if(nrow(associations[[i]]) > 0)
    {
      rr <- list()
      
      for(j in 1:nrow(associations[[i]]))
      {
        file.expressions <- "ccle_drug_both_exp_ensembl.RData"
        load(file.path(path.data,file.expressions), verbose = T)
        
        if("haematopoietic_and_lymphoid_tissue" %in% ccle.drug.tissuetype[,1])
        {
          ccle.drug.tissuetype[which(ccle.drug.tissuetype[,1] == "haematopoietic_and_lymphoid_tissue"),1] = "haematopoietic_and_lymphoid"
        }
        
        fnsensitivity.tissues.ccle <- fnSensitivity.tissue(rownames(associations[[i]])[j],as.character(associations[[i]][j,"isoform_ccle"]),paste0("drugid_",names(associations)[i]),"ccle")
        
        ccle.drug.sensitivity <- ccle.drug.sensitivity[rownames(gdsc.drug.sensitivity),colnames(gdsc.drug.sensitivity)]
        ccle.drug.microarray.exp  <- ccle.drug.microarray.exp[rownames(gdsc.drug.sensitivity),]
        ccle.genes.fpkm  <- ccle.genes.fpkm[rownames(gdsc.drug.sensitivity),]
        ccle.isoforms.fpkm  <- ccle.isoforms.fpkm[rownames(gdsc.drug.sensitivity),]
        ccle.drug.tissuetype <- subset(ccle.drug.tissuetype, rownames(ccle.drug.tissuetype) %in% rownames(gdsc.drug.sensitivity))
        
        
        fnsensitivity.tissues.gdsc <- fnSensitivity.tissue(rownames(associations[[i]])[j],as.character(associations[[i]][j,"isoform_gdsc"]),paste0("drugid_",names(associations)[i]),"gdsc")
        
        for(k in 1:length(fnsensitivity.tissues.ccle$p.values))
        {
          tissue <- names(fnsensitivity.tissues.ccle$p.values)[k]
          if(tissue %in% names(fnsensitivity.tissues.gdsc$p.values))
          {
            k2 <- which(names(fnsensitivity.tissues.gdsc$p.values)==tissue)
            if(fnsensitivity.tissues.ccle$p.values[[k]]["M0","M3B"] < cut_off) 
            {
              if(is.null(rr[[tissue]]))
              {
                rr[[tissue]] <- data.frame(matrix(NA,nrow=0,ncol=ncol(associations[[i]])))
                colnames( rr[[tissue]]) <- colnames(associations[[i]])
              }
              
              index <-  nrow(rr[[tissue]]) + 1
              rr[[tissue]][index,"Microarray.pvalue_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M0","M1"]
              rr[[tissue]][index,"Genes.pvalue_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M0","M2"]
              rr[[tissue]][index,"Isoforms.pvalue_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M0","M3B"]
              rr[[tissue]][index,"Microarray.estimate_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M1","M0"]
              rr[[tissue]][index,"Genes.estimate_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M2","M0"]        
              rr[[tissue]][index,"Isoforms.estimate_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M3B","M0"] 
              rr[[tissue]][index,sprintf("%s_ccle",stat)] <- fnsensitivity.tissues.ccle$statistics[[k]]["median","M3B"]
              rr[[tissue]][index,"n_ccle"] <- fnsensitivity.tissues.ccle$tissues[tissue]
              rr[[tissue]][index,"symbol_ccle"] <- as.character(associations[[i]][j,"symbol_ccle"])
              rr[[tissue]][index,"isoform_ccle"] <- as.character(associations[[i]][j,"isoform_ccle"])
              rr[[tissue]][index,"isoform.exp_ccle"] <- as.character(associations[[i]][j,"isoform.exp_ccle"])
              if(fnsensitivity.tissues.gdsc$p.values[[k2]]["M0","M3B"] < cut_off)
              {
                rr[[tissue]][index,"Microarray.pvalue_gdsc"] <- fnsensitivity.tissues.gdsc$p.values[[k2]]["M0","M1"]
                rr[[tissue]][index,"Genes.pvalue_gdsc"] <- fnsensitivity.tissues.gdsc$p.values[[k2]]["M0","M2"]
                rr[[tissue]][index,"Isoforms.pvalue_gdsc"] <- fnsensitivity.tissues.gdsc$p.values[[k2]]["M0","M3B"]
                rr[[tissue]][index,"Microarray.estimate_gdsc"] <- fnsensitivity.tissues.gdsc$p.values[[k2]]["M1","M0"]
                rr[[tissue]][index,"Genes.estimate_gdsc"] <- fnsensitivity.tissues.gdsc$p.values[[k2]]["M2","M0"]        
                rr[[tissue]][index,"Isoforms.estimate_gdsc"] <- fnsensitivity.tissues.gdsc$p.values[[k2]]["M3B","M0"]
                rr[[tissue]][index,sprintf("%s_gdsc",stat)] <- fnsensitivity.tissues.gdsc$statistics[[k2]]["median","M3B"]
                rr[[tissue]][index,"n_gdsc"] <- fnsensitivity.tissues.gdsc$tissues[tissue]
                rr[[tissue]][index,"isoform_gdsc"] <- as.character(associations[[i]][j,"isoform_gdsc"])
                rr[[tissue]][index,"isoform.exp_gdsc"] <- as.character(associations[[i]][j,"isoform.exp_gdsc"])
              }
              rownames(rr[[tissue]])[index] <- rownames(associations[[i]])[j]
            }
          }else{
            if(fnsensitivity.tissues.ccle$p.values[[k]]["M0","M3B"] < cut_off) 
            {
              if(is.null(rr[[tissue]]))
              {
                rr[[tissue]] <- data.frame(matrix(NA,nrow=0,ncol=ncol(associations[[i]])))
                colnames( rr[[tissue]]) <- colnames(associations[[i]])
              }
              
              index <-  nrow(rr[[tissue]]) + 1
              rr[[tissue]][index,"Microarray.pvalue_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M0","M1"]
              rr[[tissue]][index,"Genes.pvalue_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M0","M2"]
              rr[[tissue]][index,"Isoforms.pvalue_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M0","M3B"]
              rr[[tissue]][index,"Microarray.estimate_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M1","M0"]
              rr[[tissue]][index,"Genes.estimate_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M2","M0"]        
              rr[[tissue]][index,"Isoforms.estimate_ccle"] <- fnsensitivity.tissues.ccle$p.values[[k]]["M3B","M0"] 
              rr[[tissue]][index,sprintf("%s_ccle",stat)] <- fnsensitivity.tissues.ccle$statistics[[k]]["median","M3B"]
              rr[[tissue]][index,"n_ccle"] <- fnsensitivity.tissues.ccle$tissues[tissue]
              rr[[tissue]][index,"symbol_ccle"] <- as.character(associations[[i]][j,"symbol_ccle"])
              rr[[tissue]][index,"isoform_ccle"] <- as.character(associations[[i]][j,"isoform_ccle"])
              rr[[tissue]][index,"isoform.exp_ccle"] <- as.character(associations[[i]][j,"isoform.exp_ccle"])
              rownames(rr[[tissue]])[index] <- rownames(associations[[i]])[j]
            }
          }        
        }
      }
      require(WriteXLS) || stop("Library WriteXLS is not available!")
      WriteXLS::WriteXLS("rr", ExcelFileName=file.path(path.diagrams, sprintf("%s.xlsx",names(associations)[i])), row.names=TRUE)
      
      rrr <- c(rrr,list(rr))
      names(rrr)[length(rrr)] <- names(associations)[i]
    }
  }
  return(rrr)
}
associations.top.biomarkers.tissues <-function(gene, isoform, drug, specificity = c("gene","isoform","common")) {
  file.expressions <- "ccle_drug_both_exp_ensembl.RData"
  load(file.path(path.data,file.expressions), verbose = T)
  
  
  fnsensitivity.tissues.ccle <- fnSensitivity.tissue(gene, isoform, drug, "ccle")
  
  ccle.drug.sensitivity <- ccle.drug.sensitivity[rownames(gdsc.drug.sensitivity),colnames(gdsc.drug.sensitivity)]
  ccle.drug.microarray.exp  <- ccle.drug.microarray.exp[rownames(gdsc.drug.sensitivity),]
  ccle.genes.fpkm  <- ccle.genes.fpkm[rownames(gdsc.drug.sensitivity),]
  ccle.isoforms.fpkm  <- ccle.isoforms.fpkm[rownames(gdsc.drug.sensitivity),]
  ccle.drug.tissuetype <- subset(ccle.drug.tissuetype, rownames(ccle.drug.tissuetype) %in% rownames(gdsc.drug.sensitivity))
  
  fnsensitivity.tissues.gdsc <- fnSensitivity.tissue(gene, isoform, drug, "gdsc")
  
  rr <- list()
  for(k in 1:length(fnsensitivity.tissues.ccle$p.values))
  {
    tissue <- names(fnsensitivity.tissues.ccle$p.values)[k]
    if(tissue %in% names(fnsensitivity.tissues.gdsc$p.values))
    {
      k2 <- which(names(fnsensitivity.tissues.gdsc$p.values)==tissue)
      
      if(specificity == "isoform")
      {
        rr[[tissue]] <- mean(c(fnsensitivity.tissues.gdsc$statistics[[k2]]["median","M3B"],fnsensitivity.tissues.ccle$statistics[[k]]["median","M3B"])) * sign(fnsensitivity.tissues.ccle$p.values[[k]]["M3B","M0"])
      }else if(specificity == "gene"){
        rr[[tissue]] <- mean(c(fnsensitivity.tissues.gdsc$statistics[[k2]]["median","M2"],fnsensitivity.tissues.ccle$statistics[[k]]["median","M2"])) * sign(fnsensitivity.tissues.ccle$p.values[[k]]["M2","M0"])
      }else{
        rr[[tissue]] <- mean(c(fnsensitivity.tissues.gdsc$statistics[[k2]]["median","M3B"],fnsensitivity.tissues.ccle$statistics[[k]]["median","M3B"],fnsensitivity.tissues.gdsc$statistics[[k2]]["median","M2"],fnsensitivity.tissues.ccle$statistics[[k]]["median","M2"])) * sign(fnsensitivity.tissues.ccle$p.values[[k]]["M2","M0"])  
      }
    }
  }
  return(rr)
}
# BioNo = "all" or Number of biomarkers 
fnTop.significant.biomarkers <- function(associations, cut_off=.01, BioNo=50, rank.type=c("pvalue", "pvalue.adj")){
  rr <- list()
  
  for(i in 1:length(associations))
  {
    significant_table <- data.frame(matrix(NA, ncol=10, nrow=nrow(associations[[i]])))
    colnames(significant_table) <- c("symbol", "biomarker.id", "type", "specificity", "estimate", stat, adjustment.method, "pvalue", "rank", "delta.rank")
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
    significant_table[,stat] <- associations[[i]][,sprintf("gene.%s",stat)]
    significant_table[, adjustment.method] <- associations[[i]][, paste0("Genes.", adjustment.method)]
    significant_table$pvalue <- associations[[i]]$Genes.pvalue
    
    significant_table$gene.id <- rownames(associations[[i]])
    significant_table$transcript.id <- associations[[i]]$isoform
    
    significant_table$isoforms.no <- isoforms_No_List[rownames(significant_table),1]
    
    temp <- data.frame(matrix(NA, ncol=10, nrow=nrow(associations[[i]])))
    colnames(temp) <- c("symbol","biomarker.id","type","specificity","estimate",stat,adjustment.method, "pvalue","rank","delta.rank")
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
    temp[,stat] <- associations[[i]][,sprintf("isoform.%s",stat)]
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
fnTop.significant.biomarkers.all.types <- function(associations, cut_off=.01, rank.type=c("pvalue", "pvalue.adj")) {
  rr <- list()
  
  for(i in 1:length(associations))
  {
    tt <- cbind("type"=NA, "pvalue"=NA, "pvalue.adj"=NA, associations[[i]])
    
    xx <- ifelse(rank.type == "pvalue.adj", adjustment.method, "pvalue")
    gene.col <- paste0("Genes.", xx)
    isoform.col <- paste0("Isoforms.", xx)
    mut.col <- paste0("Mutations.", xx) 
    cnv.col <- paste0("CNV.", xx)
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "gene"
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "gene & isoform"
    tt[which(tt[, gene.col] < cut_off &
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "gene & mutation"
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "gene & cnv"
    tt[which(tt[, gene.col] >= cut_off &
               tt[, isoform.col] < cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "isoform"
    tt[which(tt[, gene.col] >= cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "isoform & mutation"
    tt[which(tt[, gene.col] >= cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "isoform & cnv"
    tt[which(tt[, gene.col] >= cut_off & 
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "mutation"
    tt[which(tt[, gene.col] >= cut_off & 
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "mutation & cnv"
    tt[which(tt[, gene.col] >= cut_off & 
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "cnv"
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] >= cut_off), "type"] <- "gene & isoform & mutation"
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] >= cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "gene & isoform & cnv"
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] >= cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "gene & mutation & cnv"
    tt[which(tt[, gene.col] >= cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "isoform & mutation & cnv"
    tt[which(tt[, gene.col] < cut_off & 
               tt[, isoform.col] < cut_off &
               tt[, mut.col] < cut_off &
               tt[, cnv.col] < cut_off), "type"] <- "all"
    tt <- tt[which(!is.na(tt[,"type"])), , drop=FALSE]
    cc <- grep(".pvalue$", colnames(tt))
    tt[,"pvalue"] <- apply(tt, 1, function(x){min(as.numeric(x[cc]))})
    cc <- grep(adjustment.method, colnames(tt))
    tt[,"pvalue.adj"] <- apply(tt, 1, function(x){min(as.numeric(x[cc]))})
    tt <- tt[order(tt[,rank.type]),]
    rr[[names(associations)[i]]] <- tt
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
    
    if(sum(percentage.biomarkers.type[,i]) < 100){percentage.biomarkers.type["common",i] <- percentage.biomarkers.type["common",i] + 100 - sum(percentage.biomarkers.type[,i])}
    if(sum(percentage.biomarkers.type[,i]) > 100){percentage.biomarkers.type["common",i] <- 100 - sum(percentage.biomarkers.type[c(1,3),i])}
  }
  return(percentage.biomarkers.type)
}
fnGeneIsoformCorrelation <- function(biomarkers) {
  exp <- expand.grid(c("common","gene.specific","isoform.specific"),gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)))
  cor.gene.isoform <- matrix(NA,  ncol = 3 * drugs_No, nrow = max(sapply(biomarkers,function(x){nrow(x)})))
  colnames(cor.gene.isoform) <- paste(exp[,2],exp[,1],sep="_")
  
  for (j in 1:drugs_No)
  {
    i <- DOI[j]
    drug <- gsub("drugid_","",colnames(ccle.drug.sensitivity)[i])
    index.common <- paste(drug,"common",sep ="_")
    index.gene.specific <- paste(drug,"gene.specific",sep ="_")    
    index.isoform.specific <- paste(drug,"isoform.specific",sep ="_")
    
    c <- 1;g <- 1; t <- 1;
    for( x in 1:nrow(biomarkers[[drug]]))
    {
      if(!is.na(biomarkers[[drug]][x,"specificity"]) & as.character(biomarkers[[drug]][x,"transcript.id"]) != "")
      {
        if(biomarkers[[drug]][x,"specificity"] == "common")
        {
          cor.gene.isoform[c,index.common] <- abs(cor(ccle.genes.fpkm[,as.character(biomarkers[[drug]][x,"gene.id"])],ccle.isoforms.fpkm[,as.character(biomarkers[[drug]][x,"transcript.id"])], method = "spearman"))
          c <- c + 1
        }
        if(biomarkers[[drug]][x,"specificity"] == "gene.specific")
        {
          cor.gene.isoform[g,index.gene.specific] <- abs(cor(ccle.genes.fpkm[,as.character(biomarkers[[drug]][x,"gene.id"])],ccle.isoforms.fpkm[,as.character(biomarkers[[drug]][x,"transcript.id"])], method = "spearman"))
          g <- g + 1
        }
        if(biomarkers[[drug]][x,"specificity"] == "isoform.specific")
        {
          cor.gene.isoform[t,index.isoform.specific] <- abs(cor(ccle.genes.fpkm[,as.character(biomarkers[[drug]][x,"gene.id"])],ccle.isoforms.fpkm[,as.character(biomarkers[[drug]][x,"transcript.id"])], method = "spearman"))
          t <- t + 1
        }
      }
    }
  }#drugs
  
  Filename <- file.path(path.diagrams, paste0("GeneIsoformCor.pdf"))
  pdf(file = Filename, height=7, width=14)
  boxplot(cor.gene.isoform, main="Correlation of biomarkers with their correspondant gene/isoform in terms of FPKM", xaxt="n", ylab="Correlation",cex.axis= 0.7, las=2,col = c("purple","blue","red"))
  
  text(seq(2,ncol(ccle.drug.sensitivity.ordered)*3,by=3), par("usr")[3], labels = gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7) 
  legend("topleft", legend = c("Gene specific","Isoform specific","Common genes"), fill = c("blue","red","purple"), bty="n")
  
  dev.off()
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
        tissue.stat <- -100
        tissue.stat.boot <- -100
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
              ccle.stat <- compute.stat(stat, yi=M2$ccle$dataset$drug, y.hat=fitted(M2$ccle$model))
              gdsc.stat <- compute.stat(stat, yi=M2$gdsc$dataset$drug, y.hat=fitted(M2$gdsc$model))
              
              if(nrow(summary(M2$ccle$model)$coefficients) == 2){ ccle.pvalue <- summary(M2$ccle$model)$coefficients[2, 4]} else{ ccle.pvalue <- 1}
              if(nrow(summary(M2$gdsc$model)$coefficients) == 2){ gdsc.pvalue <- summary(M2$gdsc$model)$coefficients[2, 4]} else{ gdsc.pvalue <- 1}
              
              tissue.stat <- ifelse(!is.na(ccle.stat) & !is.na(gdsc.stat), ccle.stat * weight$ccle + gdsc.stat * weight$gdsc, -100)
              tissue.pvalue <- ccle.pvalue * weight$ccle + gdsc.pvalue * weight$gdsc
              
              if(boot)
              {
                boot.models <- list("M2" = list(data = M2$ccle$dataset, formula = M2$ccle$formula, object = M2$ccle$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
            
                boot.models <- list("M2" = list(data = M2$gdsc$dataset, formula = M2$gdsc$formula, object = M2$gdsc$model))
                M_res.gdsc <- fnboot(models = boot.models, R = 100)
              
                tissue.stat.boot <- median(M_res.ccle[["M2"]]) * weight$ccle + median(M_res.gdsc[["M2"]]) * weight$gdsc
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
                ccle.stat <- compute.stat(stat, M3B$ccle$dataset$drug, fitted(M3B$ccle$model))
                gdsc.stat <- compute.stat(stat, M3B$gdsc$dataset$drug, fitted(M3B$gdsc$model))
                
                if(nrow(summary(M3B$ccle$model)$coefficients) == 2){ ccle.pvalue <- summary(M3B$ccle$model)$coefficients[2, 4]} else{ ccle.pvalue <- 1}
                if(nrow(summary(M3B$gdsc$model)$coefficients) == 2){ gdsc.pvalue <- summary(M3B$gdsc$model)$coefficients[2, 4]} else{ gdsc.pvalue <- 1}
                
                tissue.stat <- ifelse(!is.na(ccle.stat) & !is.na(gdsc.stat), ccle.stat * weight$ccle + gdsc.stat * weight$gdsc, -100)
                tissue.pvalue <- ccle.pvalue * weight$ccle + gdsc.pvalue * weight$gdsc

                if(boot)
                {
                  
                  boot.models <- list("M3B" = list(data = M3B$ccle$dataset, formula = M3B$ccle$formula, object = M3B$ccle$model))
                  M_res.ccle <- fnboot(models = boot.models, R = 100)
              
                  boot.models <- list("M3B" = list(data = M3B$gdsc$dataset, formula = M3B$gdsc$formula, object = M3B$gdsc$model))
                  M_res.gdsc <- fnboot(models = boot.models, R = 100)
                
                  tissue.stat.boot <- median(M_res.ccle[["M3B"]]) * weight$ccle + median(M_res.gdsc[["M3B"]]) * weight$gdsc
                }
              }
          }
        }, "CCLE" = {
            
          if (significant_table[j, "type"] == "gene"){
            M0 <- fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle")
            M2 <- fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0, data = ccle.genes.fpkm[,gene.id])
            if(!is.null(M2))
            {
              tissue.stat <- compute.stat(stat, M2$dataset$drug, fitted(M2$model))
              tissue.pvalue <- summary(M2$model)$coefficients[2, 4]
              if(nrow(summary(M2$model)$coefficients) == 2){ tissue.pvalue <- summary(M2$model)$coefficients[2, 4]} else{ tissue.pvalue <- 1}
              
              if(boot)
              {
                boot.models <- list("M2" = list(data = M2$dataset, formula = M2$formula, object = M2$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
                tissue.stat.boot <- median(M_res.ccle[["M2"]])
              }
              
            }
          }
          if (significant_table[j, "type"] == "isoform"){
              M0 <- fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle")
              M3B <- fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id))
              if(!is.null(M3B))
              {
                tissue.stat <- compute.stat(stat, M3B$dataset$drug, fitted(M3B$model))
                if(nrow(summary(M3B$model)$coefficients) == 2){ tissue.pvalue <- summary(M3B$model)$coefficients[2, 4]} else{ tissue.pvalue <- 1}
                
                if(boot)
                {
                  boot.models <- list("M3B" = list(data = M3B$dataset, formula = M3B$formula, object = M3B$model))
                  M_res.ccle <- fnboot(models = boot.models, R = 100)  
               
                  tissue.stat.boot <- median(M_res.ccle[["M3B"]])
                }
              }
            }
          })
          significant_table[j, tissue] <- tissue.stat
          significant_table[j, paste0(tissue,"_pvalue")] <- tissue.pvalue
          if(boot){significant_table[j, paste0(tissue,"_boot")] <- tissue.stat.boot}
      }
    }
    rr[[drug]] <- significant_table
  }
  return(rr)
}
####check if biomarkers recognised from all data is working in breast cell lines based on bootstrap 
##but for each given isoform in all tissue types it should be kept so that we would have R2 for that isoform
fnidentify.tissue.specific.biomarkers.bootstarp <- function(biomarkers, tissues.pvalue, tissues.statistics, tissue.best.isoform, tissue) {
  rr <- list()
  for(i in 1:length(biomarkers))
  {
    significant_table <- biomarkers[[i]]
    drug <- names(biomarkers)[i]
    if(!all(is.na(significant_table)))
    {
      for(j in 1:nrow(significant_table))
      {
        tissue.stat <- 0
        if (significant_table[j, "type"] == "gene"){
          if((sign(tissues.pvalue[[significant_table[j, "gene.id"]]][[drug]]["M2","M0"]) == sign(significant_table[j, "estimate"])))
          {
            tissue.stat <- tissues.statistics[[significant_table[j, "gene.id"]]][[drug]]["median","M2"]
          }
        }
        if (significant_table[j, "type"] == "isoform"){
          if(tissue.best.isoform[[significant_table[j, "gene.id"]]][[drug]] == as.character(significant_table[j, "transcript.id"])){
            if((sign(tissues.pvalue[[significant_table[j, "gene.id"]]][[drug]]["M3B","M0"]) == sign(significant_table[j, "estimate"])))
            {
              tissue.stat <- tissues.statistics[[significant_table[j, "gene.id"]]][[drug]]["median","M3B"]
            }
          }
        }
        significant_table[j, paste0(tissue,"_boot")] <- tissue.stat
      }
    }
    rr[[drug]] <- significant_table
  }
  return(rr)
}
fnTop.significant.biomarkers.heatmap <- function(top.significant.biomarkers, drug) {
  dtassociation <- top.significant.biomarkers[[drug]]
  drugid <- drug
  RowNames <- NULL
  for(i in 1:nrow(dtassociation))
  {
    temp <- associations.top.biomarkers.tissues(gene = rownames(dtassociation)[i], isoform = dtassociation[i,"isoform_ccle"], drug = drugid, specificity = dtassociation[i,"specificity"])
    if(i == 1)
    {
      df.stat.tissue <- matrix(0, ncol = length(temp) + 1,nrow = nrow(dtassociation))      
      colnames(df.stat.tissue) <- c(names(temp),"All")
    }
    if(dtassociation[i,"specificity"] == "isoform")
    {
      RowNames <- c(RowNames,dtassociation[i,"isoform_ccle"])
    }else if(dtassociation[i,"specificity"] == "gene"){
      RowNames <- c(RowNames,dtassociation[i,"symbol_ccle"])
    }else
    {
      RowNames <- c(RowNames,sprintf("%s_%s",dtassociation[i,"symbol_ccle"],dtassociation[i,"isoform_ccle"]))
    }
    for(tissue in names(temp))
    {
      df.stat.tissue[i,tissue] <- temp[[tissue]]
    }
    df.stat.tissue[i,"All"] <-  dtassociation[i,stat]
  }
  
  colnames(df.stat.tissue) <- capitalize(gsub("_", " ", colnames(df.stat.tissue)))
  rownames(df.stat.tissue) <- RowNames
  
  pdf(file.path(path.diagrams,sprintf("%s_%s_Heatmap.pdf",drug,stat)), height = 5,width = 5)
  heatmap(df.stat.tissue, Rowv=NA, Colv=NA, col = colorRampPalette(c("red","blue"))(100), scale="column", margins=c(5,10),cexRow =.7,cexCol = .7)
  dev.off()
}  
