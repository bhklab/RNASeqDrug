fnIsoforms_of_Gene <- function(Gene_Map = ensembl.map.genes.isoforms, GeneId)
{
  return(data.frame(strsplit(Gene_Map[,as.character(GeneId)],",")))
}
stacked.barplot.models.ranges <- function(model, 
                                          isoforms_No = c("1.isoform","n.isoforms","all"), 
                                          prototype, 
                                          main.title, 
                                          breakpoint, 
                                          yaxis = c("Regular", "Log")) #breakpoit Regular or the point in y axis where should be gapped
{
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
    ylabels <- mystacked.barplot(Filename = File, data = stacked.barplot.matrix, barsNo = models.NO, stackedNo = length(ranges.names), groupNo = drugs_No, group.labels = Glabels, ylab.label = "Number of Associated Genes", legend.lables= Legends, main.label = main.title, yaxis)
  }else{
    File <- file.path(path.diagrams, sprintf("Stacked_Gapped_BarPlot_ranges_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mystacked.gap.barplot(Filename = File, data = stacked.barplot.matrix, barsNo = models.NO, stackedNo = length(ranges.names), groupNo = drugs_No, group.labels = Glabels, ylab.label = "Number of Associated Genes", legend.lables= Legends, main.label = main.title, breakpoint)
  }
  return (ylabels)
}
stacked.barplot.models <- function(model, isoforms_No = c("1.isoform","n.isoforms","all"), prototype, main.title, breakpoint)
{
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
barplot.models <- function(model, isoforms_No = c("1.isoform","n.isoforms","all"), prototype, main.title, breakpoint,  yaxis = c("Regular", "Log")) #breakpoint c("Regular" ,"No", any number)
{
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
    mybarplot(Filename=File, data=barplot.matrix, barsNo=models.NO, groupNo=drugs_No, group.labels=Glabels, ylab.label="Number of Associated Genes", legend.lables=Models.names[prototype], main.label=main.title, yaxis=yaxis) 
  }else{
    File <- file.path(path.diagrams, sprintf("Gapped_BarPlot_%s_%s.pdf",isoforms_No,str_replace_all(main.title, "[^[:alnum:]]","")))    
    mybarplot.gap(Filename=File, data=barplot.matrix, barsNo=models.NO, group.labels=Glabels,ylab.label="Number of Associated Genes", legend.lables=Models.names[prototype], main.label=main.title, breakpoint) 
  }
}
fnVennDiagramTriple <- function(model, m3=c("M3","M3B","M3P"))
{
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
fnVennDiagramPair <- function(model, m3=c("M3","M3B","M3P"))
{
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
fnFancystat<- function(FDR, combinations, signed)
{
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
fnCuttOffs<-function(cutoff_statistics, combination)
{
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
fnWilcox <- function(model, signed)
{
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

associations.all.drugs <-function(model.rank, annot.ensembl.all.genes)
{
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
    table_write_dt[,"gene.R2"] <- statistics.matrix[,paste(gsub("drugid_","",drug),"M2",sep ="_")]
    table_write_dt[,"isoform.R2"] <- statistics.matrix[,paste(gsub("drugid_","",drug),"M3B",sep ="_")]
    
    table_write_dt <- cbind(table_write_dt$symbol, table_write_dt$n, table_write_dt[,models.col.names], table_write_dt$isoform, table_write_dt$isoform.exp, table_write_dt$gene.R2, table_write_dt$isoform.R2)
    colnames(table_write_dt) <- c("symbol","n",models.col.names,"isoform","isoform.exp","gene.R2","isoform.R2")
    ### Order based on pvalues of best isoform method
    table_write_dt <- table_write_dt[order(table_write_dt[,rank]),]
    
    rr <- c(rr, list(table_write_dt))
  }
  names(rr) <- gsub("drugid_", "", colnames(ccle.drug.sensitivity.ordered))
  return(rr)
  
}
associations.all.drugs.significant <-function(rank, method=c("isoform","gene"), cutoff=.01, R2_cut=.6, exp_cut=1, annot.ensembl.all.genes)
{
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
    table_write_dt[,"gene.R2"] <- statistics.matrix[,paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),"M2",sep ="_")]
    table_write_dt[,"isoform.R2"] <- statistics.matrix[,paste(gsub("drugid_","",colnames(ccle.drug.sensitivity)[i]),"M3B",sep ="_")]
    table_write_dt <- cbind(table_write_dt$symbol, table_write_dt$n, table_write_dt[,models.col.names], table_write_dt$isoform, table_write_dt$isoform.exp, table_write_dt$gene.R2, table_write_dt$isoform.R2)
    colnames(table_write_dt) <- c("symbol","n",models.col.names,"isoform","isoform.exp","gene.R2","isoform.R2")
    
    if(method == "isoform")
    {
      table_write_dt <- subset(table_write_dt,  
                                table_write_dt[,sprintf("%s.%s",Models.names[2],adjustment.method)]==1 & 
                                table_write_dt[,sprintf("%s.%s",Models.names[3],adjustment.method)] < cutoff & 
                                table_write_dt[,"isoform.R2"] > R2_cut & 
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
Check.KnownAssociations <- function(associations)
{
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
Check.KnownAssociations.PLOS1 <- function(associations)
{
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
fnComputeAssociateGenes <- function(FDR_CutOff = 0.01, signed)
{
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
fnComputeAssociateGenes.stat <- function(FDR_CutOff = 0.01, stat_CutOff = 0.6, signed)
{
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
fnComputeAssociateGenes.stat.range <- function(FDR_CutOff = 0.01, biomarkers.sign = c("positive","negative","all"))
{
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

fnMarkersPercent <- function(model, m1, m2)
{
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
fnIsoformVSGeneR2 <- function(FDR_CutOff = 0.01, m1 ,m2)
{
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
  boxplot(M1, main="Isoforms vs Genes", xaxt="n", ylab="R2",cex.axis= 0.7, las=2,col = c("blue","red","purple"))
  
  text(seq(2,ncol(ccle.drug.sensitivity.ordered)*3,by=3), par("usr")[3], labels = gsub("drugid_","",colnames(ccle.drug.sensitivity.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7) 
  legend("topleft", legend = c("Gene specific","Isoform specific","Common genes"), fill = c("blue","red","purple"), bty="n")
  
  dev.off()
  
  
}
fnIsoformVSGeneExp <- function(FDR_CutOff = 0.01, m1 ,m2)
{
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

fnSensitivity.tissue <- function (gene, isoform, drug)
{
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
associations.all.drugs.tissues.significant <-function(associations,cut_off)
{
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
              rr[[tissue]][index,"R2_ccle"] <- fnsensitivity.tissues.ccle$statistics[[k]]["median","M3B"]
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
                rr[[tissue]][index,"R2_gdsc"] <- fnsensitivity.tissues.gdsc$statistics[[k2]]["median","M3B"]
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
              rr[[tissue]][index,"R2_ccle"] <- fnsensitivity.tissues.ccle$statistics[[k]]["median","M3B"]
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
associations.top.biomarkers.tissues <-function(gene, isoform, drug, specificity = c("gene","isoform","common"))
{
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
fnTop.significant.biomarkers <- function(associations, cut_off=.01, BioNo=50, rank.type=c("pvalue", "pvalue.adj"))# BioNo = "all" or Number of biomarkers
{
  rr <- list()
  
  for(i in 1:length(associations))
  {
    significant_table <- data.frame(matrix(NA, ncol=10, nrow=nrow(associations[[i]])))
    colnames(significant_table) <- c("symbol", "biomarker.id", "type", "specificity", "estimate", "R2", adjustment.method, "pvalue", "rank", "delta.rank")
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
    significant_table$R2 <- associations[[i]]$gene.R2
    significant_table[, adjustment.method] <- associations[[i]][, paste0("Genes.", adjustment.method)]
    significant_table$pvalue <- associations[[i]]$Genes.pvalue
    
    significant_table$gene.id <- rownames(associations[[i]])
    significant_table$transcript.id <- associations[[i]]$isoform
    
    significant_table$isoforms.no <- isoforms_No_List[rownames(significant_table),1]
    
    temp <- data.frame(matrix(NA, ncol=10, nrow=nrow(associations[[i]])))
    colnames(temp) <- c("symbol","biomarker.id","type","specificity","estimate","R2",adjustment.method, "pvalue","rank","delta.rank")
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
    temp$R2 <- associations[[i]]$isoform.R2
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

fnExtractFeautures <- function(associations, cut_off = .01, BioNo = 50)# BioNo = "all" or Number of biomarkers
{
  rr <- list()
  
  for(i in 1:length(associations))
  {
    significant_table <- data.frame(matrix(NA, ncol = 9, nrow = nrow(associations[[i]])))
    colnames(significant_table) <- c("symbol","biomarker.id","type","specificity","estimate","R2",adjustment.method,"rank","delta.rank")
    rownames(significant_table) <- rownames(associations[[i]])
    significant_table$symbol <- associations[[i]]$symbol
    significant_table$biomarker.id <- associations[[i]]$symbol
    significant_table$type <- "gene"
    significant_table[rownames(subset(associations[[i]], associations[[i]][, paste0("Genes.", adjustment.method)] < cut_off & associations[[i]][, paste0("Isoforms.", adjustment.method)] >= cut_off)), "specificity"] <- "gene.specific"
    significant_table[rownames(subset(associations[[i]], associations[[i]][, paste0("Genes.", adjustment.method)] < cut_off & associations[[i]][, paste0("Isoforms.", adjustment.method)] < cut_off)), "specificity"] <- "common"
    significant_table$estimate <- associations[[i]]$Genes.estimate
    significant_table$R2 <- associations[[i]]$gene.R2
    significant_table[ , adjustment.method] <- associations[[i]][, paste0("Genes.", adjustment.method)]
    significant_table$gene.id <- rownames(associations[[i]])
    significant_table$transcript.id <- associations[[i]]$isoform
    
    significant_table$isoforms.no <- isoforms_No_List[rownames(significant_table),1]
    
    temp <- data.frame(matrix(NA, ncol = 9, nrow = nrow(associations[[i]])))
    colnames(temp) <- c("symbol","biomarker.id","type","specificity","estimate","R2",adjustment.method,"rank","delta.rank")
    rownames(temp) <- rownames(associations[[i]])
    temp$symbol <- associations[[i]]$symbol
    temp$biomarker.id <- associations[[i]]$isoform
    temp$type <- "isoform"
    temp[rownames(subset(associations[[i]], associations[[i]][, paste0("Isoforms.", adjustment.method)] < cut_off & associations[[i]][, paste0("Genes.", adjustment.method)] >= cut_off)), "specificity"] <- "isoform.specific"
    temp[rownames(subset(associations[[i]], associations[[i]][, paste0("Isoforms.", adjustment.method)] < cut_off & associations[[i]][, paste0("Genes.", adjustment.method)] < cut_off)), "specificity"] <- "common"
    temp$estimate <- associations[[i]]$Isoforms.estimate
    temp$R2 <- associations[[i]]$isoform.R2
    temp[ , adjustment.method] <- associations[[i]][, paste0("Isoforms.", adjustment.method)]
    temp$isoforms.no <- isoforms_No_List[rownames(temp),1]
    temp$gene.id <- rownames(associations[[i]])
    temp$transcript.id <- associations[[i]]$isoform
    
    
    significant_table <- rbind(significant_table,temp)
    rownames(significant_table) <- 1:nrow(significant_table)
    significant_table <- significant_table[order(significant_table[,adjustment.method]),]
    significant_table$rank <- 1:nrow(significant_table)    
    significant_table$delta.rank <- 0
    
    
    if(BioNo != "all")
    {
      for ( j in 1: BioNo)
      {
        temp <- subset(significant_table, significant_table$symbol == significant_table[j,"symbol"])
        significant_table[j,"delta.rank"] <- abs(temp[1,"rank"] - temp[2,"rank"])
      }
      significant_table <- significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE]
      last <- ifelse(nrow(significant_table) < BioNo,  nrow(significant_table), BioNo)
    }else{
      significant.NO <- nrow(significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE])
      
      for ( j in 1: significant.NO)
      {
        temp <- subset(significant_table, significant_table$symbol == significant_table[j,"symbol"])
        significant_table[j,"delta.rank"] <- abs(temp[1,"rank"] - temp[2,"rank"])
      }
      significant_table <- significant_table[which(significant_table[ , adjustment.method] < cut_off), , drop=FALSE]
      
      last <- nrow(significant_table)
    }
    
    #last <- ifelse(nrow(significant_table) < BioNo,  nrow(significant_table), BioNo)
    rr[[names(associations)[i]]] <- significant_table[1:last,]
  }
  return(rr)
}
fnPercentageBiomarkersType <- function(associations)
{
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
  }
  return(percentage.biomarkers.type)
}
fnGeneIsoformCorrelation <- function(biomarkers)
{
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

fnidentify.tissue.specific.biomarkers <- function(biomarkers, boot=FALSE)
{
  
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
        tissue.R2 <- -100
        tissue.R2.boot <- -100
        gene.id <- as.character(significant_table[j, "gene.id"])
        transcript.id <- as.character(significant_table[j, "transcript.id"])
        switch(training.type, "CCLE_GDSC"= {
          if (significant_table[j, "type"] == "gene"){
            M0 <- list(ccle = fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle"), 
                        gdsc = fnCreateNullModel(drug = Drugs_ToCheck, assay = "gdsc"))
            weight <- list(ccle = M0$ccle$n / (M0$ccle$n + M0$gdsc$n), gdsc = M0$gdsc$n / (M0$ccle$n + M0$gdsc$n))
            
            M2 <- list(ccle = fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0$ccle, data = ccle.genes.fpkm[,gene.id]), 
                        gdsc = fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0$gdsc, data = ccle.genes.fpkm[, gene.id]))
            
            if(!is.null(M2$ccle) && !is.null(M2$gdsc)){if(!is.na(M2$ccle$coefficient) && !is.na(M2$gdsc$coefficient)){if(sign(M2$ccle$coefficient) != sign(M2$gdsc$coefficient)){M2$ccle <- NULL; M2$gdsc <- NULL}}}
            if(!is.null(M2$ccle) && !is.null(M2$gdsc))
            {
              ccle.R2 <- compute.R2(M2$ccle$dataset$drug, fitted(M2$ccle$model))
              gdsc.R2 <- compute.R2(M2$gdsc$dataset$drug, fitted(M2$gdsc$model))
              
              if(nrow(summary(M2$ccle$model)$coefficients) == 2){ ccle.pvalue <- summary(M2$ccle$model)$coefficients[2, 4]} else{ ccle.pvalue <- 1}
              if(nrow(summary(M2$gdsc$model)$coefficients) == 2){ gdsc.pvalue <- summary(M2$gdsc$model)$coefficients[2, 4]} else{ gdsc.pvalue <- 1}
              
              tissue.R2 <- ccle.R2 * weight$ccle + gdsc.R2 * weight$gdsc
              tissue.pvalue <- ccle.pvalue * weight$ccle + gdsc.pvalue * weight$gdsc
              
              if(boot)
              {
                boot.models <- list("M2" = list(data = M2$ccle$dataset, formula = M2$ccle$formula, object = M2$ccle$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
            
                boot.models <- list("M2" = list(data = M2$gdsc$dataset, formula = M2$gdsc$formula, object = M2$gdsc$model))
                M_res.gdsc <- fnboot(models = boot.models, R = 100)
              
                tissue.R2.boot <- median(M_res.ccle[["M2"]]) * weight$ccle + median(M_res.gdsc[["M2"]]) * weight$gdsc
              }
            }
          }
          if (significant_table[j, "type"] == "isoform"){
              M0 <- list(ccle = fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle"), 
                         gdsc = fnCreateNullModel(drug = Drugs_ToCheck, assay = "gdsc"))
              weight <- list(ccle = M0$ccle$n / (M0$ccle$n + M0$gdsc$n), gdsc = M0$gdsc$n / (M0$ccle$n + M0$gdsc$n))
              
              M3B <- list(ccle = fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0$ccle, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id)), 
                          gdsc = fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0$gdsc, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id)))
              if(!is.null(M3B$ccle) && !is.null(M3B$gdsc)){if(!is.na(M3B$ccle$coefficient) && !is.na(M3B$gdsc$coefficient)){if(sign(M3B$ccle$coefficient) != sign(M3B$gdsc$coefficient)){M3B$ccle <- NULL; M3B$gdsc <- NULL}}}
              if(!is.null(M3B$ccle) && !is.null(M3B$gdsc))
              {
                ccle.R2 <- compute.R2(M3B$ccle$dataset$drug, fitted(M3B$ccle$model))
                gdsc.R2 <- compute.R2(M3B$gdsc$dataset$drug, fitted(M3B$gdsc$model))
                
                if(nrow(summary(M3B$ccle$model)$coefficients) == 2){ ccle.pvalue <- summary(M3B$ccle$model)$coefficients[2, 4]} else{ ccle.pvalue <- 1}
                if(nrow(summary(M3B$gdsc$model)$coefficients) == 2){ gdsc.pvalue <- summary(M3B$gdsc$model)$coefficients[2, 4]} else{ gdsc.pvalue <- 1}
                
                tissue.R2 <- ccle.R2 * weight$ccle + gdsc.R2 * weight$gdsc
                tissue.pvalue <- ccle.pvalue * weight$ccle + gdsc.pvalue * weight$gdsc

                if(boot)
                {
                  
                  boot.models <- list("M3B" = list(data = M3B$ccle$dataset, formula = M3B$ccle$formula, object = M3B$ccle$model))
                  M_res.ccle <- fnboot(models = boot.models, R = 100)
              
                  boot.models <- list("M3B" = list(data = M3B$gdsc$dataset, formula = M3B$gdsc$formula, object = M3B$gdsc$model))
                  M_res.gdsc <- fnboot(models = boot.models, R = 100)
                
                  tissue.R2.boot <- median(M_res.ccle[["M3B"]]) * weight$ccle + median(M_res.gdsc[["M3B"]]) * weight$gdsc
                }
              }
          }
        }, "CCLE" = {
            
          if (significant_table[j, "type"] == "gene"){
            M0 <- fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle")
            M2 <- fnCreateGeneModel(drug = Drugs_ToCheck, nullModel = M0, data = ccle.genes.fpkm[,gene.id])
            if(!is.null(M2))
            {
              tissue.R2 <- compute.R2(M2$dataset$drug, fitted(M2$model))
              tissue.pvalue <- summary(M2$model)$coefficients[2, 4]
              if(nrow(summary(M2$model)$coefficients) == 2){ tissue.pvalue <- summary(M2$model)$coefficients[2, 4]} else{ tissue.pvalue <- 1}
              
              if(boot)
              {
                boot.models <- list("M2" = list(data = M2$dataset, formula = M2$formula, object = M2$model))
                M_res.ccle <- fnboot(models = boot.models, R = 100)
                tissue.R2.boot <- median(M_res.ccle[["M2"]])
              }
              
            }
          }
          if (significant_table[j, "type"] == "isoform"){
              M0 <- fnCreateNullModel(drug = Drugs_ToCheck, assay = "ccle")
              M3B <- fnCreateIsoformModel(drug = Drugs_ToCheck, nullModel = M0, data = ccle.isoforms.fpkm, isoform = as.character(transcript.id))
              if(!is.null(M3B))
              {
                tissue.R2 <- compute.R2(M3B$dataset$drug, fitted(M3B$model))
                if(nrow(summary(M3B$model)$coefficients) == 2){ tissue.pvalue <- summary(M3B$model)$coefficients[2, 4]} else{ tissue.pvalue <- 1}
                
                if(boot)
                {
                  boot.models <- list("M3B" = list(data = M3B$dataset, formula = M3B$formula, object = M3B$model))
                  M_res.ccle <- fnboot(models = boot.models, R = 100)  
               
                  tissue.R2.boot <- median(M_res.ccle[["M3B"]])
                }
              }
            }
          })
          significant_table[j, tissue] <- tissue.R2
          significant_table[j, paste0(tissue,"_pvalue")] <- tissue.pvalue
          significant_table[j, paste0(tissue,"_boot")] <- tissue.R2.boot
      }
    }
    rr[[drug]] <- significant_table
  }
  return(rr)
}

####check if biomarkers recognised from all data is working in breast cell lines based on bootstrap 
##but for each give isoform in all tissue types it should be kept so that we would have R2 for that isoform
fnidentify.tissue.specific.biomarkers.bootstarp <- function(biomarkers, tissues.pvalue, tissues.statistics, tissue.best.isoform, tissue)
{
  rr <- list()
  for(i in 1:length(biomarkers))
  {
    significant_table <- biomarkers[[i]]
    drug <- names(biomarkers)[i]
    if(!all(is.na(significant_table)))
    {
      for(j in 1:nrow(significant_table))
      {
        tissue.R2 <- 0
        if (significant_table[j, "type"] == "gene"){
          if((sign(tissues.pvalue[[significant_table[j, "gene.id"]]][[drug]]["M2","M0"]) == sign(significant_table[j, "estimate"])))
          {
            tissue.R2 <- tissues.statistics[[significant_table[j, "gene.id"]]][[drug]]["median","M2"]
          }
        }
        if (significant_table[j, "type"] == "isoform"){
          if(tissue.best.isoform[[significant_table[j, "gene.id"]]][[drug]] == as.character(significant_table[j, "transcript.id"])){
            if((sign(tissues.pvalue[[significant_table[j, "gene.id"]]][[drug]]["M3B","M0"]) == sign(significant_table[j, "estimate"])))
            {
              tissue.R2 <- tissues.statistics[[significant_table[j, "gene.id"]]][[drug]]["median","M3B"]
            }
          }
        }
        significant_table[j, paste0(tissue,"_boot")] <- tissue.R2
      }
    }
    rr[[drug]] <- significant_table
  }
  return(rr)
}
fnTop.significant.biomarkers.heatmap <- function(top.significant.biomarkers, drug)
{
  dtassociation <- top.significant.biomarkers[[drug]]
  drugid <- drug
  RowNames <- NULL
  for(i in 1:nrow(dtassociation))
  {
    temp <- associations.top.biomarkers.tissues(gene = rownames(dtassociation)[i], isoform = dtassociation[i,"isoform_ccle"], drug = drugid, specificity = dtassociation[i,"specificity"])
    if(i == 1)
    {
      df.R2.tissue <- matrix(0, ncol = length(temp) + 1,nrow = nrow(dtassociation))      
      colnames(df.R2.tissue) <- c(names(temp),"All")
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
      df.R2.tissue[i,tissue] <- temp[[tissue]]
    }
    df.R2.tissue[i,"All"] <-  dtassociation[i,"R2"]
  }
  
  colnames(df.R2.tissue) <- capitalize(gsub("_", " ", colnames(df.R2.tissue)))
  rownames(df.R2.tissue) <- RowNames
  
  pdf(file.path(path.diagrams,sprintf("%s_R2_Heatmap.pdf",drug)), height = 5,width = 5)
  heatmap(df.R2.tissue, Rowv=NA, Colv=NA, col = colorRampPalette(c("red","blue"))(100), scale="column", margins=c(5,10),cexRow =.7,cexCol = .7)
  dev.off()
}


fnIsoformsExp <- function(GeneId, Isoforms_FPKM = ccle.isoforms.fpkm)
{
  linearFormula_Iso <- ""
  Isoformsdf <- data.frame(matrix(NA, nrow=nrow(Isoforms_FPKM), ncol=0, byrow=FALSE))
  rownames(Isoformsdf) <- rownames(Isoforms_FPKM)
  Isoforms_of_Genes <- fnIsoforms_of_Gene(GeneId = GeneId)
  for(j in 1: nrow(Isoforms_of_Genes))
  {
    if (as.character(Isoforms_of_Genes[j,]) %in% colnames(ccle.isoforms.fpkm))
    {
      Isoformsdf[,as.character(Isoforms_of_Genes[j,])] <- Isoforms_FPKM [,as.character(Isoforms_of_Genes[j,])]
      linearFormula_Iso <- paste(linearFormula_Iso , as.character(Isoforms_of_Genes[j,]) , "+")
    }
  }
  linearFormula_Iso <-  substr(linearFormula_Iso, 1, nchar(linearFormula_Iso)-1)
  return(list(FPKM = Isoformsdf, names = linearFormula_Iso))
}
fnSensitivity.isoforms <- function (Isoforms, drug, best, method = "bootstrap", stat = "adj.r.squared")
{
  Models <- c("M0",colnames(Isoforms$FPKM))
  drugs_No <- ncol(ccle.drug.sensitivity)
  M_Comp_cases <- list()
  M_res <- list()
  for(i in 1:length(Models))
  {
    M_res[[i]] <- numeric()
    names(M_res)[i] <- Models[i]
  }
  boot.models <- list()
  #layout(matrix(1:24,6,4))
  
  P_values <- matrix(NA,nrow=2, ncol=ncol(Isoforms$FPKM), byrow=FALSE)
  rownames(P_values) <- c("raw","boot")
  colnames(P_values) <- colnames(Isoforms$FPKM)
  
  statistics <- matrix(0,nrow=5, ncol=length(Models), byrow=FALSE)
  rownames(statistics) <- c("mean","median","min","max","var")
  colnames(statistics) <- Models
  
  Drugs_ToCheck <- drug
  Genedf0 <-data.frame(ccle.drug.sensitivity[,Drugs_ToCheck])
  colnames(Genedf0)[1] <- Drugs_ToCheck
  Genedf0[,"Tissue"] <- factor(ccle.drug.tissuetype[,"tissue.type"], ordered =FALSE)
  Genedf0 <- Genedf0[complete.cases(Genedf0),]
  Genedf0$Tissue <- factor(Genedf0$Tissue, ordered =FALSE)
  
  linearFormula0 <- paste(Drugs_ToCheck,"~ Tissue")
  
  M_Comp_cases[["M0"]] <- lm (as.formula(linearFormula0) , Genedf0)
  boot.models[["M0"]] <- list(data = Genedf0, formula = linearFormula0, object = M_Comp_cases[["M0"]])
  for(j in 1:ncol(Isoforms$FPKM))
  {
    Isoformsdf_temp <- data.frame(Genedf0, Isoforms$FPKM[rownames(Genedf0),j])
    colnames(Isoformsdf_temp)[ncol(Isoformsdf_temp)] <- colnames(Isoforms$FPKM)[j]
    linearFormula_Iso_temp <- paste(paste(Drugs_ToCheck,"~ Tissue + "), colnames(Isoforms$FPKM)[j])
    M_Comp_cases[[colnames(Isoforms$FPKM)[j]]] <- lm (linearFormula_Iso_temp, Isoformsdf_temp)
    boot.models[[colnames(Isoforms$FPKM)[j]]] <- list(data = Isoformsdf_temp, formula = linearFormula_Iso_temp, object = M_Comp_cases[[colnames(Isoforms$FPKM)[j]]])
  }
  M_res <- fnboot.tissuebased(boot.models, isoforms.no = ncol(Isoforms$FPKM), stat, R = 100)
  for(j in 1:length(Models))
  {
    statistics["mean",Models[j]] <- mean(M_res[[Models[j]]])
    statistics["median",Models[j]] <- median (M_res[[Models[j]]])
    statistics["min",Models[j]] <- min(M_res[[Models[j]]])
    statistics["max",Models[j]] <- max(M_res[[Models[j]]])
    statistics["var",Models[j]] <- var(M_res[[Models[j]]])
  }
  Gene_Coef_Index <- nrow((summary(M_Comp_cases[[best]]))$coefficients)
  
  for(j in 2:length(Models))
  {
    if((stat == "r.squared") || (stat == "adj.r.squared") || (stat == "cindex")){direction <- "greater"}else{direction <- "less"}
    P_values["boot",Models[j]] <- wilcox.test(M_res[[Models[j]]] , M_res[["M0"]], paired = TRUE, alternative = direction)$p.value
    if(is.na(P_values["boot",Models[j]])){P_values["boot",Models[j]] <- 1}
    P_values["boot",Models[j]] <- min(ncol(Isoforms$FPKM) * P_values["boot",Models[j]], 1)
    P_values["raw",Models[j]] <- ifelse(nrow(summary(M_Comp_cases[[Models[j]]])$coefficients) == Gene_Coef_Index, summary(M_Comp_cases[[Models[j]]])$coefficients[Gene_Coef_Index,4],1)
  }
  
  return (list(p.values = P_values, statistics = statistics))
}


fnCCLE.GDSC <- function(gene, isoform, drug, method=c("all","most.correlated","tissue.based","cor"), tissue.type="", celllines.col="", R2="", pvalue="", annot.ensembl.all.genes)
{
  ccle.sensitivity <- subset(ccle.drug.sensitivity, !is.na(ccle.drug.sensitivity[,drug]), select = drug)
  gdsc.sensitivity <- subset(gdsc.drug.sensitivity, !is.na(gdsc.drug.sensitivity[,drug]), select = drug)
  
  sensitivity <- NULL
  switch(method , "most.correlated" = {
    intersected.celllines <- intersect(names(ccle.sensitivity),names(gdsc.sensitivity))
    ccle.sensitivity <- ccle.sensitivity[intersected.celllines]
    gdsc.sensitivity <- gdsc.sensitivity[intersected.celllines]
    
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    
    sensitivity[,1] <- ccle.sensitivity
    rownames(sensitivity) <- names(ccle.sensitivity)
    sensitivity <- cbind(sensitivity, ccle.genes.fpkm[rownames(sensitivity),gene], ccle.isoforms.fpkm[rownames(sensitivity),isoform],ccle.drug.tissuetype[rownames(sensitivity),"tissue.type"])
    colnames(sensitivity) <- c("ccle.AUC", "gene.exp","isoform.exp","tissue.type")
    
    
    sensitivity<- cbind(sensitivity, ccle.genes.fpkm[rownames(sensitivity),gene], ccle.isoforms.fpkm[rownames(sensitivity),isoform],ccle.drug.tissuetype[rownames(sensitivity),"tissue.type"])
    colnames(sensitivity) <- c("ccle.AUC","gdsc.AUC","rank.diff", "gene.exp","isoform.exp","tissue.type")
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_cellLines.pdf",isoform,drug)), height=7, width=7)
    
    sensitivity <- cbind(sensitivity, "red")
    colnames(sensitivity)[ncol(sensitivity)] <- "col"
    my.xlim <- range(as.numeric(sensitivity[,"isoform.exp"]))
    my.ylim <- range(as.numeric(sensitivity[,"ccle.AUC"]))
    sensitivity[,"col"] <- colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity))
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="ccle AUC")
    points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"], pch = 19,col = sensitivity[,"col"])
    
    x <- c(round(my.xlim[1], digits = 1) + 0.1 , round(my.xlim[1], digits = 1) + 0.3, round(my.xlim[1], digits = 1) + 0.3,round(my.xlim[1], digits = 1) + 0.1)
    y <- c(round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1), round(my.ylim[2], digits = 1))
    legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "", limits = c("resistent","sensitive"),cex=.8)
    
    si <- (nrow(sensitivity) - 2):nrow(sensitivity)
    textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
    ri <- 1:3
    textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
    textxy(x[2], y[3] - .05, sprintf("pvalue = %s", round(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]][["M0","M3B"]], digits = 20)), cex=0.9)
    textxy(x[2]+.2, y[3] - .07, sprintf("R2 = %s", round(ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]][["median","M3B"]], digits = 2)), cex=0.9)
    title(sprintf("Correlation between %s(%s) and %s for \n the most correlated cell lines between gdsc & ccle",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene), "Symbol"], gsub("drugid_","",drug)),cex=.9)
    dev.off()
  }, 
  "tissue.based" = {
    unique.ccle.celllines <- setdiff(rownames(ccle.sensitivity),rownames(gdsc.sensitivity))
    unique <- cbind(ccle.sensitivity[unique.ccle.celllines,drug], ccle.genes.fpkm[unique.ccle.celllines,gene], ccle.isoforms.fpkm[unique.ccle.celllines,isoform],ccle.drug.tissuetype[unique.ccle.celllines,"tissue.type"])
    colnames(unique) <- c("ccle.AUC", "gene.exp","isoform.exp","tissue.type")
    unique <- subset(unique, unique[,"tissue.type"] == tissue.type)
    
    intersected.celllines <- intersect(rownames(ccle.sensitivity),rownames(gdsc.sensitivity))
    sensitivity <- cbind(ccle.sensitivity[intersected.celllines,], gdsc.sensitivity[intersected.celllines,], ccle.genes.fpkm[intersected.celllines,gene], ccle.isoforms.fpkm[intersected.celllines,isoform],ccle.drug.tissuetype[intersected.celllines,"tissue.type"])
    colnames(sensitivity) <- c("ccle.AUC","gdsc.AUC", "gene.exp","isoform.exp","tissue.type")
    sensitivity <- subset(sensitivity, sensitivity[,"tissue.type"] == tissue.type)
    #if(nrow(sensitivity) > 1){sensitivity <- sensitivity[order(sensitivity[,"gdsc.AUC"]),]}
    
    
    
    if(nrow(sensitivity) > 0) 
    {
      my.xlim <- c(min(min(as.numeric(sensitivity[,"isoform.exp"])), min(as.numeric(unique[,"isoform.exp"]))),max(max(as.numeric(sensitivity[,"isoform.exp"])), max(as.numeric(unique[,"isoform.exp"]))))
      my.ylim <- c(min(min(as.numeric(sensitivity[,"ccle.AUC"])), min(as.numeric(unique[,"ccle.AUC"]))),max(max(as.numeric(sensitivity[,"ccle.AUC"])), max(as.numeric(unique[,"ccle.AUC"]))))
      
      
      sensitivity <- cbind(sensitivity, "red")
      colnames(sensitivity)[ncol(sensitivity)] <- "col"
      sensitivity[,"col"] <- colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity)) #celllines.col[rownames(sensitivity)]#
      
      plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="ccle.AUC",cex.lab = .7, cex.axis = .7)
      points(unique[,"isoform.exp"], unique[,"ccle.AUC"], pch = 17,col = "grey")
      points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"], pch = 19,col = sensitivity[,"col"])
      
      x <- c(round(my.xlim[1], digits = 1) + 0.1 , round(my.xlim[1], digits = 1) + 0.3, round(my.xlim[1], digits = 1) + 0.3,round(my.xlim[1], digits = 1) + 0.1)
      y <- c(round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1), round(my.ylim[2], digits = 1))
      #legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "gdsc.AUC", limits = c("resistent","sensitive"),cex=.7)
      
      si <- (nrow(sensitivity) - 2):nrow(sensitivity)
      #textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
      ri <- 1:3
      #textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
      textxy(my.xlim[1] , my.ylim[2] - .04, sprintf("pvalue = %s\n    R2 = %s", round(pvalue, digits = 20), round(R2, digits = 2)), cex=0.6)
      title(sprintf("%s ",capitalize(gsub("_"," ",tissue.type))),cex.main=.8)
    }
    else if (nrow(unique) > 0) 
    {
      my.xlim <- range((as.numeric(unique[,"isoform.exp"])))
      my.ylim <- range((as.numeric(unique[,"ccle.AUC"])))
      
      
      plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="ccle.AUC",cex.lab = .7, cex.axis = .7)
      points(unique[,"isoform.exp"], unique[,"ccle.AUC"], pch = 17,col = "grey")
      
      title(sprintf("%s ",capitalize(gsub("_"," ",tissue.type))),cex.main=.8)
      
    }
  }, 
  "all"= {
    unique.ccle.celllines <- setdiff(names(ccle.sensitivity),names(gdsc.sensitivity))
    ccle.unique <- ccle.sensitivity[unique.ccle.celllines]
    unique <-  matrix(NA,nrow= length(ccle.unique))
    unique <- ccle.unique
    unique <- cbind(unique, ccle.genes.fpkm[names(unique),gene], ccle.isoforms.fpkm[names(unique),isoform],ccle.drug.tissuetype[names(unique),"tissue.type"])
    colnames(unique) <- c("ccle.AUC", "gene.exp","isoform.exp","tissue.type")
    
    intersected.celllines <- intersect(names(ccle.sensitivity),names(gdsc.sensitivity))
    ccle.sensitivity <- ccle.sensitivity[intersected.celllines]
    gdsc.sensitivity <- gdsc.sensitivity[intersected.celllines]
    
    
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    sensitivity <- ccle.sensitivity
    sensitivity <- cbind(sensitivity, gdsc.sensitivity, ccle.genes.fpkm[names(sensitivity),gene], ccle.isoforms.fpkm[names(sensitivity),isoform],ccle.drug.tissuetype[names(sensitivity),"tissue.type"])
    colnames(sensitivity) <- c("ccle.AUC","gdsc.AUC", "gene.exp","isoform.exp","tissue.type")
    
    sensitivity <- sensitivity[order(sensitivity[,"gdsc.AUC"]),]
    
    #my.xlim <- as.numeric(range(sensitivity[,"isoform.exp"]))
    #my.ylim <- as.numeric(range(sensitivity[,"ccle.AUC"]))
    my.xlim <- c(min(min(as.numeric(sensitivity[,"isoform.exp"])), min(as.numeric(unique[,"isoform.exp"]))),max(max(as.numeric(sensitivity[,"isoform.exp"])), max(as.numeric(unique[,"isoform.exp"]))))
    my.ylim <- c(min(min(as.numeric(sensitivity[,"ccle.AUC"])), min(as.numeric(unique[,"ccle.AUC"]))),max(max(as.numeric(sensitivity[,"ccle.AUC"])), max(as.numeric(unique[,"ccle.AUC"]))))
    
    
    sensitivity <- cbind(sensitivity, "red")
    colnames(sensitivity)[ncol(sensitivity)] <- "col"
    sensitivity[,"col"] <- colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity))
    
    
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="ccle.AUC",cex.lab = .7, cex.axis = .7)
    points(unique[,"isoform.exp"], unique[,"ccle.AUC"], pch = 17,col = "grey")
    points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"], pch = 19,col = sensitivity[,"col"])
    
    x <- c(round(my.xlim[1], digits = 1) + 0.1 , round(my.xlim[1], digits = 1) + 0.3, round(my.xlim[1], digits = 1) + 0.3,round(my.xlim[1], digits = 1) + 0.1)
    y <- c(round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1), round(my.ylim[2], digits = 1))
    legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "gdsc.AUC", limits = c("resistent","sensitive"),cex=.6)
    
    si <- (nrow(sensitivity) - 2):nrow(sensitivity)
    #textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
    ri <- 1:3
    #textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
    textxy(x[2], y[3] - .08, sprintf("pvalue = %s\n    R2 = %s", round(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]][["M0","M3B"]], digits = 20), round(ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]][["median","M3B"]], digits = 2)), cex=0.5)
    title(sprintf("Correlation between %s(%s) and %s \n in all tissue types",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene), "Symbol"], gsub("drugid_","",drug)),cex.main=.7)
  } , 
  "all.res"= {
    
    Genedf0 <- data.frame(ccle.isoforms.fpkm[,isoform])
    colnames(Genedf0) <- "isoform"
    Genedf0[,"Tissue"] <- factor(ccle.drug.tissuetype[,"tissue.type"], ordered =FALSE)
    Genedf0 <- Genedf0[complete.cases(Genedf0),]
    Genedf0$Tissue <- factor(Genedf0$Tissue, ordered =FALSE)
    linearFormula <- "isoform ~ Tissue"
    M3 <- lm (linearFormula, Genedf0)
    
    unique.ccle.celllines <- setdiff(names(ccle.sensitivity),names(gdsc.sensitivity))
    ccle.unique <- ccle.sensitivity[unique.ccle.celllines]
    unique <-  matrix(NA,nrow= length(ccle.unique))
    unique <- ccle.unique
    unique <- cbind(unique, ccle.genes.fpkm[names(unique),gene], ccle.isoforms.fpkm[names(unique),isoform],ccle.drug.tissuetype[names(unique),"tissue.type"],M3$residuals[names(unique)])
    colnames(unique) <- c("ccle.AUC", "gene.exp","isoform.exp","tissue.type","M3.res")
    
    intersected.celllines <- intersect(names(ccle.sensitivity),names(gdsc.sensitivity))
    ccle.sensitivity <- ccle.sensitivity[intersected.celllines]
    gdsc.sensitivity <- gdsc.sensitivity[intersected.celllines]
    
    
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    sensitivity <- ccle.sensitivity
    sensitivity <- cbind(sensitivity, gdsc.sensitivity, ccle.genes.fpkm[names(sensitivity),gene], ccle.isoforms.fpkm[names(sensitivity),isoform],ccle.drug.tissuetype[names(sensitivity),"tissue.type"],M3$residuals[names(sensitivity)])
    colnames(sensitivity) <- c("ccle.AUC","gdsc.AUC", "gene.exp","isoform.exp","tissue.type","M3.res")
    
    sensitivity <- sensitivity[order(sensitivity[,"gdsc.AUC"]),]
    
    #my.xlim <- as.numeric(range(sensitivity[,"isoform.exp"]))
    #my.ylim <- as.numeric(range(sensitivity[,"ccle.AUC"]))
    my.xlim <- c(min(min(as.numeric(sensitivity[,"M3.res"])), min(as.numeric(unique[,"M3.res"]))),max(max(as.numeric(sensitivity[,"M3.res"])), max(as.numeric(unique[,"M3.res"]))))
    my.ylim <- c(min(min(as.numeric(sensitivity[,"ccle.AUC"])), min(as.numeric(unique[,"ccle.AUC"]))),max(max(as.numeric(sensitivity[,"ccle.AUC"])), max(as.numeric(unique[,"ccle.AUC"]))))
    
    
    sensitivity <- cbind(sensitivity, "red")
    colnames(sensitivity)[ncol(sensitivity)] <- "col"
    sensitivity[,"col"] <- colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity))
    
    
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="ccle.AUC",cex.lab = .7, cex.axis = .7)
    points(unique[,"M3.res"], unique[,"ccle.AUC"], pch = 17,col = "grey")
    points(sensitivity[,"M3.res"], sensitivity[,"ccle.AUC"], pch = 19,col = sensitivity[,"col"])
    
    x <- c(round(my.xlim[1], digits = 1)  , round(my.xlim[1], digits = 1) + .3, round(my.xlim[1], digits = 1) + .3,round(my.xlim[1], digits = 1) )
    y <- c(round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1), round(my.ylim[2], digits = 1))
    legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "gdsc.AUC", limits = c("resistent","sensitive"),cex=.6)
    
    si <- (nrow(sensitivity) - 2):nrow(sensitivity)
    #textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
    ri <- 1:3
    #textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
    textxy(x[2] + 2, y[3] - .08, sprintf("pvalue = %s\n    R2 = %s", round(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]][["M0","M3B"]], digits = 20), round(ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]][["median","M3B"]], digits = 2)), cex=0.5)
    title(sprintf("Correlation between %s(%s) and %s \n in all tissue types",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene), "Symbol"], gsub("drugid_","",drug)),cex.main=.7)
  } ,
  "cor"= {   
    intersected.celllines <- intersect(names(ccle.sensitivity),names(gdsc.sensitivity))
    ccle.sensitivity <- ccle.sensitivity[intersected.celllines]
    gdsc.sensitivity <- gdsc.sensitivity[intersected.celllines]
    
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    sensitivity <- ccle.sensitivity
    sensitivity <- cbind(sensitivity, gdsc.sensitivity)
    colnames(sensitivity) <- c("ccle.AUC","gdsc.AUC")
    my.xlim <- c(min(as.numeric(sensitivity[,"ccle.AUC"])), max(as.numeric(sensitivity[,"ccle.AUC"])))
    my.ylim <- c(min(as.numeric(sensitivity[,"gdsc.AUC"])), max(as.numeric(sensitivity[,"gdsc.AUC"])))
    
    pdf(file = file.path(path.diagrams, sprintf("ccle_gdsc_%s.pdf", gsub("drugid_","",drug))), height=7, width=7)
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab="ccle.AUC", ylab="gdsc.AUC",main = sprintf("Correlation between gdsc and ccle for %s is %s", gsub("drugid_","",drug),round(cor(sensitivity[,"ccle.AUC"], sensitivity[,"gdsc.AUC"],method="spearman"),digits = 2)), cex.main = .9)
    points(sensitivity[,"ccle.AUC"], sensitivity[,"gdsc.AUC"], pch = 19,col = "blue")
    dev.off()
  })
  return(sensitivity)
}
fnCCLE.tissue <- function(gene, isoform, drug, method = c("all","tissue.based"), tissue.type= "", celllines.col = "", R2 = "", pvalue = "")
{
  ccle.sensitivity <- ccle.drug.sensitivity[complete.cases(ccle.drug.sensitivity[,drug]),drug]
  
  sensitivity <- NULL
  switch(method , "tissue.based" = {
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    sensitivity <- ccle.sensitivity
    sensitivity <- cbind(sensitivity, ccle.genes.fpkm[names(sensitivity),gene], ccle.isoforms.fpkm[names(sensitivity),isoform],ccle.drug.tissuetype[names(sensitivity),"tissue.type"])
    colnames(sensitivity) <- c("ccle.AUC","gene.exp","isoform.exp","tissue.type")
    sensitivity <- subset(sensitivity, sensitivity[,"tissue.type"] == tissue.type)
    
    if(nrow(sensitivity) > 0)
    {
      my.xlim <- range(as.numeric(sensitivity[,"isoform.exp"]))
      my.ylim <- range(as.numeric(sensitivity[,"ccle.AUC"]))
      
      sensitivity <- cbind(sensitivity, "red")
      colnames(sensitivity)[ncol(sensitivity)] <- "col"
      sensitivity[,"col"] <- celllines.col[rownames(sensitivity)]#colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity))
      
      plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="ccle.AUC",cex.lab = .7, cex.axis = .7)
      points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"], pch = 17,col = "grey")
      #points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"], pch = 19,col = sensitivity[,"col"])
      
      x <- c(round(my.xlim[1], digits = 1) + 0.1 , round(my.xlim[1], digits = 1) + 0.3, round(my.xlim[1], digits = 1) + 0.3,round(my.xlim[1], digits = 1) + 0.1)
      y <- c(round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1) - 0.1, round(my.ylim[2], digits = 1), round(my.ylim[2], digits = 1))
      #legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "gdsc.AUC", limits = c("resistent","sensitive"),cex=.7)
      
      si <- (nrow(sensitivity) - 2):nrow(sensitivity)
      #textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
      ri <- 1:3
      #textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
      textxy(my.xlim[1] , my.ylim[2] - .04, sprintf("pvalue = %s\n    R2 = %s", round(pvalue, digits = 20), round(R2, digits = 2)), cex=0.6)
      title(sprintf("%s ",capitalize(gsub("_"," ",tissue.type))),cex.main=.8)
    }
  }, 
  "all"= {
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    sensitivity <- ccle.sensitivity
    sensitivity <- cbind(sensitivity, ccle.genes.fpkm[names(sensitivity),gene], ccle.isoforms.fpkm[names(sensitivity),isoform],ccle.drug.tissuetype[names(sensitivity),"tissue.type"])
    colnames(sensitivity) <- c("ccle.AUC","gene.exp","isoform.exp","tissue.type")
    sensitivity <- sensitivity[order(sensitivity[,"ccle.AUC"]),]
    
    my.xlim <- range(as.numeric(sensitivity[,"isoform.exp"]))
    my.ylim <- range(as.numeric(sensitivity[,"ccle.AUC"]))
    
    sensitivity <- cbind(sensitivity, "red")
    colnames(sensitivity)[ncol(sensitivity)] <- "col"
    sensitivity[,"col"] <- colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity))
    
    
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="AUC",cex.lab = .7, cex.axis = .7)
    points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"], pch = 17,col = "grey")
    
    x <- c(round(my.xlim[1], digits = 1) + 0.1 , round(my.xlim[1], digits = 1) + 0.3, round(my.xlim[1], digits = 1) + 0.3,round(my.xlim[1], digits = 1) + 0.1)
    y <- c(round(my.ylim[2], digits = 1) - 0.2, round(my.ylim[2], digits = 1) - 0.2, round(my.ylim[2], digits = 1)-0.1, round(my.ylim[2], digits = 1)-.1)
    #legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "", limits = c("resistent","sensitive"),cex=.6)
    
    si <- (nrow(sensitivity) - 2):nrow(sensitivity)
    #textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
    ri = 1:3
    #textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
    textxy(my.xlim[1] , my.ylim[2] - .04, sprintf("pvalue = %s\n    R2 = %s", round(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]][["M0","M3B"]], digits = 20), round(ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]][["median","M3B"]], digits = 2)), cex=0.6)
    title(sprintf("Correlation between %s(%s) and %s \n in all tissue types",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene), "Symbol"], gsub("drugid_","",drug)),cex.main=.7)
  },
  "all.res"={
    Genedf0 <- data.frame(ccle.isoforms.fpkm[,isoform])
    colnames(Genedf0) <- "isoform"
    Genedf0[,"Tissue"] <- factor(ccle.drug.tissuetype[,"tissue.type"], ordered =FALSE)
    Genedf0 <- Genedf0[complete.cases(Genedf0),]
    Genedf0$Tissue <- factor(Genedf0$Tissue, ordered =FALSE)
    linearFormula <- "isoform ~ Tissue"
    M3 <- lm (linearFormula, Genedf0)
    
    sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
    sensitivity <- ccle.sensitivity
    sensitivity <- cbind(sensitivity, ccle.genes.fpkm[names(sensitivity),gene], ccle.isoforms.fpkm[names(sensitivity),isoform],ccle.drug.tissuetype[names(sensitivity),"tissue.type"], M3$residuals[names(sensitivity)])
    colnames(sensitivity) <- c("ccle.AUC","gene.exp","isoform.exp","tissue.type","M3.res")
    sensitivity <- sensitivity[order(sensitivity[,"M3.res"]),]
    
    my.xlim <- range(as.numeric(sensitivity[,"M3.res"]))
    my.ylim <- range(as.numeric(sensitivity[,"ccle.AUC"]))
    
    sensitivity <- cbind(sensitivity, "red")
    colnames(sensitivity)[ncol(sensitivity)] <- "col"
    sensitivity[,"col"] <- colorRampPalette(c("blue" , "light blue","red"))(nrow(sensitivity))
    
    
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="AUC",cex.lab = .7, cex.axis = .7)
    points(sensitivity[,"M3.res"], sensitivity[,"ccle.AUC"], pch = 17,col = "grey")
    
    x <- c(round(my.xlim[1], digits = 1) + 0.1 , round(my.xlim[1], digits = 1) + 0.3, round(my.xlim[1], digits = 1) + 0.3,round(my.xlim[1], digits = 1) + 0.1)
    y <- c(round(my.ylim[2], digits = 1) - 0.2, round(my.ylim[2], digits = 1) - 0.2, round(my.ylim[2], digits = 1)-0.1, round(my.ylim[2], digits = 1)-.1)
    #legend.gradient(cbind(x, y), cols = colorRampPalette(c("blue" , "light blue","red"))( nrow(sensitivity) ), title = "", limits = c("resistent","sensitive"),cex=.6)
    
    si <- (nrow(sensitivity) - 2):nrow(sensitivity)
    #textxy(as.numeric(sensitivity[si,"isoform.exp"]),as.numeric(sensitivity[si,"ccle.AUC"]),rownames(sensitivity)[si], cex=0.7)
    ri <- 1:3
    #textxy(as.numeric(sensitivity[ri,"isoform.exp"]),as.numeric(sensitivity[ri,"ccle.AUC"]),rownames(sensitivity)[ri], cex=0.7)
    textxy(my.xlim[1] + 2.5, my.ylim[2] - .04, sprintf("pvalue = %s\nR2 = %s", round(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]][["M0","M3B"]], digits = 20), round(ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]][["median","M3B"]], digits = 2)), cex=0.6)
    title(sprintf("Correlation between %s(%s) and %s \n in all tissue types",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"], gsub("drugid_","",drug)),cex.main=.7)
    
    
  })
  return(sensitivity)
}
fnCCLE <- function(gene, isoform, drug, col = c("static","tissue.based"))
{
  ccle.sensitivity <- ccle.drug.sensitivity[complete.cases(ccle.drug.sensitivity[,drug]),drug]
  sensitivity <- NULL
  
  ccle.sensitivity <- ccle.sensitivity[order(ccle.sensitivity)] 
  
  sensitivity <- matrix(NA,nrow= length(ccle.sensitivity))
  sensitivity[,1] <- ccle.sensitivity
  rownames(sensitivity) <- names(ccle.sensitivity)
  sensitivity <- cbind(sensitivity, ccle.genes.fpkm[rownames(sensitivity),gene], ccle.isoforms.fpkm[rownames(sensitivity),isoform],ccle.drug.tissuetype[rownames(sensitivity),"tissue.type"])
  colnames(sensitivity) <- c("ccle.AUC", "gene.exp","isoform.exp","tissue.type")
  col.tissuetypes <- levels(factor(sensitivity[,"tissue.type"]))
  col.tissuetypes <- cbind(col.tissuetypes, colorRampPalette(brewer.pal(n=11, name= 'Spectral'))(length(col.tissuetypes)))
  colnames(col.tissuetypes) <- c("tissue.type","col")
  rownames(col.tissuetypes) <- col.tissuetypes[,"tissue.type"]
  sensitivity <- cbind(sensitivity, "red")
  colnames(sensitivity)[ncol(sensitivity)] <- "col"
  for(i in 1:nrow(col.tissuetypes))
  {
    sensitivity[which(sensitivity[,"tissue.type"]==col.tissuetypes[i,"tissue.type"]),"col"] = col.tissuetypes[i,"col"]
  }
  my.xlim <- range(as.numeric(sensitivity[,"isoform.exp"]))
  my.ylim <- range(as.numeric(sensitivity[,"ccle.AUC"]))
  
  
  switch(col, "static" = {
    pdf(file = file.path(path.diagrams, sprintf("ccle_%s_%s_%s.pdf",gene, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"],gsub("drugid_","",drug))), height=7, width=7)
    
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="AUC", main = sprintf("Correlation between %s(%s) and %s \n in CCLE experiments",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"], gsub("drugid_","",drug)),cex=.9)
    points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"],  pch = 19,col = "blue")
    dev.off()
  }, 
  "tissue.based" = {
    
    pdf(file = file.path(path.diagrams, sprintf("ccle_%s_%s_%s.pdf",gene, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"],gsub("drugid_","",drug))), height=7, width=7)
    
    plot(NA, xlim = my.xlim, ylim = my.ylim, xlab=isoform, ylab="AUC", main = sprintf("Correlation between %s(%s) and %s \n in ccle cell lines based on tissue types",isoform, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"], gsub("drugid_","",drug)),cex=.9)
    points(sensitivity[,"isoform.exp"], sensitivity[,"ccle.AUC"],  pch = 19,col = sensitivity[,"col"])
    x <- c(round(my.xlim[1], digits = 1) , round(my.xlim[1], digits = 1) + 1, round(my.xlim[1], digits = 1) + 1,round(my.xlim[1], digits = 1) )
    y <- c(round(my.ylim[2], digits = 1) - 0.3, round(my.ylim[2], digits = 1) - 0.3, round(my.ylim[2], digits = 1), round(my.ylim[2], digits = 1))
    legend.gradient(cbind(x, y), cols =  col.tissuetypes[,"col"], title = "tissue types", limits = "",cex=.8)
    yt <- y[3]-0.009
    for(i in 1:nrow(col.tissuetypes))
    {
      textxy(x[2], yt, col.tissuetypes[i,"tissue.type"] , cex=0.5)
      yt <- yt -0.0145
    }
  })
  return(list(sensitivity = sensitivity , col.tissuetypes = col.tissuetypes))  
}
fnPlotTissuetypes <- function(drug, sensitivity, col.tissuetypes)
{
  pdf(file = file.path(path.diagrams, sprintf("tissue.type_%s.pdf", gsub("drugid_","",drug))), height=4, width=8)
  mp <- barplot(table(sensitivity[,"tissue.type"]), beside = TRUE,  space = .6, col = col.tissuetypes[,"col"], ylab = "Number of cell lines", main = sprintf("The number of cell lines with each tissue type (%s)",gsub("drugid_","",drug)),border=NA,xaxt='n',cex.main = .7, cex.lab = .7, cex.axis = .7)
  text(mp, par("usr")[3], labels = capitalize(gsub("_"," ",col.tissuetypes[,"tissue.type"])), srt = 30, adj = c(1,1), xpd = TRUE, cex=.6)
  dev.off()
  
}
fnCheckIsoform <- function(gene, isoform, drug, ccle=FALSE, annot.ensembl.all.genes)
{
  fnsensitivity.tissues <- fnSensitivity.tissue (gene, isoform, drug)
  CCLE <- fnCCLE(gene, isoform, drug, col = "static")
  sensitivity <- CCLE$sensitivity
  fnPlotTissuetypes(drug, sensitivity, CCLE$col.tissuetypes)
  
  n <- table(sensitivity[,"tissue.type"])
  tissue.types <- names(fnsensitivity.tissues$tissues)
  sensitivity.tissues <- matrix(NA,nrow= (length(tissue.types) + 1), ncol=5)
  rownames(sensitivity.tissues) <- c("all",tissue.types)
  colnames(sensitivity.tissues) <- c("log.pvalue","col","rsquared","n","significance")
  
  names.tissues.colors <- colorRampPalette(brewer.pal(n=11, name= 'Spectral'))(length(tissue.types))
  #c("blue", "green", "Blues", "RdPu", "Greens", "PuRd", "Greys", "Oranges", "PuBuGn", "Purples", "PuBu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", "BuGn", "BuPu", "GnBu","OrRd")
  ranges <- c(1,0.75,.7,.65,.60,0)
  col.tissuetypes.ranges <- matrix(NA,nrow = length(tissue.types), ncol = length(ranges))
  rownames(col.tissuetypes.ranges) <- tissue.types
  colnames(col.tissuetypes.ranges) <- ranges
  for(i in 1:nrow(col.tissuetypes.ranges))
  {
    col.tissuetypes.ranges[i,] <-  colorRampPalette(c("white",names.tissues.colors[i]))(length(ranges))
  }
  sensitivity.tissues[1,"rsquared"] <- ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]]["median","M3B"]
  sensitivity.tissues[1,"log.pvalue"] <- -log10(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]]["M0","M3B"])
  Genedf0 <- data.frame(ccle.drug.sensitivity[,which(colnames(ccle.drug.sensitivity) == drug)])
  Genedf0[,"Tissue"] <- factor(ccle.drug.tissuetype[,"tissue.type"], ordered =FALSE)
  sensitivity.tissues[1,"n"] <- nrow(Genedf0[complete.cases(Genedf0),])
  sensitivity.tissues[1,"col"] <- "grey"
  sensitivity.tissues[1,"significance"] <- " "
  for (i in 1:length(tissue.types))
  {
    sensitivity.tissues[i+1,"rsquared"] <- fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] * sign(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M3B","M0"])
    sensitivity.tissues[i+1,"log.pvalue"] <- -log10(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"])
    sensitivity.tissues[i+1,"n"] <- n[tissue.types[i]]
    if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .001)
    {sensitivity.tissues[i+1,"significance"] <- "***"} else if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .01)
    {sensitivity.tissues[i+1,"significance"] <- "**"} else if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .05)                           
    {sensitivity.tissues[i+1,"significance"] <- "*"} else if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .1)                           
    {sensitivity.tissues[i+1,"significance"] <- "-"} else  {sensitivity.tissues[i+1,"significance"] <- " "}          
    j <- ncol(col.tissuetypes.ranges)
    if(fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] <0 ){fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] <- 0}
    while(fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] < as.numeric(colnames(col.tissuetypes.ranges)[j]))
    {
      j <- j - 1
    }
    sensitivity.tissues[i+1,"col"] <- col.tissuetypes.ranges[tissue.types[i],j]
  }
  
  pdf(file = file.path(path.diagrams, sprintf("tissue.type_R2_%s_%s_%s.pdf", gene, annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"], gsub("drugid_","",drug))), height=4, width=7)
  my.ylim <- range(as.numeric(sensitivity.tissues[,"rsquared"]))
  my.ylim <- c(-1, +1)
  mp <- barplot(as.numeric(sensitivity.tissues[,"rsquared"]), beside = TRUE, ylim = my.ylim, width = .8, space = 0.6, xlim = c(0,25), col =sensitivity.tissues[,"col"], ylab = "R2", main = sprintf("Sensitivity of drug: %s for isoform: %s in different tissue types",gsub("drugid_","",drug),isoform),border=NA,xaxt='n', cex.main = .7, cex.lab = .7, cex.axis = .7)
  text(mp, 0, labels = capitalize(gsub("_", " ", rownames(sensitivity.tissues))), srt = 30, adj = c(1,1), xpd = TRUE, cex=.5)  
  text(mp + .5, as.numeric(sensitivity.tissues[,"rsquared"]) +.03, labels = sensitivity.tissues[,"significance"], adj = c(1,1),cex=1)
  text(mp +.5 , .06, labels = sprintf("n = %s", sensitivity.tissues[,"n"]), adj = c(1,1),cex=.4)
  dev.off()
  
  if(drug %in% colnames(gdsc.drug.sensitivity) & !ccle)
  {
    fnCCLE.GDSC(gene, isoform, drug, method = "cor", "","", annot.ensembl.all.genes)
    
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_%s_cellLines_all.pdf",gene,annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"],gsub("drugid_","",drug))), height=24, width=16)
    par(mfrow=c(6, 4))
    sensitivity.all <- fnCCLE.GDSC(gene, isoform, drug, method = "all.res", annot.ensembl.all.genes)
    for (tissue.type in tissue.types)
    {
      fnCCLE.GDSC(gene, isoform, drug, method="tissue.based", tissue.type, sensitivity.all[,"col"], R2=fnsensitivity.tissues$statistics[[tissue.type]]["median","M3B"], pvalue=fnsensitivity.tissues$p.values[[tissue.type]]["M0","M3B"], annot.ensembl.all.genes)
    }
    dev.off() 
    
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_%s_cellLines_significant.pdf",gene,annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"],gsub("drugid_","",drug))), height=8, width=16)
    par(mfrow=c(2, 4))
    sensitivity.all = fnCCLE.GDSC(gene, isoform, drug, method = "all.res", "", annot.ensembl.all.genes)
    for (i in 1:nrow(sensitivity.tissues))
    {
      if ((sensitivity.tissues[i,"significance"] == "***") | (sensitivity.tissues[i,"significance"] == "**"))
      {
        fnCCLE.GDSC(gene, isoform, drug, method = "tissue.based", rownames(sensitivity.tissues)[i], sensitivity.all[,"col"], R2=fnsensitivity.tissues$statistics[[rownames(sensitivity.tissues)[i]]]["median","M3B"], pvalue=fnsensitivity.tissues$p.values[[rownames(sensitivity.tissues)[i]]]["M0","M3B"], annot.ensembl.all.genes)
      }
    }
    dev.off() 
  }
  else
  {
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_%s_cellLines_all.pdf",gene,annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"],gsub("drugid_","",drug))), height=24, width=16)
    par(mfrow=c(6, 4))    
    sensitivity.all <- fnCCLE.tissue(gene, isoform, drug, method = "all.res", "","")
    for (tissue.type in tissue.types)
    {
      fnCCLE.tissue(gene, isoform, drug, method = "tissue.based", tissue.type, sensitivity.all[,"col"], R2 = fnsensitivity.tissues$statistics[[tissue.type]]["median","M3B"], pvalue = fnsensitivity.tissues$p.values[[tissue.type]]["M0","M3B"])
    }
    dev.off() 
    
    pdf(file = file.path(path.diagrams, sprintf("%s_%s_%s_cellLines_significant.pdf",gene,annot.ensembl.all.genes[which(annot.ensembl.all.genes$EnsemblGeneId == gene),"Symbol"],gsub("drugid_","",drug))), height=8, width=16)
    par(mfrow=c(2, 4))
    sensitivity.all = fnCCLE.tissue(gene, isoform, drug, method = "all.res", "")
    for (i in 1:nrow(sensitivity.tissues))
    {
      if ((sensitivity.tissues[i,"significance"] == "***") | (sensitivity.tissues[i,"significance"] == "**"))
      {
        fnCCLE.tissue(gene, isoform, drug, method = "tissue.based", rownames(sensitivity.tissues)[i], sensitivity.all[,"col"], R2 = fnsensitivity.tissues$statistics[[rownames(sensitivity.tissues)[i]]]["median","M3B"], pvalue = fnsensitivity.tissues$p.values[[rownames(sensitivity.tissues)[i]]]["M0","M3B"])
      }
    }
    dev.off() 
  }
}
fnCheckIsoform.check <- function(gene, isoform, drug, ccle=FALSE, limit=4, annot.ensembl.all.genes)
{
  fnsensitivity.tissues <- fnSensitivity.tissue (gene, isoform, drug)
  
  sensitivity <- CCLE$sensitivity
  n <- table(sensitivity[,"tissue.type"])
  tissue.types <- names(fnsensitivity.tissues$tissues)
  sensitivity.tissues <- matrix(NA,nrow= (length(tissue.types) + 1), ncol=5)
  rownames(sensitivity.tissues) <- c("all",tissue.types)
  colnames(sensitivity.tissues) <- c("log.pvalue","col","rsquared","n","significance")
  
  names.tissues.colors <- colorRampPalette(brewer.pal(n=11, name= 'Spectral'))(length(tissue.types))
  #c("blue", "green", "Blues", "RdPu", "Greens", "PuRd", "Greys", "Oranges", "PuBuGn", "Purples", "PuBu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd", "BuGn", "BuPu", "GnBu","OrRd")
  ranges <- c(1,0.75,.7,.65,.60,0)
  col.tissuetypes.ranges <- matrix(NA,nrow = length(tissue.types), ncol = length(ranges))
  rownames(col.tissuetypes.ranges) <- tissue.types
  colnames(col.tissuetypes.ranges) <- ranges
  for(i in 1:nrow(col.tissuetypes.ranges))
  {
    col.tissuetypes.ranges[i,] <-  colorRampPalette(c("white",names.tissues.colors[i]))(length(ranges))
  }
  sensitivity.tissues[1,"rsquared"] <- ccle.drug.association.statistics[[gene]][[gsub("drugid_","",drug)]]["median","M3B"]
  sensitivity.tissues[1,"log.pvalue"] <- -log10(ccle.drug.association[[gene]][[gsub("drugid_","",drug)]]["M0","M3B"])
  Genedf0 <- data.frame(ccle.drug.sensitivity[,which(colnames(ccle.drug.sensitivity) == drug)])
  Genedf0[,"Tissue"] <- factor(ccle.drug.tissuetype[,"tissue.type"], ordered =FALSE)
  sensitivity.tissues[1,"n"] <- nrow(Genedf0[complete.cases(Genedf0),])
  sensitivity.tissues[1,"col"] <- "grey"
  sensitivity.tissues[1,"significance"] <- " "
  for (i in 1:length(tissue.types))
  {
    sensitivity.tissues[i+1,"rsquared"] <- fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] * sign(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M3B","M0"])
    sensitivity.tissues[i+1,"log.pvalue"] <- -log10(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"])
    sensitivity.tissues[i+1,"n"] <- n[tissue.types[i]]
    if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .001)
    {sensitivity.tissues[i+1,"significance"] <- "***"} else if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .01)
    {sensitivity.tissues[i+1,"significance"] <- "**"} else if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .05)                           
    {sensitivity.tissues[i+1,"significance"] <- "*"} else if(fnsensitivity.tissues$p.values[[tissue.types[i]]]["M0","M3B"] < .1)                           
    {sensitivity.tissues[i+1,"significance"] <- "-"} else  {sensitivity.tissues[i+1,"significance"] <- " "}          
    j <- ncol(col.tissuetypes.ranges)
    if(fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] <0 ){fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] <- 0}
    while(fnsensitivity.tissues$statistics[[tissue.types[i]]]["median","M3B"] < as.numeric(colnames(col.tissuetypes.ranges)[j]))
    {
      j <- j - 1
    }
    sensitivity.tissues[i+1,"col"] <- col.tissuetypes.ranges[tissue.types[i],j]
  }
  
  significant.tissues <- NULL
  if(drug %in% colnames(gdsc.drug.sensitivity) & !ccle)
  {
    sensitivity.all <- fnCCLE.GDSC(gene, isoform, drug, method="all", "", annot.ensembl.all.genes)
    for (i in 2:nrow(sensitivity.tissues))
    {
      if ((sensitivity.tissues[i,"significance"] == "***") | (sensitivity.tissues[i,"significance"] == "**"))
      {
        sensitivity.tissue <- fnCCLE.GDSC(gene, isoform, drug, method="tissue.based", rownames(sensitivity.tissues)[i], sensitivity.all[,"col"], annot.ensembl.all.genes)
        ranks <- which(rownames(sensitivity.all) %in% rownames(sensitivity.tissue))
        if ((length(ranks) > (limit*2)) && all(ranks[1:limit] < nrow(sensitivity.all)%/%2) && all(ranks[(length(ranks)-limit+1):length(ranks)] > nrow(sensitivity.all)%/%2))
        {
          significant.tissues <- c(significant.tissues, rownames(sensitivity.tissues)[i])
          fnCheckIsoform(gene, isoform, drug, ccle, annot.ensembl.all.genes)
        }
      }
    }
  }
  else
  {
    sensitivity.all <- fnCCLE.tissue(gene, isoform, drug, method = "all", "")
    for (i in 1:nrow(sensitivity.tissues))
    {
      if ((sensitivity.tissues[i,"significance"] == "***") | (sensitivity.tissues[i,"significance"] == "**"))
      {
        sensitivity.tissue <- fnCCLE.tissue(gene, isoform, drug, method = "tissue.based", rownames(sensitivity.tissues)[i], sensitivity.all[,"col"])
        ranks <- which(rownames(sensitivity.all) %in% rownames(sensitivity.tissue))
        if ((length(ranks) > (limit*2)) && all(ranks[1:limit] < nrow(sensitivity.all)%/%2) && all(ranks[(length(ranks)-limit+1):length(ranks)] > nrow(sensitivity.all)%/%2))
        {
          significant.tissues <- c(significant.tissues, rownames(sensitivity.tissues)[i])
          fnCheckIsoform(gene, isoform, drug, ccle, annot.ensembl.all.genes)
        }
      }
    }
  }
  return (significant.tissues)
}
fnPlotTopBiomarkers <- function(top.biomarkers, top.no)
{
  drug.col <- c("#666666","#996666","#CC3399","#663399","#66CC00","#336600","#3333CC","#FF3399","#999999","#333333", "#0066CC","#669999","#33CC99","#339900","#CCFF99")
  top.list <- NULL;
  for(i in 1:length(top.biomarkers))
  {
    for(j in 1:top.no)
    {
      top.list <- rbind(top.list, c(names(top.biomarkers)[i],drug.col[i], top.biomarkers[[i]][j,"symbol"],top.biomarkers[[i]][j,"biomarker.id"], -log10(top.biomarkers[[i]][j,adjustment.method]), sign(top.biomarkers[[i]][j,"estimate"])* top.biomarkers[[i]][j,"R2"],top.biomarkers[[i]][j,"type"],top.biomarkers[[i]][j,"specificity"] ))
    }
  }
  colnames(top.list) <- c("drug.id","col","symbol","biomarker.id","pvalue","efficacy","type","specificity")
  
  my.ylim <- c(0,15)
  pdf(file = file.path(path.diagrams, sprintf("top_%s_list.pdf",top.no)), height=7, width=7)  
  boxplot(as.numeric(top.list[,"pvalue"]) ~top.list[,"specificity"])
  points(top.list[,"class"],jitter(as.numeric(top.list[,"pvalue"])), pch = 20, col = top.list[,"col"], ji)
  
  dev.off()
}

