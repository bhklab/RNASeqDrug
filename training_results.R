options(stringsAsFactors=FALSE)

if(!require(stringr)){install.packages("stringr");library(stringr)}
if(!require(corrplot)){install.packages("corrplot");library(corrplot)}
if(!require(VennDiagram)){install.packages("VennDiagram");library(VennDiagram)}
if(!require(SDMTools)){install.packages("SDMTools");library(SDMTools)}
if(!require(calibrate)){install.packages("calibrate");library(calibrate)}
if(!require(Hmisc)){install.packages("Hmisc");library(Hmisc)}
if(!require(grid)){install.packages("grid");library(grid)}
if(!require(lattice)){install.packages("lattice");library(lattice)}
if(!require(WriteXLS)){install.packages("WriteXLS");library(WriteXLS)}
if(!require(magicaxis)){install.packages("magicaxis");library(magicaxis)}
if(!require(np)){install.packages("np");library(np)}
if(!require(ggplot2)){install.packages("ggplot2");library(ggplot2)}
if(!require(gridBase)){install.packages("gridBase");library(gridBase)}
if(!require(gridExtra)){install.packages("gridExtra");library(gridExtra)}
if(!require(easyGgplot2)){install.packages("devtools");library(devtools);install_github("kassambara/easyGgplot2");library(easyGgplot2)}
source("https://bioconductor.org/biocLite.R")
if(!require("PharmacoGx")){biocLite("PharmacoGx");library(PharmacoGx)}
if(!require("Biobase")){biocLite("Biobase");library(Biobase)}

path.data <- file.path("../data")
path.code <- file.path(".")
path.result <- file.path("../results")

args <- commandArgs(trailingOnly=TRUE)
##for test
#args <- c("auc_recomputed_drug_association.RData", "0.55", "0.01", "breast", "glm", "gaussian", "cindex", "fdr", "auc_recomputed_ccle_gdsc", "training_ccle_gdsc.RData")

file.sensitivity <- file.path(path.data, as.character(args[1]))
effect.size.cut.off <- as.numeric(args[2])
fdr.cut.off <- as.numeric(args[3])
tissue <- as.character(args[4])
model.method <- as.character(args[5])
glm.family <- as.character(args[6]) 
effect.size <- as.character(args[7]) #c("r.squared", "cindex")
adjustment.method <- as.character(args[8])#c("bonferroni", "fdr")

breast.specific <- ifelse(regexpr("breast", file.sensitivity) < 1, FALSE, TRUE)
data.type <- ifelse(regexpr("mut", file.sensitivity) >= 1, "all", "expression")

source(file.path(path.code, "foo.R"))

path.diagrams <- file.path(path.result, as.character(args[9]))
if (!file.exists(path.diagrams)){dir.create(file.path(path.diagrams))}

path.training.data <- file.path(path.data, as.character(args[10]))
load(path.training.data, verbose=TRUE)
if("gdsc.drug.sensitivity" %in% ls()) {
  training.type <-"CCLE_GDSC"
} else {
  training.type <-"CCLE"
}

if(data.type == "all") {
  Models <- c("M3B", "M2", "MC", "MM")
  Models.names <- c("Isoforms", "Genes", "Copy Number Variations", "Mutations")
}else {
  Models <- c("M2", "M3B")
  Models.names <- c("Genes", "Isoforms")
}
names(Models.names) <- Models

ccle.tissuetype <- tissueTypes
ccle.tissuetype[,1] <- as.character(ccle.tissuetype[,1])
ccle.primary.tissuetype <- ccle.tissuetype
ccle.primary.drug.sensitivity <- ccle.drug.sensitivity

mean.ccle.isoforms.fpkm <- colMeans(ccle.isoforms.fpkm, na.rm=TRUE)
if("haematopoietic_and_lymphoid_tissue" %in% ccle.tissuetype[,1])
{
  ccle.tissuetype[which(ccle.tissuetype[,1] == "haematopoietic_and_lymphoid_tissue"), 1] <- "haematopoietic_and_lymphoid"
}

########Color GDSC
A <- NULL; for(i in 1:length(drugs)){A <- union(A, unique(gdsc.drug.sensitivity[,i]))}; A<- A[order(A)]
color.sensitivity <- matrix(NA, nrow=length(A), ncol=1)
rownames(color.sensitivity) <- A
colnames(color.sensitivity) <- "col"
color.sensitivity[1:nrow(color.sensitivity)-1, "col"] <- colorRampPalette(c("blue" , "purple", "red"))(nrow(color.sensitivity)-1)
color.sensitivity[nrow(color.sensitivity), "col"] <- "#000000"


###### drug sensitivity

objOrderDrugs <- fnOrderDrugs(data=ccle.drug.sensitivity, filename=file.path(path.diagrams, "CCLE_DrugSensitivity.pdf"), ylab="auc recomputed", main="CCLE Drug sensitivity")
invisible(fnOrderDrugs(gdsc.drug.sensitivity, file.path(path.diagrams, "GDSC_DrugSensitivity.pdf"), ylab="auc recomputed", main="GDSC Drug sensitivity"))

##gray
#objOrderDrugs <- fnOrderDrugs(ccle.drug.sensitivity, file.path(path.diagrams, "DrugSensitivity_allgenes.pdf"))
DOI <- objOrderDrugs$order
ccle.drug.sensitivity.ordered <- objOrderDrugs$ordered

load(file.path(path.data, file.sensitivity), verbose=T)
if(length(drug.association)==2) {
  drug.association <- drug.association[[effect.size]]
  drug.association.statistics <- drug.association.statistics[[effect.size]]
}

Prototype <- drug.association[[1]][[1]]

myf <- file.path(path.diagrams, "allGenes_association_matrices.RData")
if(file.exists(myf)){
  load(myf)
}else{
  models.drugs.names <- expand.grid(drugs, Models)
  FDR_List <- matrix(NA, ncol=(length(drugs) * length(Models)), nrow=length(drug.association), dimnames=list(names(drug.association), paste(models.drugs.names[, 1], models.drugs.names[, 2], sep ="_")))
  estimate_List <- Pvalues.Numinals  <-  FDR_List
  models.plus.drugs.names <- expand.grid(drugs, c("M0", Models))
  
  statistics.matrix <- matrix(NA, ncol=(length(drugs) * (length(Models) + 1)), nrow=length(drug.association), dimnames=list(names(drug.association), paste(models.plus.drugs.names[, 1], models.plus.drugs.names[, 2], sep ="_")))
  
  best.isoforms.matrix <- matrix(NA, ncol=length(drugs) , nrow=length(drug.association))
  colnames(best.isoforms.matrix) <- drugs
  rownames(best.isoforms.matrix) <- names(drug.association)
  
  Min.Pvalues <- NULL
  for(i in drugs) {
    statistics.matrix[,paste(i, "M0",sep ="_")] <- sapply(drug.association.statistics, function(x){ifelse(!is.null(x[[i]]["median", "M0"]), x[[i]]["median", "M0"], 0)})
    for(j in 1:length(Models)) {
      FDR_List[, paste(i, Models[j], sep ="_")] <- p.adjust(sapply(drug.association, function(x){ifelse(!is.null(x[[i]]["M0",Models[j]]), x[[i]]["M0", Models[j]], 1)}) ,method=adjustment.method)
      Min.Pvalues[paste(i,Models[j],sep ="_")] <- min(sapply(drug.association, function(x){ifelse(!is.null(x[[i]]["M0", Models[j]]), x[[i]]["M0", Models[j]], 1)}))
      estimate_List[,paste(i,Models[j],sep ="_")] <- sapply(drug.association, function(x){ifelse(!is.null(x[[i]][Models[j], "M0"]), x[[i]][Models[j], "M0"], 0)})
      Pvalues.Numinals[,paste(i,Models[j],sep ="_")] <- sapply(drug.association, function(x){ifelse(!is.null(x[[i]]["M0",Models[j]]), x[[i]]["M0",Models[j]], 1)})
      statistics.matrix[,paste(i,Models[j],sep ="_")] <- sapply(drug.association.statistics, function(x){ifelse(!is.null(x[[i]]["median", Models[j]]), x[[i]]["median", Models[j]], 0)})
    }
    best.isoforms.matrix[,i] <- sapply(drug.association.best.isoforms, function(x){ifelse(!is.null(x[[i]]), x[[i]], "")}) 
  }
  
  load("../data/ensembl.map.genes.isoforms.GRCh38.87.RData")
  fnNumber_Isoforms_of_Gene <- function(Gene_Map=ensembl.map.genes.isoforms, GeneId)
  {
    return(length(which(unlist(strsplit(Gene_Map[, as.character(GeneId)], ",")) %in% colnames(ccle.isoforms.fpkm))))
  }
  
  isoforms_No_List <- matrix(NA, ncol=1, nrow=nrow(annot.ensembl.all.genes))
  colnames(isoforms_No_List) <- "isoforms.NO"
  rownames(isoforms_No_List) <- rownames(annot.ensembl.all.genes)
  for( i in 1:nrow(isoforms_No_List))
  {
    isoforms_No_List[i,] <- fnNumber_Isoforms_of_Gene(GeneId=as.character(rownames(isoforms_No_List)[i]))
  }
  isoforms_No_List <- data.frame(isoforms_No_List,stringsAsFactors=FALSE)
  #  save(isoforms_No_List, file=myf)  
  #}
  isoforms_No_List <- isoforms_No_List[GeneList, , drop=FALSE]
  save(FDR_List, estimate_List, Pvalues.Numinals, statistics.matrix, best.isoforms.matrix, Min.Pvalues, isoforms_No_List, file=myf)
}


##########Analyses

result.effect.size <- fnComputeAssociateGenes.effect.size(FDR_CutOff=fdr.cut.off, effect.size_CutOff=effect.size.cut.off)
write.csv(fnWilcox(result.effect.size, TRUE)$comparison, file=file.path(path.diagrams, "comparison_test_wilcox.csv"))
###Figure 2A and 3A
### The number of significant predictive biomarkers identified in training
### biomarkerd are plotted in seperate bars for isoforms and gene models
barplot.models(model=result.effect.size, isoforms_No="all", sign="all", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=1.2)
#Figure 3B
if(data.type == "all"){
    cindexDistributions()
}else{
    barplot.models(model=result.effect.size, isoforms_No="1.isoform", sign="all", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
    barplot.models(model=result.effect.size, isoforms_No="n.isoforms", sign="all", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
    #barplot.models(model=result.effect.size, isoforms_No="all", sign="positive", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
    #barplot.models(model=result.effect.size, isoforms_No="all", sign="negative", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
    associations <- associations.all.drugs(model.rank="M3B", annot.ensembl.all.genes=annot.ensembl.all.genes)
    save(associations, file=file.path(path.diagrams, "associations.RData"))
    ##all.biomarkers would include all the significant biomarkers
    all.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=fdr.cut.off, BioNo="all", rank.type="pvalue.adj")
    ##constraint the analyses to those predicted biomarkers with effect size greater than 0.55
    for(i in names(all.biomarkers)) {
      all.biomarkers[[i]] <- all.biomarkers[[i]][which(all.biomarkers[[i]]$cindex > effect.size.cut.off),]
    }
    save(all.biomarkers, file=file.path(path.diagrams, "all.biomarkers.original.RData"))
    
    if(!breast.specific)
    {
      source(file.path(path.code, "foo_training.R"))
      all.biomarkers <- fnidentify.tissue.specific.biomarkers(biomarkers=all.biomarkers)
    }else{
      for( i in 1:length(all.biomarkers))
      {
        if(all(is.na(all.biomarkers[[i]])))
        {
          all.biomarkers[[i]][, tissue] <- NA
          all.biomarkers[[i]][, paste0(tissue, "_boot")] <- NA
        }else{
          all.biomarkers[[i]][, tissue] <- 1
          all.biomarkers[[i]][, paste0(tissue, "_boot")] <- 1
        }
      }
    }
    all.biomarkers <- lapply(all.biomarkers, function(x){if("breast" %in% colnames(x)){x[order(x[, "breast"], na.last=T, decreasing=T),]}else{x}})
    save(all.biomarkers, file=file.path(path.data, "all.biomarkers.RData"))
    #WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=file.path(path.diagrams, "all.biomarkers.xlsx"), row.names=TRUE)
    
    
    BiomarkersNo <- matrix(NA, ncol=3, nrow=length(all.biomarkers))
    rownames(BiomarkersNo) <- names(all.biomarkers)
    colnames(BiomarkersNo) <- c("gene.specific", "isoform.specific", "common")
    for( i in names(all.biomarkers))
    {
      top.significant.table <- table(all.biomarkers[[i]][, "specificity"])
      
      BiomarkersNo[i, "isoform.specific"] <- top.significant.table["isoform.specific"]
      BiomarkersNo[i, "gene.specific"] <- top.significant.table["gene.specific"]
      BiomarkersNo[i, "common"] <- top.significant.table["common"]/2
    }
    write.csv(BiomarkersNo, file=file.path(path.diagrams, "Biomarkers.csv"))
    
    
    TOP <- matrix(NA, ncol=9, nrow=length(all.biomarkers))
    rownames(TOP) <- names(all.biomarkers)
    colnames(TOP) <- c("symbol", "isoforms.no", "biomarker.id", "type", "specificity", "estimate", effect.size, "fdr", "delta.rank")
    for( i in names(all.biomarkers)) {
      xx <- which.max(all.biomarkers[[i]][,effect.size])
      #xx <- 1
      if(length(xx) > 0) {
        TOP[i, "symbol"] <- as.character(all.biomarkers[[i]][xx, "symbol"])
        TOP[i, "isoforms.no"] <- all.biomarkers[[i]][xx, "isoforms.no"]
        TOP[i, "biomarker.id"] <- as.character(all.biomarkers[[i]][xx, "biomarker.id"])
        TOP[i, "type"] <- all.biomarkers[[i]][xx, "type"]
        TOP[i, "specificity"] <- all.biomarkers[[i]][xx, "specificity"]
        TOP[i, "estimate"] <- all.biomarkers[[i]][xx, "estimate"]
        TOP[i, effect.size] <- all.biomarkers[[i]][xx, effect.size]
        TOP[i, "fdr"] <- all.biomarkers[[i]][xx, adjustment.method]
        TOP[i, "delta.rank"] <- all.biomarkers[[i]][xx, "delta.rank"]
      }
    }
    write.csv(TOP, file=file.path(path.diagrams, "TopOne.csv"))
    
    Check.KnownAssociations(associations)
    
    percentage.biomarkers.type <- fnPercentageBiomarkersType(all.biomarkers)
    ###Figure 2B
    ### The ratio of biomarkers according to their specificity in training 
    mystacked.barplot.simple(Filename=file.path(path.diagrams, "ProportionBiomarkersType.pdf"), data=percentage.biomarkers.type, main.label="Proportion of specificity of biomarkers", cex=1.2)
    
    
    ###
    ## annotation of differnet set of biomarkers based on their biotaypes
    ##update all.biomarkers with gene/transcript biotypes
    for(drug in names(all.biomarkers)) {
      ii <- which(all.biomarkers[[drug]][,"type"] == "isoform")
      gg <- which(all.biomarkers[[drug]][,"type"] == "gene")
      if(length(ii) > 0){
        all.biomarkers[[drug]][ii, "biotype"] <- annot.ensembl.all.isoforms[all.biomarkers[[drug]][ii, "transcript.id"], "TranscriptBioType"]
      }
      if(length(gg) > 0){
        all.biomarkers[[drug]][gg, "biotype"] <- annot.ensembl.all.genes[all.biomarkers[[drug]][gg, "gene.id"], "GeneBioType"]
      }
    }
    ### Supplementary Figure 10
    ## annotation of differnet set of biomarkers based on their biotaypes
    mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")[c(1, 4, 2)]
    pdf(file.path(path.result, "biomarkers_specificity_biotypes2.pdf"), width=21, height=10)
    par(mfrow=c(2,3))
    par(mar=c(7.1, 4.1, 4.1, 10.1), xpd=TRUE)
    biotypes <- c("protein_coding", "antisense", "processed_transcript", "lincRNA", "pseudogene", "other")
    rr <- list()
    for(i in 1:length(biotypes)) {
      rr[[i]] <- matrix(0, nrow=length(all.biomarkers), ncol=3, dimnames=list(names(all.biomarkers), c("isoform.specific", "common", "gene.specific")))
    }
    for(drug in names(all.biomarkers)){
      if(all(!is.na(all.biomarkers[[drug]]$specificity))){
        xx <- table(all.biomarkers[[drug]]$specificity, all.biomarkers[[drug]]$biotype)
        for(i in 1:length(biotypes)) {
          yy <- grep(biotypes[i], colnames(xx))
          if(length(yy) > 0) {
            mm <- apply(xx[ , yy, drop=FALSE], 1, sum)
            rr[[i]][drug, intersect(colnames(rr[[i]]), names(mm))] <- mm[intersect(colnames(rr[[i]]), names(mm))]/sum(mm)
            xx <- xx[ ,-yy, drop=FALSE]
          }else if (biotypes[i] == "other"){
            mm <- apply(xx, 1, sum)
            rr[[i]][drug, intersect(colnames(rr[[i]]), names(mm))] <- mm[intersect(colnames(rr[[i]]), names(mm))]/sum(mm)
            xx <- xx[ , -yy, drop=FALSE]
          }
        }
      }
    }
    
    for(i in 1:length(biotypes)) {
      barplot(t(rr[[i]]), las=2,col=mycol, main=biotypes[i], ylab="percentage")
    }
    legend("topright", inset=c(-0.2,0), legend=colnames(rr[[1]]), fill=mycol, bty="n")
    dev.off()
    ###
    
    
    ### Supplementary Figure 11
    ## annotation of differnet set of biomarkers based on the number of alternative spliced products in their corresponding gene
    mycol <- RColorBrewer::brewer.pal(n=3, name="Set1")[c(2,1)]
    isoform.specific<- gene.specific<- common <- matrix(NA, nrow=length(all.biomarkers), ncol=2, dimnames=list(names(all.biomarkers), c("1 isoform", "n isoforms")))
    pdf(file.path(path.result, "biomarkers_specificity_isoform_no.pdf"), width=21, height=5)
    par(mfrow=c(1,3))
    par(mar=c(7.1, 4.1, 4.1, 10.1), xpd=TRUE)
    rr <- list()
    for(drug in names(all.biomarkers)){
      xx <- table(all.biomarkers[[drug]]$specificity, all.biomarkers[[drug]]$isoforms.no)
      mm <- apply(xx, MARGIN=1, function(x){s <- x[which(x!=0)];c(s[1], (sum(s)-s[1]))})
      common[drug, ] <- if("common" %in% colnames(mm)){mm[,"common"]}else{c(0,0)}
      gene.specific[drug, ] <- if("gene.specific" %in% colnames(mm)){mm[,"gene.specific"]}else{c(0,0)}
      isoform.specific[drug, ] <- if("isoform.specific" %in% colnames(mm)){mm[,"isoform.specific"]}else{c(0,0)}
      rr[[drug]] <- mm
    }
    common <- common/apply(common,MARGIN = 1, function(x){ss<- sum(x);ifelse(ss!=0, ss, 1)})
    barplot(t(common), las=2,col=mycol[1:2], main="common", ylab="percentage")
    
    
    gene.specific <- gene.specific/apply(gene.specific,MARGIN = 1, function(x){ss<- sum(x);ifelse(ss!=0, ss, 1)})
    barplot(t(gene.specific), las=2,col=mycol[1:2], main="gene specific", ylab="percentage")
    
    
    isoform.specific <- isoform.specific/apply(isoform.specific,MARGIN = 1, function(x){ss<- sum(x);ifelse(ss!=0, ss, 1)})
    barplot(t(isoform.specific), las=2,col=mycol[1:2], main="isoform specific", ylab="percentage")
    legend("topright", inset=c(-0.2,0), legend = colnames(isoform.specific), fill = mycol[1:2], bty="n")
    dev.off()
    ###
    
    
    mycol3 <- RColorBrewer::brewer.pal(n=5, name="Set3")
    pdf(file.path(path.diagrams, "biomarkers_specificity_biotypes.pdf"), width=21, height=5)
    par(mfrow=c(1,3))
    par(mar=c(7.1, 4.1, 4.1, 10.1), xpd=TRUE)
    biotypes <- c("protein_coding", "antisense", "processed_transcript", "lincRNA")
    isoform.specific<- gene.specific<- common <- matrix(0, nrow=length(all.biomarkers), ncol=length(biotypes)+1, dimnames=list(names(all.biomarkers), c(biotypes, "other")))
    for(drug in names(all.biomarkers)){
      if(all(!is.na(all.biomarkers[[drug]]$specificity))){
        xx <- table(all.biomarkers[[drug]]$specificity, all.biomarkers[[drug]]$biotype)
        mm <- cbind(xx[, intersect(biotypes, colnames(xx))], "other"=apply(xx[,which(!colnames(xx) %in% biotypes), drop=FALSE], 1, sum))
        if("common" %in% rownames(mm)){
          common[drug, colnames(mm)] <- mm["common",]/apply(mm , 1, sum)["common"]
        }
        if("gene.specific" %in% rownames(mm)){
          gene.specific[drug, colnames(mm)] <- mm["gene.specific",]/apply(mm , 1, sum)["gene.specific"]
        }
        if("isoform.specific" %in% rownames(mm)){
          isoform.specific[drug, colnames(mm)] <- mm["isoform.specific",]/apply(mm , 1, sum)["isoform.specific"]
        }
      }
    }
    barplot(t(common), las=2,col=mycol3, main="common", ylab="percentage")
    barplot(t(gene.specific), las=2,col=mycol3, main="gene specific", ylab="percentage")
    barplot(t(isoform.specific), las=2,col=mycol3, main="isoform specific", ylab="percentage")
    legend("topright", inset=c(-0.2,0), legend=c(biotypes, "other"), fill = mycol3, bty="n")
    dev.off()

}
