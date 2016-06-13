file.sensitivity <- "auc_recomputed_drug_association.RData"

require(stringr) || stop("Library stringr is not available!")
require(corrplot) || stop("Library corrplot is not available!")
require(VennDiagram) || stop("Library VennDiagram is not available!")
require(RColorBrewer) || stop("Library RColorBrewer is not available!")
require(SDMTools) || stop("Library SDMTools is not available!")
require(calibrate) || stop("Library calibrate is not available!")
require(Hmisc) || stop("Library Hmisc is not available!")
require(grid) || stop("Library grid is not available!")
require(gridBase) || stop("Library gridBase is not available!")
require(lattice) || stop("Library lattice is not available!")
require(WriteXLS) || stop("Library WriteXLS is not available!")
require(SDMTools) || stop("Library SDMTools is not available!")
require(PharmacoGx) || stop("Library SDMTools is not available!")
require(Biobase) || stop("Library SDMTools is not available!")


options(stringsAsFactors=FALSE)
path.data <- "data"
path.code <- file.path("code")
path.result <- file.path("result")

tissue <- "breast"
model.method <- "glm"
glm.family <- "gaussian" 
stat <- "r.squared"


source(file.path(path.code, "foo.R"))

load(file.path(path.data, "PSets/CCLE_isoforms.RData"))
load(file.path(path.data, "PSets/GDSC.RData"))

load(file.path(path.data, file.sensitivity), verbose=T)

training.type <- ifelse(regexpr("ccle", file.sensitivity) >= 1, "CCLE", "CCLE_GDSC")
sensitivity.type <- paste0(unlist(strsplit(file.sensitivity, "_"))[1:2], collapse="_") 
breast.specific <- ifelse(regexpr("breast", file.sensitivity) < 1, FALSE, TRUE)

## It assumes the sensitive and resistent weights are 3rd and 4th strings seperated by underline in file.sensitivity name
if(sensitivity.type == "slope_recomputed")
{
  if(training.type == "CCLE_GDSC")
  {
    ss<- unlist(strsplit(file.sensitivity, "_"))
    res.weight <- as.numeric(ss[3]); sens.weight <- as.numeric(ss[4]);
    path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))
  }else
  {
    ss <- unlist(strsplit(file.sensitivity, "_"))
    res.weight <- as.numeric(ss[4]); sens.weight <- as.numeric(ss[5]);
    path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))    
  }
}else{
  ss <- unlist(strsplit(file.sensitivity, "_"))
  path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))
}
if (!file.exists(path.diagrams)){dir.create(file.path(path.diagrams))}

Models <- c("M1", "M2", "M3B")
Models.names <- c("Microarray", "Genes", "Isoforms")
names(Models.names) <- Models
adjustment.method <- "fdr"


if(training.type == "CCLE_GDSC") {
  #common <- PharmacoGx::intersectPSet(pSets=list("CCLE"=CCLE, "GDSC"=GDSC), strictIntersect=F)
  #cells <- intersect(common$CCLE@cell$cellid, pData(common$CCLE@molecularProfiles$isoforms)[, "cellid"])
  #common$CCLE@molecularProfiles$rnaseq <- CCLE@molecularProfiles$rnaseq[, rownames(pData(CCLE@molecularProfiles$isoforms))]
  #CCLE <- common$CCLE
  #GDSC <- common$GDSC
  drugs <- intersect(PharmacoGx::drugNames(CCLE), PharmacoGx::drugNames(GDSC))
}
ccle.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure=sensitivity.type, drugs=drugs))
gdsc.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC, sensitivity.measure=sensitivity.type, drugs=drugs))

if(as.character(sensitivity.type) == "slope_recomputed") {
  ##manual cutoff for sensitivity calls based on slope
  cutoff <- 0.27
  ccle.drug.sensitivity.wighted <- ccle.drug.sensitivity
  ccle.drug.sensitivity.wighted[ccle.drug.sensitivity.wighted < cutoff] <- res.weight
  ccle.drug.sensitivity.wighted[ccle.drug.sensitivity.wighted >= cutoff] <- sens.weight
  
  gdsc.drug.sensitivity.wighted <- gdsc.drug.sensitivity
  gdsc.drug.sensitivity.wighted[gdsc.drug.sensitivity.wighted < cutoff] <- res.weight
  gdsc.drug.sensitivity.wighted[gdsc.drug.sensitivity.wighted >= cutoff] <- sens.weight
  
  ccle.drug.sensitivity <- ccle.drug.sensitivity.wighted * ccle.drug.sensitivity
  gdsc.drug.sensitivity <- gdsc.drug.sensitivity.wighted * gdsc.drug.sensitivity
}

drugs <- colnames(ccle.drug.sensitivity)

ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", fill.missing=FALSE)))
ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", fill.missing=FALSE)))

ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))

ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=CCLE@cell[, "cellid"])
colnames(ccle.tissuetype) <- "tissue.type"

ccle.primary.tissuetype <- ccle.tissuetype
ccle.primary.drug.sensitivity <- ccle.drug.sensitivity

mean.ccle.isoforms.fpkm <- colMeans(ccle.isoforms.fpkm)
if("haematopoietic_and_lymphoid_tissue" %in% ccle.tissuetype[,1])
{
  ccle.tissuetype[which(ccle.tissuetype[,1] == "haematopoietic_and_lymphoid_tissue"), 1] <- "haematopoietic_and_lymphoid"
}

########Color GDSC
A <- NULL; for(i in 1:15){A <- union(A, unique(gdsc.drug.sensitivity[,i]))}; A<- A[order(A)]
color.sensitivity <- matrix(NA, nrow=length(A), ncol=1)
rownames(color.sensitivity) <- A
colnames(color.sensitivity) <- "col"
color.sensitivity[1:nrow(color.sensitivity)-1, "col"] <- colorRampPalette(c("blue" , "purple", "red"))(nrow(color.sensitivity)-1)
color.sensitivity[nrow(color.sensitivity), "col"] <- "#000000"


###### drug sensitivity

drugs_No <- ncol(ccle.drug.sensitivity)
objOrderDrugs <- fnOrderDrugs(ccle.drug.sensitivity, file.path(path.diagrams, "CCLE_DrugSensitivity.pdf"), ylab="auc recomputed", main="CCLE Drug sensitivity")
invisible(fnOrderDrugs(gdsc.drug.sensitivity, file.path(path.diagrams, "GDSC_DrugSensitivity.pdf"), ylab="auc recomputed", main="GDSC Drug sensitivity"))

##gray
#objOrderDrugs <- fnOrderDrugs(ccle.drug.sensitivity, file.path(path.diagrams, "DrugSensitivity_allgenes.pdf"))
DOI <- objOrderDrugs$order
ccle.drug.sensitivity.ordered <- objOrderDrugs$ordered

####Microarray vs RNA_sq
#fnCorrelationMicroArray_RNASeq(MicroarrayExp=ccle.drug.microarray.exp, RNA_SeqEXP=ccle.gens.fpkm, path)

Prototype <- drug.association[[1]][[1]]


FDR_List <- matrix(NA, ncol=(drugs_No * length(Models)), nrow=length(drug.association))
models.drugs.names <- expand.grid(gsub("drugid_", "", colnames(ccle.drug.sensitivity)),Models)
colnames(FDR_List) <- paste(models.drugs.names[,1],models.drugs.names[,2],sep ="_")
rownames(FDR_List) <- names(drug.association)

estimate_List <- Pvalues.Numinals <- statistics.matrix <-  FDR_List

best.isoforms.matrix <- matrix(NA, ncol=drugs_No , nrow=length(drug.association))
colnames(best.isoforms.matrix) <- gsub("drugid_", "", colnames(ccle.drug.sensitivity))
rownames(best.isoforms.matrix) <- names(drug.association)

Min.Pvalues <- NULL
Min.Pvalues.gdsc <- NULL
for(i in drugs)
{
  for(j in 1:length(Models))
  {
    FDR_List[,paste(i,Models[j],sep ="_")] <- p.adjust(sapply(drug.association, function(x){ifelse(!is.null(x[[i]]["M0",Models[j]]), x[[i]]["M0",Models[j]], 1)}) ,method= adjustment.method)
    Min.Pvalues[paste(i,Models[j],sep ="_")] <- min(sapply(drug.association, function(x){ifelse(!is.null(x[[i]]["M0",Models[j]]), x[[i]]["M0",Models[j]], 1)}))
    estimate_List[,paste(i,Models[j],sep ="_")] <- sapply(drug.association, function(x){ifelse(!is.null(x[[i]][Models[j], "M0"]), x[[i]][Models[j], "M0"], 0)})
    Pvalues.Numinals[,paste(i,Models[j],sep ="_")] <- sapply(drug.association, function(x){ifelse(!is.null(x[[i]]["M0",Models[j]]), x[[i]]["M0",Models[j]], 1)})
    statistics.matrix[,paste(i,Models[j],sep ="_")] <- sapply(drug.association.statistics, function(x){ifelse(!is.null(x[[i]]["median",Models[j]]), x[[i]]["median",Models[j]], 0)})
  }
  best.isoforms.matrix[,i] <- sapply(drug.association.best.isoforms, function(x){ifelse(!is.null(x[[i]]), x[[i]], "")}) 
}
##multiplying adjusted pvalues by 2 !!!!?
# for ( i in 1:nrow(FDR_List))
# {
#   for (j in (drugs_No + 1):ncol(FDR_List))
#   {
#     temp <- FDR_List[i,j] * 2
#     FDR_List[i,j] <- ifelse(temp > 1, 1, temp)
#   }
# }

if(file.exists("data/isoforms.no.per.gene.RData")){
  load("data/isoforms.no.per.gene.RData", verbose=TRUE)
}else{
  load("data/ensembl.map.genes.isoforms.RData")
  fnNumber_Isoforms_of_Gene <- function(Gene_Map=ensembl.map.genes.isoforms, GeneId)
  {
    return(nrow(data.frame(strsplit(Gene_Map[, as.character(GeneId)], ", "))))
  }
  
  isoforms_No_List <- matrix(NA, ncol=1, nrow=length(drug.association))
  colnames(isoforms_No_List) <- "isoforms.NO"
  rownames(isoforms_No_List) <- names(drug.association)
  for( i in 1:length(drug.association))
  {
    isoforms_No_List[i,] <- fnNumber_Isoforms_of_Gene(GeneId=as.character(names(drug.association)[i]))
  }
  isoforms_No_List <- data.frame(isoforms_No_List,stringsAsFactors=FALSE)
  save(isoforms_No_List, file="data/isoforms.no.per.gene.RData")  
}

### subsetting genes to remove pseudogenes
ff <- fData(CCLE@molecularProfiles$rnaseq)
genes <- as.character(ff[grep("pseudogene", ff$GeneBioType, invert=TRUE), "EnsemblGeneId"])
FDR_List <- FDR_List[genes, , drop=FALSE]
estimate_List <- estimate_List[genes, , drop=FALSE]
Pvalues.Numinals <- Pvalues.Numinals[genes, , drop=FALSE]
statistics.matrix <- statistics.matrix[genes, , drop=FALSE]
best.isoforms.matrix <- best.isoforms.matrix[genes, , drop=FALSE]

Min.Pvalues <- Min.Pvalues[genes]
isoforms_No_List <- isoforms_No_List[genes, , drop=FALSE]
###########################################
  
  
##########Analyses


source("code/foo_Analysis.R")

Models2 <- c("M2", "M3B")
result.unsigned <- fnComputeAssociateGenes(FDR_CutOff=0.01,signed=FALSE)
#barplot.models(model=result.unsigned, isoforms_No="all", prototype=Models2, main.title=sprintf("%s  <  1%% \n Altogether", adjustment.method), breakpoint="Regular")

#result.stat <- fnComputeAssociateGenes.stat(FDR_CutOff=0.01, stat_CutOff=0,signed=TRUE)
#barplot.models(model=result.stat, isoforms_No="all", prototype=Models2, main.title=sprintf("%s  <  1%% \n R2 > 0.60 \n Altogether", adjustment.method), breakpoint="Regular")
ranges <- 1
c <- max(statistics.matrix, na.rm=T)# c <- 0.70
for (i in 1:70)
{
  ranges <- c(ranges, c)
  c <- c - .0025
}
ranges.names <- as.character(ranges)
ranges <- c(ranges,0)
#ranges <- c(1,.68,.65,.60,.58,0)
#ranges.names <- c(">.68", "<.68", ">.65", ">.60", ">.58")
result.ranged <- fnComputeAssociateGenes.stat.range(FDR_CutOff=0.01, biomarkers.sign="all")
#stacked.barplot.models.ranges(model=result.ranged, isoforms_No="all", prototype=Models2, main.title=sprintf("%s  <  1%% ", adjustment.method), breakpoint="Regular", yaxis="Regular")


result.ranged.log <- list(length(result.ranged))
for(i in 1:drugs_No)
{
  result.ranged.log[[i]] <- list()
  names(result.ranged.log)[i] <- names(result.ranged)[i]
  
  for(j in 1:length(result.ranged[[i]]))
  {
    result.ranged.log[[i]][[j]] <- matrix(0, nrow=2, ncol=ncol(result.ranged[[i]][[j]]))
    rownames(result.ranged.log[[i]][[j]]) <- rownames(result.ranged[[i]][[j]])
    colnames(result.ranged.log[[i]][[j]]) <- colnames(result.ranged[[i]][[j]])
    if(sum(result.ranged[[i]][[j]]) == 0)
    {
      result.ranged.log[[i]][[j]] <- result.ranged[[i]][[j]]
    }else{
      for(k in 1:ncol(result.ranged[[i]][[j]]))
      {   
        result.ranged.log[[i]][[j]][1,k] <- (log10(sum(result.ranged[[i]][[j]])+1)/sum(result.ranged[[i]][[j]])) * result.ranged[[i]][[j]][1,k]
        result.ranged.log[[i]][[j]][2,k] <- (log10(sum(result.ranged[[i]][[j]])+1)/sum(result.ranged[[i]][[j]])) * result.ranged[[i]][[j]][2,k]
      }
    }
    names(result.ranged.log[[i]])[j] <- names(result.ranged[[i]])[j]
  }
}

stacked.barplot.models.ranges(model=result.ranged.log, isoforms_No="all", prototype=Models2, main.title=sprintf("%s  <  1%%", adjustment.method), breakpoint="Regular", yaxis="Log")
barplot.models(model=result.unsigned, isoforms_No="all",prototype=Models2, main.title=sprintf("%s  <  1%%", adjustment.method), breakpoint="Regular", yaxis="Log")
########In case of logarithm I should change y axis to show the real number of biomarkers instead of the log values
########I should replace the intensity marker

#fnVennDiagramTriple(model=result.unsigned, m3="M3B")
#fnVennDiagramPair(model=result.unsigned, m3="M3B")
#fnMarkersPercent(model=result.unsigned, m1="M2", m2="M3B")
#fnIsoformVSGeneR2(FDR_CutOff=0.01, m1="M2", m2="M3B")
#fnIsoformVSGeneExp(FDR_CutOff=0.01, m1="M2", m2="M3B")

result.signed <- fnComputeAssociateGenes(FDR_CutOff=0.01, signed=TRUE)
stacked.barplot.models(model=result.signed, isoforms_No="all", prototype=Models2, main.title=sprintf("%s  <  1%%", adjustment.method), breakpoint="Regular")
write.csv(fnWilcox(result.signed, TRUE)$comparison, file=file.path(path.diagrams, "comparison_test_wilcox.csv"))
write.csv(fnWilcox(result.unsigned, FALSE)$comparison, file=file.path(path.diagrams, "comparison_test_wilcox_u.csv"))


#FDR <- c(.9,.5,.2,.1,.05,.02,.01,.001)
#fancystat <- fnFancystat(FDR, combinations=c("M2_M3B"),signed=FALSE)
#fnCuttOffs(cutoff_statistics=fancystat, combination="M2_M3B")

## result csv for each drug shows the number of genes when isoform model is significant but gene models (Microarray & RNA_seq) are not (based on cutoff)
# associations.significant.isoforms <- associations.all.drugs.significant(rank="Isoforms.pvalue", method="isoform", cutoff=1, R2_cut=0, exp_cut=0, annot.ensembl.all.genes=fData(CCLE@molecularProfiles$rnaseq))
# A <- sapply(associations.significant.isoforms, function(x){nrow(x)})
# A <- cbind(A,sapply(associations.significant.isoforms, function(x){median(x[, "Isoforms.fdr"])}))
# A <- cbind(A,sapply(associations.significant.isoforms, function(x){median(x[, "isoform.R2"])}))
# A <- cbind(A,sapply(associations.significant.isoforms, function(x){median(x[, "isoform.exp"])}))
# colnames(A) <- c("bimarkers No", "median fdr", "median R2", "median isoform expression")
# write.csv(A, file=file.path(path.diagrams,sprintf("IsoformSpecificBiomarkersStat_bonf_%s_R2_%s_exp_%s.csv",1,0,0)))


# associations.significant.isoforms <- associations.all.drugs.significant("Isoforms.pvalue", method="isoform", cutoff=.01, R2_cut =.6, exp_cut=.5, annot.ensembl.all.genes=fData(CCLE@molecularProfiles$rnaseq))
# A <- sapply(associations.significant.isoforms, function(x){nrow(x)})
# A <- cbind(A,sapply(associations.significant.isoforms, function(x){median(x[, "Isoforms.fdr"])}))
# A <- cbind(A,sapply(associations.significant.isoforms, function(x){median(x[, "isoform.R2"])}))
# A <- cbind(A,sapply(associations.significant.isoforms, function(x){median(x[, "isoform.exp"])}))
# colnames(A) <- c("bimarkers No", "median fdr", "median R2", "median isoform expression")
# write.csv(A, file=file.path(path.diagrams,sprintf("IsoformSpecificBiomarkersStat_bonf_%s_R2_%s_exp_%s.csv",.01,.6,.5)))
# 
# WriteXLS::WriteXLS("associations.significant.isoforms", ExcelFileName=file.path(path.diagrams, "drug_association_significant_isoforms.xlsx"), row.names=TRUE)

#associations.significant.isoforms.tissues <- associations.all.drugs.tissues.significant(associations.significant.isoforms, cut_off=1)


associations <- associations.all.drugs(model.rank="M3B", annot.ensembl.all.genes=fData(CCLE@molecularProfiles$rnaseq))
WriteXLS::WriteXLS("associations", ExcelFileName=file.path(path.diagrams, "drug_association.xlsx"), row.names=TRUE)

#top.significant.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=0.01, BioNo=50, rank.type="pvalue.adj")
#WriteXLS::WriteXLS("top.significant.biomarkers", ExcelFileName=file.path(path.diagrams, "top.biomarkers.xlsx"), row.names=TRUE)


##all.biomarkers would include all the significant biomarkers
all.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=0.01, BioNo="all", rank.type="pvalue.adj")
if(!breast.specific)
{
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
#all.biomarkers <- sapply(all.biomarkers, function(x){if("breast" %in% colnames(x)){x[order(x[, "breast"], na.last=T, decreasing=T),]}else{x}})
WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=file.path(path.diagrams, "all.biomarkers.xlsx"), row.names=TRUE)


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
colnames(TOP) <- c("symbol", "isoforms.no", "biomarker.id", "type", "specificity", "estimate", "R2", "fdr", "delta.rank")
for( i in names(all.biomarkers))
{
  TOP[i, "symbol"] <- as.character(all.biomarkers[[i]][1, "symbol"])
  TOP[i, "isoforms.no"] <- all.biomarkers[[i]][1, "isoforms.no"]
  TOP[i, "biomarker.id"] <- as.character(all.biomarkers[[i]][1, "biomarker.id"])
  TOP[i, "type"] <- all.biomarkers[[i]][1, "type"]
  TOP[i, "specificity"] <- all.biomarkers[[i]][1, "specificity"]
  TOP[i, "estimate"] <- all.biomarkers[[i]][1, "estimate"]
  TOP[i, "R2"] <- all.biomarkers[[i]][1, "R2"]
  TOP[i, "fdr"] <- all.biomarkers[[i]][1, "fdr"]
  TOP[i, "delta.rank"] <- all.biomarkers[[i]][1, "delta.rank"]
}
write.csv(TOP, file=file.path(path.diagrams, "TopOne.csv"))

Check.KnownAssociations(associations)
Check.KnownAssociations.PLOS1(associations)

percentage.biomarkers.type <- fnPercentageBiomarkersType(all.biomarkers)
mystacked.barplot.simple(Filename=file.path(path.diagrams, "ProportionBiomarkersType.pdf"), data=percentage.biomarkers.type, main.label="Proportion of specificity of biomarkers")

#fnGeneIsoformCorrelation(all.biomarkers)
  

# sensitivity <- ccle.drug.sensitivity[complete.cases(ccle.drug.sensitivity),]
# colnames(sensitivity) <- gsub("drugid_", "", colnames(sensitivity))
# sensitivity <- sapply(sensitivity, function(x){as.numeric(x)})
# pdf(file.path(path.diagrams,sprintf("CCLE_Sensitivity_Heatmap.pdf")), height=5,width=5)
# #heatmap(sensitivity, Rowv=NA, Colv=NA, col=colorRampPalette(c("red", "blue"))(100), scale="column", margins=c(5,10), cexRow =.2, cexCol=.7)
# heatmap(sensitivity, col=colorRampPalette(c("blue", "red"))(20), scale="none", margins=c(5,10), cexRow =.1, cexCol=.7)
# 
# dev.off()




# sensitivity <- gdsc.drug.sensitivity[complete.cases(gdsc.drug.sensitivity),]
# colnames(sensitivity) <- gsub("drugid_", "", colnames(sensitivity))
# sensitivity <- sapply(sensitivity, function(x){as.numeric(x)})
# 
# pdf(file.path(path.diagrams,sprintf("GDSC_Sensitivity_Heatmap.pdf")), height=5,width=5)
# #heatmap(sensitivity, Rowv=NA, Colv=NA, col=colorRampPalette(c("red", "blue"))(100), scale="column", margins=c(5,10), cexRow =.2, cexCol=.7)
# heatmap(sensitivity, col=colorRampPalette(c("blue", "red"))(20), scale="none", margins=c(5,10), cexRow =.2, cexCol=.7)
# 
# dev.off()

## tissue
##fnTop.significant.biomarkers.heatmap(all.biomarkers, "ERLOTINIB")
save(all.biomarkers, file=file.path(path.diagrams, "all.biomarkers.RData"))
load(file.path(path.diagrams, "all.biomarkers.RData"), verbose=T)
