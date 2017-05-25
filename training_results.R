file.sensitivity <- "auc_recomputed_drug_association_ccle_gdsc.RData"

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
require(PharmacoGx) || stop("Library PharmacoGx is not available!")
require(Biobase) || stop("Library Biobase is not available!")
require(magicaxis) || stop("Library magicaxis is not available!") #add minor tick marks to the plot
require("np") || stop("Library np is not available!")


options(stringsAsFactors=FALSE)
path.data <- "data"
path.code <- file.path("code")
path.result <- file.path("result")

tissue <- "breast"
model.method <- "glm"
glm.family <- "gaussian" 
effect.size <- "cindex" #c("r.squared", "cindex")
adjustment.method <- "fdr"#c("bonferroni", "fdr")
subset.genes <- TRUE
subset.genes.method <- "wide.expression" #c("wide.expression", "biotype")
RNA_seq.normalize <- TRUE

training.type <- "CCLE_GDSC"#ifelse(regexpr("ccle", file.sensitivity) >= 1, "CCLE", "CCLE_GDSC")
sensitivity.type <- paste0(unlist(strsplit(file.sensitivity, "_"))[1:2], collapse="_") 
breast.specific <- ifelse(regexpr("breast", file.sensitivity) < 1, FALSE, TRUE)
data.type <- ifelse(regexpr("mutation", file.sensitivity) >= 1, "all", "expression")

source(file.path(path.code, "foo.R"))

load(file.path(path.data, "PSets/CCLE_hs.RData"), verbose=TRUE)
#load(file.path(path.data, "PSets/CCLE.CTRPv2.RData"), verbose=TRUE)
#CCLE <- CCLE.CTRPv2

load(file.path(path.data, "PSets/GDSC.RData"), verbose=TRUE)
#load(file.path(path.data, "PSets/GDSC1000.RData"), verbose=TRUE)
#GDSC <- GDSC1000


## It assumes the sensitive and resistent weights are 3rd and 4th strings seperated by underline in file.sensitivity name
if(sensitivity.type == "slope_recomputed") {
  if(training.type == "CCLE_GDSC") {
    ss<- unlist(strsplit(file.sensitivity, "_"))
    res.weight <- as.numeric(ss[3]); sens.weight <- as.numeric(ss[4]);
    path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))
  }else {
    ss <- unlist(strsplit(file.sensitivity, "_"))
    res.weight <- as.numeric(ss[4]); sens.weight <- as.numeric(ss[5]);
    path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))    
  }
}else{
  ss <- unlist(strsplit(file.sensitivity, "_"))
  #path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))
  path.diagrams<- "result/auc_recomputed_ccle_gdsc"
}
if (!file.exists(path.diagrams)){dir.create(file.path(path.diagrams))}

if(data.type == "all") {
  Models <- c("M2", "M3B", "MM", "MC")
  Models.names <- c("Genes", "Isoforms", "Mutations", "CNV")
}else {
  Models <- c("M2", "M3B")
  Models.names <- c("Genes", "Isoforms")
}
names(Models.names) <- Models

if(training.type == "CCLE_GDSC") {
  #common <- PharmacoGx::intersectPSet(pSets=list("CCLE"=CCLE, "GDSC"=GDSC), strictIntersect=F)
  #cells <- intersect(common$CCLE@cell$cellid, pData(common$CCLE@molecularProfiles$isoforms)[, "cellid"])
  #common$CCLE@molecularProfiles$rnaseq <- CCLE@molecularProfiles$rnaseq[, rownames(pData(CCLE@molecularProfiles$isoforms))]
  #CCLE <- common$CCLE
  #GDSC <- common$GDSC
  drugs <- intersect(PharmacoGx::drugNames(CCLE), PharmacoGx::drugNames(GDSC))
  gdsc.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC, sensitivity.measure=sensitivity.type, drugs=drugs))
}else{
  drugs <- drugNames(CCLE)
}

ccle.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure=sensitivity.type, drugs=drugs))

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

features <- fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"]
features <- intersect(fData(CCLE@molecularProfiles$isoforms)[,"EnsemblGeneId"], features)
if(data.type == "all"){
  features <- intersect(fData(CCLE@molecularProfiles$mutation)[,"EnsemblGeneId"], features)
  features <- intersect(fData(CCLE@molecularProfiles$cnv)[,"EnsemblGeneId"], features)
  mutation.features <- rownames(fData(CCLE@molecularProfiles$mutation))[match(features, fData(CCLE@molecularProfiles$mutation)[,"EnsemblGeneId"])]
  cnv.features <- rownames(fData(CCLE@molecularProfiles$cnv))[match(features, fData(CCLE@molecularProfiles$cnv)[,"EnsemblGeneId"])]
}
rnaseq.features <- rownames(fData(CCLE@molecularProfiles$rnaseq))[which(fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"] %in% features)]
isoforms.features <- rownames(fData(CCLE@molecularProfiles$isoforms))[which(fData(CCLE@molecularProfiles$isoforms)[,"EnsemblGeneId"] %in% features)]

if(data.type == "all"){
  ccle.mutation <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="mutation", features=mutation.features, fill.missing=FALSE, summary.stat="or")))
  colnames(ccle.mutation) <- fData(CCLE@molecularProfiles$mutation)[,"EnsemblGeneId"][match(colnames(ccle.mutation), rownames(fData(CCLE@molecularProfiles$mutation)))]
  ccle.cnv <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="cnv", features=cnv.features, fill.missing=FALSE)))
  colnames(ccle.cnv) <- fData(CCLE@molecularProfiles$cnv)[,"EnsemblGeneId"][match(colnames(ccle.cnv), rownames(fData(CCLE@molecularProfiles$cnv)))]
}
ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", features=rnaseq.features, fill.missing=FALSE)))
ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", features=isoforms.features, fill.missing=FALSE)))
ccle.isoforms.fpkm[which(is.na(ccle.isoforms.fpkm))] <- 0
if(RNA_seq.normalize == TRUE)
{
  ccle.genes.fpkm <- log2(ccle.genes.fpkm + 1)
  ccle.isoforms.fpkm <- log2(ccle.isoforms.fpkm + 1)
}

if(data.type == "all"){
  ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm), rownames(ccle.mutation), rownames(ccle.cnv))
  gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm), rownames(ccle.mutation), rownames(ccle.cnv))
}else{
  ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
}

ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=CCLE@cell[, "cellid"])
colnames(ccle.tissuetype) <- "tissue.type"

tissueTypes <- ccle.tissuetype
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

models.drugs.names <- expand.grid(drugs, Models)
FDR_List <- matrix(NA, ncol=(length(drugs) * length(Models)), nrow=length(drug.association), dimnames=list(names(drug.association), paste(models.drugs.names[, 1], models.drugs.names[, 2], sep ="_")))
estimate_List <- Pvalues.Numinals  <-  FDR_List
models.plus.drugs.names <- expand.grid(drugs, c("M0", Models))

statistics.matrix <- matrix(NA, ncol=(length(drugs) * (length(Models) + 1)), nrow=length(drug.association), dimnames=list(names(drug.association), paste(models.plus.drugs.names[, 1], models.plus.drugs.names[, 2], sep ="_")))

best.isoforms.matrix <- matrix(NA, ncol=length(drugs) , nrow=length(drug.association))
colnames(best.isoforms.matrix) <- drugs
rownames(best.isoforms.matrix) <- names(drug.association)

Min.Pvalues <- NULL
genes <- fData(CCLE@molecularProfiles$rnaseq)[, "EnsemblGeneId"]

if(subset.genes){
  if(subset.genes.method == "biotype"){
    ### subsetting genes to remove pseudogenes
    ff <- fData(CCLE@molecularProfiles$rnaseq)
    genes <- as.character(ff[grep("pseudogene", ff$GeneBioType, invert=TRUE), "EnsemblGeneId"])
    #genes <- as.character(ff[grep("protein_coding", ff$GeneBioType), "EnsemblGeneId"])
  }else{ #subset.genes.method == "wide.expression"
    express.cut.off <- 0.2 * nrow(ccle.genes.fpkm)
    
    xx <- apply(ccle.genes.fpkm, MARGIN=2, function(x){length(which(x > 0))})
    genes.expressed.well <- names(xx)[which(xx > express.cut.off)]
    length(genes.expressed.well) #41748
    
    tt <- apply(ccle.isoforms.fpkm, MARGIN=2, function(x){length(which(x > 0))})
    
    isoforms.expressed.well <- names(tt)[which(tt > express.cut.off)]
    length(isoforms.expressed.well) #124258
    genes.of.isoforms.expressed.well <- unique(fData(CCLE@molecularProfiles$isoforms)[match(isoforms.expressed.well, fData(CCLE@molecularProfiles$isoforms)[, "EnsemblTranscriptId"]), "EnsemblGeneId"])
    length(genes.of.isoforms.expressed.well) #30614
    genes <- intersect(genes.expressed.well, genes.of.isoforms.expressed.well)
    length(genes) #30614
  }
  FDR_List <- FDR_List[genes, , drop=FALSE]
  estimate_List <- estimate_List[genes, , drop=FALSE]
  Pvalues.Numinals <- Pvalues.Numinals[genes, , drop=FALSE]
  statistics.matrix <- statistics.matrix[genes, , drop=FALSE]
  best.isoforms.matrix <- best.isoforms.matrix[genes, , drop=FALSE]
  
  Min.Pvalues <- Min.Pvalues[genes]
  ###########################################
}
drug.association <- drug.association[genes]
drug.association.statistics <- drug.association.statistics[genes]
drug.association.best.isoforms <- drug.association.best.isoforms[genes]
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
##multiplying adjusted pvalues by 2 !!!!?
# for ( i in 1:nrow(FDR_List))
# {
#   for (j in (drugs_No + 1):ncol(FDR_List))
#   {
#     temp <- FDR_List[i,j] * 2
#     FDR_List[i,j] <- ifelse(temp > 1, 1, temp)
#   }
# }

if(file.exists("data/isoforms.no.per.gene.GRCh38.87.RData")){
  load("data/isoforms.no.per.gene.GRCh38.87.RData", verbose=TRUE)
}else{
  load("data/ensembl.map.genes.isoforms.GRCh38.87.RData")
  fnNumber_Isoforms_of_Gene <- function(Gene_Map=ensembl.map.genes.isoforms, GeneId)
  {
    return(nrow(data.frame(strsplit(Gene_Map[, as.character(GeneId)], ","))))
  }
  
  isoforms_No_List <- matrix(NA, ncol=1, nrow=nrow(fData(CCLE@molecularProfiles$rnaseq)))
  colnames(isoforms_No_List) <- "isoforms.NO"
  rownames(isoforms_No_List) <- rownames(fData(CCLE@molecularProfiles$rnaseq))
  for( i in 1:nrow(isoforms_No_List))
  {
    isoforms_No_List[i,] <- fnNumber_Isoforms_of_Gene(GeneId=as.character(rownames(isoforms_No_List)[i]))
  }
  isoforms_No_List <- data.frame(isoforms_No_List,stringsAsFactors=FALSE)
  save(isoforms_No_List, file="data/isoforms.no.per.gene.GRCh38.87.RData")  
}
isoforms_No_List <- isoforms_No_List[genes, , drop=FALSE]
save(FDR_List, estimate_List, Pvalues.Numinals, statistics.matrix, best.isoforms.matrix, Min.Pvalues, isoforms_No_List, file= file.path(path.diagrams, "allGenes_association_matrices.RData"))


##########Analyses

source("code/foo.R")

result.effect.size <- fnComputeAssociateGenes.effect.size(FDR_CutOff=0.01, effect.sizeCutOff=0.55)
write.csv(fnWilcox(result.effect.size, TRUE)$comparison, file=file.path(path.diagrams, "comparison_test_wilcox.csv"))
barplot.models(model=result.effect.size, isoforms_No="all", sign="all", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
barplot.models(model=result.effect.size, isoforms_No="1.isoform", sign="all", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
barplot.models(model=result.effect.size, isoforms_No="n.isoforms", sign="all", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
barplot.models(model=result.effect.size, isoforms_No="all", sign="positive", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
barplot.models(model=result.effect.size, isoforms_No="all", sign="negative", prototype=Models, main.title=sprintf("%s  <  1%% \n %s > 0.55", adjustment.method, effect.size), yaxis="Log", breakpoint="Regular", cex=0.7)
associations <- associations.all.drugs(model.rank="M3B", annot.ensembl.all.genes=fData(CCLE@molecularProfiles$rnaseq))
save(associations, file=file.path(path.diagrams, "associations.RData"))
##all.biomarkers would include all the significant biomarkers
all.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=0.01, BioNo="all", rank.type="pvalue.adj")
save(all.biomarkers, file=file.path(path.diagrams, "all.biomarkers.original.RData"))

if(!breast.specific)
{
  source("code/foo_training.R")
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
save(all.biomarkers, file=file.path(path.diagrams, "all.biomarkers.RData"))
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
write.csv(TOP, file=file.path(path.diagrams, "TopOne.csv"))

Check.KnownAssociations(associations)

percentage.biomarkers.type <- fnPercentageBiomarkersType(all.biomarkers)
mystacked.barplot.simple(Filename=file.path(path.diagrams, "ProportionBiomarkersType.pdf"), data=percentage.biomarkers.type, main.label="Proportion of specificity of biomarkers", cex=.7)

