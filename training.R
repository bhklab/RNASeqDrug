args <- commandArgs(trailingOnly=TRUE)
##for test
#args <- c("10841", "10842", "auc_recomputed", "1", "1", "TRUE", "ccle_gdsc", "glm", "gaussian", "all", "/mnt/work1/users/bhklab/Users/zhaleh/ccledrug/_results", "TRUE", "expression.cut.off")
##
require("np") || stop("Library np is not available!")
options(echo=TRUE) # if you want see commands in output file
print(args)
source("foo_training.R")
sensitivity.method <- as.character(args[3])# c("auc_published", "auc_recomputed", "slope_recomputed", "ic50_published", "ic50_recomputed")
RNA_seq.normalize <- as.character(args[6])# c(TRUE,FALSE)
training.method <- as.character(args[7]) # c("ccle_gdsc", "gray", "ccle", "gdsc", "gCSI")
model.method <- as.character(args[8])# c("glm", "npreg", "penalized")
glm.family <- as.character(args[9])# c("binomial", "gaussian")
tissue <- NULL
if(as.character(args[10]) != "all") {tissue <- as.character(args[10])} #c(NULL, tissue.type)
output.folder <- as.character(args[11]) #output folder
subset.genes <- as.character(args[12])# c(TRUE,FALSE)
subset.genes.method <- as.character(args[13]) #c("expression.cut.off", "biotype")

statistical.method <- "bootstrap" #c("anova", "crossvalidation", "bootstrap")
effect.size <- "r.squared & cindex" #c("adj.r.squared", "r.squared", "rmsd", "cindex","r.squared & cindex")

myf <- "/mnt/work1/users/bhklab/Users/zhaleh/ccledrug/data/training_ccle_gdsc.RData"
load("/mnt/work1/users/bhklab/Users/zhaleh/ensembl.map.genes.isoforms.GRCh38.87.RData", verbose=TRUE)
if(!file.exists(myf)){
  require(PharmacoGx) || stop("Library PharmacoGx is not available!")
  require(Biobase) || stop("Library Biobase is not available!")
  #load("/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/CCLE.CTRPv2.RData", verbose=TRUE)
  #CCLE <- CCLE.CTRPv2
  load("/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/CCLE_hs.RData", verbose=TRUE)
  #load("/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/GDSC1000.RData", verbose=TRUE)
  #GDSC <- GDSC1000
  load("/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/GDSC.RData", verbose=TRUE)
  
  ###remove noisy cases and those with not matching snp profiles
  snp.outliers <- c(
    "LC-1F",
    "HCC1937",
    "MDA-MB-468",
    "HuH-7",
    "SW403",
    "COR-L51",
    "MOG-G-CCM",
    "NB4")
  
  ccle.filter <- filterNoisyCurves(CCLE, nthread=detectCores())
  lapply(ccle.filter, length)#123 noisy#11547 ok
  
  gdsc.filter <- filterNoisyCurves(GDSC, nthread=detectCores())
  lapply(gdsc.filter, length)#$noisy 2315 #ok 77588
  
  CCLE@sensitivity$info <- CCLE@sensitivity$info[ccle.filter$ok, ]
  CCLE@sensitivity$raw <- CCLE@sensitivity$raw[ccle.filter$ok, , ]
  CCLE@sensitivity$profiles <- CCLE@sensitivity$profiles[ccle.filter$ok, ]
  
  GDSC@sensitivity$info <- GDSC@sensitivity$info[gdsc.filter$ok, ]
  GDSC@sensitivity$raw <- GDSC@sensitivity$raw[gdsc.filter$ok, , ]
  GDSC@sensitivity$profiles <- GDSC@sensitivity$profiles[gdsc.filter$ok, ]
  
  match(snp.outliers, PharmacoGx::cellNames(CCLE))
  ccle.cells <- setdiff(PharmacoGx::cellNames(CCLE), snp.outliers)
  CCLE <- PharmacoGx::subsetTo(pSet=CCLE, cells=ccle.cells)
  save(CCLE, file="/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/CCLE_hs_clarified.RData")
  
  gdsc.cells <- setdiff(PharmacoGx::cellNames(GDSC), snp.outliers)
  GDSC <- PharmacoGx::subsetTo(pSet=GDSC, cells=gdsc.cells)
  save(GDSC, file="/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/GDSC_clarified.RData")
  
  if(training.method== "ccle_gdsc") {
    ##Restrict analyses to common cells
    #   common <- PharmacoGx::intersectPSet(pSets=list("CCLE"=CCLE, "GDSC"=GDSC), strictIntersect=FALSE)
    #   cells <- intersect(common$CCLE@cell$cellid, pData(common$CCLE@molecularProfiles$isoforms)[, "cellid"])
    #   common$CCLE@molecularProfiles$rnaseq <- CCLE@molecularProfiles$rnaseq[, rownames(pData(CCLE@molecularProfiles$isoforms))]
    #   CCLE <- common$CCLE
    #   GDSC <- common$GDSC
    drugs <- intersect(PharmacoGx::drugNames(CCLE), PharmacoGx::drugNames(GDSC))
  }
  ccle.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure=as.character(sensitivity.method), drugs=drugs))
  gdsc.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GDSC, sensitivity.measure=as.character(sensitivity.method), drugs=drugs))
  
  if(as.character(sensitivity.method)== "slope_recomputed") {
    ##manual cutoff for sensitivity calls based on slope
    cutoff <- 0.27
    res.weight <- as.integer(args[4]); sens.weight <- as.integer(args[5]);
    #res.weight <- 1;sens.weight <- 10
    ccle.drug.sensitivity.wighted <- ccle.drug.sensitivity
    ccle.drug.sensitivity.wighted[ccle.drug.sensitivity.wighted < cutoff] <- res.weight
    ccle.drug.sensitivity.wighted[ccle.drug.sensitivity.wighted >= cutoff] <- sens.weight
    
    gdsc.drug.sensitivity.wighted <- gdsc.drug.sensitivity
    gdsc.drug.sensitivity.wighted[gdsc.drug.sensitivity.wighted < cutoff] <- res.weight
    gdsc.drug.sensitivity.wighted[gdsc.drug.sensitivity.wighted >= cutoff] <- sens.weight
    
    ccle.drug.sensitivity <- ccle.drug.sensitivity.wighted * ccle.drug.sensitivity
    gdsc.drug.sensitivity <- gdsc.drug.sensitivity.wighted * gdsc.drug.sensitivity
  }
  ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=CCLE@cell[, "cellid"])
  colnames(ccle.tissuetype) <- "tissue.type"
  tissueTypes <- ccle.tissuetype
  
  ###
  features <- fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"]
  features <- intersect(fData(CCLE@molecularProfiles$isoforms)[,"EnsemblGeneId"], features)
  
  rnaseq.features <- rownames(fData(CCLE@molecularProfiles$rnaseq))[which(fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"] %in% features)]
  isoforms.features <- rownames(fData(CCLE@molecularProfiles$isoforms))[which(fData(CCLE@molecularProfiles$isoforms)[,"EnsemblGeneId"] %in% features)]
  
  ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", features=rnaseq.features, fill.missing=FALSE)))
  ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", features=isoforms.features, fill.missing=FALSE)))
  #ccle.isoforms.fpkm[which(is.na(ccle.isoforms.fpkm))] <- 0
  # if(RNA_seq.normalize == FALSE)
  # {
  #   ccle.genes.fpkm <- 2 ^ ccle.genes.fpkm - 1
  #   ccle.isoforms.fpkm <- 2 ^ ccle.isoforms.fpkm -1
  # }
  if(RNA_seq.normalize == TRUE)
  {
    ccle.genes.fpkm <- log2(ccle.genes.fpkm + 1)
    ccle.isoforms.fpkm <- log2(ccle.isoforms.fpkm + 1)
  }

  ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
  
  genes <- colnames(ccle.genes.fpkm)
  if(subset.genes){
    if(subset.genes.method == "biotype"){
      ### subsetting genes to remove pseudogenes
      ff <- fData(CCLE@molecularProfiles$rnaseq)
      genes <- as.character(ff[grep("pseudogene", ff$GeneBioType, invert=TRUE), "EnsemblGeneId"])
      #genes <- as.character(ff[grep("protein_coding", ff$GeneBioType), "EnsemblGeneId"])
    }else{ #subset.genes.method == "expression.cut.off"
      exprs.cut.off <- 0.1 * nrow(ccle.genes.fpkm)
      
      xx <- apply(ccle.genes.fpkm, MARGIN=2, function(x){length(which(x > 0))})
      genes.expressed <- names(xx)[which(xx > exprs.cut.off)]
      length(genes.expressed) #46467
      
      tt <- apply(ccle.isoforms.fpkm, MARGIN=2, function(x){length(which(x > 0))})
      gg <- fData(CCLE@molecularProfiles$isoforms)[match(names(tt), fData(CCLE@molecularProfiles$isoforms)[, "EnsemblTranscriptId"]), "EnsemblGeneId"]
      aa <- table(tt, gg)
      xx <- apply(aa, MARGIN=2, function(x){max(which(x > 0)) - 1})
      isoforms.expressed <- names(xx)[which(xx > exprs.cut.off)]
      length(isoforms.expressed) #30633
      genes <- intersect(genes.expressed, isoforms.expressed)
      length(genes) #30633
      ##extra checking
      isoforms.of.remained.genes <- fData(CCLE@molecularProfiles$isoforms)[which(fData(CCLE@molecularProfiles$isoforms)[, "EnsemblGeneId"] %in% genes), "EnsemblTranscriptId"]
      tt <- tt[isoforms.of.remained.genes]
      isoforms <- names(tt)[which(tt > exprs.cut.off)]
      length(isoforms)#124314
    }
  }
  ccle.genes.fpkm <- ccle.genes.fpkm[, genes, drop=FALSE]
  ccle.isoforms.fpkm <- ccle.isoforms.fpkm[, isoforms, drop=FALSE]
  GeneList <- genes
  #dim(ccle.genes.fpkm) 925 35638
  #dim(ccle.isoforms.fpkm) 925 140331
  ## for logistic regression
  
  if(glm.family == "binomial")
  {
    for(c in colnames(ccle.drug.sensitivity)){ccle.drug.sensitivity[,c] <- factor(ccle.drug.sensitivity[, c])}
    for(c in colnames(gdsc.drug.sensitivity)){gdsc.drug.sensitivity[,c] <- factor(gdsc.drug.sensitivity[, c])}
  }
  annot.ensembl.all.genes <- Biobase::fData(CCLE@molecularProfiles$rnaseq)[genes, ]
  annot.ensembl.all.isoforms <- Biobase::fData(CCLE@molecularProfiles$isoforms)[isoforms, ]
  ccle.cell.profiles <- CCLE@cell
  
  save(drugs, ccle.drug.sensitivity, gdsc.drug.sensitivity, tissueTypes, ccle.genes.fpkm, ccle.isoforms.fpkm, ccle.cells, gdsc.cells, GeneList, annot.ensembl.all.genes, annot.ensembl.all.isoforms, ccle.cell.profiles, file=myf)
}else{
  load(myf, verbose=TRUE)
}


startIndex <- as.integer(args[1])
finishIndex <- as.integer(args[2])
pvalues <- list()
statistics <- list()
best.isoforms <- list()
pvalues.r.squared <- pvalues.cindex <- pvalues
statistics.r.squared <- statistics.cindex <- statistics

for(i in startIndex: finishIndex)
{
  #Gene <- "ENSG00000171094" ALK drugs <- c("TAE684", "Crizotinib", "PLX4720", "lapatinib", "Erlotinib") expression is supposed to be predictive to TAE684
  #Gene <- ENSG00000171094 ALK c("TAE684", "Crizotinib", "PLX4720", "lapatinib", "Erlotinib")  expression is supposed to be predictive to Crizotinib
  #Gene <- "ENSG00000157764" BRAF drugs <- c("TAE684", "Crizotinib", "PLX4720", "lapatinib", "Erlotinib")  mutation is supposed to be predictive to PLX4720
  #Gene <- "ENSG00000141736" ERBB2 drugs <- c("TAE684", "Crizotinib", "PLX4720", "lapatinib", "Erlotinib") expression is supposed to be predictive to lapatinib
  #Gene <- "ENSG00000146648" EGFR drugs <- c("TAE684", "Crizotinib", "PLX4720", "lapatinib", "Erlotinib")   expression is supposed to be predictive to Erlotinib
  #Gene <- "ENSG00000146648" EGFR drugs <- c("TAE684", "Crizotinib", "PLX4720", "lapatinib", "Erlotinib")  expression is supposed to be predictive to Erlotinib
  #Gene <-  "ENSG00000181019" "1728" #NQO1 "3480" #IGF1R
  #Gene <- "ENSG00000148426" "C10orf47"
  #MicroArrayExp <- ccle.drug.microarray.exp[,Gene]; Gene_FPKM <- ccle.genes.fpkm[,Gene]; Isoforms <- fnIsoformsExp(Isoforms_FPKM=ccle.isoforms.fpkm, GeneId=Gene); GeneID <- Gene;  model.method <- "npreg"; method <- "bootstrap"; effect.size <- "adj.r.squared";glm.family <- "gaussian"; tissue <- NULL; sample.no.threshold <- 5; assay <- "ccle"
  #which(fData(common$CCLE@molecularProfiles$rnaseq)[, "EntrezGeneId"]=="1728") #18024
  Gene <- as.character(GeneList[i])
  #print(Gene)
  if(training.method == "ccle_gdsc")
  {
    sensitivity <- fnSensitivityCompare(MicroArrayExp=NULL,
                                        Gene_FPKM=ccle.genes.fpkm[,Gene],
                                        Isoforms=fnIsoformsExp(Isoforms_FPKM=ccle.isoforms.fpkm, GeneId=Gene),
                                        GeneID=Gene,
                                        drugs=drugs,
                                        effect.size=effect.size)
  }else if(training.method == "gCSI"){
    sensitivity <- fnSensitivityOneDataSet(MicroArrayExp=NULL,
                                           Gene_FPKM= gCSI.genes.fpkm[,Gene],
                                           Isoforms=fnIsoformsExp(Isoforms_FPKM=gCSI.isoforms.fpkm, GeneId=Gene),
                                           GeneID=Gene,
                                           assay=training.method,
                                           effect.size=effect.size)
  }else{
    sensitivity <- fnSensitivityOneDataSet(MicroArrayExp=NULL,
                                           Gene_FPKM= ccle.genes.fpkm[,Gene],
                                           Isoforms=fnIsoformsExp(Isoforms_FPKM=ccle.isoforms.fpkm, GeneId=Gene),
                                           GeneID=Gene,
                                           assay=training.method,
                                           effect.size=effect.size)
  }
  if(effect.size == "r.squared & cindex") {
    pvalues.r.squared[[Gene]] <- sensitivity$p.values.r.squared
    pvalues.cindex[[Gene]] <- sensitivity$p.values.cindex
    statistics.r.squared[[Gene]] <- sensitivity$statistics.r.squared
    statistics.cindex[[Gene]] <- sensitivity$statistics.cindex
    best.isoforms[[Gene]] <- sensitivity$best.isoforms
  }else{
    pvalues[[Gene]] <- sensitivity$p.values
    statistics[[Gene]] <- sensitivity$statistics
    best.isoforms[[Gene]] <- sensitivity$best.isoforms
  }
  print(sprintf("Models of Gene %s [%d] is built at %s", Gene, i, Sys.time()))
}
if(effect.size == "r.squared & cindex") {
  both.drug.association.adj.r.squared.pvalues <- list("r.squared"=pvalues.r.squared, "cindex"=pvalues.cindex)
  both.drug.association.statistics <- list("r.squared"=statistics.r.squared, "cindex"=statistics.cindex)
}else{
  both.drug.association.adj.r.squared.pvalues <- pvalues
  both.drug.association.statistics <- statistics
}  
both.drug.association.best.isoforms <- best.isoforms
path.result <- file.path(output.folder)
if (!file.exists(path.result)){dir.create(file.path(path.result))}
save(both.drug.association.adj.r.squared.pvalues, both.drug.association.statistics,both.drug.association.best.isoforms, file=file.path(path.result ,paste0(paste(as.character(startIndex),as.character(finishIndex),sep="_"), ".RData")))

