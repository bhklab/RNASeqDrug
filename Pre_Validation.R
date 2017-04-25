#file.sensitivity <- "breast_new_drug_association.RData"
#load(file.path(path.data, file.sensitivity), verbose=T)


#both.drug.association <- both.tissue.drug.association
#both.drug.association.statistics <- both.tissue.drug.association.statistics


#path.diagrams <- file.path("result/slope0_weightes_1_1")

# load(file.diagrams, "all.biomarkers.RData", verbose=T)
#require(gdata) || stop("Library gdata is not available!")
#myfn <- file.path(path.diagrams, "all.biomarkers.xlsx")
#all.biomarkers <- gdata::read.xls(xls=myfn, sheet=12, stringsAsFactors=FALSE)

require(PharmacoGx) || stop("Library PharmacoGx is not available!")
require(Biobase) || stop("Library Biobase is not available!")
require(calibrate) || stop("Library calibrate is not available!")
require(stringr) || stop("Library stringr is not available!")
require(gdata) || stop("Library gdata is not available!")
require(genefu) || stop("Library genefu is not available!")
require(xtable) || stop("Library xtable is not available!")

options(stringsAsFactors=FALSE)
path.data <- "data"
path.code <- file.path("code")
path.result <- file.path("result")
adjustment.method <- "fdr"

source(file.path(path.code, "foo.R"))
source(file.path(path.code,"foo_PreValidation.R"))


tissue <- "breast"
model.method <- "glm"
glm.family <- "gaussian" 
stat <- "cindex" #c("r.squared", "cindex")
RNA_seq.normalize <- TRUE

#load(file.path(path.data, "annotation.RData"))
file.sensitivity <- "auc_recomputed_drug_association.RData"
training.type <- ifelse(regexpr("ccle", file.sensitivity) >= 1, "CCLE", "CCLE_GDSC")
phenotype <- sensitivity.type <- paste0(unlist(strsplit(file.sensitivity, "_"))[1:2], collapse="_") 
breast.specific <- ifelse(regexpr("breast", file.sensitivity) < 1, FALSE, TRUE)

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
load(file.path(path.data, "PSets/GRAY_hs.RData"))
annot.ensembl.all.genes <- Biobase::fData(CCLE@molecularProfiles$rnaseq)
gray.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GRAY, sensitivity.measure=sensitivity.type))

load(file.path(path.data, "PSets/CCLE_hs.RData"))
load(file.path(path.data, "PSets/GDSC.RData"))
if(training.type == "CCLE_GDSC") {
#  common <- PharmacoGx::intersectPSet(pSets=list("CCLE"=CCLE, "GDSC"=GDSC), strictIntersect=F)
#  cells <- intersect(common$CCLE@cell$cellid, pData(common$CCLE@molecularProfiles$isoforms)[,"cellid"])
#  common$CCLE@molecularProfiles$rnaseq <- CCLE@molecularProfiles$rnaseq[, rownames(pData(CCLE@molecularProfiles$isoforms))]
#  CCLE <- common$CCLE
#  GDSC <- common$GDSC
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

drugs <- intersect(colnames(ccle.drug.sensitivity), colnames(gray.drug.sensitivity))

ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", fill.missing=FALSE)))
ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", fill.missing=FALSE)))

ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=rownames(CCLE@cell), stringsAsFactors=FALSE)
colnames(ccle.tissuetype) <- "tissue.type"

gray.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="rnaseq", fill.missing=FALSE)))
gray.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="isoforms", fill.missing=FALSE)))

ccle.isoforms.fpkm[which(is.na(ccle.isoforms.fpkm))] <- 0
gray.isoforms.fpkm[which(is.na(gray.isoforms.fpkm))] <- 0
if(RNA_seq.normalize == TRUE)
{
  ccle.genes.fpkm <- log2(ccle.genes.fpkm + 1)
  ccle.isoforms.fpkm <- log2(ccle.isoforms.fpkm + 1)

  gray.genes.fpkm <- log2(gray.genes.fpkm + 1)
  gray.isoforms.fpkm <- log2(gray.isoforms.fpkm + 1)
}
#associations <- associations.all.drugs("M3B")
#top.significant.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=.5, BioNo=100, rank.type="pvalue")
#Validation.result <- fnValidation(all.biomarkers, validation.cut.off=5)

load(file.path(path.diagrams, "all.biomarkers.RData"))
max.bio.no <- max(sapply(all.biomarkers, function(x){nrow(x)}))
Validation.result <- fnValidation(top.significant.biomarkers=all.biomarkers, validation.cut.off=max.bio.no)

WriteXLS::WriteXLS("Validation.result", ExcelFileName=file.path(path.diagrams, "top.biomarkers.gray.xlsx"), row.names=TRUE)
save(Validation.result, file=file.path(path.diagrams, "top.biomarkers.gray.RData"))

for(validation.method in c(stat, "pvalue")){
  res.validated <- res.validated.breast <- res.validated.breast.boot <- NULL
  res.validated.percent <- res.validated.percent.breast <- res.validated.percent.breast.boot <- NULL
  all <- all.breast <- all.breast.boot <- NULL
  
  for(i in 1:length(Validation.result))
  {
    temp <- data.frame(lapply(Validation.result[[i]], as.character), stringsAsFactors=FALSE)
    #res.validated <- c(res.validated, ifelse("TRUE" %in% names(table(Validation.result[[i]]$all.gray.pvalue < .05)), table(Validation.result[[i]]$all.gray.pvalue < .05)["TRUE"], 0))  
    
    res.validated <-  c(res.validated, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                              sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate)), select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
    validated <- unlist(str_split(res.validated[i], pattern="_"))
    if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
    res.validated.percent <- c(res.validated.percent, round(lv/nrow(temp), digits=2))
    all <- c(all, nrow(temp))
    
    if(validation.method == stat) {
      res.validated.breast <-  c(res.validated.breast, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                                                  sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate))& 
                                                                                  #as.numeric(temp$stat) <= as.numeric(temp$breast),
                                                                                  as.numeric(temp$breast) >= .55,
                                                                                  select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
      breast <- nrow(subset(temp, as.numeric(temp$breast) >= .55))
    }else{
      res.validated.breast <-  c(res.validated.breast, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                                                sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate))& 
                                                                                #as.numeric(temp$pvalue) >= as.numeric(temp[ , paste0(tissue, "_pvalue")]),
                                                                                as.numeric(temp[ , paste0(tissue, "_pvalue")]) < .05,
                                                                                select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
      breast <- nrow(subset(temp, as.numeric(temp[ , paste0(tissue, "_pvalue")]) < .05))
    }
    validated <- unlist(str_split(res.validated.breast[i], pattern="_"))
    if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
#    breast <- nrow(subset(temp, as.numeric(temp$stat) <= as.numeric(temp$breast)))
    
    res.validated.percent.breast <- c(res.validated.percent.breast, round(lv/breast, digits=2))
    all.breast <- c(all.breast, breast)
    
    if("breast_boot" %in% colnames(temp)){
      res.validated.breast.boot <-  c(res.validated.breast.boot, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                                                      sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate))& 
                                                                                      as.numeric(temp$stat) <= as.numeric(temp$breast_boot),
                                                                                    select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
      validated <- unlist(str_split(res.validated.breast.boot[i], pattern="_"))
      if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
      breast.boot <- nrow(subset(temp, as.numeric(temp$stat) <= as.numeric(temp$breast_boot)))
      res.validated.percent.breast.boot <- c(res.validated.percent.breast.boot, round(lv/breast.boot, digits=2))
    
    
      all.breast.boot <- c(all.breast.boot, breast.boot)
    }
  }
  names(res.validated) <- names(Validation.result)
  if("breast_boot" %in% colnames(res.validated[[1]])){
    write.csv(cbind("validated"=res.validated, 
                    "ratio"=res.validated.percent, 
                    "biomarkers no"=all,
                    "validated.breast"=res.validated.breast, 
                    "ratio.breast"=res.validated.percent.breast, 
                    "biomarkers no.breast"=all.breast,
                    "validated.boot.breast"=res.validated.breast.boot, 
                    "ratio.boot.breast"=res.validated.percent.breast.boot, 
                    "biomarkers no.boot.breast"=all.breast.boot), file=file.path(path.diagrams, sprintf("validated_biomarkers_%s.csv", validation.method)))
  }else{
    write.csv(cbind("validated"=res.validated, 
                    "ratio"=res.validated.percent, 
                    "biomarkers no"=all,
                    "validated.breast"=res.validated.breast, 
                    "ratio.breast"=res.validated.percent.breast, 
                    "biomarkers no.breast"=all.breast), file=file.path(path.diagrams, sprintf("validated_biomarkers_%s.csv", validation.method)))
    
  }
  if(!"GTex.BR" %in% ls()){
    load("data/GTex_BR.RData", verbose=T)
  }
  #fnGtex(validation.method=stat)
  biomarkers <- fnGtex(validation.method=validation.method)
  save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.RData", validation.method)))
}

