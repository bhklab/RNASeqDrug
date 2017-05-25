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
effect.size <- "cindex" #c("r.squared", "cindex")
RNA_seq.normalize <- TRUE

#load(file.path(path.data, "annotation.RData"))
file.sensitivity <- "auc_recomputed_drug_association_ccle_gdsc.RData"
training.type <- "CCLE_GDSC"#ifelse(regexpr("ccle", file.sensitivity) >= 1, "CCLE", "CCLE_GDSC")
phenotype <- sensitivity.type <- paste0(unlist(strsplit(file.sensitivity, "_"))[1:2], collapse="_") 
breast.specific <- ifelse(regexpr("breast", file.sensitivity) < 1, FALSE, TRUE)

if(sensitivity.type == "slope_recomputed") {
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
  #path.diagrams <- file.path(path.result, paste(ss[1:(length(ss)-2)], collapse="_"))
  path.diagrams<- "result/auc_recomputed_ccle_gdsc"
}
load(file.path(path.data, "PSets/GRAY_hs.RData"))
gray.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GRAY, sensitivity.measure=sensitivity.type))


if(!exists("CCLE")){
  
  load(file.path(path.data, "PSets/CCLE_hs.RData"), verbose=TRUE)
  #load(file.path(path.data, "PSets/CCLE.CTRPv2.RData"), verbose=TRUE)
  #CCLE <- CCLE.CTRPv2

  load(file.path(path.data, "PSets/GDSC.RData"), verbose=TRUE)
  #load(file.path(path.data, "PSets/GDSC1000.RData"), verbose=TRUE)
  #GDSC <- GDSC1000
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


  ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", fill.missing=FALSE)))
  ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", fill.missing=FALSE)))
  
  ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=rownames(CCLE@cell), stringsAsFactors=FALSE)
  colnames(ccle.tissuetype) <- "tissue.type"
  
  ccle.isoforms.fpkm[which(is.na(ccle.isoforms.fpkm))] <- 0
  if(RNA_seq.normalize == TRUE)
  {
    ccle.genes.fpkm <- log2(ccle.genes.fpkm + 1)
    ccle.isoforms.fpkm <- log2(ccle.isoforms.fpkm + 1)
  }
}
annot.ensembl.all.genes <- Biobase::fData(CCLE@molecularProfiles$rnaseq)
drugs <- intersect(colnames(ccle.drug.sensitivity), colnames(gray.drug.sensitivity))
gray.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="rnaseq", fill.missing=FALSE)))
gray.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="isoforms", fill.missing=FALSE)))
gray.isoforms.fpkm[which(is.na(gray.isoforms.fpkm))] <- 0
if(RNA_seq.normalize == TRUE) {
  gray.genes.fpkm <- log2(gray.genes.fpkm + 1)
  gray.isoforms.fpkm <- log2(gray.isoforms.fpkm + 1)
}

load(file.path(path.diagrams, "all.biomarkers.RData"))
max.bio.no <- max(sapply(all.biomarkers, function(x){nrow(x)}))
biomarkers <- fnValidation(top.significant.biomarkers=all.biomarkers, validation.cut.off=max.bio.no, validation.method=effect.size)
save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_breast_%s_gray_%s.RData", effect.size, effect.size)))

breast <- all <- res.validated <- breast.validated.percent <- all.validated.percent <- validated.no <- NULL
for(drug in drugs) {
  xx <- which(biomarkers[[drug]][,"validation.stat"]=="validated")
  temp <- biomarkers[[drug]][xx, , drop=FALSE]
  xx <- paste(t(apply(temp, 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_")
  res.validated <-  c(res.validated, xx)
  validated <- unlist(str_split(xx, pattern="_"))
  if(all(validated == "")){lv <- 0} else{lv <- length(validated)}
  validated.no <- c(validated.no, lv)
  all.validated.percent <- c(all.validated.percent, round(lv/nrow(all.biomarkers[[drug]]), digits=2))
  all <- c(all, nrow(all.biomarkers[[drug]]))
  
  breast.validated.percent <- c(breast.validated.percent, round(lv/nrow(biomarkers[[drug]]), digits=2))
  breast <- c(breast, nrow(biomarkers[[drug]]))
}
write.csv(cbind("drug"=drugs,
                "validated"=res.validated,
                "validated no"=validated.no,
                "biomarkers no"=all,
                "ratio from all"=all.validated.percent,
                "breast significant biomarkers no"=breast,
                "ratio from breast biomarkers"=breast.validated.percent), file=file.path(path.diagrams, sprintf("validated_biomarkers_breast_%s_gray_%s.csv", effect.size, effect.size)))
##Stringtie estimated expressions are not comparable to those from cufflinks
#if(!"GTex.BR" %in% ls()){
#  load("data/GTex_BR.RData", verbose=T)
#}
#fnGtex(validation.method=effect.size)
#fnGtex(all.biomarkers=biomarkers, validation.method=validation.method)
