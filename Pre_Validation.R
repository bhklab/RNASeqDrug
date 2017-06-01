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

effect.size.cut.off <- 0.55
pvalue.cut.off <- 0.05
source(file.path("code","foo_PreValidation.R"))

if(!exists("ccle.genes.fpkm")){
  path.training.data <- "data/training_ccle_gdsc.RData"
  load(path.training.data, verbose=TRUE)
  if("gdsc.drug.sensitivity" %in% ls()) {
    training.type <-"CCLE_GDSC"
  } else {
    training.type <-"CCLE"
  }
  path.diagrams<- "result/auc_recomputed_ccle_gdsc"
  tissue <- "breast"
  model.method <- "glm"
  glm.family <- "gaussian" 
  effect.size <- "cindex" #c("r.squared", "cindex")
  path.data <- "data"
  path.code <- file.path("code")
  path.result <- file.path("result")
  adjustment.method <- "fdr"
  source(file.path(path.code, "foo.R"))
}
sensitivity.type <- phenotype <- "auc_recomputed"
RNA_seq.normalize <- TRUE
genes <- colnames(ccle.genes.fpkm)
isoforms <- colnames(ccle.isoforms.fpkm)

load(file.path(path.data, "PSets/GRAY_hs.RData"))
gray.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GRAY, sensitivity.measure=sensitivity.type))
drugs <- intersect(colnames(ccle.drug.sensitivity), colnames(gray.drug.sensitivity))

gray.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="rnaseq", features=genes, fill.missing=FALSE)))
gray.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="isoforms", features=isoforms, fill.missing=FALSE)))
gray.isoforms.fpkm[which(is.na(gray.isoforms.fpkm))] <- 0
if(RNA_seq.normalize == TRUE) {
  gray.genes.fpkm <- log2(gray.genes.fpkm + 1)
  gray.isoforms.fpkm <- log2(gray.isoforms.fpkm + 1)
}
gray.cells <- intersect(rownames(gray.drug.sensitivity), rownames(gray.genes.fpkm))
gray.drug.sensitivity <- gray.drug.sensitivity[gray.cells, , drop=FALSE]
gray.genes.fpkm <- gray.genes.fpkm[gray.cells, , drop=FALSE]
gray.isoforms.fpkm <- gray.isoforms.fpkm[gray.cells, , drop=FALSE]

load(file.path(path.diagrams, "all.biomarkers.RData"))
max.bio.no <- max(sapply(all.biomarkers, function(x){nrow(x)}))
biomarkers <- fnValidation(top.significant.biomarkers=all.biomarkers, validation.cut.off=max.bio.no, validation.method=effect.size)
save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_breast_%s_gray_%s.original.RData", effect.size, "pvalue")))

breast <- all <- res.validated <- breast.validated.percent <- all.validated.percent <- validated.no <- NULL
for(drug in drugs) {
  xx <- which(biomarkers[[drug]][,"validation.stat"]=="validated")
  temp <- biomarkers[[drug]][xx, , drop=FALSE]
  if(length(xx) > 0) {
    xx <- paste(t(apply(temp, 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_")
    res.validated <-  c(res.validated, xx)
  } else {
    res.validated <-  c(res.validated, "")
  }
  validated <- unlist(str_split(xx, pattern="_"))
  if(all(validated == "")){lv <- 0} else{lv <- length(validated)}
  validated.no <- c(validated.no, lv)
  all.validated.percent <- c(all.validated.percent, round(lv/nrow(all.biomarkers[[drug]]), digits=2))
  all <- c(all, length(which(rownames(all.biomarkers[[drug]]) != "NA")))
  
  if(!is.null(biomarkers[[drug]])){
    breast.validated.percent <- c(breast.validated.percent, round(lv/nrow(biomarkers[[drug]]), digits=2))
    breast <- c(breast, nrow(biomarkers[[drug]]))
  }else{
    breast.validated.percent <- c(breast.validated.percent, 0)
    breast <- c(breast, 0)
  }
}
write.csv(cbind("drug"=drugs,
                "validated"=res.validated,
                "validated no"=validated.no,
                "biomarkers no"=all,
                "ratio from all"=all.validated.percent,
                "breast significant biomarkers no"=breast,
                "ratio from breast biomarkers"=breast.validated.percent), file=file.path(path.diagrams, sprintf("validated_biomarkers_breast_%s_gray_%s.csv", effect.size, "pvalue")))
##Stringtie estimated expressions are not comparable to those from cufflinks
#if(!"GTex.BR" %in% ls()){
#  load("data/GTex_BR.RData", verbose=T)
#}
#fnGtex(validation.method=effect.size)
#fnGtex(all.biomarkers=biomarkers, validation.method=validation.method)
for(drug in drugs) {
  biomarkers[[drug]] <- biomarkers[[drug]][which(biomarkers[[drug]]$validation.stat == "validated"),]
}
source("code/foo_training.R")
for(drug in drugs) {
  if(nrow(biomarkers[[drug]]) > 0){
    tt <- which(biomarkers[[drug]][,"transcript.id"]!= "")
    xx <- biomarkers[[drug]][tt, ]
    biomarkers[[drug]][which(biomarkers[[drug]][,"transcript.id"]== ""), "gray.specificity"] <- "gene.specific"
    if(!is.null(xx)){
      p.values.isoform <- p.values.gene <- NULL
      for(i in 1:nrow(xx)) {
        M0 <- fnCreateNullModel(drug=drug, assay="gray")
        M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gray.genes.fpkm[ ,xx[i, "gene.id"]])
        M3B <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gray.isoforms.fpkm[ ,xx[i, "transcript.id"]])
        results <- fnRunbootstrap(models=list("M2"=M2, "M3B"=M3B), effect.size=effect.size)
        if(length(results[["M2"]])>0 && length(results[["M3B"]])>0){
          p.values.isoform <- c(p.values.isoform, wilcox.test(results[["M3B"]], results[["M2"]], paired=TRUE, alternative="greater")$p.value)
          p.values.gene <- c(p.values.gene, wilcox.test(results[["M3B"]], results[["M2"]], paired=TRUE, alternative="less")$p.value)
        }
        if(length(results[["M2"]])>0 && length(results[["M3B"]])==0){
          p.values.isoform <- c(p.values.isoform * 2, 1)
          p.values.gene <- c(p.values.gene * 2, 0.00001)
        }
        if(length(results[["M2"]])==0 && length(results[["M3B"]])>0){
          p.values.isoform <- c(p.values.isoform * 2, 0.00001)
          p.values.gene <- c(p.values.gene * 2, 1)
        }
        if(length(results[["M2"]])==0 && length(results[["M3B"]])==0){
          p.values.isoform <- c(p.values.isoform * 2, 1)
          p.values.gene <- c(p.values.gene * 2, 1)
        }
      }
      names(p.values.isoform) <- names(p.values.gene) <- rownames(xx)
      biomarkers[[drug]][tt, "gray.specificity"] <- "common"
      for( i in names(p.values.isoform)){
        if(p.values.isoform[i] < 0.05) {
          biomarkers[[drug]][i, "gray.specificity"] <- "isoform.specific"
        }else if(p.values.gene[i] < 0.05) {
          biomarkers[[drug]][i, "gray.specificity"] <- "gene.specific"
        }
        biomarkers[[drug]][i, "gray.isoform.specific.test.pvalue"] <- p.values.isoform[i]
        biomarkers[[drug]][i, "gray.gene.specific.test.pvalue"] <- p.values.gene[i]
      }
    }
  }
}

save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_validated_breast_%s_gray_%s.RData", effect.size, "pvalue")))
#load("result/auc_recomputed_ccle_gdsc/Biomarkers_validated_breast_cindex_gray_pvalue.RData")
#unlist(sapply(biomarkers, function(x){if(nrow(x)>0){length(which(x["gray.specificity"]!="gene.specific"))}}))
#for(i in names(biomarkers)){
  if(nrow(biomarkers[[i]]) > 0){
    xx <- biomarkers[[i]]$gray.isoform.specific.test.pvalue
    xx[which(xx != 1 & xx != 0.00001)] <- xx[which(xx != 1 & xx != 0.00001)] * 2
    yy <- biomarkers[[i]]$gray.gene.specific.test.pvalue
    yy[which(yy != 1 & yy != 0.00001)] <- yy[which(yy != 1 & yy != 0.00001)] * 2
    
    biomarkers[[i]]$gray.isoform.specific.test.pvalue <- xx
    biomarkers[[i]]$gray.gene.specific.test.pvalue <- yy
    
    biomarkers[[i]]$gray.specificity <- "common"
    biomarkers[[i]]$gray.specificity[which(xx < 0.05)] <- "isoform.specific"
    biomarkers[[i]]$gray.specificity[which(yy < 0.05)] <- "gene.specific"
  }
}
#unlist(sapply(biomarkers, function(x){if(nrow(x)>0){length(which(x["gray.specificity"]!="gene.specific"))}}))
# load(file.path(path.diagrams, "Biomarkers_with_validation_status_breast_cindex_gray_pvalue.original.RData"), verbose=TRUE)
# sapply(biomarkers, dim)
# for(drug in names(biomarkers)) {
#   xx <- p.adjust(biomarkers[[drug]]$gray.pvalue, method="fdr")
#   biomarkers[[drug]][,"gray.fdr"] <- xx
# }
# rr <- list()
# for(drug in names(biomarkers)) {
#   message(sprintf("%s : %s", drug, length(which(biomarkers[[drug]][,"gray.fdr"]< 0.01 & sign(biomarkers[[drug]][,"estimate"]) == sign(biomarkers[[drug]][,"gray.estimate"])))))
#   rr[[drug]] <- biomarkers[[drug]][which(biomarkers[[drug]][,"gray.fdr"]< 0.01 & sign(biomarkers[[drug]][,"estimate"]) == sign(biomarkers[[drug]][,"gray.estimate"])),]
# }
# for(drug in names(biomarkers)) {
#   message(sprintf("%s : %s", drug, length(which(biomarkers[[drug]][,"gray.pvalue"]< 0.05 & sign(biomarkers[[drug]][,"estimate"]) == sign(biomarkers[[drug]][,"gray.estimate"])))))
# }
