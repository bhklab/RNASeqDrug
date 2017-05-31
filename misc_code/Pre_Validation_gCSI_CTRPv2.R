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
require(survcomp) || stop("Library survcomp is not available!")
options(stringsAsFactors=FALSE)
path.data <- "data"
path.code <- file.path("code")
path.result <- file.path("result")
adjustment.method <- "fdr"

source(file.path(path.code, "foo.R"))
source(file.path(path.code,"foo_PreValidation_CTRPv2.R"))


tissue <- NULL
model.method <- "glm"
glm.family <- "gaussian" 
#stat <- "R2"
stat <- "cindex"

# load(file.path(path.data, "annotation.RData"))
file.sensitivity <- "auc_recomputed_drug_association_ccle_gdsc.RData"
training.type <- ifelse(regexpr("gdsc", file.sensitivity) < 1, "CCLE", "CCLE_GDSC")
phenotype <- sensitivity.type <- paste0(unlist(strsplit(file.sensitivity, "_"))[1:2], collapse="_") 
breast.specific <- ifelse(regexpr("breast", file.sensitivity) < 1, FALSE, TRUE)

if(sensitivity.type == "slope_recomputed"){
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
load(file.path(path.data, "PSets/gCSI_hs.RData"))

tissueCounts <- table(cellInfo(gCSI)$tissue)
keep.tissues <- names(tissueCounts)[tissueCounts>2]
keep.cells <- cellNames(gCSI)[cellInfo(gCSI)$tissueid%in%keep.tissues]
keep.cells2 <- unique(gCSI@sensitivity$info$cellid)
keep.cells <- intersect(keep.cells, keep.cells2)

gCSI <- subsetTo(gCSI, cells=keep.cells)

gCSI.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=gCSI, sensitivity.measure=sensitivity.type))

load(file.path(path.data, "PSets/CCLE.GDSC.RData"))
load(file.path(path.data, "PSets/CCLE.CTRPv2.RData"))
if(training.type == "CCLE_GDSC") {
#  common <- PharmacoGx::intersectPSet(pSets=list("CCLE"=CCLE, "GDSC"=GDSC), strictIntersect=F)
#  cells <- intersect(common$CCLE@cell$cellid, pData(common$CCLE@molecularProfiles$isoforms)[,"cellid"])
#  common$CCLE@molecularProfiles$rnaseq <- CCLE@molecularProfiles$rnaseq[, rownames(pData(CCLE@molecularProfiles$isoforms))]
#  CCLE <- common$CCLE
#  GDSC <- common$GDSC
  drugs <- intersect(PharmacoGx::drugNames(CCLE.GDSC), PharmacoGx::drugNames(CCLE.CTRPv2))
}

CCLE <- CCLE.CTRPv2
GDSC <- CCLE.GDSC

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

drugs <- intersect(colnames(ccle.drug.sensitivity), colnames(gCSI.drug.sensitivity))

ccle.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="rnaseq", fill.missing=FALSE)))
ccle.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=CCLE, mDataType="isoforms", fill.missing=FALSE)))

ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=CCLE@cell[, "cellid"], stringsAsFactors=FALSE)
colnames(ccle.tissuetype) <- "tissue.type"


gCSI.tissuetype <- as.data.frame(gCSI@cell[, "tissueid"], row.names=rownames(gCSI@cell), stringsAsFactors=FALSE)
colnames(gCSI.tissuetype) <- "tissue.type"

tissueTypes <- gCSI.tissuetype 

gCSI.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=gCSI, mDataType="rnaseq", fill.missing=FALSE)))
gCSI.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=gCSI, mDataType="isoforms", fill.missing=FALSE)))
gCSI.genes.fpkm[is.na(gCSI.genes.fpkm)] <- 0
gCSI.isoforms.fpkm[is.na(gCSI.isoforms.fpkm)] <- 0

#associations <- associations.all.drugs("M3B")
#top.significant.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=.5, BioNo=100, rank.type="pvalue")
#Validation.result <- fnValidation(all.biomarkers, validation.cut.off=5)

load(file.path(path.diagrams, "all.biomarkers.ctrpv2.RData"))
max.bio.no <- max(sapply(all.biomarkers, function(x){nrow(x)}))
annot.ensembl.all.genes <- fData(CCLE@molecularProfiles$rnaseq)
annot.ensembl.all.genes$gene_id <- annot.ensembl.all.genes$EnsemblGeneId
annot.ensembl.all.genes$gene_name <- annot.ensembl.all.genes$Symbol
annot.ensembl.all.genes$gene_biotype <- annot.ensembl.all.genes$GeneBioType

common.cells.ccle <- intersect(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm))
common.cells.gCSI <- intersect(rownames(gCSI.drug.sensitivity), rownames(gCSI.genes.fpkm))

ccle.genes.fpkm <- ccle.genes.fpkm[common.cells.ccle,]
ccle.isoforms.fpkm <- ccle.isoforms.fpkm[common.cells.ccle,]
ccle.drug.sensitivity <- ccle.drug.sensitivity[common.cells.ccle,]

gCSI.genes.fpkm <- gCSI.genes.fpkm[common.cells.gCSI,]
gCSI.isoforms.fpkm <- gCSI.isoforms.fpkm[common.cells.gCSI,]
gCSI.drug.sensitivity <- gCSI.drug.sensitivity[common.cells.gCSI,]


myx <- apply(gCSI.genes.fpkm, 2, function(x) return(!all(is.na(x)&&!all(x==0))))
gCSI.genes.fpkm <- gCSI.genes.fpkm[,myx]
myx <- apply(gCSI.isoforms.fpkm, 2, function(x) return(!all(is.na(x)&&!all(x==0))))
gCSI.isoforms.fpkm <- gCSI.isoforms.fpkm[,myx]

tissue <- "all"

Validation.result <- fnValidation(top.significant.biomarkers=all.biomarkers, validation.cut.off=max.bio.no, stat)
tissue <- NULL

WriteXLS::WriteXLS("Validation.result", ExcelFileName=file.path(path.diagrams, "CTRPv2.top.biomarkers.gCSI.xlsx"), row.names=TRUE)
# "pvalue" validation maybe broken?
for(validation.method in c(stat)){ 
  # res.validated <- res.validated.breast <- res.validated.breast.boot <- NULL
  # res.validated.percent <- res.validated.percent.breast <- res.validated.percent.breast.boot <- NULL
  # all <- all.breast <- all.breast.boot <- NULL
  
  # for(i in 1:length(Validation.result))
  # {
  #   temp <- data.frame(lapply(Validation.result[[i]], as.character), stringsAsFactors=FALSE)
  #   #res.validated <- c(res.validated, ifelse("TRUE" %in% names(table(Validation.result[[i]]$all.gCSI.pvalue < .05)), table(Validation.result[[i]]$all.gCSI.pvalue < .05)["TRUE"], 0))  
    
  #   res.validated <-  c(res.validated, paste(t(apply(subset(temp, temp$all.gCSI.pvalue < .05 & 
  #                                                             sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gCSI.estimate)), select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
  #   validated <- unlist(str_split(res.validated[i], pattern="_"))
  #   if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
  #   res.validated.percent <- c(res.validated.percent, round(lv/nrow(temp), digits=2))
  #   all <- c(all, nrow(temp))
  #   # what is going on here????
  #   if(validation.method == stat) {
  #     res.validated.breast <-  c(res.validated.breast, paste(t(apply(subset(temp, temp$all.gCSI.pvalue < .05 & 
  #                                                                                 sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gCSI.estimate))& 
  #                                                                                 as.numeric(temp[,stat]) <= as.numeric(temp$breast), 
  #                                                                                 select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
  #   }else{ 
  #     res.validated.breast <-  c(res.validated.breast, paste(t(apply(subset(temp, temp$all.gCSI.pvalue < .05 & 
  #                                                                               sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gCSI.estimate))& 
  #                                                                               as.numeric(temp[ , paste0(tissue, "_pvalue")]) <= as.numeric(temp$breast),
  #                                                                               select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
  #   }
  #   validated <- unlist(str_split(res.validated.breast[i], pattern="_"))
  #   if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
  #   breast <- nrow(subset(temp, as.numeric(temp[,stat]) <= as.numeric(temp$breast)))
    
  #   res.validated.percent.breast <- c(res.validated.percent.breast, round(lv/breast, digits=2))
  #   all.breast <- c(all.breast, breast)
    
  #   res.validated.breast.boot <-  c(res.validated.breast.boot, paste(t(apply(subset(temp, temp$all.gCSI.pvalue < .05 & 
  #                                                                                     sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gCSI.estimate))& 
  #                                                                                     as.numeric(temp[,stat]) <= as.numeric(temp$breast_boot),
  #                                                                                   select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
  #   validated <- unlist(str_split(res.validated.breast.boot[i], pattern="_"))
  #   if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
  #   breast.boot <- nrow(subset(temp, as.numeric(temp[,stat]) <= as.numeric(temp$breast_boot)))
  #   res.validated.percent.breast.boot <- c(res.validated.percent.breast.boot, round(lv/breast.boot, digits=2))
    
    
  #   all.breast.boot <- c(all.breast.boot, breast.boot)
  # }
  # names(res.validated) <- names(Validation.result)
  # write.csv(cbind("validated"=res.validated, 
  #                 "ratio"=res.validated.percent, 
  #                 "biomarkers no"=all,
  #                 "validated.breast"=res.validated.breast, 
  #                 "ratio.breast"=res.validated.percent.breast, 
  #                 "biomarkers no.breast"=all.breast,
  #                 "validated.boot.breast"=res.validated.breast.boot, 
  #                 "ratio.boot.breast"=res.validated.percent.breast.boot, 
  #                 "biomarkers no.boot.breast"=all.breast.boot), file=file.path(path.diagrams, sprintf("validated_biomarkers_%s.csv", validation.method)))
  # load("data/GTex_BR.RData", verbose=T)
  #fnPanCancer(validation.method=stat)
  biomarkers <- fnPanCancer(validation.method=validation.method)
  save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.ctrpv2.RData", validation.method)))
}


## Heatmaps, effect sizes and sensitivity plots for tarining and pre validation
biomarkers.category <- "all"
validation.method <- "cindex"
tissue.type <- "all"

mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
red <- mycol[1]  
blue <- mycol[2]

model.method <- "glm"
glm.family <- "gaussian"
stat <- "cindex"

load(file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.ctrpv2.RData", validation.method)), verbose=TRUE)
biomarkers.temp <- biomarkers
bb <- list()
for(drug in names(biomarkers)) {
  drug.name <- drug
  print(paste0("Now processing drug: ", drug))
  tt <- biomarkers[[drug]]
  if(!"gtex.specificity" %in% colnames(tt)){
    tt <- cbind(tt, "gCSI.specificity"=NA)
  }
  stat.colname <- paste("gCSI", stat, sep=".")
  delta.stat.colname <- paste("gCSI", "delta", stat, sep=".")

  # tt <- cbind(tt, stat.colname=NA)
  #tt <- subset(tt,  sign(tt$estimate) >=0)
  #biomarkers.toPlot <- c("ITGA6", "S100A11", "SH3TC2", "MBP")
  ## select those expressed in CCLE 
  # tt <- subset(tt,  tt$ccle > 0)
  # vtt <- subset(tt, tt$validation.stat == "validated")
  vtt <- tt
  tt[,"bootstrap.validation.stat"] <- "unvalidated"
  if(!is.null(vtt) && nrow(vtt) > 0){
    biomarkers <- fnFetchBiomarkers(top.significant.biomarkers=vtt, drug=drug, indices=1:nrow(vtt))
    if(biomarkers.category == "isoforms") {
      temp <- do.call(rbind, biomarkers)
      temp <- temp[which(temp[,"type"] == "isoform"), , drop=FALSE]
      if(nrow(temp) > 0)
      {
        temp[,"short.label"] <- gsub(".ISO$","",temp[,"short.label"])
        biomarkers <- apply(temp, 1, function(x){x})
      }else {
        biomarkers <- NULL
      }
    }
    if(!is.null(biomarkers)){
      p.values <- NULL
      M0 <- fnCreateNullModel(drug=drug, assay="gCSI")
      for(i in 1:length(biomarkers)){
        if((i%%25)==0){
           print(paste0("Completed: ",i/length(biomarkers), "%"))
        }
        # M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gCSI.genes.fpkm[ ,biomarkers[[i]]$gene.id])
        # M3B <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gCSI.isoforms.fpkm[ ,biomarkers[[i]]$isoform.id])
        # if(biomarkers[[i]]$type =="isoform"){
        #   results <- fnRunbootstrap(models=list("M0"=M0, "M2"=M2, "M3B"=M3B))
        #   p.values <- c(p.values, wilcox.test(results[["M3B"]], results[["M2"]], paired=TRUE, alternative="greater")$p.value)
        # }
        p.values <- c(p.values, 1)
        if(biomarkers[[i]]$type == "gene"){
          M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gCSI.genes.fpkm[ ,biomarkers[[i]]$gene.id])
          if(!is.null(M2$model)){
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,stat.colname] <- compute.stat(stat, yi=M2$dataset$drug, y.hat=fitted(M2$model))
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,delta.stat.colname] <- tt[tt$biomarker.id==biomarkers[[i]]$symbol,stat.colname] - compute.stat(stat, yi=M0$dataset$drug, y.hat=fitted(M0$model))
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,"gCSI.R2"] <- compute.stat("r.squared", yi=M2$dataset$drug, y.hat=fitted(M2$model))
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,"delta.gCSI.R2"] <- tt[tt$biomarker.id==biomarkers[[i]]$symbol,"gCSI.R2"] - compute.stat("r.squared", yi=M0$dataset$drug, y.hat=fitted(M0$model))
            results <- fnRunbootstrap(models=list("M0"=M0, "M2"=M2))
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,"gCSI.bootstrap.pvalue"] <- wilcox.test(results[["M2"]], results[["M0"]], paired=TRUE, alternative="greater")$p.value
            if(tt[tt$biomarker.id==biomarkers[[i]]$symbol,"gCSI.bootstrap.pvalue"]<=0.05){
              tt[tt$biomarker.id==biomarkers[[i]]$symbol,"bootstrap.validation.stat"] <- "validated"
            }
          }
        }
        if(biomarkers[[i]]$type == "isoform"){
          M3B <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gCSI.isoforms.fpkm[ ,biomarkers[[i]]$isoform.id])
          if(!is.null(M3B$model)){
            tt[tt$biomarker.id==biomarkers[[i]]$id,stat.colname] <- compute.stat(stat, yi=M3B$dataset$drug, y.hat=fitted(M3B$model))
            tt[tt$biomarker.id==biomarkers[[i]]$id,delta.stat.colname] <- tt[tt$biomarker.id==biomarkers[[i]]$id,stat.colname] - compute.stat(stat, yi=M0$dataset$drug, y.hat=fitted(M0$model))        
            tt[tt$biomarker.id==biomarkers[[i]]$id,"gCSI.R2"] <- compute.stat("r.squared", yi=M3B$dataset$drug, y.hat=fitted(M3B$model))
            tt[tt$biomarker.id==biomarkers[[i]]$id,"delta.gCSI.R2"] <- tt[tt$biomarker.id==biomarkers[[i]]$id,"gCSI.R2"] - compute.stat("r.squared", yi=M0$dataset$drug, y.hat=fitted(M0$model))
            results <- fnRunbootstrap(models=list("M0"=M0, "M3B"=M3B))
            tt[tt$biomarker.id==biomarkers[[i]]$id,"gCSI.bootstrap.pvalue"] <- wilcox.test(results[["M3B"]], results[["M0"]], paired=TRUE, alternative="greater")$p.value        
            if(tt[tt$biomarker.id==biomarkers[[i]]$id,"gCSI.bootstrap.pvalue"]<=0.05){
              tt[tt$biomarker.id==biomarkers[[i]]$id,"bootstrap.validation.stat"] <- "validated"
            }
          }
        }
      }
      # names(p.values) <- sapply(biomarkers, function(x){x[["id"]]})
      ## Plot Prevalidated biomarkers heatmap and sensitivity plot, cell lines are ordered by sensitivity measurement and similar biomarkers are clustered together
      # rr <- fnPlotAUCoverCellLinesgCSI(drug=drug, tissue.type="all", biomarkers=biomarkers, p.values=p.values, suffix="labeled")#, biomarkers.toPlot)
      ##plot effect sizes in gCSI dataset
      # biomarkers.order <- NULL
      # if(all(!is.na(rr))){
      #   biomarkers.order <- rr$hv$rowInd
      #   tt[match(names(rr$label.col), tt$biomarker.id), "gCSI.specificity"] <- "isoform.specific"
      # }
      # fnPlotEffectSize(drug=drug, biomarkers=biomarkers, effect.size="gCSI.estimate", biomarkers.order=biomarkers.order)
      # fnPlotLogPvalue(drug=drug, biomarkers=biomarkers, pvalue.col.name="gCSI.pvalue", effect.size.col="gCSI.estimate", biomarkers.order=biomarkers.order)
      
      # ## Plot Prevalidated biomarkers heatmap and sensitivity plot in training set
      # fnPlotAUCoverCellLinesCCLE.GDSC(drug, tissue.type="breast", biomarkers, biomarkers.order=biomarkers.order)#, biomarkers.toPlot)
      # fnPlotAUCoverCellLinesCCLE.GDSC.union.cells(drug, tissue.type="breast", biomarkers, biomarkers.order=biomarkers.order)#, biomarkers.toPlot)
      # ##plot effect sizes in gCSI dataset
      # fnPlotEffectSize(drug=drug, biomarkers=biomarkers, effect.size="estimate", biomarkers.order=biomarkers.order)
      # fnPlotLogPvalue(drug=drug, biomarkers=biomarkers, pvalue.col.name="pvalue", effect.size.col="estimate", biomarkers.order=biomarkers.order)
      
      # if(!is.null(biomarkers.order)) {
      #   delta.ranks <-  do.call(rbind, biomarkers)[rev(rr$hv$rowInd), c("isoform.no","delta.rank")]
      #   temp <- do.call(rbind, biomarkers)[rev(rr$hv$rowInd), c("short.label","gtex")]
      #   #rownames(delta.ranks) <- paste(temp[,"short.label"],ifelse(temp[,"gtex"] == "tumor.specific","*",""),sep="")
      #   #delta.ranks[,"delta.rank"] <- gsub("-", "{\\\color{red}-}", delta.ranks[,"delta.rank"])
      #   mode(delta.ranks) <- "character"
      #   xtable::print.xtable(xtable::xtable(delta.ranks, digits=0), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, sprintf("%s_delta.rank.tex", drug)), append=FALSE)
      # }
    }
    #biomarkers.toPlot <-  biomarkers.toPlot[hv$rowInd]
    #for (i in 1:length(biomarkers)){
    #bb[[i]] <- biomarkers[[hv$rowInd[i]]]
    #}
    biomarkers <- biomarkers.temp
    bb[[drug]] <- tt
  }
}
biomarkers <- bb
save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_2_%s.ctrpv2.RData", validation.method)))
####################################################

validated.biomarkers <- lapply(biomarkers, function(x) return(subset(x, bootstrap.validation.stat=="validated")))

validated.isoforms <- lapply(validated.biomarkers, function(x) return(subset(x, type=="isoform")))
validated.genes <- lapply(validated.biomarkers, function(x) return(subset(x, type=="gene")))

save(validated.biomarkers, validated.genes, validated.isoforms, file=file.path(path.diagrams, "validated.biomarkers.gCSI.ctrpv2.RData"))

WriteXLS::WriteXLS("validated.isoforms", ExcelFileName=file.path(path.diagrams, "validated.isoforms.gCSI.ctrpv2.xlsx"), row.names=TRUE)
WriteXLS::WriteXLS("validated.biomarkers", ExcelFileName=file.path(path.diagrams, "validated.biomarkers.gCSI.ctrpv2.xlsx"), row.names=TRUE)

#############

biomarkers.fdr <- lapply(biomarkers, function(x) return(cbind(x, "all.gCSI.fdr"=p.adjust(x[,"all.gCSI.pvalue"], method="fdr"))))


validated.biomarkers.fdr <- lapply(biomarkers.fdr, function(x) return(subset(x, validation.stat=="validated")))

validated.isoforms.fdr <- lapply(validated.biomarkers.fdr, function(x) return(subset(x, type=="isoform")))
validated.genes.fdr <- lapply(validated.biomarkers.fdr, function(x) return(subset(x, type=="gene")))

save(validated.biomarkers.fdr, validated.genes.fdr, validated.isoforms.fdr, file=file.path(path.diagrams, "validated.biomarkers.fdr.gCSI.RData"))

WriteXLS::WriteXLS("validated.isoforms.fdr", ExcelFileName=file.path(path.diagrams, "validated.isoforms.fdr.gCSI.xlsx"), row.names=TRUE)
WriteXLS::WriteXLS("validated.biomarkers.fdr", ExcelFileName=file.path(path.diagrams, "validated.biomarkers.fdr.gCSI.xlsx"), row.names=TRUE)

# validated.isoforms.fdr <- 

#############

fnValidatedPerCat <- function(biomarkers){

  valid.gene.n <- nrow(subset(biomarkers, specificity=="gene.specific"))
  valid.gene.percent <- sum(subset(biomarkers, specificity=="gene.specific")$validation.stat=="validated")/valid.gene.n

  valid.isoform.n <- nrow(subset(biomarkers, specificity=="isoform.specific"))
  valid.isoform.percent <- sum(subset(biomarkers, specificity=="isoform.specific")$validation.stat=="validated")/valid.isoform.n

  valid.common.n <- nrow(subset(biomarkers, specificity=="common"))
  valid.common.percent <- sum(subset(biomarkers, specificity=="common")$validation.stat=="validated")/valid.common.n

  return(c(gene.specific.percent = valid.gene.percent, gene.specific.n=valid.gene.n, 
           isoform.specific.percent = valid.isoform.percent, isoform.specific.n=valid.isoform.n,
           common.percent = valid.common.percent, common.n=valid.common.n))
}


fnValidatedPerType <- function(biomarkers){

  valid.gene.n <- nrow(subset(biomarkers, type=="gene"))
  valid.gene.percent <- sum(subset(biomarkers, type=="gene")$validation.stat=="validated")/valid.gene.n

  valid.isoform.n <- nrow(subset(biomarkers, type=="isoform"))
  valid.isoform.percent <- sum(subset(biomarkers, type=="isoform")$validation.stat=="validated")/valid.isoform.n

  return(c(gene.percent = valid.gene.percent, gene.n=valid.gene.n, 
           isoform.percent = valid.isoform.percent, isoform.n=valid.isoform.n))
}

data.frame(lapply(biomarkers, fnValidatedPerCat))

data.frame(lapply(biomarkers, fnValidatedPerType))


#############


if (FALSE) {
## histogram of expression of all the isoforms of a given gene to investigate why a given gene isoform is better correlated to drug response than the other mates in pre validation set
gene.name <- "S100A11"
drug <- "AZD6244"
xx <- rownames(gCSI.drug.sensitivity)[which(!is.na(gCSI.drug.sensitivity[ , drug]))]
trans <- rownames(annot.ensembl.all.transcripts)[which(annot.ensembl.all.transcripts$gene_name == gene.name)]
pdf(file=file.path(path.diagrams, sprintf("%s_%s.pdf", drug, gene.name)), height=4, width=4 * length(trans))
par(mfrow=c(1, length(trans)))
for( i in 1: length(trans)) {
  hist(gCSI.isoforms.fpkm[xx, trans[i]], main=trans[i], xlab="expression")
}
dev.off()
####################################################

## Plot expression and linear models tumor specific biomarkers for a certain drug
load(file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.RData", validation.method)), verbose=TRUE)
drug <- "AZD6244"
tt <- biomarkers[[drug]]
## select those which are expressed in CCLE with average expression greater than zero
tt <- subset(tt,  tt$ccle > 0)
## select vaidated biomarkers
vtt <- subset(tt, tt$validation.stat == "validated")
## tumor specific biomarkers
vtt <- subset(tt, tt$Gtex.train.sit < 0.05)
biomarkers <- fnFetchBiomarkers(top.significant.biomarkers=vtt, drug=drug, indices=1:nrow(vtt))

pdf(file=file.path(path.diagrams, sprintf("%s_tumour_specific.pdf", drug)), height=4 * (length(biomarkers) + 1), width=16)
par(mfrow=c((length(biomarkers) + 1), 4))
for(i in 1:length(biomarkers)) {
  biomarker <- biomarkers[[i]]
  if(biomarker$type == "isoform") {f <- biomarker$id %in% colnames(gCSI.isoforms.fpkm)}
  if(biomarker$type == "gene") {f <- biomarker$id %in% colnames(gCSI.genes.fpkm)}
  if (f) {
    fnPlotEXP(first=paste0("ccle.",biomarker$type), second=paste0("gCSI.",biomarker$type), biomarker)    
    fnPlotCCLEGDSC(biomarker, tissue.type=tissue, drug=drug)
    gCSI <- fnValidateWithgCSI(biomarker, method="all", color="#66CC00", drug=drug)
    gCSI <- fnValidateWithgCSI(biomarker, method="common", color="#336600", drug=drug)
  }
}
dev.off()
####################################################

## plot a gradiant legend for GDSC colors 
pdf(file=file.path(path.diagrams,"legend.pdf"), height=4, width=4)
barplot(NA, ylim=c(0,2) , yaxt="n", border=NA,axisnames=FALSE)
legend.gradient(cbind(c(0.2,.5,.5,0.2), c(0.1,0.1,1.3,1.3)), cols=colorRampPalette(c("blue","white","red"))(50), title="GDSC sensitivity",limits =c("resistent","sensitive"),cex=.7)
#legend.gradient("topleft", cols=exp.col, title="GDSC sensitivity", limits=c("R2<0.55","R2>0.7"),cex=.9)
dev.off()
####################################################


## plot a gradiant legend for expression heatmaps
pdf(file=file.path(path.diagrams,"legend_exp.pdf"), height=4, width=4)
barplot(NA, ylim=c(0,2) , yaxt="n", border=NA,axisnames=FALSE)
legend.gradient(cbind(c(0.2,.5,.5,0.2),c(0.1,0.1,1.3,1.3)), cols=colorRampPalette(c(blue,"white",red))(50), title="Expression",limits =c("low","high"),cex=.7)
#legend.gradient("topleft", cols=exp.col, title="GDSC sensitivity", limits=c("R2<0.55","R2>0.7"),cex=.9)
dev.off()
####################################################

#Distribution of breast samples gene expressions in GTex, gCSI and CCLE
pdf(file=file.path(path.diagrams, "Expression_breast.pdf"), height=14, width=14)
par(mfrow=c(2, 1))
fnPlotDensity(gtex.breast.genes, "GTEX")
ccle.cells <- intersect(rownames(ccle.tissuetype)[which(ccle.tissuetype == "breast")], rownames(ccle.genes.fpkm))
fnPlotDensity(ccle.genes.fpkm[ccle.cells, ], "CCLE")
#fnPlotDensity(gCSI.breast.genes, "gCSI")
dev.off()
####################################################

###box plots of isoforms expressions in gCSI dataset for selected biomarkers
cases <- c("MBP","SLC16A3","TOP1MT","CDCA4", "C10orf54","C10orf47","SLC12A4","TFAP2C","C19orf24")
locs <- (c(46, 10, 6, 1, 1, 1, 14, 2, 5))
library("genefu")
mycol <- RColorBrewer::brewer.pal(n=8, name="Set2")[c(1,4)]
i=1
for(gene in cases)
{
  loc <- locs[i]
  trans <- annot.ensembl.all.transcripts[which(annot.ensembl.all.transcripts$gene_name == gene),"transcript_id"]
  pdf(file=file.path(path.diagrams, sprintf("%s_isoforms_exp.pdf", gene)), height=7, width=10)
  par(mar=c(9,5,5,2))
  par(oma=c(2,2,2,2))
  
  invisible(genefu::boxplotplus2(t(gCSI.isoforms.fpkm[,trans]),
                                 pt.col=c(rep(mycol[1],times=loc - 1), mycol[2], rep(mycol[1],length(trans) - loc)),
                                 #names=c(paste0("none (",length(res.nosignifseeds),")"), paste0("all (",length(res.allsignifseeds),")")),
                                 ylab="expression", 
                                 #xlab="signif seeds (# of targets)",
                                 .las=2,
                                 pt.cex=.5))
  dev.off()
  i <- i + 1
}
####################################################

###box plots of isoforms expressions in CCLE dataset for selected biomarkers
i=1
for(gene in cases)
{
  gene <- cases[i]
  loc <- locs[i]
  trans <- annot.ensembl.all.transcripts[which(annot.ensembl.all.transcripts$gene_name == gene),"transcript_id"]
  pdf(file=file.path(path.diagrams, sprintf("%s_isoforms_exp_ccle.pdf", gene)), height=7, width=10)
  par(mar=c(9,5,5,2))
  par(oma=c(2,2,2,2))
  ccle.cells <- intersect(rownames(ccle.tissuetype)[which(ccle.tissuetype == tissue)], rownames(ccle.isoforms.fpkm))
  
  invisible(genefu::boxplotplus2(t(ccle.isoforms.fpkm[ccle.cells, trans]),
                                 pt.col=c(rep(mycol[1],times=loc - 1), mycol[2], rep(mycol[1],length(trans) - loc)),
                                 #names=c(paste0("none (",length(res.nosignifseeds),")"), paste0("all (",length(res.allsignifseeds),")")),
                                 ylab="expression", 
                                 #xlab="signif seeds (# of targets)",
                                 .las=2,
                                 pt.cex=.5))
  dev.off()
  i <- i + 1
}
####################################################



}
