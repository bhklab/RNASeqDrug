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
stat <- "r.squared"

load(file.path(path.data, "annotation.RData"))
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
load(file.path(path.data, "PSets/GRAY_isoforms.RData"))
gray.drug.sensitivity <- t(PharmacoGx::summarizeSensitivityProfiles(pSet=GRAY, sensitivity.measure=sensitivity.type))

load(file.path(path.data, "PSets/CCLE_isoforms.RData"))
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

ccle.tissuetype <- as.data.frame(CCLE@cell[, "tissueid"], row.names=CCLE@cell[, "cellid"], stringsAsFactors=FALSE)
colnames(ccle.tissuetype) <- "tissue.type"

gray.genes.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="rnaseq", fill.missing=FALSE)))
gray.isoforms.fpkm <- t(Biobase::exprs(PharmacoGx::summarizeMolecularProfiles(pSet=GRAY, mDataType="isoforms", fill.missing=FALSE)))

#associations <- associations.all.drugs("M3B")
#top.significant.biomarkers <- fnTop.significant.biomarkers(associations, cut_off=.5, BioNo=100, rank.type="pvalue")
#Validation.result <- fnValidation(all.biomarkers, validation.cut.off=5)

load(file.path(path.diagrams, "all.biomarkers.RData"))
max.bio.no <- max(sapply(all.biomarkers, function(x){nrow(x)}))
Validation.result <- fnValidation(top.significant.biomarkers=all.biomarkers, validation.cut.off=max.bio.no)

WriteXLS::WriteXLS("Validation.result", ExcelFileName=file.path(path.diagrams, "top.biomarkers.gray.xlsx"), row.names=TRUE)

for(validation.method in c("R2", "pvalue")){
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
    
    if(validation.method == "R2") {
      res.validated.breast <-  c(res.validated.breast, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                                                  sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate))& 
                                                                                  as.numeric(temp$R2) <= as.numeric(temp$breast),
                                                                                  select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
    }else{
      res.validated.breast <-  c(res.validated.breast, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                                                sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate))& 
                                                                                as.numeric(temp[ , paste0(tissue, "_pvalue")]) <= as.numeric(temp$breast),
                                                                                select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
    }
    validated <- unlist(str_split(res.validated.breast[i], pattern="_"))
    if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
    breast <- nrow(subset(temp, as.numeric(temp$R2) <= as.numeric(temp$breast)))
    
    res.validated.percent.breast <- c(res.validated.percent.breast, round(lv/breast, digits=2))
    all.breast <- c(all.breast, breast)
    
    res.validated.breast.boot <-  c(res.validated.breast.boot, paste(t(apply(subset(temp, temp$all.gray.pvalue < .05 & 
                                                                                      sign(as.numeric(temp$estimate)) == sign(as.numeric(temp$all.gray.estimate))& 
                                                                                      as.numeric(temp$R2) <= as.numeric(temp$breast_boot),
                                                                                    select=c("symbol", "type")), 1, function(x){sprintf("%s(%s)", x["symbol"], x["type"])}) ), collapse="_"))
    validated <- unlist(str_split(res.validated.breast.boot[i], pattern="_"))
    if(all(validated == "")){lv <- 0} else{lv <-length(validated)}
    breast.boot <- nrow(subset(temp, as.numeric(temp$R2) <= as.numeric(temp$breast_boot)))
    res.validated.percent.breast.boot <- c(res.validated.percent.breast.boot, round(lv/breast.boot, digits=2))
    
    
    all.breast.boot <- c(all.breast.boot, breast.boot)
  }
  names(res.validated) <- names(Validation.result)
  write.csv(cbind("validated"=res.validated, 
                  "ratio"=res.validated.percent, 
                  "biomarkers no"=all,
                  "validated.breast"=res.validated.breast, 
                  "ratio.breast"=res.validated.percent.breast, 
                  "biomarkers no.breast"=all.breast,
                  "validated.boot.breast"=res.validated.breast.boot, 
                  "ratio.boot.breast"=res.validated.percent.breast.boot, 
                  "biomarkers no.boot.breast"=all.breast.boot), file=file.path(path.diagrams, sprintf("validated_biomarkers_%s.csv", validation.method)))
  load("data/GTex_BR.RData", verbose=T)
  #fnGtex(validation.method="R2")
  biomarkers <- fnGtex(validation.method=validation.method)
  save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.RData", validation.method)))
}

## Heatmaps, effect sizes and sensitivity plots for tarining and pre validation
biomarkers.category <- "isoforms"
validation.method <- "R2"
tissue.type <- "all"

mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
red <- mycol[1]  
blue <- mycol[2]

model.method <- "glm"
glm.family <- "gaussian"
stat <- "r.squared"

load(file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.RData", validation.method)), verbose=TRUE)
biomarkers.temp <- biomarkers
bb <- list()
for(drug in names(biomarkers)) {
  drug.name <- drug
  tt <- biomarkers[[drug]]
  if(!"gtex.specificity" %in% colnames(tt)){
    tt <- cbind(tt, "gray.specificity"=NA)
  }
  #tt <- subset(tt,  sign(tt$estimate) >=0)
  #biomarkers.toPlot <- c("ITGA6", "S100A11", "SH3TC2", "MBP")
  ## select those expressed in CCLE 
  tt <- subset(tt,  tt$ccle > 0)
  vtt <- subset(tt, tt$validation.stat == "validated")
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
      for(i in 1:length(biomarkers)) {
        M0 <- fnCreateNullModel(drug=drug, assay="gray")
        M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gray.genes.fpkm[ ,biomarkers[[i]]$gene.id])
        M3B <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gray.isoforms.fpkm[ ,biomarkers[[i]]$isoform.id])
        results <- fnRunbootstrap(models=list("M0"=M0, "M2"=M2, "M3B"=M3B))
        p.values <- c(p.values, wilcox.test(results[["M3B"]], results[["M2"]], paired=TRUE, alternative="greater")$p.value)
      }
      names(p.values) <- sapply(biomarkers, function(x){x[["id"]]})
      ## Plot Prevalidated biomarkers heatmap and sensitivity plot, cell lines are ordered by sensitivity measurement and similar biomarkers are clustered together
      rr <- fnPlotAUCoverCellLinesGray(drug=drug, tissue.type="all", biomarkers=biomarkers, p.values=p.values, suffix="labeled")#, biomarkers.toPlot)
      ##plot effect sizes in gray dataset
      biomarkers.order <- NULL
      if(all(!is.na(rr))){
        biomarkers.order <- rr$hv$rowInd
        tt[match(names(rr$label.col), tt$biomarker.id), "gray.specificity"] <- "isoform.specific"
      }
      fnPlotEffectSize(drug=drug, biomarkers=biomarkers, effect.size="gray.estimate", biomarkers.order=biomarkers.order)
      fnPlotLogPvalue(drug=drug, biomarkers=biomarkers, pvalue.col.name="gray.pvalue", effect.size.col="gray.estimate", biomarkers.order=biomarkers.order)
      
      ## Plot Prevalidated biomarkers heatmap and sensitivity plot in training set
      fnPlotAUCoverCellLinesCCLE.GDSC(drug, tissue.type="breast", biomarkers, biomarkers.order=biomarkers.order)#, biomarkers.toPlot)
      fnPlotAUCoverCellLinesCCLE.GDSC.union.cells(drug, tissue.type="breast", biomarkers, biomarkers.order=biomarkers.order)#, biomarkers.toPlot)
      ##plot effect sizes in gray dataset
      fnPlotEffectSize(drug=drug, biomarkers=biomarkers, effect.size="estimate", biomarkers.order=biomarkers.order)
      fnPlotLogPvalue(drug=drug, biomarkers=biomarkers, pvalue.col.name="pvalue", effect.size.col="estimate", biomarkers.order=biomarkers.order)
      
      if(!is.null(biomarkers.order)) {
        delta.ranks <-  do.call(rbind, biomarkers)[rev(rr$hv$rowInd), c("isoform.no","delta.rank")]
        temp <- do.call(rbind, biomarkers)[rev(rr$hv$rowInd), c("short.label","gtex")]
        #rownames(delta.ranks) <- paste(temp[,"short.label"],ifelse(temp[,"gtex"] == "tumor.specific","*",""),sep="")
        #delta.ranks[,"delta.rank"] <- gsub("-", "{\\\color{red}-}", delta.ranks[,"delta.rank"])
        xtable::print.xtable(xtable::xtable(delta.ranks, digits=0), include.rownames=TRUE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, sprintf("%s_delta.rank.tex", drug)), append=FALSE)
      }
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
save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_2_%s.RData", validation.method)))
####################################################

## histogram of expression of all the isoforms of a given gene to investigate why a given gene isoform is better correlated to drug response than the other mates in pre validation set
gene.name <- "S100A11"
drug <- "AZD6244"
xx <- rownames(gray.drug.sensitivity)[which(!is.na(gray.drug.sensitivity[ , drug]))]
trans <- rownames(annot.ensembl.all.transcripts)[which(annot.ensembl.all.transcripts$gene_name == gene.name)]
pdf(file=file.path(path.diagrams, sprintf("%s_%s.pdf", drug, gene.name)), height=4, width=4 * length(trans))
par(mfrow=c(1, length(trans)))
for( i in 1: length(trans)) {
  hist(gray.isoforms.fpkm[xx, trans[i]], main=trans[i], xlab="expression")
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
  if(biomarker$type == "isoform") {f <- biomarker$id %in% colnames(gray.isoforms.fpkm)}
  if(biomarker$type == "gene") {f <- biomarker$id %in% colnames(gray.genes.fpkm)}
  if (f) {
    fnPlotEXP(first=paste0("ccle.",biomarker$type), second=paste0("gray.",biomarker$type), biomarker)    
    fnPlotCCLEGDSC(biomarker, tissue.type=tissue, drug=drug)
    gray <- fnValidateWithGray(biomarker, method="all", color="#66CC00", drug=drug)
    gray <- fnValidateWithGray(biomarker, method="common", color="#336600", drug=drug)
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

#Distribution of breast samples gene expressions in GTex, GRAY and CCLE
pdf(file=file.path(path.diagrams, "Expression_breast.pdf"), height=14, width=14)
par(mfrow=c(2, 1))
fnPlotDensity(gtex.breast.genes, "GTEX")
ccle.cells <- intersect(rownames(ccle.tissuetype)[which(ccle.tissuetype == "breast")], rownames(ccle.genes.fpkm))
fnPlotDensity(ccle.genes.fpkm[ccle.cells, ], "CCLE")
#fnPlotDensity(gray.breast.genes, "GRAY")
dev.off()
####################################################

###box plots of isoforms expressions in GRAY dataset for selected biomarkers
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
  
  invisible(genefu::boxplotplus2(t(gray.isoforms.fpkm[,trans]),
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




