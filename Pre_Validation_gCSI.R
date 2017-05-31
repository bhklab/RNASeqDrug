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

source(file.path(path.code, "foo_gCSI.R"))
source(file.path(path.code,"foo_PreValidation_gCSI.R"))
set.seed(12345)

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

load(file.path(path.data, "PSets/CCLE_hs.RData"))
load(file.path(path.data, "PSets/GDSC.RData"))
if(training.type == "CCLE_GDSC") {

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


load(file.path(path.diagrams, "all.biomarkers.RData"))
max.bio.no <- max(sapply(all.biomarkers, function(x){nrow(x)}))
annot.ensembl.all.genes=fData(CCLE@molecularProfiles$rnaseq)
annot.ensembl.all.genes$gene_id <- annot.ensembl.all.genes$EnsemblGeneId
annot.ensembl.all.genes$gene_name <- annot.ensembl.all.genes$Symbol
annot.ensembl.all.genes$gene_biotype <- annot.ensembl.all.genes$GeneBioType

common.cells.ccle <- intersect(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm))
common.cells.gCSI <- intersect(rownames(gCSI.drug.sensitivity), rownames(gCSI.genes.fpkm))


ccle.cells <- intersectList(rownames(ccle.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))
gdsc.cells <- intersectList(rownames(gdsc.drug.sensitivity), rownames(ccle.genes.fpkm), rownames(ccle.isoforms.fpkm))


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

tissue <- NULL


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

biomarkers <- all.biomarkers

# load(file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_%s.RData", validation.method)), verbose=TRUE)
biomarkers.temp <- biomarkers
bb <- list()
for(drug in drugs) {
  drug.name <- drug
  print(paste0("Now processing drug: ", drug))
  tt <- biomarkers[[drug]]
  stat.colname <- paste("gCSI", stat, sep=".")
  delta.stat.colname <- paste("gCSI", "delta", stat, sep=".")

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

        p.values <- c(p.values, 1)
        if(biomarkers[[i]]$type == "gene"){
          M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gCSI.genes.fpkm[ ,biomarkers[[i]]$gene.id])
          if(!is.null(M2$model)){
            results <- fnRunbootstrap(models=list("M0"=M0, "M2"=M2))
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,stat.colname] <- median(results[["M2"]])
            tt[tt$biomarker.id==biomarkers[[i]]$symbol,"gCSI.bootstrap.pvalue"] <- wilcox.test(results[["M2"]], results[["M0"]], paired=TRUE, alternative="greater")$p.value
            if(tt[tt$biomarker.id==biomarkers[[i]]$symbol,"gCSI.bootstrap.pvalue"]<=0.05 && tt[tt$biomarker.id==biomarkers[[i]]$symbol,stat.colname] >= 0.55){
              tt[tt$biomarker.id==biomarkers[[i]]$symbol,"bootstrap.validation.stat"] <- "validated"
            }
          }
        }
        if(biomarkers[[i]]$type == "isoform"){
          M3B <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gCSI.isoforms.fpkm[ ,biomarkers[[i]]$isoform.id])
          if(!is.null(M3B$model)){
            results <- fnRunbootstrap(models=list("M0"=M0, "M3B"=M3B))
            tt[tt$biomarker.id==biomarkers[[i]]$id,stat.colname] <- median(results[["M3B"]])
            tt[tt$biomarker.id==biomarkers[[i]]$id,"gCSI.bootstrap.pvalue"] <- wilcox.test(results[["M3B"]], results[["M0"]], paired=TRUE, alternative="greater")$p.value        
            if(tt[tt$biomarker.id==biomarkers[[i]]$id,"gCSI.bootstrap.pvalue"]<=0.05 && tt[tt$biomarker.id==biomarkers[[i]]$id,stat.colname] >= 0.55){
              tt[tt$biomarker.id==biomarkers[[i]]$id,"bootstrap.validation.stat"] <- "validated"
            }
          }
        }
      }
    }
    biomarkers <- biomarkers.temp
    bb[[drug]] <- tt
  }
}
biomarkers <- bb
save(biomarkers, file=file.path(path.diagrams, sprintf("Biomarkers_with_validation_status_2_%s.RData", validation.method)))
####################################################

validated.biomarkers <- lapply(biomarkers, function(x) return(subset(x, bootstrap.validation.stat=="validated")))

validated.isoforms <- lapply(validated.biomarkers, function(x) return(subset(x, type=="isoform")))
validated.genes <- lapply(validated.biomarkers, function(x) return(subset(x, type=="gene")))

save(validated.biomarkers, validated.genes, validated.isoforms, file=file.path(path.diagrams, "validated.biomarkers.gCSI.RData"))

WriteXLS::WriteXLS("validated.isoforms", ExcelFileName=file.path(path.diagrams, "validated.isoforms.gCSI.xlsx"), row.names=TRUE)
WriteXLS::WriteXLS("validated.biomarkers", ExcelFileName=file.path(path.diagrams, "validated.biomarkers.gCSI.xlsx"), row.names=TRUE)

# #############

# biomarkers.fdr <- lapply(biomarkers, function(x) return(cbind(x, "all.gCSI.fdr"=p.adjust(x[,"all.gCSI.pvalue"], method="fdr"))))


# validated.biomarkers.fdr <- lapply(biomarkers.fdr, function(x) return(subset(x, validation.stat=="validated")))

# validated.isoforms.fdr <- lapply(validated.biomarkers.fdr, function(x) return(subset(x, type=="isoform")))
# validated.genes.fdr <- lapply(validated.biomarkers.fdr, function(x) return(subset(x, type=="gene")))

# save(validated.biomarkers.fdr, validated.genes.fdr, validated.isoforms.fdr, file=file.path(path.diagrams, "validated.biomarkers.fdr.gCSI.RData"))

# WriteXLS::WriteXLS("validated.isoforms.fdr", ExcelFileName=file.path(path.diagrams, "validated.isoforms.fdr.gCSI.xlsx"), row.names=TRUE)
# WriteXLS::WriteXLS("validated.biomarkers.fdr", ExcelFileName=file.path(path.diagrams, "validated.biomarkers.fdr.gCSI.xlsx"), row.names=TRUE)

# validated.isoforms.fdr <- 

#############

fnValidatedPerCat <- function(biomarkers){

  gene.n <- nrow(subset(biomarkers, specificity=="gene.specific"))
  valid.gene.percent <- sum(subset(biomarkers, specificity=="gene.specific")$bootstrap.validation.stat=="validated")#/gene.n

  valid.isoform.n <- nrow(subset(biomarkers, specificity=="isoform.specific"))
  valid.isoform.percent <- sum(subset(biomarkers, specificity=="isoform.specific")$bootstrap.validation.stat=="validated")#/valid.isoform.n

  valid.common.n <- nrow(subset(biomarkers, specificity=="common"))
  valid.common.percent <- sum(subset(biomarkers, specificity=="common")$bootstrap.validation.stat=="validated")#/valid.common.n

  return(c( "Gene Specific Not Validated"=gene.n, "Gene Specific Validated" = valid.gene.percent, 
            "Isoform Specific Not Validated"=valid.isoform.n , "Isoform Specific Validated" = valid.isoform.percent,
            "Common Not Validated"=valid.common.n , "Common Validated" = valid.common.percent))
}


fnValidatedPerType <- function(biomarkers){

  valid.gene.n <- nrow(subset(biomarkers, type=="gene"))
  valid.gene.percent <- sum(subset(biomarkers, type=="gene")$bootstrap.validation.stat=="validated")/valid.gene.n

  valid.isoform.n <- nrow(subset(biomarkers, type=="isoform"))
  valid.isoform.percent <- sum(subset(biomarkers, type=="isoform")$bootstrap.validation.stat=="validated")/valid.isoform.n

  return(c(gene.percent = valid.gene.percent, gene.n=valid.gene.n, 
           isoform.percent = valid.isoform.percent, isoform.n=valid.isoform.n))
}

valid.per.cat <- data.frame(lapply(biomarkers, fnValidatedPerCat), check.names=FALSE)
valid.per.cat

data.frame(lapply(biomarkers, fnValidatedPerType))

#############
require(grid)
require(gridExtra)
require(ggthemes)
require(reshape2)
require(RColorBrewer)


lay <- rbind(c(rep(c(1,1,1,1),3),NA, rep(2,3),2))

# valid.per.cat <- valid.per.cat[c(1,3,5),] * valid.per.cat[c(1,3,5)+1,]
# valid.percent <- valid.per.cat[c(1,3,5),]

# valid.percent2 <- apply(valid.percent, c(1,2), function(x) return(ifelse(is.finite(x), x, 0)))

pdf("barplot_gCSI_valid.pdf", height=10, width=12.5)
{

valid.per.cat <- valid.per.cat[,c("Crizotinib", "lapatinib", "Erlotinib", "paclitaxel", "PD-0325901")]
grobs <- list()
toPlot <- melt(as.matrix(valid.per.cat), varnames=c("Status", "Drug"), value.name="Number of Biomarkers")

toPlot$Category <- gsub(toPlot$Status, pattern="( Not)?.Validated", rep="")

toPlot$Category <- factor(toPlot$Category, levels=c("Isoform Specific", "Common", "Gene Specific"))

toPlot$Status <- gsub(toPlot$Status, pattern=".+Not Validated", rep="Not Validated")

toPlot$Status <- factor(toPlot$Status, levels=c("Isoform Specific Validated", "Common Validated", "Gene Specific Validated",  "Not Validated"))
toPlot1 <- subset(toPlot, Drug!="PD-0325901")
# toPlot <- cbind(toPlot, "Status"=c("Validated", "Not Validated"))
#subset(toPlot, Drug==drug)
myPal <- RColorBrewer::brewer.pal(n=9, name="Set1")

mycols <- c(myPal[1], myPal[4], myPal[2], "#b3b3b3")

g <- ggplot(toPlot1, aes(x=NA, y=`Number of Biomarkers`, group=Category, fill=Status)) + 
            geom_col(width=0.8, position = position_dodge(width=0.9)) +
            # expand_limits(y=(max(toPlot[,"Validation Rate"]) + 0.05)) +
            #coord_cartesian(ylim=c(-0.001,1)) + 
            theme_minimal()  + 
            #scale_y_reverse(expand=c(0,0), labels=function(x){return(gsub(sprintf("%1.2f",as.numeric(x)), pattern="0.00", replace="", fixed=TRUE))}) + 
            theme(legend.position = c(0.3,0.8), legend.text=element_text(size=16)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            facet_wrap(~Drug, nrow=1, strip.position ="bottom") +
            # theme(axis.title.x=element_blank(),axis.text.y = element_blank(), legend.background = element_rect(fill="white", size=0), plot.title = element_text(size=16,face="bold", hjust = 0.5), axis.title.y = element_text(size=12,face="bold")) + 
            theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face="bold"), legend.background = element_rect(fill="white", size=0), plot.title = element_text(size=24,face="bold", hjust = 0.5), axis.title = element_text(size=30)) + 
            xlab("") +  #coord_flip() +
            theme(plot.margin=unit(c(0.005,0,0.005,0.005),"npc"), strip.text.x = element_text(size=30)) + 
            scale_fill_manual(values=mycols) #+ 
            #guides(fill=guide_legend(reverse=TRUE))


grobs[["barplot1"]] <-  g

toPlot2 <- subset(toPlot, Drug=="PD-0325901")

g <- ggplot(toPlot2, aes(x=NA, y=`Number of Biomarkers`, group=Category, fill=Status)) + 
            geom_col(width=0.8, position = position_dodge(width=0.9)) +
            # expand_limits(y=(max(toPlot[,"Validation Rate"]) + 0.05)) +
            #coord_cartesian(ylim=c(-0.001,1)) + 
            theme_minimal()  + 
            #scale_y_reverse(expand=c(0,0), labels=function(x){return(gsub(sprintf("%1.2f",as.numeric(x)), pattern="0.00", replace="", fixed=TRUE))}) + 
            theme(legend.position = "none") + scale_y_continuous(position="right") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            facet_wrap(~Drug, nrow=1, strip.position ="bottom") +
            # theme(axis.title.x=element_blank(),axis.text.y = element_blank(), legend.background = element_rect(fill="white", size=0), plot.title = element_text(size=20,face="bold", hjust = 0.5), axis.title.y = element_text(size=12,face="bold")) + 
            theme(axis.text.x = element_blank(), axis.text.y = element_text(size=20, face="bold"), legend.background = element_rect(fill="white", size=0), plot.title = element_text(size=16,face="bold", hjust = 0.5), axis.title = element_text(size=30)) + 
            xlab("") +  ylab("") +  #coord_flip() +
            theme(plot.margin=unit(c(0.005,0,0.005,0.005),"npc"), strip.text.x = element_text(size=30)) + 
            #ggtitle("Pre-Validation Rate of Biomarkers in gCSI") + 
            scale_fill_manual(values=mycols)

grobs[["barplot2"]] <-  g

title1=textGrob("Pre-Validation Rate of Biomarkers in gCSI", gp=gpar(fontface="bold", cex=3))
title2=textGrob("Drug", gp=gpar(cex=2.5))


grid.arrange(grobs = grobs, layout_matrix=lay, bottom=title2)


}
dev.off()

pdf("AUC_dist_gCSI_Valid.pdf", height=7, width=11)
{
drugs <- names(biomarkers)
toPlot <- cbind(Dataset="gCSI", melt(gCSI.drug.sensitivity[,drugs],varnames=c("Cell", "Drug"),value.name = "AAC"))

toPlot2 <- cbind(Dataset="CCLE", melt(ccle.drug.sensitivity[,drugs],varnames=c("Cell", "Drug"),value.name = "AAC"))

toPlot3 <- cbind(Dataset="GDSC", melt(gdsc.drug.sensitivity[,drugs],varnames=c("Cell", "Drug"),value.name = "AAC"))

toPlot <- rbind(toPlot, toPlot2, toPlot3)

toPlot$Dataset <- factor(toPlot$Dataset, levels=c("CCLE", "GDSC", "gCSI"))
toPlot$Drug <- factor(toPlot$Drug, levels=c("Crizotinib", "lapatinib", "Erlotinib", "paclitaxel", "PD-0325901"))
# toPlot1 <- subset(toPlot, Dataset != "gCSI")
toPlot$DatasetPos <- as.numeric(toPlot$Dataset)
toPlot$DatasetPos <- as.numeric(gsub(toPlot$DatasetPos, pat="3", rep="3.3"))
# g <- 
# grobs[["violin"]] <-  g


ggplot(toPlot, aes(x=DatasetPos, y=AAC, fill=Dataset)) + 
            geom_violin(show.legend=TRUE) + scale_fill_grey() + theme(legend.position = "none") +
            theme_minimal() + #theme(legend.position = "bottom") +
            #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            facet_wrap(~Drug, nrow=1) + scale_x_continuous(breaks=c(1,2,3.3), labels=c("CCLE", "GDSC", "gCSI")) +
            theme(plot.title = element_text(size=16,face="bold", hjust = 0.5)) +
            # theme(axis.title = element_blank(), axis.text.y = element_blank()) + 
            xlab("") + theme(axis.text.y = element_blank(), axis.text.x=element_text(size=16, angle=45)) +
            coord_flip() + geom_vline(xintercept=2.65, linetype = "1313") + 
            theme(plot.margin = unit(c(0.005,0.03,0.005,0.01-0.01),"npc") , strip.text.x = element_text(size=22)) + 
            scale_y_continuous(expand=c(0,0), position="left") 
            #+ ggtitle("Distribution of Sensitivity Measures for each Drug by Dataset")
}
dev.off()

# valid.percent <- valid.per.cat[c(1,3,5)+1,]


# {
# grobs <- list()
# toPlot <- melt(as.matrix(valid.percent), varnames=c("Category", "Drug"), value.name="Validation Rate")
# #subset(toPlot, Drug==drug)
# g <- ggplot(toPlot, aes(x=NA, y=`Validation Rate`, group=Category, fill=Category)) + geom_col(position = "dodge", alpha = 4/5) +
#             # expand_limits(y=(max(toPlot[,"Validation Rate"]) + 0.05)) +
#             coord_cartesian(ylim=c(-0.001,1)) + 
#             theme_minimal()  + 
#             scale_y_reverse(expand=c(0,0), labels=function(x){return(gsub(sprintf("%1.2f",as.numeric(x)), pattern="0.00", replace="", fixed=TRUE))}) + 
#             theme(legend.position = c(0.2,0.6)) +
#             #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#             facet_wrap(~Drug, ncol=1) +
#             theme(axis.title.x=element_blank(),axis.text.y = element_blank(), legend.background = element_rect(fill="white", size=0), plot.title = element_text(size=16,face="bold", hjust = 0.5), axis.title.y = element_text(size=12,face="bold")) + 
#             xlab("Drug") +  coord_flip() +
#             theme(plot.margin=unit(c(0.005,0,0.005,0.005),"npc"), strip.text.x = element_text(face="bold")) + 
#             ggtitle("Number of training biomarkers") + 
#             scale_fill_brewer(palette="Set1") + 
#             guides(fill=guide_legend(reverse=TRUE))

# grobs[["boxplot"]] <-  g

# drugs <- names(biomarkers)
# # toPlot <- cbind(Dataset="gCSI", melt(gCSI.drug.sensitivity[,drugs],varnames=c("Cell", "Drug"),value.name = "AAC"))

# toPlot2 <- cbind(Dataset="CCLE", melt(ccle.drug.sensitivity[,drugs],varnames=c("Cell", "Drug"),value.name = "AAC"))

# toPlot3 <- cbind(Dataset="GDSC", melt(gdsc.drug.sensitivity[,drugs],varnames=c("Cell", "Drug"),value.name = "AAC"))

# toPlot <- rbind(toPlot2, toPlot3)

# g <- 
# ggplot(toPlot, aes(x=Dataset, y=AAC, fill=Dataset)) + 
#             geom_violin() + scale_fill_grey() + theme(legend.position = c(0.2,0.8)) +
#             theme_minimal() + #theme(legend.position = "bottom") +
#             #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#             facet_wrap(~Drug, ncol=1) + 
#             theme(plot.title = element_text(size=16,face="bold", hjust = 0.5)) +
#             # theme(axis.title = element_blank(), axis.text.y = element_blank()) + 
#             # xlab(NULL) + coord_flip() +
#             theme(plot.margin = unit(c(0.005,0.03,0.005,0.01-0.01),"npc") , strip.text.x = element_text(face="bold")) + scale_y_continuous(expand=c(0,0), position="right") + ggtitle("AAC")
# grobs[["violin"]] <-  g
# grid.arrange(grobs = grobs, layout_matrix=lay)
# }


# grid.arrange(grobs = grobs, layout_matrix=lay)
#############

#############
## gCSI AAC vs CCLE



myScatterPlot <- function(Name, x, y, method=c("plain", "annotated", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), cex=1.5, annot=NA, ...){
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  
  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch,cex=cex, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             # pdf(file = Name, height=7, width=7)
             par(mar=c(4.5,5,4,4))
             plot(x=x, y=y, col=col, pch=pch, cex.lab=cex, cex.main=cex, cex.axis=cex/1.25, ...)
             abline(0, 1, lty=2, col="gray")
             legend("topright", legend=sprintf("r=%1.2g", cor(as.numeric(x),as.numeric(y), method="spearman")), cex=cex, bty="n")
             # dev.off()
           },
            "annotated"={
             # pdf(file = Name, height=7, width=7)
             par(mar=c(4.5,5,4,4))
             plot(x=x, y=y, col=col, pch=pch, cex.lab=cex, cex.main=cex, cex.axis=cex/1.25, ...)
             abline(0, 1, lty=2, col="gray")
             if(!is.na(annot)){
                legend("topright", legend=annot, cex=cex, bty="n")
             }
             # dev.off()
           },
           "transparent"={
             pdf(file = Name, height=7, width=7)
             myrgb <- grDevices::col2rgb(col, alpha=FALSE) / 255
             plot(x=x, y=y, pch=pch, col=rgb(red=myrgb[1], green=myrgb[2], blue=myrgb[3], alpha=transparency, maxColorValue=1), ...)
             dev.off()
           },
           "smooth"={
             pdf(file = Name, height=7, width=7)
             smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
             #abline(0,max(y)/max(x),col="red")
             abline(lm(y~x), col = "red")
             #lines(x=c(0,0), y=c(max(x),max(y)), lty="solid", col="red")
             dev.off()
           }
    )
  }
}

common.cells.gCSI.ccle <- intersect(common.cells.gCSI, common.cells.ccle)
common.cells.gCSI.gdsc <- intersectList(common.cells.gCSI, common.cells.ccle, cellNames(GDSC))


for(drug in sort(names(biomarkers))){

  toPlot1 <- data.frame("gCSI"=gCSI.drug.sensitivity[common.cells.gCSI.ccle,drug], 
                        "CCLE" = ccle.drug.sensitivity[common.cells.gCSI.ccle,drug], check.names=FALSE)
  toPlot2 <- data.frame("gCSI"=gCSI.drug.sensitivity[common.cells.gCSI.gdsc,drug], 
                        "GDSC" = gdsc.drug.sensitivity[common.cells.gCSI.gdsc,drug], check.names=FALSE)
  pdf(file = sprintf("gCSI vs train AAC for %s.pdf", drug), height=14, width=7)
  par(mfrow=c(2,1))
  # ScatterHist(toPlot1, "gCSI", "CCLE", title=sprintf("%s", drug), smoothmethod="lm", se=FALSE)
  myScatterPlot("", x=toPlot1$gCSI, y=toPlot1$CCLE, ylab="CCLE AAC", xlab="gCSI AAC", main=sprintf("%s", drug), cex=2, xlim=c(0,1), ylim=c(0,1))
  myScatterPlot("", x=toPlot2$GDSC, y=toPlot2$gCSI, ylab="GDSC AAC", xlab="gCSI AAC", cex=2, xlim=c(0,1), ylim=c(0,1))
  dev.off()
}



expression.cors.valid <- lapply(biomarkers, function(x){
  tt <- subset(x, type=="gene" & bootstrap.validation.stat=="validated")
  xx <- ccle.genes.fpkm[common.cells.gCSI.ccle,tt$gene.id]
  yy <- gCSI.genes.fpkm[common.cells.gCSI.ccle,tt$gene.id]
  tt <- subset(x, type=="isoform" & bootstrap.validation.stat=="validated")
  xx <- cbind(xx,ccle.isoforms.fpkm[common.cells.gCSI.ccle,tt$transcript.id])
  yy <- cbind(yy,gCSI.isoforms.fpkm[common.cells.gCSI.ccle,tt$transcript.id])
  if(NROW(tt)==0){
    return(NA)
  }
  return(diag(cor(xx, yy, method="pearson", use="pairwise.complete.obs")))
})


expression.cors.unvalid <- lapply(biomarkers, function(x){
  tt <- subset(x, type=="gene" & bootstrap.validation.stat=="unvalidated")
  xx <- ccle.genes.fpkm[common.cells.gCSI.ccle,tt$gene.id]
  yy <- gCSI.genes.fpkm[common.cells.gCSI.ccle,tt$gene.id]
  tt <- subset(x, type=="isoform" & bootstrap.validation.stat=="unvalidated")
  xx <- cbind(xx,ccle.isoforms.fpkm[common.cells.gCSI.ccle,tt$transcript.id])
  yy <- cbind(yy,gCSI.isoforms.fpkm[common.cells.gCSI.ccle,tt$transcript.id])
  if(NROW(tt)==0){
    return(NA)
  }
  return(diag(cor(xx, yy, method="pearson", use="pairwise.complete.obs")))
})

expression.cors.all <- lapply(biomarkers, function(x){
  tt <- subset(x, type=="gene")
  xx <- ccle.genes.fpkm[common.cells.gCSI.ccle,tt$gene.id]
  yy <- gCSI.genes.fpkm[common.cells.gCSI.ccle,tt$gene.id]
  tt <- subset(x, type=="isoform")
  xx <- cbind(xx,ccle.isoforms.fpkm[common.cells.gCSI.ccle,tt$transcript.id])
  yy <- cbind(yy,gCSI.isoforms.fpkm[common.cells.gCSI.ccle,tt$transcript.id])
  if(NROW(tt)==0){
    return(NA)
  }
  return(diag(cor(xx, yy, method="pearson", use="pairwise.complete.obs")))
})

toPlot1 <- cbind(melt(expression.cors.valid), status="Validated")
toPlot2 <- cbind(melt(expression.cors.unvalid), status="Unvalidated")
toPlot <- rbind(toPlot1,toPlot2)
colnames(toPlot) <- c("Cor", "Drug", "Status")

# myPal <- RColorBrewer::brewer.pal(n=8, name="Set2")
pdf("distribution_expression_biomarker.pdf", height=7, width=14)
ggplot(toPlot, aes(x=Cor)) + geom_density(fill="grey") + facet_wrap(~Drug, nrow=1) + theme_minimal() +
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
       #scale_fill_manual(values=myPal[c(1,2)]) + 
       theme(axis.text.x=element_text(size=16, angle=45), axis.text.y=element_text(size=16)) +
       xlab("") + ylab("")
dev.off()


# require(genefu)
# myPal <- RColorBrewer::brewer.pal(n=9, name="Set3")
# for (drug in names(gene.expression.cors)){
#   pdf(sprintf("%s_biomarker_cors.pdf", drug))
#   par(mfrow=c(1,1))
#   boxplot(list("Validated" = expression.cors.valid[[drug]], "Not Validated"=expression.cors.unvalid[[drug]]),
#           cex.axis = 0.85, border = "black",
#           col = c(myPal[1],myPal[4]), boxwex = 0.5, ylab="Pearson Correlation of Expression", ylim=c(-.2,1))
#   dev.off()

# }

require(VennDiagram)
myPal <- RColorBrewer::brewer.pal(n=9, name="Set1")
pdf(file.path(file.path(path.diagrams, "gCSI_ccle_gdsc_celllines.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(CCLE@cell$cellid), 
                                          area2 = length(GDSC@cell$cellid),
                                          area3 = length(cellNames(gCSI)),
                                          n12 = length(intersect(CCLE@cell$cellid, GDSC@cell$cellid)),
                                          n13 = length(intersect(CCLE@cell$cellid, cellNames(gCSI))),
                                          n23 = length(intersect(GDSC@cell$cellid, cellNames(gCSI))),
                                          n123 = length(PharmacoGx::intersectList(CCLE@cell$cellid, GDSC@cell$cellid, cellNames(gCSI))),
                                          category = c("CCLE", "GDSC", "gCSI"),
                                          col = mycol[c(1:2,4)],
                                          fill = mycol[c(1:2,4)],
                                          margin=0.10,
                                          cat.dist=.1,
                                          cat.default.pos="text",
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()

pdf(file.path(file.path(path.diagrams, "gCSI_ccle_gdsc_drugs.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(drugNames(CCLE)), 
                                          area2 = length(drugNames(GDSC)),
                                          area3 = length(drugNames(gCSI)),
                                          n12 = length(intersect(drugNames(CCLE), drugNames(GDSC))),
                                          n13 = length(intersect(drugNames(CCLE), drugNames(gCSI))),
                                          n23 = length(intersect(drugNames(GDSC), drugNames(gCSI))),
                                          n123 = length(PharmacoGx::intersectList(drugNames(CCLE), drugNames(GDSC), drugNames(gCSI))),
                                          category = c("CCLE", "GDSC", "gCSI"),
                                          col = mycol[c(1:2,4)],
                                          fill = mycol[c(1:2,4)],
                                          margin=0.10,
                                          cat.dist=.2,
                                          cat.default.pos="text",
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()

#lapply(validated.biomarkers, function(x)



fnValidatedPerCat <- function(biomarkers){

  gene.n <- nrow(subset(biomarkers, specificity=="gene.specific"))
  valid.gene.percent <- sum(subset(biomarkers, specificity=="gene.specific")$bootstrap.validation.stat=="validated")#/gene.n

  valid.isoform.n <- nrow(subset(biomarkers, specificity=="isoform.specific"))
  valid.isoform.percent <- sum(subset(biomarkers, specificity=="isoform.specific")$bootstrap.validation.stat=="validated")#/valid.isoform.n

  valid.common.n <- nrow(subset(biomarkers, specificity=="common"))
  valid.common.percent <- sum(subset(biomarkers, specificity=="common")$bootstrap.validation.stat=="validated")#/valid.common.n

  return(c( "Gene Specific Total Num"=gene.n, "Gene Specific percent" = valid.gene.percent/gene.n, 
            "Isoform Specific Total Num"=valid.isoform.n , "Isoform Specific percent" = valid.isoform.percent/valid.isoform.n ,
            "Common Total Num"=valid.common.n , "Common Validated percent" = valid.common.percent/valid.common.n))
}
valid.per.cat <- data.frame(lapply(biomarkers, fnValidatedPerCat), check.names=FALSE)


pdf("supplementary_figure_cindex.pdf", height=15, width=9)
par(mfrow=c(5,3))
for (drug in names(validated.biomarkers)){

  tt <- validated.biomarkers[[drug]]

  top <- head(tt[order(tt$gCSI.cindex, decreasing = TRUE),], n=1)

  for (jj in seq_len(3)){
    dataset <- c("ccle", "gdsc", "gCSI")[jj]
    if(dataset!="gCSI"){
      tissueTypes <- ccle.tissuetype 
      gene.data <- ccle.genes.fpkm
      isoform.data <- ccle.isoforms.fpkm
    } else {
      tissueTypes <- gCSI.tissuetype 
      gene.data <- gCSI.genes.fpkm
      isoform.data <- gCSI.isoforms.fpkm
    }

    M0 <- fnCreateNullModel(drug=drug, assay=dataset)
    mTitle <- paste0(c("CCLE", "GDSC", "gCSI")[jj], "-", top[,"biomarker.id"])
    if(top[,"type"]=="gene"){
      M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=gene.data[ ,top[,"gene.id"]])
      an <- NA#sprintf("Cindex: %1.2g", compute.stat("cindex", M2$dataset[,"drug"],fitted(M2$model)))
      tissueOnlyValues <- fitted(M2$model)-M2$model$coefficients[["Gene"]]*M2$data$Gene
      myScatterPlot("", x=abs(M2$model$coefficients[["Gene"]]*M2$data$Gene), y=abs(M2$data$drug - tissueOnlyValues), ylab="Actual Residual", xlab="Predicted Effect on AAC", cex=1.5, main=mTitle, annot=an, method="annotated")

    } else if(top[,"type"]=="isoform"){
      M3B <- fnCreateGeneModel(drug=drug, nullModel=M0, data=isoform.data[ ,top[,"transcript.id"]])
      an <- NA# sprintf("Cindex: %1.2g", compute.stat("cindex", M3B$dataset[,"drug"],fitted(M3B$model)))
      tissueOnlyValues <- fitted(M3B$model)-M3B$model$coefficients[["Gene"]]*M3B$data$Gene
      myScatterPlot("", x=abs(M3B$model$coefficients[["Gene"]]*M3B$data$Gene), y=abs(M3B$data$drug - tissueOnlyValues), ylab="Actual Residual", xlab="Predicted Effect on AAC", cex=1.5, main=mTitle, annot=an, method="annotated")

    }


  }

}
dev.off()

common.cells.gdsc.ccle <- intersect(common.cells.ccle, cellNames(GDSC))

pdf(file = sprintf("gdsc_ccle_inconsistencies.pdf"), height=15, width=9)
par(mfrow=c(5,3))
for(drug in sort(names(all.biomarkers))){

  toPlot1 <- data.frame("GDSC"=gdsc.drug.sensitivity[common.cells.gdsc.ccle,drug], 
                        "CCLE" = ccle.drug.sensitivity[common.cells.gdsc.ccle,drug], check.names=FALSE)
  # ScatterHist(toPlot1, "gCSI", "CCLE", title=sprintf("%s", drug), smoothmethod="lm", se=FALSE)
  myScatterPlot("", x=toPlot1$GDSC, y=toPlot1$CCLE, ylab="CCLE AAC", xlab="GDSC AAC", main=sprintf("%s", drug), cex=2, xlim=c(0,1), ylim=c(0,1))

}
dev.off()








myx <- apply(ccle.genes.fpkm, 2, function(x) return(!all(is.na(x))&&!all(x==0, na.rm=TRUE)))
ccle.genes.fpkm <- ccle.genes.fpkm[,myx]
myx <- apply(ccle.isoforms.fpkm, 2, function(x) return(!all(is.na(x))&&!all(x==0, na.rm=TRUE)))
ccle.isoforms.fpkm <- ccle.isoforms.fpkm[,myx]

test <- apply(ccle.genes.fpkm, 2,function(x) shapiro.test(x)$statistic)
test2 <- apply(ccle.genes.fpkm, 2,function(x) shapiro.test(2^x-1)$statistic)
test3 <- apply(ccle.genes.fpkm, 2,function(x) shapiro.test(log2(2^x-1 + .Machine$double.eps))$statistic)

summary(test)
summary(test2)

wilcox.test(test, test2, paired=TRUE, alternative="greater")

wilcox.test(test3, test, paired=TRUE, alternative="greater")

boxplot(list(test,test2))
