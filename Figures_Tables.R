load("data/PSets/CCLE_hs.RData")
load("data/PSets/GDSC.RData")
load("data/PSets/GRAY_hs.RData")
load("data/PSets/gCSI_hs.RData")
load("data/PSets/UHN_hs.RData")
require(PharmacoGx)
source("code/foo.R")

###Figures
###Figure 1 is created in isoforms.graffle
###Figure 2 is created in training_results.R
###Figure 3 ??
###Figure 4 is created in Final_Validation.R
###Figure 5 is created in Final_Validation.R
###

##Supplementary Tables
##Supplementary Table 1
xx <- as.data.frame(matrix(NA, ncol=6, nrow=5, dimnames=(list(c("CCLE", "GDSC", "gCSI", "GRAY", "UHN"), c("Dataset", "Compounds", "Cell lines", "Tissue types", "RNA-seq availability", "Assay")))))
xx["CCLE", ] <- c(pSetName(CCLE), length(drugNames(CCLE)), length(cellNames(CCLE)), length(table(cellInfo(CCLE)$tissueid)), "rnaseq" %in% names(CCLE@molecularProfiles), "CellTiter-Glo")
xx["GDSC", ] <- c(pSetName(GDSC), length(drugNames(GDSC)), length(cellNames(GDSC)), length(table(cellInfo(GDSC)$tissueid)), "rnaseq" %in% names(GDSC@molecularProfiles), "Syto60")
xx["gCSI", ] <- c(pSetName(gCSI), length(drugNames(gCSI)), length(cellNames(gCSI)), length(table(cellInfo(gCSI)$tissueid)), "rnaseq" %in% names(gCSI@molecularProfiles), "CellTiter-Glo")
xx["GRAY", ] <- c(pSetName(GRAY), length(drugNames(GRAY)), length(cellNames(GRAY)), length(table(cellInfo(GRAY)$tissueid)), "rnaseq" %in% names(GRAY@molecularProfiles), "CellTiter-Glo")
xx["UHN", ] <- c(pSetName(UHN), length(drugNames(UHN)), length(cellNames(UHN)), length(table(cellInfo(UHN)$tissueid)), "rnaseq" %in% names(UHN@molecularProfiles), "CellTiter-Glo")
xtable::print.xtable(xtable::xtable(xx), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "datasets.tex"), append=FALSE)
###

##Supplementary table 2
##validation rate of isformic biomarkers 
drugs <- unionList(intersectList(drugNames(GRAY),
                                 drugNames(CCLE),
                                 drugNames(GDSC)), 
                   intersectList(drugNames(gCSI),
                                 drugNames(CCLE),
                                 drugNames(GDSC)))
drugs <- sort(drugs)
xx <- matrix("-", ncol=6, nrow=length(drugs), dimnames=list(drugs, c("Compound", "Training", "Pre-validation Pan-Cancer","Training Breast Specific", "Pre-validation Breast", "Final validation")))
xx[drugs, "Compound"] <- drugs
load(file.path(path.diagrams, "all.biomarkers.RData"), verbose=TRUE)
xx[drugs, "Training"] <- sapply(all.biomarkers, function(x){length(which(x["type"]=="isoform"))})[drugs]
load(file.path(path.diagrams, "validated.biomarkers.gCSI.RData"), verbose=TRUE)
xx[intersect(drugs, names(validated.biomarkers)), "Pre-validation Pan-Cancer"] <- sapply(validated.biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})[intersect(drugs, names(validated.biomarkers))]
load(file.path(path.diagrams, "breast.biomarkers.RData"), verbose=TRUE)
xx[intersect(drugs, names(breast.biomarkers)), "Training Breast Specific"] <- sapply(breast.biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})[intersect(drugs, names(breast.biomarkers))]
load(file.path(path.diagrams, "validated.biomarkers.gray.RData"), verbose=TRUE)
xx[intersect(drugs, names(validated.biomarkers)), "Pre-validation Breast"] <- sapply(validated.biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})[intersect(drugs, names(validated.biomarkers))]
load(file.path(path.diagrams, "validated.biomarkers.uhn.RData"), verbose=TRUE)
xx[intersect(drugs, names(final.validated.biomarkers)), "Final validation"] <- sapply(final.validated.biomarkers, function(x){if(nrow(x) > 0){length(which(x["type"]=="isoform"))}else{0}})
xtable::print.xtable(xtable::xtable(xx, digits=0, align=c("l","l|","l","l||","l","l","l")), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "validation_rate.tex"), append=FALSE)
###

##Supplementary table 3
## doubling time of UHN cell lines
xx <- read.csv(file.path(path.data, "GrowthCurve_dataCombined.csv"))
dec.batch <- xx[3:10, 2:6]
colnames(dec.batch) <- xx[2, 2:6]
rownames(dec.batch) <- xx[3:10, 1]

nov.batch.1 <- xx[14:18, 2:17]
colnames(nov.batch.1) <- xx[13, 2:17]
rownames(nov.batch.1) <- xx[14:18, 1]

nov.batch.2 <- xx[c(22:26, 28:29), 2:15]
colnames(nov.batch.2) <- xx[21, 2:15]
rownames(nov.batch.2) <- xx[c(22:26, 28:29),1]

nov.batch.3 <- xx[c(33:35, 38:41), 2:16]
colnames(nov.batch.3) <- xx[32, 2:16]
rownames(nov.batch.3) <- xx[c(33:35, 38:41),1]

nov.batch.4 <- xx[44:50, 2:14]
colnames(nov.batch.4) <- xx[43, 2:14]
rownames(nov.batch.4) <- xx[44:50, 1]

doubling.time <- function(initial.conc, final.conc, culture.duration) {
  return((culture.duration*log10(2))/(log10(final.conc)- log10(initial.conc)))
}
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[\xa0]|[ ]"

dec.batch.dt <- sapply(1:ncol(dec.batch), function(x){doubling.time(as.numeric(dec.batch[1,x]), as.numeric(dec.batch[nrow(dec.batch), x]), 8)})
names(dec.batch.dt) <- colnames(dec.batch)
names(dec.batch.dt) <- cellNames(UHN)[match(tolower(gsub(badchars, "", names(dec.batch.dt))),  tolower(gsub(badchars, "", cellNames(UHN))))]
dec.batch.dt <- dec.batch.dt[which(!is.nan(dec.batch.dt))]

nov.batch.1.dt <- sapply(1:ncol(nov.batch.1), function(x){doubling.time(as.numeric(nov.batch.1[1, x]), as.numeric(nov.batch.1[nrow(nov.batch.1), x]), 5)})
names(nov.batch.1.dt) <- colnames(nov.batch.1)
names(nov.batch.1.dt)[which(names(nov.batch.1.dt)=="MPE-600")] <- "600MPE"
names(nov.batch.1.dt)[which(names(nov.batch.1.dt)=="MDA231")] <- "MDA-MB-231"
names(nov.batch.1.dt)[which(names(nov.batch.1.dt)=="MDA436")] <- "MDA-MB-436"
names(nov.batch.1.dt)[which(names(nov.batch.1.dt)=="MDA468")] <- "MDA-MB-468"
names(nov.batch.1.dt)[which(names(nov.batch.1.dt)=="OCUB1")] <- "OCUB-M"
names(nov.batch.1.dt) <- cellNames(UHN)[match(tolower(gsub(badchars, "", names(nov.batch.1.dt))),  tolower(gsub(badchars, "", cellNames(UHN))))]
nov.batch.1.dt <- nov.batch.1.dt[which(!is.nan(nov.batch.1.dt))]

nov.batch.2.dt <- sapply(1:ncol(nov.batch.2), function(x){doubling.time(as.numeric(nov.batch.2[1,x]), as.numeric(nov.batch.2[nrow(nov.batch.2), x]), 8)})
names(nov.batch.2.dt) <- colnames(nov.batch.2)
names(nov.batch.2.dt)[which(names(nov.batch.2.dt) == "MFN223")] <- "MFM223"
names(nov.batch.2.dt)[which(names(nov.batch.2.dt) == "MDA157")] <- "MDA-MB-157"
names(nov.batch.2.dt)[which(names(nov.batch.2.dt)=="OCUB1")] <- "OCUB-M"
names(nov.batch.2.dt) <- cellNames(UHN)[match(tolower(gsub(badchars, "", names(nov.batch.2.dt))),  tolower(gsub(badchars, "", cellNames(UHN))))]
nov.batch.2.dt <- nov.batch.2.dt[which(!is.nan(nov.batch.2.dt))]

nov.batch.3.dt <- sapply(1:ncol(nov.batch.3), function(x){doubling.time(as.numeric(nov.batch.3[1,x]), as.numeric(nov.batch.3[nrow(nov.batch.3), x]), 8)})
names(nov.batch.3.dt) <- colnames(nov.batch.3)
names(nov.batch.3.dt)[which(names(nov.batch.3.dt)=="MDA231")] <- "MDA-MB-231"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt)=="MDA436")] <- "MDA-MB-436"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt)=="MDA468")] <- "MDA-MB-468"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "MDA157")] <- "MDA-MB-157"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "HS578T")] <- "Hs578T"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "JMT1")] <- "JIMT1"
names(nov.batch.3.dt)[which(names(nov.batch.3.dt) == "OCUB")] <- "OCUB-M"
names(nov.batch.3.dt) <- cellNames(UHN)[match(tolower(gsub(badchars, "", names(nov.batch.3.dt))),  tolower(gsub(badchars, "", cellNames(UHN))))]
nov.batch.3.dt <- nov.batch.3.dt[which(!is.nan(nov.batch.3.dt))]

nov.batch.4.dt <- sapply(1:ncol(nov.batch.4), function(x){doubling.time(as.numeric(nov.batch.4[1,x]), as.numeric(nov.batch.4[nrow(nov.batch.4), x]), 7)})
names(nov.batch.4.dt) <- colnames(nov.batch.4)
names(nov.batch.4.dt)[which(names(nov.batch.4.dt)=="MDA175VII")] <- "MDA-MB-175-VII"
names(nov.batch.4.dt)[which(names(nov.batch.4.dt)=="MDA361")] <- "MDA-MB-361"
names(nov.batch.4.dt)[which(names(nov.batch.4.dt)=="MDA330")] <- "MDA-MB-330"
names(nov.batch.4.dt) <- cellNames(UHN)[match(tolower(gsub(badchars, "", names(nov.batch.4.dt))),  tolower(gsub(badchars, "", cellNames(UHN))))]
nov.batch.4.dt <- nov.batch.4.dt[which(!is.nan(nov.batch.4.dt))]


xx <- sort(unionList(names(dec.batch.dt), names(nov.batch.1.dt), names(nov.batch.2.dt), names(nov.batch.3.dt), names(nov.batch.4.dt)))
tt <- NULL
for(cell in xx){
  t <- NA
  if(cell %in% names(dec.batch.dt)) {
    x <- dec.batch.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.1.dt)) {
    x <- nov.batch.1.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.2.dt)) {
    x <- nov.batch.2.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.3.dt)) {
    x <- nov.batch.3.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  if(cell %in% names(nov.batch.4.dt)) {
    x <- nov.batch.4.dt[cell]
    if(!is.nan(x)){
      t <- ifelse(is.na(t), x, median(t, x))
    }
  }
  names(t) <- cell
  tt <- c(tt, round(t, digits=2))
}

dd <- cbind("Cell line"=names(tt), "Doubling time"=tt)
xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "doubling_time.tex"), append=FALSE)
###

###
##Supplementary table 4
##known biomarkers
file.sensitivity.all.types <- "auc_recomputed_drug_association_mut_cnv.RData"
file.sensitivity <- "auc_recomputed_drug_association.RData"
load("data/training_ccle_gdsc_mut_cnv.RData")
annot.ensembl.all.genes.mut <- annot.ensembl.all.genes
annot.ensembl.all.isoforms.mut <- annot.ensembl.all.isoforms
load("data/training_ccle_gdsc.RData")
load(file.path(path.data, file.sensitivity.all.types), verbose=T)
drug.association.mut <- drug.association
drug.association.statistics.mut <- drug.association.statistics
drug.association.best.isoforms.mut <- drug.association.best.isoforms
load(file.path(path.data, file.sensitivity), verbose=T)
knownBiomarkersCheck <-
  function()
  {
    known.biomarkers <- read.csv(file="data/known.biomarkers.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, na.strings=c("", " ", "NA"))
    known.biomarkers <- known.biomarkers[which(!is.na(known.biomarkers$type) & known.biomarkers$type != "fusion"),1:3]
    #known.biomarkers <- known.biomarkers[which(known.biomarkers$type == "expression"),1:2]
    known.biomarkers <- cbind(known.biomarkers,  "gene cindex"=NA, "gene pvalue"=NA, "isoform"=NA, "isoform cindex"=NA, "isoform pvalue"=NA, "null model cindex"=NA)
    
    for(i in 1:nrow(known.biomarkers)) {
      feature <- rownames(annot.ensembl.all.genes)[which(annot.ensembl.all.genes$Symbol == known.biomarkers[i ,"gene"])]
      if(known.biomarkers[i,"type"]== "expression" & length(feature) > 0){
        known.biomarkers[i ,"gene cindex"] <- drug.association.statistics[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M2"]
        known.biomarkers[i ,"gene pvalue"] <- drug.association[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M2"]
        #known.biomarkers[i ,"gene R2"] <- drug.association.statistics[["r.squared"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M2"]
        
        known.biomarkers[i ,"isoform"] <- drug.association.best.isoforms[[feature]][[known.biomarkers[i,"drug"]]]
        
        known.biomarkers[i ,"isoform cindex"] <- drug.association.statistics[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M3B"]
        known.biomarkers[i ,"isoform pvalue"] <- drug.association[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M3B"]
        #known.biomarkers[i ,"isoform R2"] <- drug.association.statistics[["r.squared"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M3B"]
        
        known.biomarkers[i ,"null model cindex"] <- drug.association.statistics[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M0"]
        
        #xx <- rownames(ccle.drug.sensitivity)[which(!is.na(ccle.drug.sensitivity[,known.biomarkers[i,"drug"]]))]
        #xx <- intersect(xx, rownames(ccle.genes.fpkm))
        #known.biomarkers[i ,"mean expr"] <- mean(ccle.genes.fpkm[xx, feature], na.rm=TRUE)
        #breast.cells <- rownames(ccle.cell.profiles)[which(ccle.cell.profiles == "breast")]
        #breast.cells <- intersect(breast.cells, xx)
        #known.biomarkers[i ,"breast expr"] <- mean(ccle.genes.fpkm[breast.cells, feature], na.rm=TRUE)
        #known.biomarkers[i ,"mean isoform expr"] <- mean(ccle.isoforms.fpkm[xx, known.biomarkers[i ,"isoform"]], na.rm=TRUE)
        #known.biomarkers[i ,"breast isoform expr"] <- mean(ccle.isoforms.fpkm[breast.cells, known.biomarkers[i ,"isoform"]], na.rm=TRUE)
        
      }else{
        feature <- rownames(annot.ensembl.all.genes.mut)[which(annot.ensembl.all.genes.mut$Symbol == known.biomarkers[i ,"gene"])]
        
        known.biomarkers[i ,"null model cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M0"]
        if(known.biomarkers[i,"type"]== "expression"){
          known.biomarkers[i ,"gene cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M2"]
          known.biomarkers[i ,"gene pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M2"]
          
          known.biomarkers[i ,"isoform"] <- drug.association.best.isoforms.mut[[feature]][[known.biomarkers[i,"drug"]]]
          
          known.biomarkers[i ,"isoform cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "M3B"]
          known.biomarkers[i ,"isoform pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "M3B"]
        }else if(known.biomarkers[i,"type"]== "mutation"){
          known.biomarkers[i ,"gene cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "MM"]
          known.biomarkers[i ,"gene pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "MM"]
        }else if(known.biomarkers[i,"type"]== "amplification"){
          known.biomarkers[i ,"gene cindex"] <- drug.association.statistics.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]]["mean", "MC"]
          known.biomarkers[i ,"gene pvalue"] <- drug.association.mut[["cindex"]][[feature]][[known.biomarkers[i,"drug"]]][1, "MC"]
        }
      }
    }
    known.biomarkers <- known.biomarkers[which(known.biomarkers$`gene cindex` != 0), ]
    #colnames(known.biomarkers)[1:3] <- capitalize(colnames(known.biomarkers)[1:3])
    xtable::print.xtable(xtable::xtable(known.biomarkers, digits=c(0, 0, 0, 0, 2, -1, 0, 2, -1, 2)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "known_biomarkers.tex"), append=FALSE)
    #write.csv(known.biomarkers, file=file.path(path.diagrams, "known_biomarkers.csv"), row.names = F)
  }
knownBiomarkersCheck()
###

#######
##Supplementary table 5
##Novel isoformic biomarkers
load(file.path(path.diagrams, "validated.biomarkers.uhn.RData"), verbose=TRUE)
rr <- c("drug", "gene", "isoform", "Training cindex", "Training corrected pvalue", "Final validation cindex")
xx <- as.data.frame(matrix(NA, ncol=length(rr), nrow=sum(sapply(final.validated.biomarkers,dim)[1,]), dimnames=(list(unlist(sapply(final.validated.biomarkers,function(x){x["biomarker.id"]})), rr))))
i <- 0
drugs <- sort(names(final.validated.biomarkers))
for(drug in drugs){
  rr <- order(final.validated.biomarkers[[drug]]$UHN.cindex, decreasing=T)
  for(j in rr){
    i <- i + 1
    xx[i ,"drug"] <- drug
    isoform <- final.validated.biomarkers[[drug]][j,"biomarker.id"]
    #feature <- final.validated.biomarkers[[drug]][i, "gene.id"]
    feature <- annot.ensembl.all.isoforms[which(annot.ensembl.all.isoforms$EnsemblTranscriptId == isoform), "EnsemblGeneId"]
    
    xx[i ,"gene"] <- annot.ensembl.all.genes[feature, "Symbol"]
    #xx[i ,"gene cindex"] <- drug.association.statistics[["cindex"]][[feature]][[drug]]["mean", "M2"]
    #xx[i ,"gene pvalue"] <- drug.association[["cindex"]][[feature]][[drug]][1, "M2"]
    xx[i ,"isoform"] <- isoform
    xx[i ,"Training cindex"] <- final.validated.biomarkers[[drug]][j,"cindex"]
    xx[i ,"Training corrected pvalue"] <- final.validated.biomarkers[[drug]][j,"fdr"]
    xx[i ,"Final validation cindex"] <- round(as.numeric(final.validated.biomarkers[[drug]][j,"UHN.cindex"]), digits=2)
    
  }
}
xtable::print.xtable(xtable::xtable(xx, digits=c(0, 0, 0, 0, 2, -1, 2)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "novel_biomarkers.tex"), append=FALSE)
###

##Supplementary Figures
##Supplementary Figure 1 is created in isoforms.graffle
##Supplementary Figure 2 is created in isoforms.graffle
##Supplementary Figure 3
##
require(VennDiagram) || stop("Library gdata is not available!")
mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
common <- intersectPSet(pSets = list("CCLE"=CCLE, "GDSC"=GDSC), intersectOn = c("cell.lines", "drugs"), strictIntersect = TRUE)

### venn diagram of common cell lines between studies
pdf(file.path(path.diagrams, "celllines.CCLE.GDSC.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@cell), area2=nrow(GDSC@cell), cross.area=nrow(common$CCLE@cell), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()

### venn diagram of common drugs between studies
pdf(file.path(path.diagrams, "drugs.CCLE.GDSC.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@drug), area2=nrow(GDSC@drug), cross.area=nrow(common$CCLE@drug), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()                 

##pie chart of common tissues between studies
#mycol <- RColorBrewer::brewer.pal(n=7, name="Set3")
mycol <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f",
           "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")#,"#b15928")
pdf(file = file.path(path.diagrams, "tissues.pdf"), height=8, width=10)
temp <- c(table(common$GDSC@cell[,"tissueid"])[15], table(common$GDSC@cell[,"tissueid"])[1:14],   table(common$GDSC@cell[,"tissueid"])[16] ,table(common$GDSC@cell[,"tissueid"])[18:23], table(common$GDSC@cell[,"tissueid"])[17])               
pie(temp, labels = gsub("_", " ", capitalize(names(temp))), col=mycol, main="", radius=1, cex=0.8)#,theta=pi/4, labelcex = .8, labelrad = 1.1)
dev.off()
###

##Supplementary Figure 4
##rnaseq/microarray concordance
if(file.exists(file.path(path.data, "rnaseq_microarray.RData"))){
  load(file.path(path.data, "rnaseq_microarray.RData"))
}else{
  cells <- intersect(pData(CCLE@molecularProfiles$rnaseq)[,"cellid"], pData(CCLE@molecularProfiles$rna)[,"cellid"])
  rnaseq.samples <- rownames(pData(CCLE@molecularProfiles$rnaseq))[match(cells, pData(CCLE@molecularProfiles$rnaseq)[,"cellid"])]
  rna.samples <- rownames(pData(CCLE@molecularProfiles$rna))[match(cells, pData(CCLE@molecularProfiles$rna)[,"cellid"])]
  
  
  features <- intersect(fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"], fData(CCLE@molecularProfiles$rna)[,"EnsemblGeneId"])
  rnaseq.features <- rownames(fData(CCLE@molecularProfiles$rnaseq))[match(features, fData(CCLE@molecularProfiles$rnaseq)[,"EnsemblGeneId"])]
  rna.features <- rownames(fData(CCLE@molecularProfiles$rna))[match(features, fData(CCLE@molecularProfiles$rna)[,"EnsemblGeneId"])]
  ccle.rnaseq.microarray.cor <- cor(exprs(CCLE@molecularProfiles$rna)[rna.features, rna.samples], exprs(CCLE@molecularProfiles$rnaseq)[rnaseq.features, rnaseq.samples], use="pairwise.complete.obs", method="spearman")
  save(ccle.rnaseq.microarray.cor, file=file.path(path.data, "rnaseq_microarray.RData"))
}

pdf(file.path(path.diagrams, "rnaseq_microarray.pdf"), height=7, width=7)
par(mar=c(9,5,5,2))
boxplot(cbind("Identical"=diag(ccle.rnaseq.microarray.cor), "Different"=c(ccle.rnaseq.microarray.cor[upper.tri(ccle.rnaseq.microarray.cor)], ccle.rnaseq.microarray.cor[lower.tri(ccle.rnaseq.microarray.cor)])), las = 2, col = "gray", cex.lab=1, cex.axis=1, pch=19, ylab="spearman correlation", outpch=20, outcex=0.5)
dev.off()
###

###supplementary Figure 5 is created in training_results by training on mutation_cnv data

###supplementary Figure 6
##intersection of cells/drugs between ccle, gdsc, gray
mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
library(VennDiagram)
pdf(file.path(file.path(path.diagrams, "gray_ccle_gdsc_celllines.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(CCLE@cell$cellid), 
                                          area2 = length(GDSC@cell$cellid),
                                          area3 = length(GRAY@cell$cellid),
                                          n12 = length(intersect(CCLE@cell$cellid, GDSC@cell$cellid)),
                                          n13 = length(intersect(CCLE@cell$cellid, GRAY@cell$cellid)),
                                          n23 = length(intersect(GDSC@cell$cellid, GRAY@cell$cellid)),
                                          n123 = length(PharmacoGx::intersectList(CCLE@cell$cellid, GDSC@cell$cellid, GRAY@cell$cellid)),
                                          category = c("CCLE", "GDSC", "GRAY"),
                                          col = mycol[1:3],
                                          fill = mycol[1:3],
                                          margin=0.10,
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()

pdf(file.path(file.path(path.diagrams, "gray_ccle_gdsc_breast_celllines.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")]), 
                                          area2 = length(GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")]),
                                          area3 = length(GRAY@cell$cellid),
                                          n12 = length(intersect(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")], GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")])),
                                          n13 = length(intersect(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")], GRAY@cell$cellid)),
                                          n23 = length(intersect(GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")], GRAY@cell$cellid)),
                                          n123 = length(PharmacoGx::intersectList(CCLE@cell$cellid[which(CCLE@cell$tissueid=="breast")], GDSC@cell$cellid[which(GDSC@cell$tissueid=="breast")], GRAY@cell$cellid)),
                                          category = c("CCLE", "GDSC", "GRAY"),
                                          col = mycol[1:3],
                                          fill = mycol[1:3],
                                          margin=0.10,
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()


pdf(file.path(file.path(path.diagrams, "gray_ccle_gdsc_drugs.pdf")), height=4, width=4)
venn.plot <-VennDiagram::draw.triple.venn(area1 = length(rownames(CCLE@drug)), 
                                          area2 = length(rownames(GDSC@drug)),
                                          area3 = length(rownames(GRAY@drug)),
                                          n12 = length(intersect(rownames(CCLE@drug), rownames(GDSC@drug))),
                                          n13 = length(intersect(rownames(CCLE@drug), rownames(GRAY@drug))),
                                          n23 = length(intersect(rownames(GDSC@drug), rownames(GRAY@drug))),
                                          n123 = length(PharmacoGx::intersectList(rownames(CCLE@drug), rownames(GDSC@drug), rownames(GRAY@drug))),
                                          category = c("CCLE", "GDSC", "GRAY"),
                                          col = mycol[1:3],
                                          fill = mycol[1:3],
                                          margin=0.10,
                                          lty = "blank",cex = 1,cat.cex = 1,cat.col = c("black", "black","black"))
dev.off()
###

###Supplementary Figure 7
###intersection of cells/drugs between ccle, gdsc, gCSI
### ???
###

###Supplementary Figure 8
###inconsistency of pharmacological profiles between CCLE and GDSC
### ???
###

##Supplementary Figure 9
##inconsistency of pharmacological profiles for sorafenib and critozenib between 3 studies
gray.ccle.cells <- intersect(cellNames(GRAY), cellNames(CCLE))
gray.gdsc.cells <- intersect(cellNames(GRAY), cellNames(GDSC))
gdsc.ccle.cells <- intersect(cellNames(GDSC), cellNames(CCLE))

gray.sensitivity <- PharmacoGx::summarizeSensitivityProfiles(GRAY, sensitivity.measure = "auc_recomputed", drugs = c("Crizotinib", "Sorafenib"))
ccle.sensitivity <- PharmacoGx::summarizeSensitivityProfiles(CCLE, sensitivity.measure = "auc_recomputed", drugs = c("Crizotinib", "Sorafenib"))
gdsc.sensitivity <- PharmacoGx::summarizeSensitivityProfiles(GDSC, sensitivity.measure = "auc_recomputed", drugs = c("Crizotinib", "Sorafenib"))
myScatterPlot(Name=file.path(path.diagrams, "gray_ccle_sorafenib.pdf"), 
              x=ccle.sensitivity["Sorafenib", gray.ccle.cells], 
              y=gray.sensitivity["Sorafenib", gray.ccle.cells], 
              method="plain", 
              minp=10,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "gray_gdsc_sorafenib.pdf"), 
              x=gdsc.sensitivity["Sorafenib", gray.gdsc.cells], 
              y=gray.sensitivity["Sorafenib", gray.gdsc.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="GDSC",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "ccle_gdsc_sorafenib.pdf"), 
              x=ccle.sensitivity["Sorafenib", gdsc.ccle.cells], 
              y=gdsc.sensitivity["Sorafenib", gdsc.ccle.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GDSC")


myScatterPlot(Name=file.path(path.diagrams, "gray_ccle_crizotinib.pdf"), 
              x=ccle.sensitivity["Crizotinib", gray.ccle.cells], 
              y=gray.sensitivity["Crizotinib", gray.ccle.cells], 
              method="plain", 
              minp=10,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "gray_gdsc_crizotinib.pdf"), 
              x=gdsc.sensitivity["Crizotinib", gray.gdsc.cells], 
              y=gray.sensitivity["Crizotinib", gray.gdsc.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="GDSC",
              ylab="GRAY")

myScatterPlot(Name=file.path(path.diagrams, "ccle_gdsc_crizotinib.pdf"), 
              x=ccle.sensitivity["Crizotinib", gdsc.ccle.cells], 
              y=gdsc.sensitivity["Crizotinib", gdsc.ccle.cells], 
              method="plain", 
              minp=1,
              xlim=c(0, 1), 
              ylim=c(0, 1),
              xlab="CCLE",
              ylab="GDSC")

###
##Supplementary Figure 10
###predicted AAC versus actual values for top predicted biomarkers in for 5 drugsin common between CCLE, GDSC and gCSI
###??? drug names are missed in the plot
###

###Supplementary Figure 11 is created in training_results.R
###Supplementary Figure 12 is created in training_results.R

###Supplementary Figure 13
###the heatmaps for the other drugs not included in uhn (not in Fig 4)
source("code/foo.R")
source("code/foo_PreValidation.R")
source("code/foo_FinalValidation.R")
mycol <- RColorBrewer::brewer.pal(n=4, name="Set1")
red <- mycol[1]  
blue <- mycol[2]

load(file.path(path.diagrams, "validated.biomarkers.gray.RData"), verbose=TRUE)
drugs <- setdiff(names(validated.biomarkers), drugNames(UHN))
for(drug in drugs) {
  if(nrow(validated.biomarkers[[drug]]) > 0){
    vtt <- validated.biomarkers[[drug]]
    #vtt <- vtt[which(vtt[,"type"] == "isoform"), , drop=F]
    #vtt <- vtt[which(vtt[, "isoforms.no"] > 1), , drop=FALSE]
    gray.specificity <- vtt[ ,"gray.specificity"]
    xx <- fnFetchBiomarkers(top.significant.biomarkers=vtt, drug=drug, indices=1:nrow(vtt))
    xx <- do.call(rbind, xx)
    xx[,"short.label"] <- gsub(".ISO$","",xx[,"short.label"])
    xx <- apply(xx, 1, function(x){x})
    
    rr <- fnPlotAUCoverCellLinesGray(drug=drug, tissue.type="all", biomarkers=xx, suffix="isoform.specific", gray.specificity=gray.specificity)#, biomarkers.toPlot)
    if(all(!is.na(rr))){
      biomarkers.order <- rr$hv$rowInd
      fnPlotEffectSize(drug, biomarkers=xx, effect.size=effect.size, biomarkers.order) 
      xx <- do.call(rbind, xx)
      names(biomarkers.order) <- xx[biomarkers.order, "isoform.id"]
    }
    
  }
}
###
###Supplementary Figure 14 is created in Final_Validation.R
###Supplementary Figure 15 is created in Final_Validation.R
###Supplementary Figure 16 ???

###supplementary Files 
###supplementary File 1
### list of cell lines in all studies
library(gdata)
dd <- gdata::cbindX(data.frame(cellNames(CCLE)), 
             data.frame(cellNames(GDSC)) , 
             data.frame(intersect(cellNames(CCLE), cellNames(GDSC))), 
             data.frame(cellNames(gCSI)) , 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(gCSI))), 
             data.frame(cellNames(GRAY)), 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(GRAY))), 
             data.frame(cellNames(UHN)), 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(GRAY), cellNames(UHN))))
colnames(dd) <- c("CCLE", "GDSC", "CCLE_GDSC", "gCSI", "CCLE_GDSC_gCSI", "GRAY", "CCLE_GDSC_GRAY", "UHN", "CCLE_GDSC_GRAY_UHN")
#xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "cells.tex"), append=FALSE)
write.csv(dd, file=file.path(path.diagrams, "cells.csv"), na="", row.names=FALSE)
###

###supplementary File 2
### list of drugs in all studies
dd <- cbindX(data.frame(drugNames(CCLE)), 
             data.frame(drugNames(GDSC)) , 
             data.frame(intersect(drugNames(CCLE), drugNames(GDSC))), 
             data.frame(drugNames(gCSI)), 
             data.frame(intersectList(drugNames(CCLE), drugNames(GDSC), drugNames(gCSI))), 
             data.frame(drugNames(GRAY)), 
             data.frame(intersectList(drugNames(CCLE), drugNames(GDSC), drugNames(GRAY))), 
             data.frame(drugNames(UHN)), 
             data.frame(intersectList(drugNames(CCLE), drugNames(GDSC), drugNames(GRAY), drugNames(UHN))))
colnames(dd) <- c("CCLE", "GDSC", "CCLE_GDSC", "gCSI", "CCLE_GDSC_gCSI", "GRAY", "CCLE_GDSC_GRAY", "UHN", "CCLE_GDSC_GRAY_UHN")
#xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "drugs.tex"), append=FALSE)
write.csv(dd, file=file.path(path.diagrams, "drugs.csv"), na="", row.names=FALSE)

###
###supplementary File 3
##list of predicted biomarkers in training phase
load(file.path(path.diagrams, "all.biomarkers.RData"))
WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=file.path(path.diagrams, "predicted_biomarkers_training.xlsx"), row.names=TRUE)
###
###supplementary File 4
##list of validated biomarkers in gCSI
load(file.path(path.diagrams, "validated.biomarkers.gCSI.RData"), verbose=TRUE)
WriteXLS::WriteXLS("validated.biomarkers", ExcelFileName=file.path(path.diagrams, "validated_biomarkers_gCSI.xlsx"), row.names=TRUE)
###
###supplementary File 4
##list of pre-validated breast biomarkers in GRAY
load(file.path(path.diagrams, "validated.biomarkers.gray.RData"), verbose=TRUE)
WriteXLS::WriteXLS("validated.biomarkers", ExcelFileName=file.path(path.diagrams, "prevalidated_biomarkers_gray.xlsx"), row.names=TRUE)
###
###supplementary File 4
##list of validated breast biomarkers in UHN
load(file.path(path.diagrams, "validated.biomarkers.uhn.RData"), verbose=TRUE)
WriteXLS::WriteXLS("biomarkers", ExcelFileName=file.path(path.diagrams, "validated_biomarkers_uhn.xlsx"), row.names=TRUE)
###
