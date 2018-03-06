args <- commandArgs(trailingOnly=TRUE)
##for test
#args <- c("../results/training_results", ../data/training_ccle_gdsc.RData", "../data/auc_recomputed_drug_association.RData")

if(!require("PharmacoGx")){biocLite("PharmacoGx");library(PharmacoGx)}
if(!require("Biobase")){biocLite("Biobase");library(Biobase)}
stat <- "r.squared & cindex"

path.training.results <- file.path(as.character(args[1]))
pvalues.files <- list.files(file.path(path.training.results))

path.training.data <- as.character(args[2])
load(path.training.data)
GeneList <- colnames(ccle.genes.fpkm)
drugs <- colnames(ccle.drug.sensitivity)

pvalues = list()
best.isoforms = list()
statistics = list()

length(pvalues) <- length(GeneList)
length(best.isoforms) <-length(GeneList)
length(statistics) <-length(GeneList)

###
if(stat == "r.squared & cindex"){
  pvalues <- list("r.squared"=pvalues, "cindex"=pvalues)
  statistics <- list("r.squared"=statistics, "cindex"=statistics)
}
for ( i in 1: length(pvalues.files))
{
  load(file.path(path.training.results, pvalues.files[i]))
  Index <- unlist(strsplit(pvalues.files[i],"[.,_]"))
  for(j in Index[1]:Index[2])
  {
    t <- j-as.numeric(Index[1])+1
    for(k in 1:length(drugs)) #24 ccle #15 ccle & gdsc
    {   
      pvalues[["r.squared"]][[j]][[k]] <- both.drug.association.adj.r.squared.pvalues[["r.squared"]][[t]][[k]]
      pvalues[["cindex"]][[j]][[k]] <- both.drug.association.adj.r.squared.pvalues[["cindex"]][[t]][[k]]
      statistics[["r.squared"]][[j]][[k]] <- both.drug.association.statistics[["r.squared"]][[t]][[k]]
      statistics[["cindex"]][[j]][[k]] <- both.drug.association.statistics[["cindex"]][[t]][[k]]
      names(pvalues[["r.squared"]][[j]])[k] <- drugs[k]
      names(pvalues[["cindex"]][[j]])[k] <- drugs[k]
      names(statistics[["r.squared"]][[j]])[k] <- drugs[k]
      names(statistics[["cindex"]][[j]])[k] <- drugs[k]
    }
    best.isoforms[[j]] <- both.drug.association.best.isoforms[[t]]
    names(best.isoforms)[j] <-names(both.drug.association.best.isoforms)[t]
    names(pvalues[["r.squared"]])[j] <- names(both.drug.association.adj.r.squared.pvalues[["r.squared"]])[t]
    names(pvalues[["cindex"]])[j] <- names(both.drug.association.adj.r.squared.pvalues[["cindex"]])[t]
    names(statistics[["r.squared"]])[j] <- names(both.drug.association.statistics[["r.squared"]])[t]
    names(statistics[["cindex"]])[j] <- names(both.drug.association.statistics[["cindex"]])[t]
  }
}

drug.association <- pvalues
drug.association.statistics <- statistics
drug.association.best.isoforms <- best.isoforms
save(drug.association,drug.association.statistics, drug.association.best.isoforms, file=as.character(args[3]))
