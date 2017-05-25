R2gauss<- function(y,model){
  moy<-mean(y)
  N<- length(y)
  p<-length(model$coefficients)-1
  SSres<- sum((y-predict(model))^2)
  SStot<-sum((y-moy)^2)
  R2<-1-(SSres/SStot)
  Rajust<-1-(((1-R2)*(N-1))/(N-p-1))
  return(data.frame(R2,Rajust,SSres,SStot))
}
rsq <- function(model, indices) {
  
  indices <- unique(indices)
  data <- model$data
  formula <- model$formula
  object <- model$object
  
  VARS <- colnames(data)
  LHS <- VARS[1]
  RHS <- VARS[-1]
  PRED.pos <- match(RHS, colnames(data))
  PRED.name <- RHS[which(!is.na(PRED.pos))]
  PRED.pos <- as.numeric(na.omit(PRED.pos))
  RESP.pos <- match(LHS, colnames(data))
  
  newDATA <- as.data.frame(data[indices,]) # allows boot to select sample
  newDATA <- newDATA[which(!is.na(newDATA[,ncol(newDATA)])), ,drop=FALSE]
  if(is.null(tissue))
  {
    newDATA <- subset(newDATA, !(newDATA$Tissue %in% row.names(table(newDATA$Tissue))[which(table(newDATA$Tissue)== 1)] ))
    newDATA$Tissue <- factor(newDATA$Tissue)
  }
  colnames(newDATA) <- colnames(data)
  
  if(is.null(tissue))
  {
    newPRED <- as.data.frame(data[-indices, PRED.pos])
    colnames(newPRED) <- PRED.name
    rownames(newPRED) <- rownames(data[-indices,])
    newPRED <- newPRED[which(!is.na(newPRED[,ncol(newPRED)])), ,drop=FALSE]
    
    newPRED <- subset(newPRED, newPRED$Tissue %in% levels(newDATA$Tissue))
    newPRED$Tissue <- factor(newPRED$Tissue)
    
    # newDATA <- subset(newDATA, newDATA$Tissue %in% levels(newPRED$Tissue))
    # newDATA$Tissue <- factor(newDATA$Tissue)
  }else
  {
    newPRED <- as.data.frame(data[-indices,])
    colnames(newPRED) <- colnames(data)
    rownames(newPRED) <- rownames(data)[-indices]
  }
  
  switch(model.method, "glm"={
    if(length(levels(newDATA$Tissue)) == 1 || length(levels(as.factor(newDATA$Gene))) == 1 || length(levels(as.factor(newDATA$drug))) == 1){
      return(list(rmsd=0, cindex=0, r.squared=-100, adj.r.squared=-100))
    }
    newMOD <- glm(formula=formula, data=newDATA, family=glm.family)
    if(glm.family== "binomial"){
      y.hat <- predict(newMOD, newdata=newPRED,  type="response")
    }else{
      y.hat <- predict(newMOD, newdata=newPRED)
    }
  }, "penalized"={
    newMOD <- penalized(newDATA[,1], penalized=newDATA[,3:ncol(newDATA)],unpenalized=~0+Tissue, data=newDATA, lambda1=1, lambda2=0)
    y.hat <- as.numeric(predict(newMOD, newPRED[,2:ncol(newPRED)], data=newPRED)[, "mu"])
  },
  "npreg"={
    ff <- paste(paste0("newDATA$", colnames(data)[1]), paste( paste0("newDATA$", colnames(data)[2:ncol(data)]), collapse=" + "),sep=" ~ ")
    newMOD <- npreg(as.formula(formula),gradients=TRUE, data=newDATA)
    y.hat <- predict(newMOD, newdata=newPRED)
  })
  if(glm.family== "binomial"){
    test <- as.numeric(data[rownames(newPRED),RESP.pos]) - 1 - 1
    test[test== -1] <- 1
    boot.res <- test - y.hat#1 - y.hat
  }else{
    boot.res <- data[rownames(newPRED),RESP.pos] - y.hat
  }
  
  RSS <- sum(boot.res^2)
  RMSD <- sqrt(RSS/length(RSS))
  
  if(is.null(tissue))
  {
    CIndex <- survcomp::concordance.index(x=-y.hat, surv.time=data[rownames(newPRED),RESP.pos], surv.event=rep(1, length(y.hat)), na.rm=TRUE, outx=TRUE)[[1]]
  }else
  {
    CIndex <- 0
  }
  
  Yi <- data[rownames(newPRED),RESP.pos] #fitted(newMOD) + residuals(newMOD) # fitted(object) - residuals(object)
  TSS <- sum((Yi - mean(Yi))^2)
  
  R.square <- 1- (RSS/TSS)
  
  p <- ncol(data) - 1
  n <- nrow(data)
  if (p== (n - 1)){ p <- p - 1}
  R.square <- ifelse(R.square== -Inf, -100, R.square)
  Adjusted.R.square <- 1 - (1 - R.square) * ((n - 1)/(n - p - 1))
  Adjusted.R.square <- ifelse(Adjusted.R.square == -Inf, -100, Adjusted.R.square)
  xx<- R2gauss(y=newMOD$y, model=newMOD)
  return(list( rmsd=RMSD, cindex=CIndex, r.squared=R.square, adj.r.squared=Adjusted.R.square, r.squared.train=xx[1,"R2"], adj.r.squared.train=xx[1,"Rajust"]))
}
fnboot <- function(models, R) {
  boot.res <- list()
  for(i in 1:length(models))
  {
    boot.res[[i]] <- numeric()
    names(boot.res)[i] <- names(models)[i]
  }
  boot.res2 <- boot.res
  idx <- 1:nrow(models[[1]]$data)
  set.seed(1)
  for(i in 1:R)
  {
    indices <- sample(idx,replace=TRUE)
    for(j in 1:length(models))
    {
      M <- names(models)[j]
      if(!is.null(models[[M]]$object))
      {
        switch(stat, "r.squared"={
          boot.res[[M]] <- c(boot.res[[M]], rsq(model=models[[M]], indices)$r.squared)
        }, "rmsd"= {
          boot.res[[M]] <- c(boot.res[[M]], rsq(model=models[[M]], indices)$rmsd)
        }, "cindex"={
          boot.res[[M]] <- c(boot.res[[M]], rsq(model=models[[M]], indices)$cindex)
        }, "adj.r.squared"={
          boot.res[[M]] <- c(boot.res[[M]], rsq(model=models[[M]], indices)$adj.r.squared)
        }, "r.squared & cindex"={
          xx <- rsq(model=models[[M]], indices)
          boot.res[[M]] <- c(boot.res[[M]], xx$r.squared)
          boot.res2[[M]] <- c(boot.res2[[M]], xx$cindex)
        })
      }
    }
  }
  if(stat == "r.squared & cindex") {
    return(list(cindex=boot.res2, r.squared=boot.res))
  }
  return(boot.res)
}
fncounter <- function (i) {
  if (i%%10== 0)
    cat(i)
  else cat(".")
  if (i%%50== 0)
    cat("\n")
  flush.console()
}
fnfetchData  <- function (data) {
  DATA <- data
  VARS <- colnames(DATA)
  LHS <- VARS[1]
  RHS <- VARS[-1]
  PRED.pos <- match(RHS, colnames(DATA))
  PRED.name <- RHS[which(!is.na(PRED.pos))]
  PRED.pos <- as.numeric(na.omit(PRED.pos))
  RESP.pos <- match(LHS, colnames(DATA))
  return(list(data=DATA, pred.pos=PRED.pos, resp.pos=RESP.pos,
              pred.name=PRED.name))
}
fnPRESS <- function (object, data, formula, verbose=TRUE) {
  fetchDATA <- fnfetchData(data)
  DATA <- fetchDATA$data
  
  PRED.pos <- fetchDATA$pred.pos
  RESP.pos <- fetchDATA$resp.pos
  PRED.name <- fetchDATA$pred.name
  PRESS.res <- vector("numeric", nrow(DATA))
  for (i in 1:nrow(DATA)) {
    if (verbose) {
      fncounter(i)
      flush.console()
    }
    newDATA <- DATA[-i, ]
    #if (class(object)== "pcrfit")
    #newMOD <- pcrfit(newDATA, cyc=1, fluo=2, model=object$MODEL,
    #verbose=FALSE)
    #else newMOD <- update(object, data=newDATA)
    switch(model.method, "glm"={
      newMOD <- glm(formula=formula, data=newDATA, family=glm.family)
      newPRED <- as.data.frame(DATA[i, PRED.pos])
      colnames(newPRED) <- PRED.name
      y.hat <- as.numeric(predict(newMOD, newdata=newPRED))
      PRESS.res[i] <- DATA[i, RESP.pos] - y.hat
    }, "Penalized"={
      newDATA[, "Tissue"] <- factor(newDATA[, "Tissue"])
      newMOD <- penalized(newDATA[,1], penalized=newDATA[,3:ncol(newDATA)],unpenalized=~0+Tissue, data=newDATA, lambda1=1, lambda2=0)
      newPRED <- as.data.frame(DATA[i, PRED.pos])
      colnames(newPRED) <- PRED.name
      y.hat <- as.numeric(predict(newMOD, newPRED[,2:ncol(newPRED)], data=newPRED)["mu"])
      PRESS.res[i] <- DATA[i, RESP.pos] - y.hat
    })
  }
  if (verbose)
    cat("\n")
  Yi <- residuals(object) - fitted(object)
  TSS <- sum((Yi - mean(Yi))^2)
  RSS <- sum(PRESS.res^2)
  P.square <- 1 - (RSS/TSS)
  return(list(stat=sum(PRESS.res^2), residuals=PRESS.res,
              P.square=P.square))
}
fnIsoforms_of_Gene <- function(Gene_Map=ensembl.map.genes.isoforms, GeneID) {
  return(data.frame(strsplit(Gene_Map[,as.character(GeneID)], split=","), stringsAsFactors=FALSE))
}
fnIsoformsExp <- function(Isoforms_FPKM, GeneId) {
  linearFormula_Iso <- ""
  Isoformsdf <- data.frame(matrix(NA, nrow=nrow(Isoforms_FPKM), ncol=0, byrow=FALSE))
  rownames(Isoformsdf) <- rownames(Isoforms_FPKM)
  Isoforms_of_Genes <- fnIsoforms_of_Gene(GeneID=GeneId)
  for(j in 1: nrow(Isoforms_of_Genes))
  {
    if (as.character(Isoforms_of_Genes[j,]) %in% colnames(Isoforms_FPKM))
    {
      Isoformsdf[,as.character(Isoforms_of_Genes[j,])] <- Isoforms_FPKM [,as.character(Isoforms_of_Genes[j,])]
      linearFormula_Iso <- paste(linearFormula_Iso , as.character(Isoforms_of_Genes[j,]) , "+")
    }
  }
  linearFormula_Iso <-  substr(linearFormula_Iso, 1, nchar(linearFormula_Iso)-1)
  return(list(FPKM=Isoformsdf, names=linearFormula_Iso))
}
fnCreateNullModel <- function(drug, assay=c("ccle", "gdsc", "gray", "gCSI")) {
  tissues.no <- 1
  if(is.null(tissue))
  {
    linearFormula0 <- "drug ~ Tissue"
  }else{
    linearFormula0 <- "drug ~ 1"
  }
  switch(assay, 
         "ccle"={Sensitivity=ccle.drug.sensitivity[ccle.cells, ]}, 
         "gdsc"={Sensitivity=gdsc.drug.sensitivity[gdsc.cells, ]}, 
         "gray"={Sensitivity=gray.drug.sensitivity},
         "gCSI"={Sensitivity=gCSI.drug.sensitivity})
  Genedf0 <- data.frame(Sensitivity[, drug])
  rownames(Genedf0) <- rownames(Sensitivity)
  colnames(Genedf0)[1] <- "drug"
  if(is.null(tissue))
  {
    Genedf0[, "Tissue"] <- factor(tissueTypes[rownames(Genedf0), "tissue.type"], ordered=FALSE)
    Genedf0 <- Genedf0[complete.cases(Genedf0),]
    Genedf0$Tissue <- factor(Genedf0$Tissue, ordered=FALSE)
    tissues.no <- nrow(table(Genedf0$Tissue))
  }else{
    CellLines <- rownames(subset(tissueTypes, tissueTypes$tissue.type== tissue))
    Genedf0 <- subset(Genedf0, rownames(Genedf0) %in% CellLines)
    Genedf0 <- subset(Genedf0, !is.na(Genedf0[,1]))
  }
  switch(model.method,
         "npreg"={
           M0 <- npreg(as.formula(linearFormula0), data=Genedf0, gradients=TRUE)
         },
         "glm"={
           if(glm.family== "linear")
           {
             M0 <- lm (formula=as.formula(linearFormula0) , data=Genedf0)
           }else{
             M0 <- glm (formula=as.formula(linearFormula0) , data=Genedf0, family=glm.family)
           }
         })
  return(list(n=nrow(Genedf0), model=M0, dataset=Genedf0, formula=linearFormula0, tissues.no=tissues.no))
}
fnCreateGeneModel <- function(drug, nullModel, data) {
  if(is.null(tissue))
  {
    linearFormula <- "drug ~ Tissue + Gene"
  }else{
    linearFormula <- "drug ~ Gene"
  }
  
  GenedfM <- nullModel$dataset
  GenedfM[, "Gene"] <- data[rownames(nullModel$dataset)]
  if(length(table(GenedfM$Tissue)) == 1 || length(table(GenedfM$Gene)) == 1 || length(table(GenedfM$drug)) == 1){
    return(list(model=NULL, dataset=GenedfM, formula=linearFormula, coefficient=0))
  }
  
  switch(model.method, "glm"={
    if(glm.family== "linear")
    {
      M1 <- lm (as.formula(linearFormula) , GenedfM)
    }else{
      M1 <- glm (as.formula(linearFormula) , GenedfM, family=glm.family)
    }
  }, "npreg"={
    M1 <- npreg(as.formula(linearFormula), data=GenedfM, gradients=TRUE)
  })
  marker.index <- grep("Gene", names(coefficients(M1)))
  if(length(marker.index) == 0) {
    return(list(model=NULL, dataset=GenedfM, formula=linearFormula, coefficient=0))
  }
  return(list(model=M1, dataset=GenedfM, formula=linearFormula, coefficient=coefficients(M1)[marker.index]))
}
fnCreateIsoformsModel <- function(drug, nullModel, isoforms, correlation.threshold=0.51) {
  Isoformsdf <- data.frame(nullModel$dataset, isoforms$FPKM[rownames(nullModel$dataset),])
  colnames(Isoformsdf)[(ncol(nullModel$dataset) + 1):ncol(Isoformsdf)] <- colnames(isoforms$FPKM)
  isoforms.col <- grep("ENST", colnames(Isoformsdf))
  test.cor <- cor(Isoformsdf[,isoforms.col],Isoformsdf[,isoforms.col])
  accepted.isoforms <- colnames(Isoformsdf)[isoforms.col]
  if(!is.na(test.cor))
  {
    if(!is.null(dim(test.cor)))
    {
      correlated.col <- NULL
      for(i in 1:(nrow(test.cor) - 1))
      {
        for(j in (i+1):(ncol(test.cor)))
        {
          if(!is.na(test.cor[i,j]))
          {
            if(test.cor[i,j] > correlation.threshold)
            {
              correlated.col <- c(correlated.col, colnames(test.cor)[j])
            }
          }
        }
      }
      accepted.isoforms <- setdiff(colnames(Isoformsdf)[isoforms.col], correlated.col)
    }
  }
  linear.isoforms.comb <- paste(accepted.isoforms, collapse=" + ")
  if(is.null(tissue))
  {
    linearFormula_Iso <- paste(paste(drug, "~ Tissue + "), linear.isoforms.comb)
  }else{
    linearFormula_Iso <- paste(paste(drug, "~ "), linear.isoforms.comb)
  }
  drug.tissue.col <- colnames(Isoformsdf)[1:(isoforms.col[1]-1)]
  M3 <- glm (formula=linearFormula_Iso, data=Isoformsdf[,c(drug.tissue.col, accepted.isoforms)], family=glm.family)
  
  M3_sum <- summary(M3)$coefficients
  coefficient <- 0
  marker.index <- grep("ENST", rownames(M3_sum))[1]
  if(nrow(M3_sum) >= marker.index)
  {
    isoforms_indices <- marker.index:nrow(M3_sum)
    if(!is.nan(min(M3_sum[isoforms_indices,4])))
    {
      significant <- which(M3_sum[isoforms_indices,4]== min(M3_sum[isoforms_indices,4]))[1] + marker.index - 1
      coefficient <- M3_sum[significant,1]
    }
  }
  return(list(model=M3, dataset=Isoformsdf, formula=linearFormula_Iso, coefficient=coefficient))
}
fnCreateIsoformModel <- function(drug, nullModel, data, isoform) {
  Isoformsdf_temp <- data.frame(nullModel$dataset, data[rownames(nullModel$dataset),isoform])
  colnames(Isoformsdf_temp)[ncol(Isoformsdf_temp)] <- isoform
  if(is.null(tissue))
  {
    linearFormula_Iso_temp <- paste(paste("drug ~ Tissue + "), isoform)
  }else{
    linearFormula_Iso_temp <- paste(paste("drug ~ "), isoform)
  }
  if(length(table(Isoformsdf_temp$Tissue)) == 1 || length(table(Isoformsdf_temp$isoform)) == 1 || length(table(Isoformsdf_temp$drug)) == 1){
    return(list(model=NULL, dataset=Isoformsdf_temp, formula=linearFormula_Iso_temp, coefficient=0, pvalue=1))
  }
  
  switch(model.method, "glm"={
    if(glm.family== "linear")
    {
      M3_temp <- lm (formula=linearFormula_Iso_temp, data=Isoformsdf_temp)
    }else{
      M3_temp <- glm (formula=linearFormula_Iso_temp, data=Isoformsdf_temp, family=glm.family)
    }
    M3_temp_summary <- summary(M3_temp)$coefficients
    coefficient <- 0;
    pvalue <- 1;
    if (isoform %in% rownames(M3_temp_summary))
    {
      coefficient <- M3_temp_summary[isoform,1]
      pvalue <- M3_temp_summary[isoform,4]
    }
  }, "npreg"={
    M3_temp <- npreg(as.formula(linearFormula_Iso_temp), data=Isoformsdf_temp, gradients=TRUE)
    test <- npsigtest(M3_temp, boot.num=100)
    coefficient <- test$bws$bw[2]
    pvalue <- test$P[2]
  })
  if (coefficient == 0){
    M3_temp <- NULL
  }
  
  return(list(model=M3_temp, dataset=Isoformsdf_temp, formula=linearFormula_Iso_temp, coefficient=coefficient, pvalue=pvalue))
}
fnSelectBestIsoform <- function(drug, nullModel, weight, isoforms) {
  min.pvalue <- 1
  best <- ""
  ccle <- NULL; gdsc <- NULL
  for(j in 1:ncol(isoforms$FPKM))
  {
    isoform <- colnames(isoforms$FPKM)[j]
    
    ccle.M3B <- fnCreateIsoformModel(drug, nullModel=nullModel$ccle, data=isoforms$FPKM, isoform)
    gdsc.M3B <- fnCreateIsoformModel(drug, nullModel=nullModel$gdsc, data=isoforms$FPKM, isoform)
    
    if(ccle.M3B$coefficient != 0 && gdsc.M3B$coefficient != 0)
    {
      if(sign(ccle.M3B$coefficient)== sign(gdsc.M3B$coefficient))
      {
        pvalue <- weight$ccle * ccle.M3B$pvalue  + weight$gdsc * gdsc.M3B$pvalue
        if(pvalue <= min.pvalue)
        {
          min.pvalue <- pvalue
          ccle <- list(model=ccle.M3B$model, dataset=ccle.M3B$dataset, formula=ccle.M3B$formula, coefficient=ccle.M3B$coefficient)
          gdsc <- list(model=gdsc.M3B$model, dataset=gdsc.M3B$dataset, formula=gdsc.M3B$formula, coefficient=gdsc.M3B$coefficient)
          best <-  isoform
        }
      }
    }
  }
  return(list(ccle=ccle, gdsc=gdsc, best.isoform=best))
}
fnSelectBestIsoformOneDataSet <- function(drug, nullModel, weight, isoforms, assay) {
  min.pvalue <- 1
  best <- ""
  best.model <- NULL
  for(j in 1:ncol(isoforms$FPKM))
  {
    isoform <- colnames(isoforms$FPKM)[j]
    M3B <- fnCreateIsoformModel(drug, nullModel=nullModel, data=isoforms$FPKM, isoform)
    pvalue <- M3B$pvalue
    if(pvalue <= min.pvalue)
    {
      min.pvalue <- pvalue
      best.model <- list(model=M3B$model, dataset=M3B$dataset, formula=M3B$formula, coefficient=M3B$coefficient)
      best <-  isoform
    }
  }
  return(list(best.model=best.model, best.isoform=best))
}
fnRunbootstrap <- function (models) {
  boot.models <- list()
  for (i in 1:length(models)){
    boot.models[[i]] <- list(data=models[[i]]$dataset, formula=models[[i]]$formula, object=models[[i]]$model)
    names(boot.models)[i] <- names(models)[i]
  }
  M_res <- fnboot(models=boot.models, R=100)
  return(M_res)
}
fnSensitivityCompare <- function (MicroArrayExp, Gene_FPKM, Isoforms, GeneID, sample.no.threshold=5, drugs) {
  Models <- c("M0", "M1", "M2", "M3", "M3B")
  if(is.null(drugs)){drugs <- colnames(ccle.drug.sensitivity)}
  P_values <- list()
  best.isoforms <- character(length=length(drugs))
  statistics <- list()
  if(stat == "r.squared & cindex"){statistics.cindex <- statistics; P_values.cindex <- P_values}
  #layout(matrix(1:24,6,4))
  i <- 0
  for(drug in drugs)
  {
    i <- i + 1
    names(best.isoforms)[i] <- drug
    best.isoforms[i] <- "-"
    
    P_values[[i]] <- matrix(NA,nrow=length(Models), ncol=length(Models), byrow=FALSE)
    rownames(P_values[[i]]) <- Models
    colnames(P_values[[i]]) <- Models
    names(P_values)[i] <- drug
    
    statistics[[i]] <- matrix(0,nrow=5, ncol=length(Models), byrow=FALSE)
    rownames(statistics[[i]]) <- c("mean", "median", "min", "max", "var")
    colnames(statistics[[i]]) <- Models
    names(statistics)[i] <- drug
    
    P_values[[i]][2:nrow(P_values[[i]]),1] <- 0; for(ip in 1:(nrow(P_values[[i]]) -1)){P_values[[i]][ip,(ip+1):ncol(P_values[[i]])] <- 1}
    if(stat == "r.squared & cindex"){
      statistics.cindex[[i]] <- statistics[[i]]
      P_values.cindex[[i]] <- P_values[[i]]
      names(P_values.cindex)[i] <- drug
      names(statistics.cindex)[i] <- drug
    }
    if(!is.null(tissue))
    {
      ccle.sample.no <- 0
      gdsc.sample.no <- 0
      sensitivity.ccle <- ccle.drug.sensitivity[,drug]
      names(sensitivity.ccle) <- rownames(ccle.drug.sensitivity)
      sensitivity.ccle <- sensitivity.ccle[complete.cases(sensitivity.ccle)]
      tissue.types.ccle <- table(tissueTypes[names(sensitivity.ccle),])
      if (tissue %in% names(tissue.types.ccle))
      {
        ccle.sample.no <- tissue.types.ccle[tissue]
      }
      sensitivity.gdsc <- gdsc.drug.sensitivity[,drug]
      names(sensitivity.gdsc) <- rownames(gdsc.drug.sensitivity)
      sensitivity.gdsc <- sensitivity.gdsc[complete.cases(sensitivity.gdsc)]
      tissue.types.gdsc <- table(tissueTypes[names(sensitivity.gdsc),])
      if (tissue %in% names(tissue.types.gdsc))
      {
        gdsc.sample.no <- tissue.types.gdsc[tissue]
      }
    }else{
      ccle.sample.no <- table(complete.cases(ccle.drug.sensitivity[,drug]))["TRUE"]
      gdsc.sample.no <- table(complete.cases(gdsc.drug.sensitivity[,drug]))["TRUE"]
    }
    
    if((ccle.sample.no >= sample.no.threshold) & (gdsc.sample.no >= sample.no.threshold))
    {
      M0 <- list(ccle=fnCreateNullModel(drug=drug, assay="ccle"), gdsc=fnCreateNullModel(drug=drug, assay="gdsc"))
      weight <- list(ccle=M0$ccle$n / (M0$ccle$n + M0$gdsc$n), gdsc=M0$gdsc$n / (M0$ccle$n + M0$gdsc$n))
      
      if(is.null(MicroArrayExp)) {
        M1 <- list(ccle=NULL, gdsc=NULL)
      }else {
        M1 <- list(ccle=fnCreateGeneModel(drug=drug, nullModel=M0$ccle, data=MicroArrayExp), gdsc=fnCreateGeneModel(drug=drug, nullModel=M0$gdsc, data=MicroArrayExp))
        if(!is.null(M1$ccle) && !is.null(M1$gdsc)){if(!is.na(M1$ccle$coefficient) && !is.na(M1$gdsc$coefficient)){if(sign(M1$ccle$coefficient) != sign(M1$gdsc$coefficient)){M1$ccle <- NULL; M1$gdsc <- NULL}}}
      }
      
      if(is.null(Gene_FPKM)) {
        M2 <- list(ccle=NULL, gdsc=NULL)
      }else {
        M2 <- list(ccle=fnCreateGeneModel(drug=drug, nullModel=M0$ccle, data=Gene_FPKM), gdsc=fnCreateGeneModel(drug=drug, nullModel=M0$gdsc, data=Gene_FPKM))
        if(!is.null(M2$ccle) && !is.null(M2$gdsc)){if(!is.na(M2$ccle$coefficient) && !is.na(M2$gdsc$coefficient)){if(sign(M2$ccle$coefficient) != sign(M2$gdsc$coefficient)){M2$ccle <- NULL; M2$gdsc <- NULL}}}
      }
      
      #M3 <- list(ccle=fnCreateIsoformsModel(drug=drug, nullModel=M0$ccle, isoforms=Isoforms), gdsc=fnCreateIsoformsModel(drug=drug, nullModel=M0$gdsc, isoforms=Isoforms))
      #if(!is.null(M3$ccle) && !is.null(M3$gdsc)){if(!is.na(M3$ccle$coefficient) && !is.na(M3$gdsc$coefficient)){if(sign(M3$ccle$coefficient) != sign(M3$gdsc$coefficient)){M3$ccle <- NULL; M3$gdsc <- NULL}}}
      M3 <- list(ccle=NULL, gdsc=NULL)
      
      if(is.null(Isoforms)) {
        M3B <- list(ccle=NULL, gdsc=NULL)
        best.isoforms[i] <- ""
      }else {
        M3B <- fnSelectBestIsoform(drug=drug, nullModel=M0, weight, isoforms=Isoforms)
        best.isoforms[i] <- M3B$best.isoform
      }
      
      #       if(M3B$best.isoform== "")
      #       {
      #         M3B$ccle <- M3$ccle
      #         M3B$gdsc <- M3$gdsc
      #       }
      
      predictors.ccle <- list(M0=M0$ccle, M1=M1$ccle, M2=M2$ccle, M3=M3$ccle, M3B=M3B$ccle)
      predictors.gdsc <- list(M0=M0$gdsc, M1=M1$gdsc, M2=M2$gdsc, M3=M3$gdsc, M3B=M3B$gdsc)
      
      
      switch(statistical.method,
             "bootstrap"={
               
               results <- list(ccle=fnRunbootstrap(models=predictors.ccle),
                               gdsc=fnRunbootstrap(models=predictors.gdsc))
               
               if((stat == "r.squared") || (stat == "adj.r.squared") || (stat == "cindex") || (stat == "r.squared & cindex")){direction="greater"}else{direction="less"}
               if(stat == "r.squared & cindex"){
                 tt <- results
                 results <- list(ccle=results[["ccle"]][["cindex"]], 
                                 gdsc=results[["gdsc"]][["cindex"]])
                 for(j in 1:length(Models))
                 {
                   if(length(results$ccle[[Models[j]]]) > 0 & length(results$gdsc[[Models[j]]]) > 0)
                   {
                     statistics.cindex[[i]]["mean",Models[j]] <- mean(results$ccle[[Models[j]]]) * weight$ccle + mean(results$gdsc[[Models[j]]]) * weight$gdsc
                     statistics.cindex[[i]]["median",Models[j]] <- median(results$ccle[[Models[j]]]) * weight$ccle + median(results$gdsc[[Models[j]]]) * weight$gdsc
                     statistics.cindex[[i]]["min",Models[j]] <- min(results$ccle[[Models[j]]]) * weight$ccle + min(results$gdsc[[Models[j]]]) * weight$gdsc
                     statistics.cindex[[i]]["max",Models[j]] <- max(results$ccle[[Models[j]]]) * weight$ccle + max(results$gdsc[[Models[j]]]) * weight$gdsc
                     statistics.cindex[[i]]["var",Models[j]] <- var(results$ccle[[Models[j]]]) * weight$ccle + var(results$gdsc[[Models[j]]]) * weight$gdsc
                   }
                 }
                 
                 for(j in 1:(length(Models) - 1))
                 {
                   for(k in (j + 1):length(Models))
                   {
                     if((length(results$ccle[[Models[j]]]) > 0) & 
                        (length(results$ccle[[Models[k]]]) > 0) & 
                        (length(results$gdsc[[Models[j]]]) > 0) & 
                        (length(results$gdsc[[Models[k]]]) > 0)){
                       
                       p.value.ccle <- wilcox.test(results$ccle[[Models[k]]] , results$ccle[[Models[j]]], paired=TRUE, alternative=direction)$p.value
                       p.value.gdsc <- wilcox.test(results$gdsc[[Models[k]]] , results$gdsc[[Models[j]]], paired=TRUE, alternative=direction)$p.value
                       P_values.cindex[[i]][Models[j],Models[k]] <- p.value.ccle * weight$ccle + p.value.gdsc * weight$gdsc
                     }
                     else{
                       P_values.cindex[[i]][Models[j],Models[k]] <- 1
                     }
                     
                     if(is.na(P_values.cindex[[i]][Models[j],Models[k]])){P_values.cindex[[i]][Models[j],Models[k]]=1}
                   }
                   ### For Select Based
                   P_values.cindex[[i]][Models[j], "M3B"] <- min(ncol(Isoforms$FPKM) * P_values.cindex[[i]][Models[j], "M3B"], 1)
                 }
                 results <- tt
                 results <- list(ccle=results[["ccle"]][["r.squared"]], 
                                 gdsc=results[["gdsc"]][["r.squared"]])
               }
               for(j in 1:length(Models))
               {
                 if(length(results$ccle[[Models[j]]]) > 0 & length(results$gdsc[[Models[j]]]) > 0)
                 {
                   statistics[[i]]["mean",Models[j]] <- mean(results$ccle[[Models[j]]]) * weight$ccle + mean(results$gdsc[[Models[j]]]) * weight$gdsc
                   statistics[[i]]["median",Models[j]] <- median(results$ccle[[Models[j]]]) * weight$ccle + median(results$gdsc[[Models[j]]]) * weight$gdsc
                   statistics[[i]]["min",Models[j]] <- min(results$ccle[[Models[j]]]) * weight$ccle + min(results$gdsc[[Models[j]]]) * weight$gdsc
                   statistics[[i]]["max",Models[j]] <- max(results$ccle[[Models[j]]]) * weight$ccle + max(results$gdsc[[Models[j]]]) * weight$gdsc
                   statistics[[i]]["var",Models[j]] <- var(results$ccle[[Models[j]]]) * weight$ccle + var(results$gdsc[[Models[j]]]) * weight$gdsc
                 }
               }
               
               for(j in 1:(length(Models) - 1))
               {
                 for(k in (j + 1):length(Models))
                 {
                   if((length(results$ccle[[Models[j]]]) > 0) & 
                      (length(results$ccle[[Models[k]]]) > 0) & 
                      (length(results$gdsc[[Models[j]]]) > 0) & 
                      (length(results$gdsc[[Models[k]]]) > 0)){
                     
                     p.value.ccle <- wilcox.test(results$ccle[[Models[k]]] , results$ccle[[Models[j]]], paired=TRUE, alternative=direction)$p.value
                     p.value.gdsc <- wilcox.test(results$gdsc[[Models[k]]] , results$gdsc[[Models[j]]], paired=TRUE, alternative=direction)$p.value
                     P_values[[i]][Models[j],Models[k]] <- p.value.ccle * weight$ccle + p.value.gdsc * weight$gdsc
                   }
                   else{
                     P_values[[i]][Models[j],Models[k]] <- 1
                   }
                   
                   if(is.na(P_values[[i]][Models[j],Models[k]])){P_values[[i]][Models[j],Models[k]]=1}
                 }
                 ### For Select Based
                 P_values[[i]][Models[j], "M3B"] <- min(ncol(Isoforms$FPKM) * P_values[[i]][Models[j], "M3B"], 1)
               }
             }, "anova"={
               for(j in 1:(length(Models) - 1))
               {
                 for(k in (j + 1):length(Models))
                 {
                   if(!is.null(predictors.ccle[[Models[k]]]$model) & !is.null(predictors.ccle[[Models[j]]]$model)){
                     
                     p.value.ccle <- anova(predictors.ccle[[Models[j]]]$model, predictors.ccle[[Models[k]]]$model, test="Chisq")$"Pr(>Chi)"[2]
                     p.value.gdsc <- anova(predictors.gdsc[[Models[j]]]$model, predictors.gdsc[[Models[k]]]$model, test="Chisq")$"Pr(>Chi)"[2]
                     if(!is.null(p.value.ccle) & !is.null(p.value.gdsc))
                     {
                       P_values[[i]][Models[j],Models[k]] <- p.value.ccle * weight$ccle + p.value.gdsc * weight$gdsc
                     }else{
                       P_values[[i]][Models[j],Models[k]] <- 1
                     }
                   }
                   else{
                     P_values[[i]][Models[j],Models[k]] <- 1
                   }
                   if(is.na(P_values[[i]][Models[j],Models[k]])){P_values[[i]][Models[j],Models[k]] <- 1}
                 }
                 ### For Select Based
                 P_values[[i]][Models[j], "M3B"] <- min(ncol(Isoforms$FPKM) * P_values[[i]][Models[j], "M3B"], 1)
               }
             })
      ### save the coefficients of the models for the M1, M2 and M3 in the forst column of the p_value matrix
      ### I assume that for penalized package there will exist just one isoform after penalization and its coef is reported
      for(j in 2:length(Models))
      {
        P_values[[i]][j, "M0"] <- 0
        if(!is.null(predictors.ccle[[Models[j]]]) & !is.null(predictors.gdsc[[Models[j]]]))
        {
          if(!is.na(predictors.ccle[[Models[j]]]$coefficient) & !is.na(predictors.gdsc[[Models[j]]]$coefficient))
          {
            
            if(sign(predictors.ccle[[Models[j]]]$coefficient) == sign(predictors.gdsc[[Models[j]]]$coefficient)){
              P_values[[i]][j, "M0"] <-  predictors.ccle[[Models[j]]]$coefficient * weight$ccle + predictors.gdsc[[Models[j]]]$coefficient * weight$gdsc
            }
          }
        }
      }
      if(stat == "r.squared & cindex"){
        for(j in 2:length(Models))
        {
          P_values.cindex[[i]][j, "M0"] <- 0
          if(!is.null(predictors.ccle[[Models[j]]]) & !is.null(predictors.gdsc[[Models[j]]]))
          {
            if(!is.na(predictors.ccle[[Models[j]]]$coefficient) & !is.na(predictors.gdsc[[Models[j]]]$coefficient))
            {
              
              if(sign(predictors.ccle[[Models[j]]]$coefficient) == sign(predictors.gdsc[[Models[j]]]$coefficient)){
                P_values.cindex[[i]][j, "M0"] <-  predictors.ccle[[Models[j]]]$coefficient * weight$ccle + predictors.gdsc[[Models[j]]]$coefficient * weight$gdsc
              }
            }
          }
        }
      }
    }
    
  }
  if(stat == "r.squared & cindex"){
    return (list(p.values.r.squared=P_values, p.values.cindex=P_values.cindex, statistics.r.squared=statistics, statistics.cindex=statistics.cindex, best.isoforms=best.isoforms))
  }
  return (list(p.values=P_values, statistics=statistics, best.isoforms=best.isoforms))
}

fnSensitivityOneDataSet <- function (MicroArrayExp, Gene_FPKM, Isoforms, GeneID, assay=c("gray", "gCSI",  "ccle", "gdsc"), sample.no.threshold=5, drugs) {
  Models <- c("M0", "M1", "M2", "M3", "M3B")
  if(is.null(drugs)){drugs <- colnames(ccle.drug.sensitivity)}
  P_values <- list()
  statistics <- list()
  #layout(matrix(1:24,6,4))
  if(assay == "ccle"){Sensitivity <- ccle.drug.sensitivity}else if(assay == "gdsc"){Sensitivity <- gdsc.drug.sensitivity}else if(assay == "gray"){Sensitivity <- gray.drug.sensitivity}else if(assay == "gCSI"){Sensitivity <- gCSI.drug.sensitivity}
  drugs_No <- ncol(Sensitivity)
  best.isoforms <- character(length=drugs_No)
  
  i <- 0
  for(drug in drugs)
  {
    i <- i + 1
    names(best.isoforms)[i] <- drug
    best.isoforms[i] <- "-"
    
    P_values[[i]] <- matrix(NA,nrow=length(Models), ncol=length(Models), byrow=FALSE)
    rownames(P_values[[i]]) <- Models
    colnames(P_values[[i]]) <- Models
    names(P_values)[i] <- drug
    
    statistics[[i]] <- matrix(0,nrow=5, ncol=length(Models), byrow=FALSE)
    rownames(statistics[[i]]) <- c("mean", "median", "min", "max", "var")
    colnames(statistics[[i]]) <- Models
    names(statistics)[i] <- drug
    
    P_values[[i]][2:nrow(P_values[[i]]),1] <- 0; for(ip in 1:(nrow(P_values[[i]]) -1)){P_values[[i]][ip,(ip+1):ncol(P_values[[i]])] <- 1}
    
    if(!is.null(tissue) & assay != "gray")
    {
      sample.no <- 0
      sensitivity.drug <- Sensitivity[,drug]
      names(sensitivity.drug) <- rownames(Sensitivity)
      sensitivity.drug <- sensitivity.drug[complete.cases(sensitivity.drug)]
      tissue.types.ccle <- table(tissueTypes[names(sensitivity.drug),])
      if (tissue %in% names(tissue.types.ccle))
      {
        sample.no <- tissue.types.ccle[tissue]
      }
    }else{
      sample.no <- table(complete.cases(Sensitivity[,drug]))["TRUE"]
    }
    
    
    if(sample.no >= sample.no.threshold)
    {
      M0 <- fnCreateNullModel(drug=drug, assay=assay)
      if(is.null(MicroArrayExp)) {
        M1 <- NULL
      }else{
        M1 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=MicroArrayExp)
      }
      if(is.null(Gene_FPKM)) {
        M2 <- NULL
      }else{
        M2 <- fnCreateGeneModel(drug=drug, nullModel=M0, data=Gene_FPKM)
      }
      #M3 <- fnCreateIsoformsModel(drug=drug, nullModel=M0, isoforms=Isoforms)
      M3 <- NULL
      if(is.null(Isoforms)) {
        M3B <- NULL
        best.isoforms[i] <- ""
      }else{
        M3B <- fnSelectBestIsoformOneDataSet(drug=drug, nullModel=M0, weight, isoforms=Isoforms, assay=assay)
        best.isoforms[i] <- M3B$best.isoform
      }
      
      #        if(M3B$best.isoform == "")
      #        {
      #            M3B <- M3
      #        }
      
      predictors <- list(M0=M0, M1=M1, M2=M2, M3=M3, M3B=M3B$best.model)
      
      
      switch(statistical.method,
             "bootstrap"={
               
               results <- fnRunbootstrap(models=predictors)
               
               for(j in 1:length(Models))
               {
                 if(length(results[[Models[j]]]) > 0)
                 {
                   statistics[[i]]["mean",Models[j]] <- mean(results[[Models[j]]])
                   statistics[[i]]["median",Models[j]] <- median(results[[Models[j]]])
                   statistics[[i]]["min",Models[j]] <- min(results[[Models[j]]])
                   statistics[[i]]["max",Models[j]] <- max(results[[Models[j]]])
                   statistics[[i]]["var",Models[j]] <- var(results[[Models[j]]])
                 }
               }
               
               if((stat== "r.squared") || (stat== "adj.r.squared") || (stat== "cindex")){direction="greater"}else{direction="less"}
               
               for(j in 1:(length(Models) - 1))
               {
                 for(k in (j + 1):length(Models))
                 {
                   if((length(results[[Models[j]]]) > 0) &(length(results[[Models[k]]]) > 0)){
                     
                     p.value <- wilcox.test(results[[Models[k]]] , results[[Models[j]]], paired=TRUE, alternative=direction)$p.value
                     P_values[[i]][Models[j],Models[k]] <- p.value
                   }
                   else{
                     P_values[[i]][Models[j],Models[k]] <- 1
                   }
                   
                   if(is.na(P_values[[i]][Models[j],Models[k]])){P_values[[i]][Models[j],Models[k]] <- 1}
                 }
                 ### For Select Based
                 P_values[[i]][Models[j], "M3B"] <- min(ncol(Isoforms$FPKM) * P_values[[i]][Models[j], "M3B"], 1)
               }
             }, "anova"={
               for(j in 1:(length(Models) - 1))
               {
                 for(k in (j + 1):length(Models))
                 {
                   if(!is.null(predictors[[Models[k]]]$model) & !is.null(predictors[[Models[j]]]$model)){
                     
                     p.value <- anova(predictors[[Models[j]]]$model, predictors[[Models[k]]]$model, test="Chisq")$"Pr(>Chi)"[2]
                     
                     P_values[[i]][Models[j],Models[k]] <- p.value
                   }
                   else{
                     P_values[[i]][Models[j],Models[k]] <- 1
                   }
                   if(is.na(P_values[[i]][Models[j],Models[k]])){P_values[[i]][Models[j],Models[k]] <- 1}
                 }
                 ### For Select Based
                 P_values[[i]][Models[j], "M3B"] <- min(ncol(Isoforms$FPKM) * P_values[[i]][Models[j], "M3B"], 1)
               }
             })
      ### save the coefficients of the models for the M1, M2 and M3 in the forst column of the p_value matrix
      ### I assume that for penalized package there will exist just one isoform after penalization and its coef is reported
      for(j in 2:length(Models))
      {
        P_values[[i]][j, "M0"] <- 0
        if(!is.null(predictors[[Models[j]]]))
        {
          if(!is.na(predictors[[Models[j]]]$coefficient))
          {
            P_values[[i]][j, "M0"] <-  predictors[[Models[j]]]$coefficient
          }
        }
      }
    }
    
  }
  return (list(p.values=P_values, statistics=statistics, best.isoforms=best.isoforms))
}
