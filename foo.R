require(magicaxis) || stop("Library magicaxis is not available!") #add minor tick marks to the plot
require("np") || stop("Library np is not available!")

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
compute.stat <- function(stat, yi, y.hat){
  if(stat == "r.squared"){
    RSS <- sum((yi - y.hat)^2)   
    TSS <- sum((yi - mean(yi))^2)
    return(1-(RSS/TSS))
  }else if(stat == "cindex"){
    return(survcomp::concordance.index(x=-y.hat, surv.time=yi, surv.event=rep(1, length(yi)), na.rm=TRUE, outx=TRUE)[[1]])
  }
}
plot.multi <- function(x,y){
  cl = colors(TRUE)  
  xr <- range(x)
  yr <- range(y)
  plot(y[[1]], xlim = xr, ylim = yr, xlab="Cutoff",ylab = "Statistics")
  for(i in 1:length(y)) {
    lines(x, y[[i]], xlim = xr, ylim = yr, col = cl[i*5 + i])
  }
}
plot.multi.dens <- function(s,t){
  junk.x = NULL
  junk.y = NULL
  for(i in 1:length(s)) {
    junk.x = c(junk.x, density(s[[i]])$x)
    junk.y = c(junk.y, density(s[[i]])$y)
  }
  xr <- range(junk.x)
  yr <- range(junk.y)
  plot(density(s[[1]]), xlim = xr, ylim = yr, xlab="expression", main = t)
  for(i in 1:length(s)) {
    lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
  }
}
fnPlotDensity <- function(Expression,t){  
  L= list()
  for(i in 1:nrow(Expression))
  {
    if(!all(is.na(Expression[i,]))) {
      L[[i]] <- Expression[i, which(!is.na(Expression[i, ]))]
    }
  }
  plot.multi.dens(L,t)
}
myScatterPlot <- function(Name, x, y, method=c("plain", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...){
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  
  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             pdf(file = Name, height=7, width=7)
             plot(x=x, y=y, col=col, pch=pch, ...)
             abline(0, 1, lty=2, col="gray")
             legend("topright", legend=sprintf("r=%.1e", cor(as.numeric(x),as.numeric(y), method="spearman")), cex=1.5, bty="n")
             dev.off()
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
mybarplot <- function(Filename, data, barsNo, groupNo, group.labels, ylab.label, legend.lables, main.label, yaxis = c("Regular", "Log"), cex=1.1){
  #compare all
  barplot.colors= RColorBrewer::brewer.pal(n=7, name="Set1")#c("steelblue3","palegreen3", "mediumpurple3","darkorange3", "indianred3")
  Sp <- c(.75,.25,.25,.25,.25)
  yr <- c(0, max(data,na.rm=TRUE)*1.02)
  idx <- NULL
  k <- 0
  for(i in 1:round(ncol(data)/barsNo))
  {
      idx <- c(idx,k + round(barsNo/2))
      k <- k + barsNo
  }
  
  Groups.labels <- matrix(NA, nrow = barsNo, ncol = groupNo)
  Sp <- NULL
  legend.colors <- NULL
  for(i in 1:barsNo)
  {
    Groups.labels[i,] <- ""
    Sp <- c(Sp,.85)
  }
  Groups.labels[round(barsNo/2+1),] <- group.labels
  Sp[1] <- .15
  
  pdf(file = sprintf("%s.pdf", Filename) , height=7, width=14)
  par(mar=c(9,5,5,8), xpd=T)
  par(oma=c(2,2,2,2))
  temp <- NULL
  for(i in 1:barsNo) {
    temp <- rbind(temp, data[seq(i, length(data), by=barsNo)])
  }
  if(yaxis == "Log") {
    if(0 %in% temp){
      temp <- temp + 1
    }
    #mp <- barplot(log10(temp), beside=TRUE, col=barplot.colors[1:barsNo], ylim=log10(c(0, max(temp)) + 1), ylab=ylab.label, space=Sp, border=NA, axes=FALSE)
    mp <- barplot(log10(temp), beside=TRUE, col=barplot.colors[1:barsNo], ylim=log10(c(0, max(temp)) + 1), ylab=ylab.label, space=c(.15,.85), border=NA, axes=FALSE)
    text(mp, par("usr")[3], labels=Groups.labels, srt=50, adj=c(1.1,1.1), xpd=TRUE, cex=cex)
    magicaxis::magaxis(unlog = 'y', side=2, tcl=-.3)
    legend("topright", inset=c(-.05,0), legend=legend.lables, fill=barplot.colors[1:barsNo], bty="n")
  }else {
    
#    mp <- barplot(temp, beside=TRUE, col=barplot.colors[1:barsNo], ylim=c(0, max(temp)*1.02), ylab=ylab.label, space=Sp, border=NA, axes=FALSE)
    mp <- barplot(temp, beside=TRUE, col=barplot.colors[1:barsNo], ylim=c(0, max(temp)*1.02), ylab=ylab.label, space=c(.15,.85), border=NA, axes=FALSE)
    
    text(mp, par("usr")[3], labels=Groups.labels, srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=1.1)
    magicaxis::magaxis(side=2, tcl=-.3)
    legend("topright", inset=c(-.05,0), legend=legend.lables, fill=barplot.colors[1:barsNo], bty="n")
  }

  dev.off()
}
mybarplot.gap <- function(Filename, data, barsNo, group.labels, ylab.label, legend.lables, main.label, breakpoint){
  #compare all
  barplot.colors= c("steelblue3","palegreen3", "mediumpurple3","darkorange3", "indianred3")
  Sp = c(.75,.25,.25,.25,.25)
  Groups.labels <- matrix("", nrow = 1, ncol = ncol(data))
  Groups.labels[1, seq(from=1,to=ncol(Groups.labels),by=barsNo)] = group.labels
  
  Max = max(data,na.rm=TRUE)
  yr = c(0, Max*1.02)
  if(breakpoint == "No")
  {
    lower = c(0,round(Max/10))
  }
  else
  {
    lower = c(0,breakpoint)
  }
  
  upper = c(round(Max/2),(Max*1.02))
  #gr = c(100, round(max(Sum,na.rm=TRUE)/2))
  y_outer=21
  lowspan=c(0,11)
  topspan=c(lowspan[2]+1,21)
  
  pdf(file = Filename, height=7, width=14)
  
  plot(c(0,1),c(0,y_outer),type='n',axes=FALSE,ylab=ylab.label,xlab='', cex=0.7)
  subplot(barplot(t(data), beside = TRUE, col = barplot.colors[1:barsNo], ylim=lower,xpd=FALSE, space=Sp[1:barsNo], names.arg = Groups.labels,las =3, cex.names=0.7,cex.axis=0.7),x=c(0,1),y=lowspan)
  subplot(barplot(t(data), beside = TRUE, col = barplot.colors[1:barsNo], ylim=upper,xpd=FALSE, space=Sp[1:barsNo], main = main.label,cex.axis=0.7), x=c(0,1),y=topspan)
  
  legend("topright", legend = legend.lables, fill = barplot.colors[1:barsNo], bty="n")
  
  
  dev.off()
}
mystacked.barplot.destiny <- function(Filename, data, barsNo, stackedNo, groupNo, group.labels, ylab.label, legend.lables, main.label){
  #compare all
  stacked.barplot.colors= c("steelblue4","palegreen4", "mediumpurple4","darkorange4", "indianred4")
  #compare rna-seq
  #stacked.barplot.colors= c("royalblue4","red4","seagreen4","slategray4")
  stacked.barplot.density= c(100,40,20,1)
  pdf(file = Filename, height=7, width=14)
  Groups.labels <- matrix(NA, nrow = barsNo, ncol = groupNo)
  Sp <- NULL
  legend.colors <- NULL
  legend.density <- NULL
  for(i in 1:barsNo)
  {
    Groups.labels[i,] = ""
    Sp = c(Sp,.25)
    for(k in 1:stackedNo)
    {
      legend.colors = c(legend.colors, stacked.barplot.colors[i])
      legend.density = c(legend.density, stacked.barplot.density[k])
    }
  }
  Groups.labels[round(barsNo/2+1),] <- group.labels
  Sp[1] = .75
  
  Dest <- NULL
  Sum = 0
  for(i in 1:stackedNo)
  {
    Dest = c(Dest,stacked.barplot.density[i]) 
    Sum = Sum + data[i,]
  }
  Max = max(Sum,na.rm=TRUE)
  yr = c(0, Max*1.02)
  for(i in 1:barsNo)
  {
    start = (groupNo * (i - 1)) + 1
    finish = start + groupNo - 1
    stacked.barplot.matrix <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    stacked.barplot.matrix[,seq(from=i,to=ncol(data),by=barsNo)] <- data[,start:finish]
    C <- NULL
    for(j in 1:stackedNo)
    {
      C = c(C,stacked.barplot.colors[i])
    }
    if(i == 1)
    {
      mp <- barplot(stacked.barplot.matrix, col = C, density = Dest, ylim = yr, ylab = ylab.label, space=Sp, main = main.label)
      text(mp, par("usr")[3], labels = Groups.labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7)
      #axis(1, at=seq(1,length(xx), by=1),labels=FALSE)
      axis(2)
      legend("topleft", legend = legend.lables, fill = legend.colors, density = legend.density, bty="n")
    }else{
      barplot(stacked.barplot.matrix, col = C,density = Dest, ylim = yr, space=Sp, add=TRUE)
    }
  }
  dev.off()
}
mystacked.barplot <- function(Filename, data, barsNo, stackedNo, groupNo, group.labels, ylab.label, legend.lables, main.label, yaxis = c("Regular", "Log"), cex=1.1){ 
  colors.brewer= c("Blues", "Reds","Greens", "Purples")
  colors.legend= c("deepskyblue4", "firebrick3","seagreen4", "purple3")
#  colors.start= c("steelblue1", "mediumpurple1","palegreen1")
#  colors.finish= c("steelblue4", "mediumpurple4","palegreen4")
  pdf(file = Filename, height=7, width=14)
  Groups.labels <- matrix(NA, nrow = barsNo, ncol = groupNo)
  Sp <- NULL
  legend.colors <- NULL
  for(i in 1:barsNo)
  {
    Groups.labels[i,] = ""
    Sp = c(Sp,.25)
  }
  Groups.labels[round(barsNo/2+1),] <- group.labels
  Sp[1] = .75
  
  Sum = 0
  for(i in 1:stackedNo)
  {
    Sum = Sum + data[i,]
  }
  Max = max(Sum,na.rm=TRUE)

  yr = c(0, Max*1.02)
  for(i in 1:barsNo)
  {
    start = (groupNo * (i - 1)) + 1
    finish = start + groupNo - 1
    stacked.barplot.matrix <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    stacked.barplot.matrix[,seq(from=i,to=ncol(data),by=barsNo)] <- data[,start:finish]
    #C <- colorRampPalette(c(colors.start[i],colors.finish[i]))(stackedNo)
    C = colorRampPalette(brewer.pal(n=9, name= colors.brewer[i]))(stackedNo)
    if(i == 1)
    {
      if(yaxis == "Log")
      {
        mp <- barplot(stacked.barplot.matrix, col = C, ylim = yr, ylab = ylab.label, space=Sp, border=NA, axes =FALSE)
        text(mp, par("usr")[3], labels = Groups.labels, srt = 50, adj = c(1.1,1.1), xpd = TRUE, cex=cex)
        #magicaxis::magaxis(side=2,tcl=-.3,majorn=c(5,3),minorn=c(5,2))
        magicaxis::magaxis(unlog = 'y', side = 2)
      }else{
        
        mp <- barplot(stacked.barplot.matrix, col = C, ylim = yr, ylab = ylab.label, space=Sp, main = main.label,border=NA)
        #text(mp, par("usr")[3], labels = Groups.labels, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
        #axis(1, at=seq(1,length(xx), by=1),labels=FALSE)
        axis(2)
      }
    }else{
      barplot(stacked.barplot.matrix, col = C, ylim = yr, space=Sp, add=TRUE,border=NA,axes =FALSE)
    }
  }
  legend("topright", legend = legend.lables, fill = colors.legend,bty="n", cex = 1.2)
  # x = c(round(mp[1], digits = 1) -3 , round(mp[1], digits = 1)-1, round(mp[1], digits = 1)-1 ,round(mp[1], digits = 1) -3)
  #  y = c(round(yr[2], digits = 1) - 1200, round(yr[2], digits = 1) - 1200, round(yr[2], digits = 1)-700, round(yr[2], digits = 1)-700)
  if(yaxis == "Log")
  {
    x = c(round(mp[1], digits = 1) -1.5 , round(mp[1], digits = 1), round(mp[1], digits = 1) ,round(mp[1], digits = 1) -1.5)
    y = c(round(yr[2], digits = 1) - .5, round(yr[2], digits = 1) - .5, round(yr[2], digits = 1)-.8, round(yr[2], digits = 1)-.8)
  }else{
    x = c(round(mp[1], digits = 1) -2 , round(mp[1], digits = 1), round(mp[1], digits = 1) ,round(mp[1], digits = 1) -2)
    y = c(round(yr[2], digits = 1) - 450, round(yr[2], digits = 1) - 450, round(yr[2], digits = 1)-200, round(yr[2], digits = 1)-200)
  }

  legend.gradient(cbind(x, y), cols = colorRampPalette(brewer.pal(n=9, name= 'Greys'))(100), title = "", limits = c("R2<0.55","R2>0.7"),cex=1.2)
  #colorRampPalette(c("gray88","black"))(stackedNo)
  dev.off()

}
mystacked.barplot.simple <- function(Filename, data, main.label, cex=1.3){
  mycol <- c(RColorBrewer::brewer.pal(n=9, name="Set1"), RColorBrewer::brewer.pal(n=5, name="Set3"), "black")
  
  C = mycol[1:nrow(data)]#c(mycol[1],mycol[4],mycol[2])
  #C = c("purple","blue","red")
  
  pdf(file = Filename, height=9, width=14)
  #par(mar=c(8,5,2,2))
  par(mar=c(5.1, 4.1, 4.1, 12.1), xpd=TRUE)
  yr = c(0, 100)
 # mp <- barplot(data, col = C, ylim = yr, ylab = "Percentage", main = main.label,border=NA,axisnames = FALSE, space = 0.5)
  mp <- barplot(data, col = C, ylim = yr, ylab = "Percentage",border=NA,axisnames = FALSE, space = 0.5)
  
  text(mp, par("usr")[3], labels = colnames(data), srt = 50, adj = c(1.1,1.1), xpd = TRUE, cex=cex)
  #axis(1, at=seq(1,length(xx), by=1),labels=FALSE)
  #axis(2)
  legend("topright", inset=c(-0.2,0), legend = rownames(data), fill = C, bty="n")
  dev.off()
  #pdf(file = paste0(Filename,"_legend.pdf"), height=2, width=2)
  #barplot(NA, col = C, ylim = c(0,1) , yaxt = "n", border=NA,axisnames = FALSE, space = 0.5)
  #legend("topright", legend = rownames(data), fill = C, bty="n") 
  #dev.off()
}
mystacked.gap.barplot <- function(Filename, data, barsNo, stackedNo, groupNo, group.labels, ylab.label, legend.lables, main.label, breakpoint){
  stacked.barplot.colors= c("steelblue4","palegreen4", "mediumpurple4","darkorange4", "indianred4")
  stacked.barplot.density= c(100,40,20,1)
  pdf(file = Filename, height=7, width=14)
  
  Groups.labels <- matrix("", nrow = 1, ncol = (groupNo*barsNo))
  Groups.labels[1, seq(from=1,to=ncol(Groups.labels),by=barsNo)] = group.labels
  Sp <- NULL
  legend.colors <- NULL
  legend.density <- NULL
  for(i in 1:barsNo)
  {
    Sp = c(Sp,.25)
    for(k in 1:stackedNo)
    {
      legend.colors = c(legend.colors, stacked.barplot.colors[i])
      legend.density = c(legend.density, stacked.barplot.density[k])
    }
  }
  Sp[1] = .75
  
  Dest <- NULL
  Sum = 0
  for(i in 1:stackedNo)
  {
    Dest = c(Dest,stacked.barplot.density[i]) 
    Sum = Sum + data[i,]
  }
  Max = max(Sum,na.rm=TRUE)
  yr = c(0, Max)
  if(breakpoint == "No")
  {
    lower = c(0,round(Max/10))
  }
  else
  {
    lower = c(0,breakpoint)
  }
  
  upper = c(round(Max/2),(Max*1.02))
  #gr = c(100, round(max(Sum,na.rm=TRUE)/2))
  y_outer=21
  lowspan=c(0,11)
  topspan=c(lowspan[2]+1,21)
  plot(c(0,1),c(0,y_outer),type='n',axes=FALSE,ylab=ylab.label,xlab='', cex=0.7)
  
  for(i in 1:barsNo)
  {
    start = (groupNo * (i - 1)) + 1
    finish = start + groupNo - 1
    stacked.barplot.matrix <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    stacked.barplot.matrix[,seq(from=i,to=ncol(data),by=barsNo)] <- data[,start:finish]
    C <- NULL
    for(j in 1:stackedNo)
    {
      C = c(C,stacked.barplot.colors[i])
    }
    if(i == 1)
    {
      subplot(barplot(stacked.barplot.matrix,col=C, density = Dest,ylim=lower,xpd=FALSE, space=Sp,names.arg = Groups.labels,las =3, cex.names=0.7,cex.axis=0.7),x=c(0,1),y=lowspan)
      subplot(barplot(stacked.barplot.matrix,col=C, density = Dest,ylim=upper,xpd=FALSE, space=Sp, main = main.label,cex.axis=0.7), x=c(0,1),y=topspan)
      
      legend("topright", legend = legend.lables, fill = legend.colors, density = legend.density, bty="n")
    }else{
      subplot(barplot(stacked.barplot.matrix,col=C, density = Dest,ylim=lower,xpd=FALSE, space=Sp,cex.axis=0.7),x=c(0,1),y=lowspan)
      subplot(barplot(stacked.barplot.matrix,col=C, density = Dest,ylim=upper,xpd=FALSE, space=Sp,cex.axis=0.7), x=c(0,1),y=topspan)
    }
  }
  dev.off()
}
cnvrt.coords <-function(x,y=NULL){
  # Stolen from the teachingDemos library, simplified for this use case
  xy <- xy.coords(x,y, recycle=TRUE)
  cusr <- par('usr')
  cplt <- par('plt')	
  plt <- list()
  plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
  plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
  fig <- list()
  fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
  fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
  return( list(fig=fig) )
}
subplot <- function(fun, x, y=NULL){
  # Stolen from the teachingDemos library, simplified for this use case
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  xy <- xy.coords(x,y)
  xy <- cnvrt.coords(xy)$fig
  par(plt=c(xy$x,xy$y), new=TRUE)
  fun
  tmp.par <- par(no.readonly=TRUE)
  return(invisible(tmp.par))
}
fnOrderDrugs<- function(data,filename, ylab, main){
  mean.data = double()
  for(i in 1:ncol(data)){ mean.data[i] = mean(data[complete.cases(data[,i]),i])}
  data.ordered <- data[,order(mean.data)]
  data.order <- order(mean.data)
  
  pdf(file = filename, height=8, width=16)
  #A <- gsub("drugid_","",colnames(data.ordered))
  mp = boxplot(data.ordered, main = main, xaxt="n", ylab = ylab,cex.axis= 0.7, las=2,col = rainbow(ncol(data)))
  text(1:ncol(data), par("usr")[3], labels = gsub("drugid_","",colnames(data.ordered)), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.7)
  dev.off()
  return(list(order = data.order,ordered = data.ordered))
}
fnPie_Isoforms<- function(data, filename){
  Isoforms_No = numeric()
  lbls = character()
  max = 0
  GeneID = 0
  for(i in 1:8)
  {
    Isoforms_No[i] = 0
    lbls[i] = i
  }
  lbls[7] = ">=7"
  lbls[8] = ">=10"
  
  for(i in 1:length(data))
  {
    No = data[i]
    if(No < 7)
    {
      Isoforms_No[No] = Isoforms_No[No] + 1
    }
    else if(No < 10)
    {
      Isoforms_No[7] = Isoforms_No[7] + 1
    }
    else
    {
      Isoforms_No[8] = Isoforms_No[8] + 1
    }
  }
  require("plotrix")
  pdf(file = filename, height=7, width=7)
  pie3D(Isoforms_No, labels = lbls, col=rainbow(length(lbls)), main = "Distribution of genes based on the number of isoforms",theta=pi/4)
  #pie(Isoforms_No, lbls, col=rainbow(length(lbls)), main = "Distribution of Genes based on the number of isoforms")
  dev.off()
}
fncorr_two_approach<- function(drug, app1, app2, Filename, app1.model, app2.model, app1.name, app2.name, drug.name){
  Pvalues.Numerical.app1 = sapply(app1,function(x){x[[drug]]["M0",app1.model]})
  Pvalues.Numerical.app2 = sapply(app2,function(x){x[[drug]]["M0",app2.model]})
  
  
  Pvalues.Numerical.order.app1 = order(Pvalues.Numerical.app1)
  Pvalues.Numerical.order.app2 = order(Pvalues.Numerical.app2)
  overlapped.first100 = length(intersect(Pvalues.Numerical.order.app1[1:100],Pvalues.Numerical.order.app2[1:100]))  
  message(sprintf("The number of overlapped genes between two approched in the set of first 100 genes with least pvalues is %s !",overlapped.first100))

  cor_Sp = cor(Pvalues.Numerical.app1,Pvalues.Numerical.app2,method="spearman")
  message(sprintf("correlation between two approaches is %s !",cor_Sp))
  pdf(file = Filename, height=7, width=7)
  plot( -log10(Pvalues.Numerical.app1),-log10(Pvalues.Numerical.app2),xlab=app1.name, ylab=app2.name, cex.main = .9, main=sprintf("Correlation between %s.%s & %s.%s for %s is %s \n and the Number of overpalled genes in first best 100 is %s", app1.name, app1.model, app2.name, app2.model, drug.name, round(cor_Sp,digits=2), overlapped.first100))
  dev.off()
  
}
fnCorrelationMicroArray_RNASeq <- function(MicroarrayExp, RNA_SeqEXP, path.result){
  pdf(file = file.path(path.result, "Identical_VS_different_cv_spearman.pdf"), height=7, width=7)
  
  MicroarrayExp <- MicroarrayExp[,intersect(colnames(MicroarrayExp),colnames(RNA_SeqEXP))]
  RNA_SeqEXP <- RNA_SeqEXP[,intersect(colnames(MicroarrayExp),colnames(RNA_SeqEXP))]
  cor_Sp = cor(t(MicroarrayExp),t(RNA_SeqEXP),method="spearman")
  boxplot(cbind(diag(cor_Sp), c(cor_Sp[upper.tri(cor_Sp)],cor_Sp[lower.tri(cor_Sp)])), names= c("Identical", "Different"), ylab="Correlation", col = "grey", ylim = c(0.5,1))
  cellline.min.corr = which.max(diag(cor_Sp))
  cellline.max.corr = which.min(diag(cor_Sp))
  #title(sprintf("minimum correlation between two technologies for %s and maximum correlation for %s", rownames(MicroarrayExp)[cellline.min.corr], rownames(MicroarrayExp)[cellline.max.corr]), cex.main=.9)
  dev.off()
  
  ###Transparent
  #myScatterPlot(Name = "Worst_cv.pdf", x=MicroarrayExp[cellline.min.corr,], y=RNA_SeqEXP[cellline.min.corr,], xlab="Microarray", ylab="RNA-seq", main=sprintf("The worst case (%s), Spearman Corr = %f ", rownames(MicroarrayExp)[cellline.min.corr],min(diag(cor_Sp))), pch=16, method="transparent", transparency=0.75
  ###Smooth
  myScatterPlot(Name = file.path(path.result, "Worst_cv.pdf"), x=MicroarrayExp[cellline.min.corr,], y=RNA_SeqEXP[cellline.min.corr,], xlab="Microarray", ylab="RNA-seq", main=sprintf("The worst case (%s), Spearman Corr = %f ", rownames(MicroarrayExp)[cellline.min.corr],min(diag(cor_Sp))), pch=16, method="smooth")
  myScatterPlot(Name = file.path(path.result,"Best_cv.pdf"),x=MicroarrayExp[cellline.max.corr,], y=RNA_SeqEXP[cellline.max.corr,],xlab="Microarray", ylab="RNA-seq", main=sprintf("The Best case (%s), Spearman Corr = %f ", rownames(MicroarrayExp)[cellline.max.corr], max(diag(cor_Sp))), pch=16, method="smooth")
}


