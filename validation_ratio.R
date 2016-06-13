load("result/auc_recomputed/Biomarkers_uhn_status.RData")

isoforms <-matrix(NA, ncol=2, nrow=length(biomarkers))
colnames(isoforms) <- c("all", "validated")
rownames(isoforms) <- names(biomarkers)
for(i in 1:length(biomarkers)) {
  isoforms[i,] <- c(length(which(biomarkers[[i]][,"type"] == "isoform")),
                length(which(biomarkers[[i]][,"type"] == "isoform" &
                                                        biomarkers[[i]][,"validation.stat"] == "validated")))
}
mycol <- RColorBrewer::brewer.pal(n=8, name="Set2")[c(3,4,6)]

pdf(file = file.path(path.diagrams, "pre_validation_ratio.pdf"), width = 5, height = 5)
par(mar=c(7,3,1,1))
barplot(t(isoforms), col=mycol, las=2)
dev.off()