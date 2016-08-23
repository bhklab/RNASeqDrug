load(file.path(path.data, "PSets/CCLE_isoforms.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/GDSC.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/GRAY_isoforms.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/UHN.RData"), verbose=TRUE)
library(gdata)
dd <- cbindX(data.frame(cellNames(CCLE)), 
             data.frame(cellNames(GDSC)) , 
             data.frame(intersect(cellNames(CCLE), cellNames(GDSC))), 
             data.frame(cellNames(GRAY)), 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(GRAY))), 
             data.frame(cellNames(UHN)), 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(GRAY), cellNames(UHN))))
colnames(dd) <- c("CCLE", "GDSC", "CCLE_GDSC", "GRAY", "CCLE_GDSC_GRAY", "UHN", "CCLE_GDSC_GRAY_UHN")
xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "cells.tex"), append=FALSE)
write.csv(dd, file=file.path(path.diagrams, "cells.csv"), na="", row.names=FALSE)

dd <- cbindX(data.frame(drugNames(CCLE)), 
             data.frame(drugNames(GDSC)) , 
             data.frame(intersect(drugNames(CCLE), drugNames(GDSC))), 
             data.frame(drugNames(GRAY)), 
             data.frame(intersectList(drugNames(CCLE), drugNames(GDSC), drugNames(GRAY))), 
             data.frame(drugNames(UHN)), 
             data.frame(intersectList(drugNames(CCLE), drugNames(GDSC), drugNames(GRAY), drugNames(UHN))))
colnames(dd) <- c("CCLE", "GDSC", "CCLE_GDSC", "GRAY", "CCLE_GDSC_GRAY", "UHN", "CCLE_GDSC_GRAY_UHN")
xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "drugs.tex"), append=FALSE)
write.csv(dd, file=file.path(path.diagrams, "drugs.csv"), na="", row.names=FALSE)
