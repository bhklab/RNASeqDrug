load(file.path(path.data, "PSets/CCLE_hs.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/GDSC.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/GRAY_hs.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/UHN_hs.RData"), verbose=TRUE)
load(file.path(path.data, "PSets/gCSI_hs.RData"), verbose=TRUE)
library(gdata)
dd <- cbindX(data.frame(cellNames(CCLE)), 
             data.frame(cellNames(GDSC)) , 
             data.frame(intersect(cellNames(CCLE), cellNames(GDSC))), 
             data.frame(cellNames(gCSI)) , 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(gCSI))), 
             data.frame(cellNames(GRAY)), 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(GRAY))), 
             data.frame(cellNames(UHN)), 
             data.frame(intersectList(cellNames(CCLE), cellNames(GDSC), cellNames(GRAY), cellNames(UHN))))
colnames(dd) <- c("CCLE", "GDSC", "CCLE_GDSC", "gCSI", "CCLE_GDSC_gCSI", "GRAY", "CCLE_GDSC_GRAY", "UHN", "CCLE_GDSC_GRAY_UHN")
xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "cells.tex"), append=FALSE)
write.csv(dd, file=file.path(path.diagrams, "cells.csv"), na="", row.names=FALSE)

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
xtable::print.xtable(xtable::xtable(dd, digits=0), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=file.path(path.diagrams, "drugs.tex"), append=FALSE)
write.csv(dd, file=file.path(path.diagrams, "drugs.csv"), na="", row.names=FALSE)
