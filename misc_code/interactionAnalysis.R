
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


load(file.path(path.data, "PSets/CCLE_hs.RData"))
load(file.path(path.data, "PSets/GDSC.RData"))

drugs <- intersect(drugNames(GDSC), drugNames(CCLE))


source(file.path(path.code, "mergePSets.R"))
GDSC <- mergePSets(mDataPSet = CCLE, sensDataPSet=GDSC)



source(file.path(path.code, "geneDrugSensitivity.R"))

source(file.path(path.code, "rankGeneDrugSensitivity.R"))

source(file.path(path.code, "drugSensitivitySig.R"))



# undebug(drugSensitivitySig)


GDSC.NAs <- drugSensitivitySig(GDSC, mDataType = "rnaseq", sensitivity.measure="auc_recomputed", drugs=drugs, nthread=16)
save(GDSC.NAs, file="GDSCmissing.RData")

CCLE.NAs <- drugSensitivitySig(CCLE, mDataType = "rnaseq", sensitivity.measure="auc_recomputed", drugs=drugs, nthread=3)

save(CCLE.NAs, file="CCLEmissing.RData")

load(file.path(path.data, "PSets/gCSI_hs.RData"))

drugs <- intersect(drugs, drugNames(gCSI))

gCSI.NAs <- drugSensitivitySig(gCSI, mDataType = "rnaseq", sensitivity.measure="auc_recomputed", drugs=drugs, nthread=3)
save(gCSI.NAs, file="gCSImissing.RData")

