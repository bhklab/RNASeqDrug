Gene isoforms as expression-based biomarkers predictive of drug response in vitro
=================================================================


Abstract
--------

Background. Advances in genome-wide molecular profiling and high-throughput drug screening technologies offer a unique opportunity to identify new biomarkers predictive of response to anticancer therapies. If successfully validated in preclinical studies, these biomarkers could be advanced for testing in clinical trials to provide clinicians with tools to select the best therapy for each individual cancer patient, one of the main challenges in precision medicine. In this context, the vast majority of predictive biomarkers for targeted therapies are based on genetic aberrations or protein expressions, as opposed to transcriptomic biomarkers. However, the recent adoption of next-generation sequencing technologies enables accurate profiling of not only gene expression but also alternatively-spliced transcripts in large-scale pharmacogenomic studies. Given that altered mRNA splicing has been shown to be prominent in cancers, linking this feature to drug response will open new avenues of research in biomarker discovery. 
Methods. To identify transcriptomic biomarkers for drug response while addressing the lack of reproducibility of drug sensitivity measurements across studies, we developed a meta-analytical framework combining the pharmacological data generated within the Cancer Cell Line Encyclopedia (CCLE) and the Genomics of Drug Sensitivity in Cancer (GDSC). Our framework consists of fitting predictive models with CCLE RNA-seq data as predictor variables, controlled for tissue type, and combined GDSC and CCLE drug sensitivity values as dependent variables. The accuracy and significance of the fitted models have been assessed using cross-validation, embedding both feature selection and model fitting.
Results. We applied a two step validation model to assess the predictive value of the most promising isoform-based biomarkers in breast cancer. First, an existing independent pharmacogenomic dataset has been exploited to pre-validate the biomarkers that predict response to 8 drugs in a set of 70 breast cancer cell lines. Second, we selected four drugs with the most promising biomarkers (AZD6244, lapatinib, erlotinib and paclitaxel) to test whether their predictive value is robust to change in pharmacological assay. We successfully validated 10 isoform-based biomarkers predictive of drug response in breast cancer, including TNKS1BP1-005 as a biomarker for the EGFR/HER2 tyrosine kinase inhibitor (TKI) Lapatinib and TGFA-001 for the MEK TKI AZD6244. 
Conclusion. The results of our meta-analysis of pharmacogenomic data suggest that isoforms represent a rich resource for biomarkers predictive of response to chemo- and targeted therapies. Our study also showed that the validation rate for this type of biomarkers is low (<50%) for most drugs, supporting the requirements for independent datasets to identify reproducible predictors of response to anticancer drugs.
 

Citation
--------

To cite this work in publication, use

??


Reproducibility of the Analysis Results
--------------------------------------------

1.  Set up the software environment

2.   Download pSets

3.   Run the R scripts


Set up the software environment
-------------------------------

Analysis of this paper has been done in the following session environment. So to reproduce all the paper results the same session should be prepared by installing all the mentioned packages.

#sessionInfo()

R version 3.1.0 Patched (2014-06-08 r65888)
Platform: x86_64-apple-darwin13.1.0 (64-bit)

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] PharmacoGx_0.9.1          RCurl_1.95-4.8            bitops_1.0-6             
 [4] hgu133plus2.db_2.14.0     hgu133plus2frmavecs_1.3.0 hgu133afrmavecs_1.3.0    
 [7] hgu133a.db_2.14.0         hthgu133a.db_2.14.0       hthgu133acdf_2.14.0      
[10] hthgu133afrmavecs_1.1.0   frma_1.16.0               affyio_1.32.0            
[13] affxparser_1.36.0         affy_1.42.3               Hmisc_3.17-2             
[16] ggplot2_2.1.0             Formula_1.2-1             lattice_0.20-29          
[19] R.utils_2.2.0             R.oo_1.20.0               R.methodsS3_1.7.1        
[22] gdata_2.17.0              sva_3.10.0                mgcv_1.7-29              
[25] nlme_3.1-124              corpcor_1.6.8             MetaGx_0.9.10            
[28] lsa_0.73.1                SnowballC_0.5.1           mRMRe_2.0.5              
[31] igraph_1.0.1              genefu_1.14.0             biomaRt_2.20.0           
[34] mclust_5.1                survcomp_1.14.0           prodlim_1.5.7            
[37] survival_2.38-3           jetset_3.1.3              org.Hs.eg.db_2.14.0      
[40] RSQLite_0.11.4            DBI_0.2-7                 AnnotationDbi_1.26.1     
[43] GenomeInfoDb_1.0.2        rcdk_3.3.2                fingerprint_3.5.2        
[46] WriteXLS_4.0.0            magicaxis_1.9.4           sm_2.2-5.4               
[49] plotrix_3.6-1             MASS_7.3-33               RColorBrewer_1.1-2       
[52] downloader_0.4            piano_1.4.2               caTools_1.17.1           
[55] Biobase_2.24.0            BiocGenerics_0.10.0   

Download pSets
-------------------------------

These pSets are required to be accessible for the scripts in a directory named data/PSets:

CCLE_isoforms.RData, GRAY_isoforms.RData, UHN.RData, GDSC.RData

#You can contact the authors to get access to these datasets
benjamin.haibe.kains@utoronto.ca, zhaleh.safikhani@utoronto.ca

Run the R scripts
-------------------------------

#training_phase.R: 

Script is computing signature for the effect of each gene and its isoforms expression on the molecular profile of cell lines for each drug. It returns the estimated coefficient, the the p-values, statistics(mean, median, min, max and variance) of adjusted r squared for the association of each gene and its isoforms to each drug in a dataset named auc_recomputed_drug_association.RData.

*This script takes a long time to be completed. You can ask authors to provide you with the auc_recomputed_drug_association.RData to skip this long step

#training_results.R: 

Script is computing the false discovery rate of the results and report the strongest breast cancer treatment biomarkers for all 15 drugs in common between CCLE and GDSC in a dataset named all.biomarkers.RData.

#Pre_Validation.R: 

Script is validating biomarkers against GRAY dataset and put the validation results in Biomarkers_with_validation_status_2_R2.RData

#Final_Validation.R: 

Script is validating pre validated biomarkers against UHN dataset and pute the final validation results in Biomarkers_uhn_status.RData and creates the main figures of paper.

#Supplementary_Files.R: 

Script is creating supplemntary figures, files and tables

