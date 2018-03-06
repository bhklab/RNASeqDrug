#!/bin/bash

# check to see if all required arguments were provided
if [ $# -eq 26 ]; then
    # assign the provided arguments to variables
    path_diagrams=$1
    effect_size_cut_off=$2
    fdr_cut_off=$3
    tissue=$4
    model_method=$5
    glm_family=$6
    effect_size=$7
    adjustment_method=$8
    sensitivity_method=$9
    RNA_seq_normalize=$10
    training_input=$11
    training_phase=$12
    training_input_all_type=$13
    training_phase_all_type=$14
    fdr_cut_off_training=$15
    subset_genes=$16
    subset_genes_method=$17
    training=$18
    taining_all_type=$19
    training_analysis=$20
    training_all_type_analysis=$21
    pre_validation_gcsi=$22
    pre_validation_gray=$23
    final_validation=$24
    figures=$25
    parallel=$26
else
    # assign the default values to variables
    path_diagrams="auc_recomputed_ccle_gdsc"
    effect_size_cut_off="0.55"
    fdr_cut_off="0.05"
    tissue="breast"
    model_method="glm"
    glm_family="gaussian"
    effect_size="cindex"
    adjustment_method="fdr"
    sensitivity_method="auc_recomputed"
    RNA_seq_normalize="TRUE"
    training_input="training_ccle_gdsc.RData"
    training_phase="auc_recomputed_drug_association.RData"
    training_input_all_type="training_ccle_gdsc_mut_cnv.RData"
    training_phase_all_type="auc_recomputed_drug_association_mut_cnv.RData"
    fdr_cut_off_training="0.01"
    subset_genes="TRUE"
    subset_genes_method="expression.cut.off"
    training=false
    taining_all_type=false
    training_analysis=true
    training_all_type_analysis=true
    pre_validation_gcsi=true
    pre_validation_gray=true
    final_validation=true
    figures=true
    parallel=false
fi
##TRAINING PHASE##
if $training && $parallel; then
echo "Running training_batches.sh with arguments 100 35638 $sensitivity_method 1 1 $RNA_seq_normalize ccle_gdsc $model_method $glm_family all ../data/training_results $subset_genes $subset_genes_method ../data/jobs $training_input"
training_batches.sh "100" "35638" $sensitivity_method "1" "1" $RNA_seq_normalize "ccle_gdsc" $model_method $glm_family "all" "../data/training_results" $subset_genes $subset_genes_method "../data/jobs" $training_input

for i in $(find -L "../data/jobs/0" -name '*.sh')
do
    #if there is cluster system available to do the training in parallel
    #qsub $i
    $i
done

echo "Running merge_training_batches.R ../data/training_results ../data/$training_input ../data/$training_phase"
Rscript merge_training_batches.R ../data/training_results ../data/$training_input ../data/$training_phase
fi
##TRAINING PHASE for genes with molecular data available for mutation, cnv and expression##
if $taining_all_type; then
echo "Running training.R with arguments 1 1550 $sensitivity_method 1 1 $RNA_seq_normalize $training_method $model_method $glm_family $tissue ../data/training_mut_cnv $subset_genes $subset_genes_method $training_input_all_type"
Rscript training.R 1 1550 $sensitivity_method 1 1 $RNA_seq_normalize $training_method $model_method $glm_family $tissue ../data/training_mut_cnv $subset_genes $subset_genes_method $training_input_all_type
echo "Running merge_training_batches.R ../data/training_mut_cnv ../data/$training_input_all_type ../data/$training_phase_all_type"
Rscript merge_training_batches.R ../data/training_mut_cnv ../data/$training_input_all_type ../data/$training_phase_all_type
fi

##ANALYIS TRAINING RESULTS##
if $training_analysis; then
echo "Running training_results.R with args $training_phase $effect_size_cut_off $fdr_cut_off_training $tissue $model_method $glm_family $effect_size $adjustment_method $path_diagrams $training_input"
Rscript training_results.R $training_phase $effect_size_cut_off $fdr_cut_off_training $tissue $model_method $glm_family $effect_size $adjustment_method $path_diagrams $training_input
fi

##ANALYIS TRAINING RESULTS FOR 1550 genes with ALL TYPES##
if $training_all_type_analysis; then
echo "Running training_results.R with args $training_phase_all_type $effect_size_cut_off $fdr_cut_off_training $tissue $model_method $glm_family $effect_size $adjustment_method ${path_diagrams}_all_types $training_input_all_type"
Rscript training_results.R $training_phase_all_type $effect_size_cut_off $fdr_cut_off_training $tissue $model_method $glm_family $effect_size $adjustment_method ${path_diagrams}_all_types $training_input_all_type
fi

##PRE VALIDATION ON gCSI##
if [$pre_validation_gcsi eq "TRUE"]; then
echo "Running Pre_Validation_gCSI.R with args $path_diagrams $effect_size_cut_off $fdr_cut_off $model_method $glm_family $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input"
Rscript Pre_Validation_gCSI.R $path_diagrams $effect_size_cut_off $fdr_cut_off $model_method $glm_family $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input
fi

##PRE VALIDATION ON GRAY##
if $pre_validation_gray; then
echo "Running Pre_Validation.R with args $path_diagrams $effect_size_cut_off $fdr_cut_off $tissue $model_method $glm_family $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input"
Rscript Pre_Validation.R $path_diagrams $effect_size_cut_off $fdr_cut_off $tissue $model_method $glm_family $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input
fi

##FINAL VALIDATION ON UHN##
if $final_validation; then
echo "Running Final_Validation.R with args $path_diagrams $effect_size_cut_off $fdr_cut_off $tissue $model_method $glm_family $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input"
Rscript Final_Validation.R $path_diagrams $effect_size_cut_off $fdr_cut_off $tissue $model_method $glm_family $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input
fi

##CREATE SUPLEMENTAL FIGURES AND TABLES##
if $figures; then
echo "Running Figures_Tables.R with args $path_diagrams $training_phase $training_phase_all_type $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input"
Rscript Figures_Tables.R $path_diagrams $training_phase $training_phase_all_type $effect_size $adjustment_method $sensitivity_method $RNA_seq_normalize $training_input
fi


