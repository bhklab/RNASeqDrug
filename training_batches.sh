#!/bin/bash

# check to see if all required arguments were provided
if [ $# -eq 15 ]; then
    # assign the provided arguments to variables
    batch_no=$1
    last_gene_index=$2
    sensitivity_method=$3
    resistants_weight=$4
    sensitives_weight=$5
    RNA_seq_normalize=$6
    training_method=$7
    model_method=$8
    glm_family=$9
    tissue=$10
    output_folder=$11
    subset_genes=$12
    subset_genes_method=$13
    jobs_path=$14
    training_input=$15
else
    # assign the default values to variables
    batch_no=200
    last_gene_index=35638
    sensitivity_method="auc_recomputed"
    resistants_weight="1"
    sensitives_weight="1"
    RNA_seq_normalize="TRUE"
    training_method="ccle_gdsc"
    model_method="glm"
    glm_family="gaussian"
    tissue="all"
    output_folder="../data/training_results"
    subset_genes="TRUE"
    subset_genes_method="expression.cut.off"
    jobs_path="../data/jobs"
    training_input="training_ccle_gdsc.RData"
fi

Interval=$((last_gene_index/batch_no))
Interval=$((Interval+1))
batch_no=$((last_gene_index/Interval))
echo $Interval
c=0
j=0
TEMP=$jobs_path$c

f=0

for (( i=1; i<=batch_no+1; i++ ))
do
if (( j > 100)); then
((c++))
TEMP=$jobs_path$c
j=0
fi

if [ ! -d "$TEMP" ]; then
mkdir $TEMP
j=0
fi
((j++))

s=$((f+1))

if ((i <= batch_no));
then
f=$((f+Interval))
else
f=$last_gene_index
fi

SAMPLE="training$i.sh"

echo '#!/bin/bash'>$TEMP/$SAMPLE
cmd="Rscript training.R $s $f $sensitivity_method $resistants_weight $sensitives_weight $RNA_seq_normalize $training_method $model_method $glm_family $tissue $output_folder $subset_genes $subset_genes_method $training_input"
echo $cmd >>$TEMP/$SAMPLE
echo "date" >>$TEMP/$SAMPLE
done

#echo "Running training.R with arguments $firstGene $lastGene $sensitivity_method $resistants_weight $sensitives_weight $RNA_seq_normalize $training_method $model_method $glm_family $tissue $output_folder $subset_genes $subset_genes_method"
#Rscript training.R $firstGene $lastGene $sensitivity_method $resistants_weight $sensitives_weight $RNA_seq_normalize $training_method $model_method $glm_family $tissue $output_folder $subset_genes $subset_genes_method


