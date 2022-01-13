# exosomePurity
A tumour purity deconvolution model to estimate tumour purity in serum exosomes of cancer patients based on miRNA signatures

## Prerequisites
exosomePurity is written in R(version 4.0.3), requiring R packages DESeq2, preprocessCore, quadprog and edgeR.

Preparing expresssion miRNA-Seq datas of cancer cell line-derived exosomes, healthy cell-derived exosomes and serum exosomes of cancer patients.


## Usage

The steps of exosomePurity are as follows:

### Step1
Performing the differential analysiss between cancer cell line-derived exosomes and healthy cell-derived exosomes using miRNA-Seq data.

Applying R script 01signature_DEGS.r for this step. 
#### Input data
The miRNA-Seq data of cancer cell line-derived exosomes and healthy cell-derived exosomes(rawcount).
#### Output data
Differential miRNAs between cancer cell line-derived exosomes and healthy cell-derived exosomes.

##
### Step2
Selecting miRNA signatures which are differentially expressed between groups and stably expressed within groups.
Applying R script 02signature_Var_Exp.r for this step. 
##### Input data
The miRNA-Seq data of cancer cell line-derived exosomes and healthy cell-derived exosomes(CPM).

The Differential gene list obtained from step1.
#### Output data
Signatures between cancer cell line-derived exosomes and healthy cell-derived exosomes.

##
### Step3
Builting tumour purity deconvolution model to quantify the proportions of cancer exosomes and normal exosomes in serum

Applying R script 03mode_0_1.R for this step. 
#### Input data
The miRNA-Seq data of cancer cell line-derived exosomes and healthy cell-derived exosomes(CPM).

The miRNA-Seq data of serum exosomes of cancer patients(CPM).

The signature list obtained from step2.
#### Output data
The tumor purity of cancer patients' serum exosomes.

##
### Step4
Utilizing the exosome purity calculated by the model to correct differentially expressed miRNAs

Applying R script 04gene_correct.R for this step.
#### Input data
The miRNA-Seq data of serum exosomes of cancer patients(CPM).

The tumor purity of cancer patients' serum exosomes obtained from step3.
#### Output data
Differentially expressed miRNAs corrected by tumour purity


## Running the tests

### Step1
Performing the differential analysiss between cancer cell line-derived exosomes and healthy cell-derived exosomes using miRNA-Seq data.

Applying R script 01signature_DEGS.r for this step. 
#### Input data
[cancer_cellline_derived_exosomes.csv]
[healthy_cellline_derived_exosomes.csv]
#### Output data
[CRC_0.01_200.Rdata]

###
### Step2
Selecting miRNA signatures which are differentially expressed between groups and stably expressed within groups.

Applying R script 02signature_Var_Exp.r for this step. 
##### Input data
[cell_line_cpm_miExpr_var.csv]
[normal_cpm_miExpr_var.csv]
[CRC_0.01_200.Rdata]
#### Output data
[CRC_var2_h20_51refgenes.Rdata]

### Step3
Builting tumour purity deconvolution model to quantify the proportions of cancer exosomes and normal exosomes in serum

Applying R script 03mode_0_1.R for this step. 
#### Input data
[cell_line_cpm_miExpr_var.csv]
[normal_cpm_miExpr_var.csv]
[CRC_var2_h20_51refgenes.Rdata]
#### Output data
[normal_purity_0_1.csv]
[celline_purity_0_1.csv]

### Step4
Utilizing the exosome purity calculated by the model to correct differentially expressed miRNAs

Applying R script 04gene_correct.R for this step.
#### Input data
[serum_cancer_cpm.csv]
[serum_normal_cpm.csv]
[crc_tumor_putity.csv]
#### Output data
[crc_correct_113.Rdata]
