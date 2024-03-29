# ExosomePurity
A tumour purity deconvolution model to estimate tumour purity in serum exosomes of cancer patients based on miRNA signatures.

## Prerequisites
exosomePurity is written in R(version 4.0.3), requiring R packages DESeq2, preprocessCore, quadprog and edgeR.

Preparing expresssion miRNA-Seq datasets of cancer cell line-derived exosomes, healthy cell-derived exosomes and serum exosomes of cancer patients.


## Usage 

The steps of exosomePurity are as follows:

### Step1
Performing the differential analysiss between cancer cell line-derived exosomes and healthy cell-derived exosomes using miRNA-Seq data.

Applying R script **01signature_DEGS.r** for this step. 
#### Input data
The miRNA-Seq data of cancer cell line-derived exosomes and healthy cell-derived exosomes(rawcount).
#### Output data
Differential miRNAs between cancer cell line-derived exosomes and healthy cell-derived exosomes.

### Step2
Selecting miRNA signatures which are differentially expressed between groups and stably expressed within groups.

Applying R script **02signature_Var_Exp.r** for this step. 
##### Input data
The miRNA-Seq data of cancer cell line-derived exosomes and healthy cell-derived exosomes(CPM).

The Differential gene list obtained from step1.
#### Output data
Signatures between cancer cell line-derived exosomes and healthy cell-derived exosomes.

### Step3
Building tumour purity deconvolution model to quantify the proportions of cancer exosomes and normal exosomes in serum.

Applying R script **03purity_mode.R** for this step. 
#### Input data
The miRNA-Seq data of cancer cell line-derived exosomes and healthy cell-derived exosomes(CPM).

The miRNA-Seq data of serum exosomes of cancer patients(CPM).

The signature list obtained from step2.
#### Output data
The tumor purity of cancer patients' serum exosomes.

### Step4
Utilizing the exosome purity calculated by the model to correct differentially expressed miRNAs.

Applying R script **04gene_correct.R** for this step.
#### Input data
The miRNA-Seq data of serum exosomes of cancer patients(CPM).
The tumor purity of cancer patients' serum exosomes obtained from step3.
#### Output data
Differentially expressed miRNAs corrected by tumour purity.



## Running the tests
### Step1
Performing the differential analysiss between cancer cell line-derived exosomes and healthy cell-derived exosomes using miRNA-Seq data.
Applying R script **01signature_DEGS.r** for this step.
#### Input data
./testdata/01cancer_healthy_cellline_derived_count_exosomes/healthy_cellline_derived_exosomes.csv

./testdata/01cancer_healthy_cellline_derived_count_exosomes/cancer_cellline_derived_exosomes1.csv

./testdata/01cancer_healthy_cellline_derived_count_exosomes/cancer_cellline_derived_exosomes2.csv

./testdata/01cancer_healthy_cellline_derived_count_exosomes/cancer_cellline_derived_exosomes3.csv

#### Output data
./resultdata/01signature_DEGS_data/CRC_0.01_160.Rdata


### Step2
Selecting miRNA signatures which are differentially expressed between groups and stably expressed within groups.
Applying R script **02signature_Var_Exp.r** for this step.
#### Input data
./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/healthy_cellline_derived_cpm_exosomes.csv

./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes1.csv

./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes2.csv

./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes3.csv

./resultdata/01signature_DEGS_data/CRC_0.01_160.Rdata

#### Output data
./resultdata/02signature_Var_data/CRC_var2_h20_48refgenes.Rdata


### Step3
Building tumour purity deconvolution model to quantify the proportions of cancer exosomes and normal exosomes in serum.
Applying R script **03purity_mode.R** for this step.
#### Input data
./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/healthy_cellline_derived_cpm_exosomes.csv

./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes1.csv

./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes2.csv

./testdata/02cancer_healthy_cellline_derived_cpm_exosomes/cancer_cellline_derived_cpm_exosomes3.csv

./testdata/03cancer_serum_cpm_exosomes/cancer_serum_cpm_exosomes.csv

./resultdata/02signature_Var_data/CRC_var2_h20_48refgenes.Rdata

#### Output data
./resultdata/03cancer_serum_purity/serum_putity.csv


### Step4
Utilizing the exosome purity calculated by the model to correct differentially expressed miRNAs.
Applying R script **04gene_correct.R** for this step.
#### Input data
./testdata/04correct_DEG_gene/healthy_cpm_exosomes.csv

./testdata/04correct_DEG_gene/cancer_cpm_exosomes.csv

./testdata/04correct_DEG_gene/01before_crc_0.05_188.Rdata

./resultdata/03cancer_serum_purity/serum_putity.csv

#### Output data
./resultdata/04correct_DEG_gene/01after_crc_p0.05_71.Rdata

./resultdata/04correct_DEG_gene/01new_crc_p0.05_27.Rdata
