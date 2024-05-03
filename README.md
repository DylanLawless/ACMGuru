![Alt Text](./images/DALLE_guru.jpg)

# ACMGuru

[Manual: ACMGuru_0.0.0.9000.pdf](https://github.com/DylanLawless/ACMGuru/blob/main/inst/doc/ACMGuru_0.0.0.9000.pdf)

ACMGuru is an R package designed to facilitate the interpretation of genetic determinants of disease in genomic data. It applies extensive annotation and filtering based on the ACMG (American College of Medical Genetics and Genomics) and AMP (Association for Molecular Pathology) guideline standards. 

Richards et al.'s joint consensus recommendation on the standards and guidelines for the interpretation of sequence variants is a foundational document in this field (2015, [doi: 10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)). ACMGuru is designed to implement the initial steps of our filtering protocol for the addition of ACMG-standardised labels to candidate causal variants, according to these standards.

Accepts VCF format.

## gitignore

Untracked test files in `./output` (automatic user ouput) and in `./data/` due to file sizes:

```
data/samples.tsv
data/dataset_v1_chr21_40411318_41411317.csv
data/dataset_v1_chr21_41411318_42411317.csv
data/dataset_v1_chr21_42411318_43411317.csv
data/dataset_v1_chr21_43411318_44411317.csv
data/dataset_v1_chr21_44411318_45411317.csv
data/dataset_v1_chr21_45411318_46411317.csv
data/dataset_v1_chr21_46411318_47411317.csv
data/dataset_v1_chr21_47411318_48129895.csv
data/uniprot/README.md
data/uniprot/humsavar.txt
data/uniprot/humsvar_corrected.txt
data/uniprot/humsvar_corrected_withID.txt
data/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.gff.zip
data/uniprot/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22.tab.zip
data/uniprot/uniprot_naIDs_convert.tsv
```

## Features

- Automatic annotation and filtering of genomic data based on the ACMG/AMP guideline standards.
- Generation of standardised labels for candidate causal variants.
- Customisation options for filtering and annotation criteria.
- Comprehensive visualization tools for interpreting filtering and annotation results.

## NOTE: THIS PROJECT IS NOT COMPLETE

The following directions are for development and will be updated for users shortly.

## Installation

To install the latest development version of ACMGuru from GitHub:

```r
# install.packages("devtools")
devtools::install_github("your_github_username/ACMGuru")
```

## Usage

```
library(ACMGuru)

# Load your genomic data
genomic_data <- ...

# Apply ACMGuru annotation and filtering
result <- ACMGuru::acmg_filter(genomic_data)

# Visualize results
ACMGuru::plot_result(result)
```

## Documentation

For more detailed information on using ACMGuru, see the package vignette or the reference manual.

## Rule columns

Where available, AMGuru will use annotation features to apply rules. 
The minimal recommended features are typically useful regardless of ACMGuru.

The features listed here are applied using the VEP plugins and their associated databases
(<https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html>)
and/or the dbNSFP plugin and database (SNVs only) 
(<http://database.liulab.science/dbNSFP>).

### Minimal recommended features:
```
sample.id,
SYMBOL,
HGVSp,
HGVSc,
Consequence,
IMPACT,
genotype,
gnomAD_AF,
AC,
AF,
(AF.x, AF.y)*,
AN,
start,
end,
width,
strand,
Allele,
Gene,
Feature_type,
Feature,
BIOTYPE,
EXON,
INTRON,
cDNA_position,
CDS_position,
Protein_position,
Amino_acids,
Codons,
CANONICAL,
```

### Secondary recommended features:
```
ACMG_count,
ACMG_highest,
ACMG_total_score,
comp_het_flag,
Inheritance,
CLIN_SIG,
seqnames,
Existing_variation,
DISTANCE,
STRAND,
FLAGS,
VARIANT_CLASS,
SYMBOL_SOURCE,
HGNC_ID,
CCDS,
ENSP,
SWISSPROT,
TREMBL,
HGVS_OFFSET,
CADD_PHRED
```

### QC recommended features:
```
BaseQRankSum,
DB,
DP,
DS,
END,
ExcessHet,
FS,
InbreedingCoeff,
MLEAC,
MLEAF,
MQ,
MQRankSum,
QD,
RAW_MQ,
ReadPosRankSum,
SOR,
OLD_MULTIALLELIC,
```

### Where these are used
```
#### Apply ACMG PM2 Criterion
gnomAD_AF
gnomad_max

#### Apply ACMG PM3 Criterion
ACMG_PS1
ACMG_PS5
ACMG_PM3
ACMG_PP3
CADD_PHRED
comp_het_flag
MetaLR_pred
MutationAssessor_pred
PolyPhen
REVEL_rankscore
sample
SIFT
SYMBOL

#### Apply ACMG PS1 Criterion
ACMG_PS1
CLIN_SIG

#### Apply ACMG PS5 Criterion
ACMG_PS5
comp_het_flag
IMPACT
sample
SYMBOL

#### Apply ACMG PVS1 Criterion
ACMG_PVS1
comp_het_flag
genotype
IMPACT
Inheritance
sample
SYMBOL

#### Apply Varsome In-silico Predictions for ACMG PP3 Criterion
Engine
pathogenic_type
Moderate_pathogenic_GE
Strong_pathogenic_GE
Supporting_pathogenic_GE

##### Renaming logic: Match engine names with df columns
BayesDel_addAF, BayesDel_addAF_score
BayesDel_noAF, BayesDel_noAF_score
CADD, CADD_PHRED
DANN, DANN_score
EIGEN, Eigen.raw_coding
EIGEN-PC, Eigen.PC.phred_coding
FATHMM, FATHMM_score
FATHMM-MKL, fathmm.MKL_coding_score
FATHMM-XF, fathmm.XF_coding_score
LRT, LRT_score
M-CAP, M.CAP_score
MetaLR, MetaLR_score
MetaSVM, MetaSVM_score
MetaRNN, MetaRNN_score
MutPred, MutPred_score
MutationAssessor, MutationAssessor_score
MutationTaster, MutationTaster_score
phastCons100way_vertebrate, phastCons100way_vertebrate
Polyphen2-HDIV, Polyphen2_HDIV_score
Polyphen2-HVAR, Polyphen2_HVAR_score
PROVEAN, PROVEAN_score
REVEL, REVEL_score
SIFT, SIFT_score
```

## Contributing

I welcome contributions of all types! Please see Contributing Guide for more details.

## Contact

Please report any issues or questions on the GitHub Issues page.

"link_to_vignette", "link_to_reference_manual", "link_to_contributing_guide", "license",
