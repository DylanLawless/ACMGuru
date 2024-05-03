# Rule columns

## Minimal recommended features:
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

## Secondary recommended features:
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

## QC recommended features:
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


## Where these are used
### Apply ACMG PM2 Criterion
gnomAD_AF
gnomad_max

### Apply ACMG PM3 Criterion
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

### Apply ACMG PS1 Criterion
ACMG_PS1
CLIN_SIG

### Apply ACMG PS5 Criterion
ACMG_PS5
comp_het_flag
IMPACT
sample
SYMBOL

### Apply ACMG PVS1 Criterion
ACMG_PVS1
comp_het_flag
genotype
IMPACT
Inheritance
sample
SYMBOL

### Apply Varsome In-silico Predictions for ACMG PP3 Criterion
Engine
pathogenic_type
Moderate_pathogenic_GE
Strong_pathogenic_GE
Supporting_pathogenic_GE

#### Renaming logic: Match engine names with df columns
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
