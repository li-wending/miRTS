# Introduction to miRTS

## Overview
The **`miRTS`** package is designed to infer biological signals originating from different tissue types using circulating microRNA (miRNA) expression data.

The main function, ***`miRTS_score()`***, estimates *miR-TS (miRNA-based Tissue Signal) scores* from bulk extracellular miRNA expression profiles. For each individual, the function returns a tissue-specific score for each modeled tissue type.

*miR-TS score* has been tested to reflect tissue-specific health and disease status across different tissue types in three population cohorts and 11 external validation datasets. The CIBERSORT (default) deconvolution method is recommended for optimal performance.

## Quick start
```{r setup}
# install.packages("remotes")
remotes::install_github("li-wending/miRTS", build_vignettes = TRUE)

# Alternatively, download the latest .tar.gz from the GitHub Releases page and install:
# install.packages("C:/path/to/miRTS_1.0.0.tar.gz", repos = NULL, type = "source")

library(miRTS)
# vignette(topic = "Intro_to_miRTS", package = "miRTS")
# ?miRTS_score
# ?CIBERSORT_download
# miR_TS.output <- miRTS_score(Input_df = example_counts)
```
### Input & Output

- **Input:** A matrix with **miRNAs as rows** and **samples as columns**
- **Output:** A matrix with **samples as rows** and **tissue types as columns** (miR-TS scores)

## Example: Hepatitis C dataset

This example demonstrates how to compute *miR-TS scores* using an included example expression matrix and visualizes liver-associated scores by hepatitis C disease status. The original data are publicly available from GEO under accession [GSE74872](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74872).

All example outputs shown here are generated using data and reference signatures bundled with the package.
Note that in the `example_counts` data, the original dataset was in log2 scale, and has thus been transformed back to linear scale for compatibility with the default count input (see `?example_counts` for details).

### Construct miR-TS scores
```{r hepatitisC-deconvolute}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("This vignette requires ggplot2. Please install it to build the vignette.")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("This vignette requires scales Please install it to build the vignette.")
}

data(example_counts)
data(miRTS_signature_v1)

# Deconvolute miRNA expression data to obtain miR-TS scores.
# CIBERSORT is essential for optimal miR-TS performance (see below):
#
# source('C:/path/to/CIBERSORT.R') # path to your downloaded file

# miR_TS.output <- miRTS_score(
#   Input_df = example_counts,
#   signature_matrix = miRTS_signature_v1,
#   method = "cibersort"  #"cibersort", "xCell2", "MCP-counter"
# )
# all.equal(miR_TS.output , miR_TS.output.bkup)

# Without obtaining the CIBERSORT script,
#  user can alternatively use the saved miR_TS.output.bkup object built-in with the package:
#
miR_TS.output <- miR_TS.output.bkup

```

### Inspect returned objects
```{r hepatitisC-check}
# miR_TS.output is a list containing:
#  (1) a samples × tissues matrix of miR-TS scores,
#  (2) scaled mixture matrix (the `Input_df`),
#  (3) scaled signature matrix (the `signature_matrix`)
str(miR_TS.output)

library(ggplot2)
plt_df <- as.data.frame(miR_TS.output$proportions)
if ("Correlation" %in% colnames(plt_df)){plt_df <- dplyr::select(plt_df, -"P-value",-"Correlation",-"RMSE") }
plt_df <- reshape::melt(plt_df, variable_name = "tissue")
ggplot(plt_df) + aes(tissue, value) + geom_boxplot() +
  xlab("") + ylab("score") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),n.breaks = 5) +
  theme(axis.text.x = element_text(angle = 45,size = 12, hjust=1))

```

### Plot liver scores by disease group
```{r hepatitisC-example}
data(hepatitis_C.meta)

hepatitis_C.output <- Add_metadata(miR_TS.output, hepatitis_C.meta)

hepatitis_C.output$sample_type <- factor(
  hepatitis_C.output$sample_type,
  levels = c(
    "healthy control",
    "Hepatitis C, \nno fibrosis",
    "Hepatitis C, \nwith fibrosis"
  )
)

ggplot(hepatitis_C.output) + aes(sample_type, liver) + geom_boxplot() +
  xlab("") + ylab("Liver score") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2), n.breaks = 5) +
  theme(axis.text = element_text(size = 12))

```

## Other Demo Data

In addition to the above example, miR-TS scores were tested across the following tissue types and health conditions using publicly available datasets:

| Tissue     | Condition                             | GEO Accession |
|------------|----------------------------------------|---------------|
| Adipocyte  | Adipose inflammation                   | [GSE240273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240273) |
| Bone       | Osteoporosis                           | [GSE201543](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201543) |
| Bowel      | Ulcerative colitis                     | [GSE32273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32273)   |
| Brain      | Traumatic brain injuries               | [GSE131695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131695) |
| Heart      | Fulminant myocarditis                  | [GSE148153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148153) |
| Kidney     | T2D + diabetic kidney disease          | [GSE262414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262414) |
| Liver      | Acetaminophen overdose                 | [GSE59565](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59565)   |
| Liver      | Liver allograft rejection              | [GSE69579](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69579)   |
| Lung       | COVID-19                               | [GSE178246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178246) |
| Skin       | Medicamentosa-like dermatitis          | [GSE247297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247297) |



## Note: CIBERSORT prerequisite

The default `method = "cibersort"` requires access to the CIBERSORT source code (**CIBERSORT.R**).

CIBERSORT is licensed but free of charge for non-commercial use only.
Following registration, the **CIBERSORT.R** file is available via the
[cibersort website](https://cibersortx.stanford.edu/).

## Reference
Matsuura, Kentaro et al. “Circulating let-7 levels in plasma and extracellular vesicles correlate with hepatic fibrosis progression in chronic hepatitis C.” Hepatology (Baltimore, Md.) vol. 64,3 (2016): 732-45. doi:10.1002/hep.28660

Keller, Andreas et al., “miRNATissueAtlas2: an update to the human miRNA tissue atlas.” Nucleic acids research vol. 50,D1 (2022): D211-D221. doi:10.1093/nar/gkab808

Li et al., "Circulating extracellular microRNAs as tissue-specific biomarkers of human health and disease" (in revision).
