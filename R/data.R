#' miR-TS tissue signature matrix,  v1
#'
#' A tissue-specific miRNA signature matrix used to compute miR-TS scores.
#' Rows correspond to miRNAs and columns correspond to tissue types.
#' To achieve best performance, use in conjunction with CIBERSORT deconvolution.
#'
#' @format A numeric matrix with 257 miRNAs (rows) and 23 tissues (columns).
#'   - Row names: miRNA IDs (e.g., "hsa-miR-122-5p")
#'   - Column names: tissue types (e.g., "liver", "heart", "brain", ...)
#'
#' @details
#' The miRTS_signature_v1 object contains the final tissue signature optimized
#' from miRNATissueAtlas2 data (PMID: 34570238) as described in the manuscript.
#' It is used as the default `signature_matrix` in `miRTS_score()`.
#'
#' @source
#' Keller, Andreas et al., “miRNATissueAtlas2: an update to the human miRNA tissue atlas.” Nucleic acids research vol. 50,D1 (2022): D211-D221. doi:10.1093/nar/gkab808
#'
#' Li et al., "Circulating extracellular microRNAs as tissue-specific
#'   biomarkers of human health and disease" (in revision).
#'
#' @examples
#' data(miRTS_signature_v1)
#' dim(miRTS_signature_v1)
#' head(rownames(miRTS_signature_v1))
#' head(colnames(miRTS_signature_v1))
"miRTS_signature_v1"


#' Example extracellular miRNA counts
#'
#' A small example dataset of circulating extracellular miRNA expression
#' measured in human plasma/serum samples. This object is intended for
#' demonstrating the use of the `miRTS` package.
#'
#' @format A numeric matrix (or data frame) with m miRNAs (rows) and n samples
#'   (columns).
#'   - Row names: miRNA IDs (e.g., "hsa-miR-122-5p", "hsa-let-7a-5p").
#'   - Column names: sample IDs (e.g., "sample_01", "sample_02", ...).
#'
#'   The values represent normalized expression / intensity units derived from
#'   microarray data (see manuscript), and are suitable as input to
#'   `miRTS_score()`. They are not raw sequencing counts.
#'
#' @details
#' This dataset is provided to illustrate a typical input to `miRTS_score()`.
#' It can be used together with the `miRTS_signature_v1` object to compute
#' tissue-specific miR-TS scores:
#'
#' \preformatted{
#'   data(example_counts)
#'   data(miRTS_signature_v1)
#'
#'   scores <- miRTS_score(
#'     counts           = example_counts,
#'     signature_matrix = miRTS_signature_v1
#'   )
#' }
#'
#' The original data were derived from PMID: 27227815,
#' which quantified miRNAs using Affymetrix GeneChip miRNA 4.0 Array (Affymetrix, Santa Clara, CA).
#' Note that miR-TS was originally designed for sequencing data,
#' but here microArray data was intentionally chosen to demonstrate that
#' miR-TS can be easily extended to miRNA data generated from other platforms.
#'
#' The authors provided "counts" data in log2-transformed, normalized format:
#' "Raw microarray data (CEL files) were imported into Partek Genomics Suite (Partek Inc., St. Louis, MO), and probe set summaries were computed using Robust Multi-Array algorithm. To adjust for differences in labeling intensities and hybridization, a global normalization was made by aligning signal intensities of data arrays across the medians."
#' The code to re-generate the example_counts from published data:
#' \preformatted{
#' GSE74872 <- read.csv("./Example_data/GSE74872.csv", header = TRUE)
#' GSE74872.name <- read.table("./Example_data/GPL19117-74051.txt", sep = "\t")
#' library(dplyr)
#' GSE74872.name <- GSE74872.name %>% select(V1, V4)
#' GSE74872 <- GSE74872 %>%
#'   inner_join(., GSE74872.name, by="V1")
#' GSE74872 <- GSE74872 %>%
#'  rename(miRNA=V4) %>%
#'  filter(grepl("hsa", miRNA)) %>%
#'  select(-V1)
#' rownames(GSE74872) <- GSE74872$miRNA
#' GSE74872$miRNA <- NULL
#' GSE74872 <- 2^GSE74872 # Convert back to the count scale
#' GSE74872[is.na(GSE74872)] <- 0
#' dim(GSE74872)
#' example_counts <- GSE74872
#' }

#'
#' The samples and values have been anonymized and reduced for the purpose of
#' package examples and do not correspond to identifiable individuals.
#'
#' @source Matsuura, Kentaro et al. “Circulating let-7 levels in plasma and extracellular vesicles correlate with hepatic fibrosis progression in chronic hepatitis C.” Hepatology (Baltimore, Md.) vol. 64,3 (2016): 732-45. doi:10.1002/hep.28660
#'
#' @examples
#' data(example_counts)
#' dim(example_counts)
#' head(rownames(example_counts))
#' head(colnames(example_counts))
"example_counts"


#' Output object saved for testing the miRTS_score function
#'
#' The miR_TS.output.bkup object was archived in this package to
#' enable the testing of the `miRTS_score` function. The
#' output of the miRTS_score(Input_df = example_counts) should be this object.
#'
#' @details
#' This list object includes the following:
#' 1) the "proportions" is the estimated miR_TS score;
#' 2) the "mix" is the z-scaled input mixture data (`Input_df `);
#' 3) the "signature" is the z-scaled signature matrix (`miRTS_signature_v1`)
#' @examples
#' data(miR_TS.output.bkup)
#' dim(miR_TS.output.bkup)
"miR_TS.output.bkup"


#' Hepatitis C example metadata (GSE74872)
#'
#' Sample-level metadata for the example Hepatitis C dataset (GSE74872) used in
#' the `miRTS` package vignettes and tests. This table can be merged with
#' miR-TS output (`miR_TS.output$proportions`) to enable group comparisons
#' by disease/fibrosis status.
#'
#' @format A data frame with one row per sample and columns including:
#' \describe{
#'   \item{Col_names}{Sample names used in the Input_df, which are GEO sample accession ID (e.g., "GSM...").}
#'   \item{sample_type}{Derived grouping variable with
#'   levels "healthy control", "Hepatitis C, no fibrosis", and "Hepatitis C, with fibrosis".}
#' }
#'
#' Additional columns from the source metadata file may be present as from the original publication, such as the Ishak fibrosis score.
#'
#' @details
#' This metadata is intended to be used with helper functions (e.g., `Add_metadata()`)
#' to join phenotype information to miR-TS scores.
#'
#' @source
#' GEO series GSE74872 (metadata file `GSE74872_meta.csv`) included in `Example_data/`.
#'
#' Matsuura, Kentaro et al. “Circulating let-7 levels in plasma and extracellular vesicles correlate with hepatic fibrosis progression in chronic hepatitis C.” Hepatology (Baltimore, Md.) vol. 64,3 (2016): 732-45. doi:10.1002/hep.28660
#'
#' @examples
#' data(hepatitis_C.meta)
#' head(hepatitis_C.meta)
#' table(hepatitis_C.meta$sample_type)
"hepatitis_C.meta"


#' Tissue-enriched miRNA marker map used for MCP-counter-style scoring
#'
#' A two-column mapping of tissue-enriched miRNAs to their enriched tissue/organ
#' labels. This object is used to construct `ts_df` (the marker table) when
#' running `miRTS_score(method = "MCP-counter")`, and is also compatible with
#' `compute_mcp_mirna_scores()`.
#'
#' This mapping relationship is consistent with the tissue-specific miRNAs in the `miRTS_signature_v1`
#'
#' @format A data frame with 2 columns:
#' \describe{
#'   \item{miRNA}{Character. miRNA identifier matching rownames of the input
#'   expression matrix (e.g., `"hsa-miR-122-5p"`).}
#'   \item{enriched_organ}{Character. Tissue/organ label for which the miRNA is
#'   enriched (e.g., `"liver"`, `"heart"`, etc.).}
#' }
#'
#' @details
#' When `method = "MCP-counter"` is selected in `miRTS_score()`, a marker list is
#' constructed by grouping miRNAs by `enriched_organ`. Tissue scores are then
#' computed by aggregating expression across marker miRNAs within each tissue
#' (median by default; see `compute_mcp_mirna_scores()`).
#'
#'
#' @seealso
#' \code{\link{miRTS_score}}, \code{\link{compute_mcp_mirna_scores}}
#'
#' @source
#' Curated tissue-specific miRNA list distributed with the package
#' (file: `Signature_matrix/tissue_specific_miRNAs.csv`).
#'
#' @examples
#' data(ts_miRNA)
#' head(ts_miRNA)
#'
#' # Example: build marker list (one vector of miRNAs per tissue)
#' marker_list <- split(ts_miRNA$miRNA, ts_miRNA$enriched_organ)
#' length(marker_list)
"ts_miRNA"


#' The "xCell2Object" for `miRTS_score()` when applying xCell2
#'
#' This xCell2Object was created by applying the xCell2::xCell2Train() function
#'  to the miRTS_signature_v1 with default settings.
#'
#' @details
#' When `method = "xCell2"` is selected in `miRTS_score()`, this xCell2Object is needed as input.
#'
#'
"xCell2Object.miRTS_signature"


#' 17 tissue types with miRTS signals detectable in plasma/serum
#'
#' @details
#' A vector of characters. The 17 tissue types with miRTS signals detectable in plasma/serum.
#'
#'
"Tissues.blood_detect_17"
