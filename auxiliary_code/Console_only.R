library(devtools)

use_r("miRTS_score")

use_vignette("Intro_to_miRTS", title = "Introduction to miRTS")

# I need to load the data needed for this package:
load("../../Rpackage/miRTS/Data_file_for_miRTS_package.RData")
usethis::use_data(miRTS_signature_v1,overwrite = TRUE)
usethis::use_data(miR_TS.output.bkup,overwrite = TRUE)
usethis::use_data(example_counts,overwrite = TRUE)
usethis::use_data(Tissues.blood_detect_17,overwrite = TRUE)
usethis::use_data(xCell2Object.miRTS_signature,overwrite = TRUE)
hepatitis_C.meta <- hepatitis_C.meta[, c("Col_names", "sample_type"), drop = FALSE]
usethis::use_data(hepatitis_C.meta, overwrite = TRUE)
usethis::use_data(ts_miRNA,overwrite = TRUE)

dim(miRTS::miR_TS.output.bkup$proportions)

#
usethis::use_package("e1071")
usethis::use_package("edgeR")
usethis::use_package("parallel")
usethis::use_package("preprocessCore")

# export the vignette:

dir.create("inst/doc", recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  input       = "vignettes/Intro_to_miRTS.Rmd",
  output_dir  = "inst/doc",
  clean       = TRUE
)

# build the package!!!
pkgbuild::build(clean_doc = FALSE)
usethis::use_mit_license("Wending Li; Haotian Wu")

devtools::check()

devtools::build_vignettes()

# Windows (CRAN-required-ish sanity):
devtools::check_win_devel()

# Initialize git + GitHub repo:
usethis::use_git() # need to commit
gitcreds::gitcreds_set()
# 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
usethis::use_github()
usethis::use_git_config(
  user.name  = "li-wending",
  user.email = "wdli.wending@gmail.com"
)
# Add GitHub Actions CI for 'R CMD check' on major OSes (good CRAN prep):
usethis::use_github_action_check_standard()


# Submit to CRAN
#
# Use the “release” helper (it walks you through common CRAN checks):
  # Or manually: build tar.gz, then submit via CRAN web form.
  devtools::release()
