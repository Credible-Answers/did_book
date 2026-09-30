# Packages from CRAN
cran <- c(
  "fixest", "haven", "dplyr", "tidyr", "car", "ggplot2",   # data & regressions
  "sandwich", "lmtest",                                     # HC2 SEs (ch07)
  "panelView",                                              # panel view plots
  "TwoWayFEWeights",                                        # ch05-ch08
  "DIDmultiplegtDYN",                                       # ch04, ch06, ch08
  "did",                                                    # Callaway & Sant'Anna (ch06)
  "didimputation",                                          # Borusyak et al. (ch06)
  "fect",                                                   # IFE (ch04)
  "nprobust"                                                # needed by DIDHAD (ch07)
)

# Packages from GitHub
github <- c(
  pretrends = "jonathandroth/pretrends",                    # ch03
  HonestDiD = "asheshrambachan/HonestDiD",                  # ch04
  synthdid  = "synth-inference/synthdid",                   # ch04
  DIDHAD    = "chaisemartinPackages/did_had/R"              # ch07
)

# polars (required by DIDmultiplegtDYN) is not on CRAN
if (!requireNamespace("polars", quietly = TRUE)) {
  install.packages("polars", repos = "https://community.r-multiverse.org")
}

# Install missing CRAN packages
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)

# Install missing GitHub packages
for (p in names(github)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github(github[[p]])
  }
}

# Create output folder for figures
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
# call from GitHub: source("https://raw.githubusercontent.com/Credible-Answers/did_book/Version2/did_book_final/requirements/setup.R")
