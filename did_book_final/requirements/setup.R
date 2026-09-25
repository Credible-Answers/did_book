cran   <- c("fixest", "haven", "dplyr", "car")
github <- c(pretrends = "jonathandroth/pretrends")

faltan <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan) > 0) install.packages(faltan)

for (p in names(github)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github(github[[p]])
  }
}

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
# call from Github source("https://raw.githubusercontent.com/Credible-Answers/did_book/Version2/did_book_final/requirements/setup.R")
