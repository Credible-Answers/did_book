# ============================================================
# did_book - setup.R
# Nuvolos / renv / R 4.6
# ============================================================


# ------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------

options(
  repos = c(
    CRAN = "https://cloud.r-project.org"
  )
)

os <- Sys.info()[["sysname"]]

message("Operating system: ", os)
message("R version: ", R.version.string)
message("R platform: ", R.version$platform)


# Use precompiled binaries on Windows/macOS
if (
  .Platform$OS.type == "windows" ||
  os == "Darwin"
) {
  
  options(
    pkgType = "binary",
    install.packages.compile.from.source = "never"
  )
}


# ------------------------------------------------------------
# 2. CRAN PACKAGES
# ------------------------------------------------------------

cran <- c(
  
  # Data and regressions
  "fixest",
  "haven",
  "dplyr",
  "tidyr",
  "car",
  "ggplot2",
  
  # HC2 standard errors
  "sandwich",
  "lmtest",
  
  # Panel plots
  "panelView",
  
  # TWFE weights
  "TwoWayFEWeights",
  
  # DID
  "DIDmultiplegtDYN",
  "did",
  
  # Other DID dependencies
  "zoo",
  "didimputation",
  "fect",
  "nprobust",
  
  # HonestDiD
  "CVXR",
  "ECOSolveR",
  "Rglpk",
  "lpSolveAPI",
  "TruncatedNormal",
  "matrixStats",
  "pracma",
  "latex2exp",
  "mvtnorm",
  "foreach",
  "purrr",
  "tibble"
)


# ------------------------------------------------------------
# 3. GITHUB PACKAGES
# ------------------------------------------------------------

github <- c(
  
  pretrends = "jonathandroth/pretrends",
  
  HonestDiD = "asheshrambachan/HonestDiD",
  
  synthdid = "synth-inference/synthdid",
  
  DIDHAD = "chaisemartinPackages/did_had/R"
)


# ------------------------------------------------------------
# 4. MAKE SURE renv IS AVAILABLE
# ------------------------------------------------------------

if (!requireNamespace("renv", quietly = TRUE)) {
  
  install.packages(
    "renv",
    repos = "https://cloud.r-project.org"
  )
}


# ------------------------------------------------------------
# 5. GLPK
# Required by Rglpk / HonestDiD
# ------------------------------------------------------------

if (!requireNamespace("Rglpk", quietly = TRUE)) {
  
  if (os == "Linux") {
    
    has_glpk <-
      nzchar(Sys.which("glpsol")) ||
      length(
        Sys.glob(
          c(
            "/usr/lib/*/libglpk.so*",
            "/usr/lib/libglpk.so*"
          )
        )
      ) > 0
    
    if (!has_glpk) {
      
      can_sudo <- system(
        "sudo -n true",
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
      ) == 0
      
      if (can_sudo) {
        
        message(
          "Installing GLPK system library..."
        )
        
        system(
          "sudo apt-get update && ",
          "sudo apt-get install -y libglpk-dev"
        )
        
      } else {
        
        message(
          "GLPK system library not found and no sudo access."
        )
        
        message(
          "Administrator command:"
        )
        
        message(
          "sudo apt-get install libglpk-dev"
        )
      }
    }
    
  } else if (os == "Darwin") {
    
    if (
      nzchar(Sys.which("brew")) &&
      !nzchar(Sys.which("glpsol"))
    ) {
      
      message(
        "Installing GLPK via Homebrew..."
      )
      
      system("brew install glpk")
    }
  }
}


# ------------------------------------------------------------
# 6. S7
#
# Required by polars
#
# IMPORTANT:
# Nuvolos uses renv.
# Therefore S7 must be installed through renv.
# ------------------------------------------------------------

message("")
message("============================================================")
message("Installing/checking S7")
message("============================================================")


if (!requireNamespace("S7", quietly = TRUE)) {
  
  message(
    "Installing S7 into the active renv project..."
  )
  
  renv::install(
    "S7"
  )
}


if (!requireNamespace("S7", quietly = TRUE)) {
  
  stop(
    "S7 could not be installed in the active renv project."
  )
  
} else {
  
  message(
    "[OK] S7 ",
    as.character(packageVersion("S7"))
  )
}


# ------------------------------------------------------------
# 7. POLARS
#
# Required by DIDmultiplegtDYN
#
# polars is NOT installed from standard CRAN.
#
# Nuvolos:
#   R 4.6
#   Ubuntu Noble
#   x86_64
#   renv
#
# rpolars.r-universe.dev provides the package.
# The package downloads the precompiled Rust library.
# ------------------------------------------------------------

message("")
message("============================================================")
message("Installing/checking polars")
message("============================================================")


if (!requireNamespace("polars", quietly = TRUE)) {
  
  message(
    "Installing polars..."
  )
  
  # Tell polars that this is not a CRAN build
  Sys.setenv(
    NOT_CRAN = "true"
  )
  
  tryCatch(
    
    {
      
      renv::install(
        "polars",
        repos = c(
          polars = "https://rpolars.r-universe.dev",
          CRAN = "https://cloud.r-project.org"
        ),
        rebuild = TRUE
      )
      
    },
    
    error = function(e) {
      
      message("")
      message(
        "ERROR installing polars:"
      )
      
      message(
        conditionMessage(e)
      )
    }
  )
}


if (requireNamespace("polars", quietly = TRUE)) {
  
  message(
    "[OK] polars ",
    as.character(packageVersion("polars"))
  )
  
} else {
  
  stop(
    "polars could not be installed in the active renv project."
  )
}


# ------------------------------------------------------------
# 8. CRAN PACKAGES
# ------------------------------------------------------------

message("")
message("============================================================")
message("Installing/checking CRAN packages")
message("============================================================")


for (p in cran) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    
    message(
      "Installing ",
      p,
      "..."
    )
    
    tryCatch(
      
      {
        
        install.packages(
          p,
          repos = "https://cloud.r-project.org"
        )
        
      },
      
      error = function(e) {
        
        message(
          "Could not install ",
          p,
          ": ",
          conditionMessage(e)
        )
      }
    )
    
  } else {
    
    message(
      "[OK] ",
      p
    )
  }
}


# ------------------------------------------------------------
# 9. GITHUB PACKAGES
# ------------------------------------------------------------

message("")
message("============================================================")
message("Installing/checking GitHub packages")
message("============================================================")


if (!requireNamespace("remotes", quietly = TRUE)) {
  
  install.packages(
    "remotes",
    repos = "https://cloud.r-project.org"
  )
}


for (p in names(github)) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    
    message(
      "Installing GitHub package: ",
      p
    )
    
    tryCatch(
      
      {
        
        remotes::install_github(
          github[[p]],
          dependencies = TRUE
        )
        
      },
      
      error = function(e) {
        
        message(
          "Could not install ",
          p,
          ": ",
          conditionMessage(e)
        )
      }
    )
    
  } else {
    
    message(
      "[OK] ",
      p
    )
  }
}


# ------------------------------------------------------------
# 10. FINAL CHECK
# ------------------------------------------------------------

all_pkgs <- c(
  "polars",
  cran,
  names(github)
)


installed <- vapply(
  all_pkgs,
  requireNamespace,
  logical(1),
  quietly = TRUE
)


still_missing <- all_pkgs[
  !installed
]


message("")
message("============================================================")
message("FINAL PACKAGE CHECK")
message("============================================================")


if (length(still_missing) == 0) {
  
  message("")
  message(
    "All packages installed successfully."
  )
  
} else {
  
  message("")
  message(
    "These packages could not be installed:"
  )
  
  for (p in still_missing) {
    message(
      "  - ",
      p
    )
  }
}


# ------------------------------------------------------------
# 11. FIGURES DIRECTORY
# ------------------------------------------------------------

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)


message("")
message("============================================================")
message("SETUP FINISHED")
message("============================================================")


# Original GitHub command:
#
# source(
#   "https://raw.githubusercontent.com/Credible-Answers/did_book/Version2/did_book_final/requirements/setup.R"
# )
