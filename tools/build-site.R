# Run from the repository root: Rscript tools/build-site.R
# Build from a temporary copy so README assets and Quarto output never alter
# the package source. README.md at the repository root is the canonical home.
build_site <- function() {
  stopifnot(file.exists("README.md"), file.exists("cosTapered/DESCRIPTION"))
  if (!requireNamespace("pkgdown", quietly = TRUE) ||
      utils::packageVersion("pkgdown") < "2.2.1") {
    stop("Install pkgdown >= 2.2.1 before building the site.")
  }

  repo <- normalizePath(".", winslash = "/")
  stage <- tempfile("cosTapered-site-")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  pkg <- file.path(stage, "cosTapered")
  fs::dir_copy(file.path(repo, "cosTapered"), pkg)

  homepage <- readLines(file.path(repo, "README.md"), warn = FALSE)
  homepage <- gsub("https://finleya.github.io/cosTapered/", "", homepage,
                   fixed = TRUE)
  homepage <- gsub("cosTapered/vignettes/cosTapered.pdf",
                   "articles/cosTapered.pdf", homepage, fixed = TRUE)
  homepage <- gsub("README_files/", "man/figures/", homepage, fixed = TRUE)
  writeLines(homepage, file.path(pkg, "README.md"))
  fs::file_copy(fs::dir_ls(file.path(repo, "README_files"), type = "file"),
                file.path(pkg, "man", "figures"), overwrite = TRUE)

  # pkgdown copies these assets into the finished site, including the PDF
  # linked from the model article and the homepage.
  assets <- file.path(pkg, "pkgdown", "assets", "articles")
  dir.create(assets, recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(pkg, "vignettes", "cosTapered.pdf"), assets,
                     overwrite = TRUE))

  destination <- file.path(stage, "docs")
  pkgdown::build_site(
    pkg = pkg,
    override = list(destination = destination),
    preview = FALSE,
    new_process = TRUE,
    install = TRUE
  )
  file.create(file.path(destination, ".nojekyll"))
  pkgdown::check_pkgdown(pkgdown::as_pkgdown(
    pkg, override = list(destination = destination)
  ))
  fs::dir_copy(destination, file.path(repo, "docs"), overwrite = TRUE)
  message("Site built at ", file.path(repo, "docs", "index.html"))
}

build_site()
