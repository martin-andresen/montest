## The version stamp stored in the output of montest(), seqtest() and the seed wrappers.
## Keep this file free of `:=` and of `[` with two index slots: data-raw/reserved_names.R scans every function of
## the package, and new hits would change which user columns montest() aliases.

## `stamp` is filled when the package is loaded, so that installing another version or committing during a session
## does not change what is reported for the code in memory. `code_md5` is added at first use.
.mt_version_state <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  .mt_version_state$stamp <- mtv_read(getNamespaceInfo(pkgname, "path"))
}

#' Which version of montest is running
#'
#' Identifies the code of the loaded package. \code{\link{montest}}, \code{\link{seqtest}},
#' \code{\link{montest_seeds}} and \code{\link{seqtest_seeds}} store it in their output as \code{$version}.
#'
#' The commit is known if the package was installed from a git repository by a tool that records it, such as
#' \code{remotes::install_github()}, or if it was loaded from a git clone with \code{devtools::load_all()}. It is
#' not known after installing from a local folder with \code{devtools::install()} or \code{R CMD INSTALL}.
#' \code{version} and \code{sha} describe the package as it was when loaded.
#'
#' @return A list with
#'   \item{version}{The version number in \code{DESCRIPTION}.}
#'   \item{sha}{The git commit, or \code{NA} if not known.}
#'   \item{dirty}{\code{FALSE} if the code is that of \code{sha}, \code{TRUE} if it was loaded from a git clone
#'     with uncommitted changes to the package files, \code{NA} if not known.}
#'   \item{code_md5}{Checksum of the functions and constants of the loaded package. It is available however the
#'     package was installed: equal values mean equal code. Different versions of R may give different values for
#'     the same code.}
#'
#' @examples
#' montest_version()
#' @export
montest_version <- function() {
  if (is.null(.mt_version_state$code_md5)) .mt_version_state$code_md5 <- mtv_code_md5()
  c(.mt_version_state$stamp, list(code_md5 = .mt_version_state$code_md5))
}

## Version, commit and state of the package in `path`, an installed package or a source tree.
mtv_read <- function(path) {
  desc <- drop(read.dcf(file.path(path, "DESCRIPTION"), fields = c("Version", "RemoteSha")))
  ## RemoteSha is written by the installer. For sources other than git it holds something else than a commit.
  sha <- desc[["RemoteSha"]]
  dirty <- FALSE
  if (!grepl("^[0-9a-f]{40}$", sha)) {
    sha <- NA_character_
    dirty <- NA
    ## `.git` must be in `path` itself: git would otherwise report a repository that merely contains the library.
    if (file.exists(file.path(path, ".git"))) {
      head <- mtv_git(path, c("rev-parse", "HEAD"))
      if (length(head) == 1L) {
        sha <- head
        changed <- mtv_git(path, c("status", "--porcelain", "--", "DESCRIPTION", "NAMESPACE", "R", "src", "data", "inst"))
        if (!is.null(changed)) dirty <- length(changed) > 0L
      }
    }
  }
  list(version = desc[["Version"]], sha = sha, dirty = dirty)
}

## The output of a git command run in `path`, or NULL if git is missing or the command fails.
mtv_git <- function(path, args) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("--no-optional-locks", "-C", shQuote(path), args), stdout = TRUE, stderr = FALSE)),
    error = function(e) NULL)
  if (!is.null(attr(out, "status"))) return(NULL)
  out
}

## Checksum of the text of the code, which unlike an rds file does not record the R version that wrote it.
mtv_code_md5 <- function() {
  f <- tempfile()
  on.exit(unlink(f))
  con <- file(f, "wb")
  writeLines(deparse(mts_code(skip = "^([.]__|[.]packageName$)")), con, useBytes = TRUE)
  close(con)
  unname(tools::md5sum(f))
}
