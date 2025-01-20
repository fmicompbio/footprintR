# This script is provided as a utility via the swissknife package
# (https://github.com/fmicompbio/swissknife). This script is provided under
# the MIT license, and package authors are permitted to
# include the code as-is in other packages, as long as this note and the
# information provided below crediting the authors of the respective
# functions is retained.

#' Utility function to check validity of scalar variable values.
#'
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria.
#'
#' @param x The variable to be checked.
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x},
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x},
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' @param allowNULL Logical, whether or not \code{NULL} is an acceptable
#'     value for \code{x}.
#'
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom methods is
.assertScalar <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          allowNULL = FALSE) {

    .assertVector(x = x, type = type, rngIncl = rngIncl,
                  rngExcl = rngExcl, validValues = validValues,
                  len = 1, rngLen = NULL, allowNULL = allowNULL)

}

#' Utility function to check validity of vector variable values.
#'
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria.
#'
#' @param x The variable to be checked
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x},
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x},
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' @param len The required length of \code{x}.
#' @param rngLen The allowed range for the length of \code{x}.
#' @param allowNULL Logical, whether or not \code{NULL} is an acceptable
#'     value for \code{x}.
#'
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom methods is
.assertVector <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          len = NULL,
                          rngLen = NULL,
                          allowNULL = FALSE) {
    sc <- sys.calls()
    mycall <- sc[[length(sc)]]
    if (length(sc) >= 2 &&
        identical(as.character(sc[[length(sc) - 1]])[1], ".assertScalar")) {
        mycall <- sc[[length(sc) - 1]]
    }
    args <- lapply(mycall, as.character)[-1]
    xname <- if ("x" %in% names(args)) args$x else "argument"

    ## Check arguments
    stopifnot(is.null(type) || (length(type) == 1L && is.character(type)))
    stopifnot(is.null(rngIncl) || (length(rngIncl) == 2L && is.numeric(rngIncl)))
    stopifnot(is.null(rngExcl) || (length(rngExcl) == 2L && is.numeric(rngExcl)))
    stopifnot(is.null(len) || (length(len) == 1L && is.numeric(len)))
    stopifnot(is.null(rngLen) || (length(rngLen) == 2L && is.numeric(rngLen)))
    stopifnot(is.logical(allowNULL) && length(allowNULL) == 1L)
    if (!is.null(rngIncl) && !is.null(rngExcl)) {
        stop("'rngIncl' and 'rngExcl' can not both be specified")
    }

    ## If there are too many valid values, print only the first 15
    if (length(validValues) > 15) {
        vvPrint <- paste(c(validValues[seq_len(15)],
                           "...(truncated)"),
                         collapse = ", ")
    } else {
        vvPrint <- paste(validValues, collapse = ", ")
    }

    if (is.null(x)) {
        if (allowNULL) {
            return(invisible(TRUE))
        } else {
            stop("'", xname, "' must not be NULL", call. = FALSE)
        }
    }

    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }

    if (!is.null(type) && !is(x, type)) {
        stop("'", xname, "' must be of class '", type, "'", call. = FALSE)
    }

    if (!is.null(rngIncl)) {
        if (!is.null(validValues)) {
            if (any((x < rngIncl[1] | x > rngIncl[2]) & !(x %in% validValues))) {
                stop("'", xname, "' must be within [", rngIncl[1], ",",
                     rngIncl[2], "] (inclusive), or one of: ", vvPrint,
                     call. = FALSE)
            }
        } else {
            if (any(x < rngIncl[1] | x > rngIncl[2])) {
                stop("'", xname, "' must be within [", rngIncl[1], ",",
                     rngIncl[2], "] (inclusive)", call. = FALSE)
            }
        }
    } else if (!is.null(rngExcl)) {
        if (!is.null(validValues)) {
            if (any((x <= rngExcl[1] | x >= rngExcl[2]) & !(x %in% validValues))) {
                stop("'", xname, "' must be within (", rngExcl[1], ",",
                     rngExcl[2], ") (exclusive), or one of: ", vvPrint,
                     call. = FALSE)
            }
        } else {
            if (any(x <= rngExcl[1] | x >= rngExcl[2])) {
                stop("'", xname, "' must be within (", rngExcl[1], ",",
                     rngExcl[2], ") (exclusive)", call. = FALSE)
            }
        }
    } else {
        if (!is.null(validValues) && !all(x %in% validValues)) {
            stop("All values in '", xname, "' must be one of: ", vvPrint,
                 call. = FALSE)
        }
    }


    if (!is.null(len) && length(x) != len) {
        stop("'", xname, "' must have length ", len, call. = FALSE)
    }

    if (!is.null(rngLen) && (length(x) < rngLen[1] || length(x) > rngLen[2])) {
        stop("length of '", xname, "' must be within [", rngLen[1], ",",
             rngLen[2], "] (inclusive)", call. = FALSE)
    }

    return(invisible(TRUE))
}

#' Utility function that makes sure that packages are available
#'
#' The function tries loading the namespaces of the packages given in
#' \code{pkgs}, and throws an exception with an informative error message if
#' that is not the case.
#'
#' @param pkgs Character vector with package names. Can be either just a
#'   package name or a string of the form \code{"githubuser/packagename"} for
#'   packages hosted on GitHub.
#' @param suggestInstallation Logical scalar. If \code{TRUE}, include an
#'   expression to install the missing package(s) as part of the generated
#'   error message.
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @noRd
#' @keywords internal
.assertPackagesAvailable <- function(pkgs, suggestInstallation = TRUE) {
    stopifnot(exprs = {
        is.character(pkgs)
        is.logical(suggestInstallation)
        length(suggestInstallation) == 1L
    })

    avail <- unlist(lapply(sub("^[^/]+/", "", pkgs),
                           function(pkg) {
                               requireNamespace(pkg, quietly = TRUE)
                           }))

    if (any(!avail)) {
        caller <- deparse(sys.calls()[[sys.nframe() - 1]])
        callerfunc <- sub("\\(.+$", "", caller)
        haveBioc <- requireNamespace("BiocManager", quietly = TRUE)
        msg <- paste0("The package", ifelse(sum(!avail) > 1, "s '", " '"),
                      paste(sub("^[^/]+/", "", pkgs[!avail]), collapse = "', '"),
                      "' ",
                      ifelse(sum(!avail) > 1, "are", "is"), " required for ",
                      callerfunc, "(), but not installed.\n")
        if (suggestInstallation) {
            msg <- paste0(msg,
                          "Install ", ifelse(sum(!avail) > 1, "them", "it"), " using:\n",
                          ifelse(haveBioc, "", "install.packages(\"BiocManager\")\n"),
                          "BiocManager::install(c(\"",
                          paste(pkgs[!avail], collapse = "\", \""), "\"))")
        }
        stop(msg, call. = FALSE)
    }

    invisible(TRUE)
}

#' Linearly interpolate NA-value gaps in columns of a NaArray
#'
#' @param assaydat A \code{\link[SparseArray]{NaArray}} object
#'     with read-level footprinting data (positions in rows and reads in
#'     columns).
#' @param pos A numerical vector giving the positions of rows in \code{assaydat}.
#' @param maxgap A numeric scalar giving the maximum number of consecutive
#'     \code{NA}s to fill. Any longer gaps will be left unchanged.
#'
#' @returns A dense matrix with \code{diff(range(pos)) + 1} rows (corresponding
#'     to all positions \code{seq(min(pos), max(pos))}) and \code{ncol(assaydat)}
#'     columns. In each column, runs of \code{NA} values flanked by non-\code{NA}
#'     values have been linearly interpolated.
#'
#' @importFrom BiocGenerics colnames
#' @importFrom SparseArray is_nonna nnavals
#' @importFrom zoo na.approx
#'
#' @noRd
#' @keywords internal
.interpolateColumns <- function(assaydat, pos, maxgap = Inf) {
    stopifnot(exprs = {
        !is.null(colnames(assaydat))
    })
    idx <- is_nonna(assaydat)
    pos_filled <- seq(min(pos), max(pos))
    npos <- length(pos_filled)
    res <- do.call(
        cbind,
        lapply(
            structure(colnames(assaydat), names = colnames(assaydat)),
            function(nm) {
                x <- rep(NA, npos)
                x[match(pos[which(idx[, nm], useNames = FALSE)], pos_filled)] <-
                    nnavals(assaydat[, nm])
                nna <- which(!is.na(x))
                if (length(nna) > 0) {
                    irng <- range(nna)
                    ii <- seq(irng[1], irng[2])
                    x[ii] <- na.approx(object = x, maxgap = maxgap)
                }
                return(x)
            })
    )
    attr(res, "pos") <- pos_filled
    return(res)
}

#' Generate console output messages
#'
#' This is a drop-in replacement for \code{base::message}, which only
#' creates a message if \code{verbose} exists in the calling environment and
#' is set to \code{TRUE}. It also supports inline-markup via the \code{cli}
#' package.
#'
#' @param message The message to be written to the console. It will be
#'     forwarded to \code{\link[cli]{cli_progress_step}} and thus supports
#'     inline markup (see \code{\link[cli]{inline-markup}}).
#' @param noTimer Logical scalar. If \code{FALSE} (the default), the message is
#'     generated using \code{\link[cli]{cli_progress_step}}, which will first
#'     show it with an "info" icon and then again with a "check" icon and timing
#'     information when the next message is generated or the function terminates.
#'     If \code{TRUE}, the message is generated using
#'     \code{\link[cli]{cli_alert_info}}, which means it will be show only once
#'     with an "info" icon.
#' @param ... Additional arguments passed to \code{\link[cli]{cli_progress_step}}
#'     (ignored if \code{noTimer = TRUE}).
#'
#' @importFrom cli cli_progress_step cli_alert_info
#'
#' @noRd
#' @keywords internal
.message <- function(message, noTimer = FALSE, ...) {
    # Try to get 'verbose' from the calling environment
    env <- parent.frame()
    verbose <- tryCatch(get("verbose", envir = env),
                        error = function(e) FALSE)
    if (verbose) {
        if (noTimer) {
            cli_alert_info(text = message, .envir = env)
        } else {
            cli_progress_step(msg = message, .envir = env, ...)
        }
    }
}

#' Convert character region(s) to a \code{GRanges} object
#'
#' This function takes a character vector with one or several strings
#' corresponding to genomic regions and converts them to a
#' \code{\link[GenomicRanges]{GRanges}} object. The supported forms
#' of strings are the same as the ones supported by the \code{htslib} C API
#' that is for example used in \code{samtools}:
#' \describe{
#'      \item{"REF" or "REF:"}{: All of seqname REF}
#'      \item{"REF:START"}{: Seqname REF from START to end of REF}
#'      \item{"REF:-END"}{: Seqname REF from 1 to END}
#'      \item{"REF:START-END"}{: Seqname REF from START to END}
#'      \item{"."}{: All seqnames from 1 to their ends}
#' }
#' Please note that the returned \code{GRanges} is generally parallel
#' to the elements of \code{regions}, with the exception of \code{"."}:
#' This special region is only allowed as a length-one \code{regions}
#' argument and will possibly result several returned regions
#' corresponding to the sequences in \code{seqinfo}.
#'
#' @param regions Character vector with region strings. See details for
#'     supported forms.
#' @param seqinfo One of \code{NULL}, an object for which a
#'     \code{\link[GenomeInfoDb]{seqlengths}} method is available,
#'     such as a \code{BSgenome}, \code{Seqinfo} or \code{SummarizedExperiment}
#'     object or a named numeric vector with genomic sequence names and
#'     lengths. If not \code{NULL}, it will be used to obtain "END" for
#'     \code{regions} that do not specify it. Otherwise, \code{maxend} will
#'     be used for "END".
#' @param maxend Integer scalar giving the maximal value for END in cases
#'     where it is not given in \code{regions} (for example "REF:START")
#'     and was also not provided through other parameters.
#'
#' @author Michael Stadler
#'
#' @importFrom cli cli_abort cli_warn
#' @importFrom GenomeInfoDb seqlengths seqlevels
#' @importFrom GenomicRanges GRanges trim
#' @importFrom IRanges IRanges start end
#'
#' @noRd
#' @keywords internal
.regionStringToGRanges <- function(regions,
                                   seqinfo = NULL,
                                   maxend = .Machine$integer.max) {
    # check arguments
    .assertVector(x = regions, type = "character")
    if (is.null(seqinfo)) {
        reflens <- FALSE
    } else {
        if (is.numeric(seqinfo) && !is.null(names(seqinfo))) {
            reflens <- seqinfo
        } else {
            supportsSeqlengths <- tryCatch({
                seqlengths(seqinfo)
                TRUE
            }, error = function(e) {
                FALSE
            })
            if (supportsSeqlengths) {
                reflens <- seqlengths(seqinfo)
            } else {
                cli_abort(
                    paste0(
                        "`seqinfo` must be `NULL`, an object supporting `seqlengths` or a",
                        " named numeric vector with genomic sequence lengths."))
            }
        }
    }
    .assertScalar(x = maxend, type = "numeric", rngIncl = c(1, .Machine$integer.max))

    # check and convert regions
    if (any(regions == ".")) {
        if (length(regions) != 1L) {
            cli_abort("regions='.' can only be given as a single region")
        } else if (!identical(reflens, FALSE)) {
            gr <- GRanges(seqnames = names(reflens),
                          ranges = IRanges(start = 1, width = unname(reflens)),
                          seqlengths = reflens)
        } else {
            cli_abort(
                paste0("For regions='.' a `seqinfo` argument is required",
                       " that supports `seqlengths(seqinfo)`"))
        }
    } else {
        # regions must be one of:
        # - "REF" or "REF:"
        # - "REF:START"
        # - "REF:-END"
        # - "REF:START-END"
        pat <- "^([^:]+)(:|(:([0-9]+)?-?([0-9]+)?))?$"
        if (any(i <- which(!grepl(pattern = pat, x = regions)))) {
            cli_abort(
                paste0("unrecognized format in {length(i)} region{?s}: '",
                       paste(regions[i[seq.int(min(3, length(i)))]],
                             collapse = "', '"), "'")
            )
        } else {
            df <- strcapture(
                pattern = pat,
                x = regions,
                proto = data.frame(seqnames = character(0),
                                   null1 = character(0),
                                   null2 = character(0),
                                   start = integer(0),
                                   end = integer(0)))
            df$start[is.na(df$start)] <- 1L
            df$end[is.na(df$end)] <- ifelse(df$seqnames[is.na(df$end)] %in% names(reflens),
                                            reflens[df$seqnames[is.na(df$end)]],
                                            rep(maxend, sum(is.na(df$end))))
            ir <- IRanges(start = df$start, end = df$end)
            suppressWarnings({ # avoid out-of-range warning (will trim anyway)
                gr <- GRanges(seqnames = df$seqnames, ranges = ir)
                missingchrs <- setdiff(df$seqnames, names(reflens))
                seqlengths(gr) <- c(reflens,
                                    structure(rep(maxend, length(missingchrs)),
                                              names = missingchrs))[seqlevels(gr)]
                gr <- trim(gr)
            })
            if (any(neq <- start(ir) != start(gr) | end(ir) != end(gr))) {
                cli_warn(
                    paste0("'regions' contained {sum(neq)} out-of-bound ",
                           "range{?s} that were trimmed to the sequence bounds"))
            }
        }
    }

    return(gr)
}

