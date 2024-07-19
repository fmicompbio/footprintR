#' Call `modkit extract` to extract read-level base modifications from modBAM files
#'
#' @description
#' This function is a wrapper around the `modkit extract` sub-command to extract
#' read-level base modification information from modBAM files into tab-separated
#' values table(s).
#' For more information on available `modkit extract` arguments and output
#' tables specification see https://nanoporetech.github.io/modkit/intro_extract.html
#'
#' @param modkit_bin Character scalar specifying the path to the \code{modkit}
#'     binary.
#' @param bamfile Character scalar specifying the path to a \code{modBAM}
#'     file. An indexed BAM file can significantly speed up certain extract
#'     operations.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object specifying which
#'     genomic regions to extract the reads from. Note that the reads are not
#'     trimmed to the boundaries of the specified ranges. As a result, returned
#'     positions will typically extend out of the specified regions.
#' @param num_reads Number of reads to extract per specified genomic region.
#'     When \code{N} genomic ranges are specified, the total number of extracted
#'     reads will be at most \code{num_reads * N}.
#' @param out_extract_table Character scalar specifying the path for the
#'     `extract table` output. Can be \code{NULL} if only a `read-calls`
#'     table is needed.
#' @param out_read_calls Character scalar specifying the path for the
#'     `read-calls table` output. Can be \code{NULL} if only an `extract table`
#'     is needed.
#' @param out_log_file Character scalar specifying the path for the the command
#'     run log output. Can be \code{NULL} if no log file is needed (not
#'     recommended).
#' @param modkit_args Character vector with additional `modkit extract`
#'     arguments. Please refer to the `modkit extract` documentation for a
#'     complete list of possible arguments.
#' @param tempdir_base Character scalar specifying the path to create the
#'     `modkit_temp` directory. This temporary directory is only created when
#'     multiple genomic regions are specified.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A character vector with the elements "extract-table", "read-calls"
#'     and "run-log", specifying the paths to the generated table, call and log
#'     files, respectively. An \code{NA} value indicates that a given output
#'     was not generated.
#'
#' @author Panagiotis Papapasaikas
#'
#' @examples
#' \dontrun{
#' GR <- as(c("chr1:12340678-12345678", "chr2:12340678-12345678"), "GRanges")
#'
#' further_args <- c(
#'    '-t 8', '--force',
#'    '--mapped-only', '--edge-filter 100',
#'    '-p  0.2',
#'    '--mod-threshold m:0.8',
#'    '--mod-threshold a:0.5'
#' )
#'
#' # Produce both `extract table` and `read calls` files for multiple GRanges
#' modkitExtract(modkit = modkit_bin_PATH, bamfile = BAMF, num_reads = 10,
#'               regions = GR[1:2], out_extract_table = "test.etbl",
#'               out_read_calls = "test.rdcl", modkit_args = further_args)
#' # Produce only `extract table` file
#' modkitExtract(modkit = modkit_bin_PATH, bamfile = BAMF, num_reads = 10,
#'               regions = GR[1], out_extract_table = "test.etbl",
#'               out_read_calls = NULL, modkit_args = further_args)
#' # Produce only `read calls` file
#' modkitExtract(modkit = modkit_bin_PATH, bamfile = BAMF, num_reads = 10,
#'               regions = GR[1], out_extract_table = NULL,
#'               out_read_calls = "test.rdcl", modkit_args = further_args)
#' }
#'
#' @seealso [`modkit` software](https://nanoporetech.github.io/modkit),
#'     [`modkit extract` documentation and tabulated output formats specification](https://nanoporetech.github.io/modkit/intro_extract.html),
#'     \code{\link[GenomicRanges]{GRanges}} for the object used to specify
#'     genomic regions.
#'
#' @import GenomicRanges
#'
#' @export
modkitExtract <- function(modkit_bin,
                          bamfile,
                          regions = NULL,
                          num_reads = NULL,
                          out_extract_table = NULL,
                          out_read_calls = NULL,
                          out_log_file = NULL,
                          modkit_args = NULL,
                          tempdir_base = tempdir(),
                          verbose = TRUE) {

    # digest arguments
    # --------------------------------------------------------------------------
    .assertScalar(x = modkit_bin, type = "character")
    modkitFAIL <- suppressWarnings(
        system(paste0(modkit_bin, " --version"), intern = FALSE,
               ignore.stdout = TRUE, ignore.stderr = TRUE))
    if (modkitFAIL != 0) {
        stop("A valid path to a modkit executable  has not been provided")
    } else {
        if (verbose) {
            message("Using ", system(paste0(modkit_bin, " --version"), intern=TRUE))
        }
    }

    .assertScalar(x = bamfile, type = "character")
    if (!file.exists(bamfile)) {
        stop("BAM file not found at: ", normalizePath(bamfile, mustWork = FALSE))
    }

    .assertVector(x = regions, type = "GRanges", allowNULL = TRUE)

    .assertScalar(x = num_reads, type = "numeric", allowNULL = TRUE)

    .assertScalar(x = out_read_calls, type = "character", allowNULL = TRUE)

    .assertScalar(x = out_extract_table, type = "character", allowNULL = TRUE)

    .assertScalar(x = out_log_file, type = "character", allowNULL = TRUE)

    .assertScalar(x = tempdir_base, type = "character")

    .assertScalar(x = verbose, type = "logical")

    # Convert GRanges to `chr:start-end` format
    if (!is.null(regions)){
        regions_char <- as.character(reduce(GRanges(regions, strand="*")))
    }


    # Generate character vector of arguments to be passed to modkit.
    # --------------------------------------------------------------------------
    # first argument should always be the modkit sub-command (here: 'extract')
    pass_ARGS <- c('extract', '--suppress-progress')

    # Prepare --num-reads argument
    if (!is.null(num_reads)) {
        pass_ARGS <- c(pass_ARGS, paste("--num-reads", num_reads, sep=" "))
        if (is.null(regions)) {
            # in absence of a specified region, --num-reads is only respected
            #    if bam index is ignored
            pass_ARGS <- c(pass_ARGS, "--ignore-index")
        }
    }

    # Prepare --log-filepath argument
    if(!is.null(out_log_file)){
        pass_ARGS <- c(pass_ARGS, paste("--log-filepath",out_log_file, sep=" ") )
        print(c( "Specified path to run log:", normalizePath(out_log_file, mustWork=FALSE) ))
        if (file.exists(out_log_file)){
            print("Warning: Specified `out_log_file` already exists. The log will be appened to the existing file")
        }
    }

    # combine function exposed args with user-passed modkit args:
    pass_ARGS <- c(pass_ARGS, modkit_args)

    # Prepare separately --read-calls-path argument
    if (!is.null(out_read_calls)) {
        pass_out_read_calls <- paste('--read-calls-path', out_read_calls,
                                     sep = " ")
        if (verbose) {
            message("Specified path to read-calls table: ",
                    normalizePath(out_read_calls, mustWork = FALSE))
        }
    } else {
        pass_out_read_calls <- NULL
    }

    # Prepare <OUT_PATH> argument
    if(is.null(out_extract_table) ){
        pass_out_extract_table <-  'null'
    } else{
        pass_out_extract_table <- out_extract_table
        if (verbose) {
            message("Specified path to extract table: ",
                    normalizePath(out_extract_table, mustWork = FALSE))
        }
    }


    # Execute system modkit command(s):
    # --------------------------------------------------------------------------
    if (is.null(regions)) {
        system2(
            modkit_bin,
            args=c(
                pass_ARGS,
                pass_out_read_calls,
                bamfile,
                pass_out_extract_table
            ),
            wait=TRUE
        )
    }

    if (!is.null(regions) && length(regions) == 1) {
        system2(
            modkit_bin,
            args = c(
                pass_ARGS,
                paste('--region', regions_char, sep=" "),
                pass_out_read_calls,
                bamfile,
                pass_out_extract_table
            ),
            wait = TRUE
        )
    }

    if (!is.null(regions) && length(regions) > 1) {
        # Prepare tempdir:
        tempdir <- file.path(tempdir_base, "modkit_temp")
        dir.create(tempdir, showWarnings = TRUE)
        file.remove(list.files(tempdir, full.names = TRUE))

        for (region in regions_char) {

            temp_out_extract_table <- ifelse(is.null(out_extract_table),
                                             'null',
                                             tempfile(pattern = region,
                                                      fileext = '.etbl',
                                                      tmpdir = tempdir))
            temp_out_read_calls <- tempfile(pattern = region,
                                            fileext = '.rdcl',
                                            tmpdir = tempdir)
            temp_pass_out_read_calls <- if (is.null(out_read_calls)) NULL else {
                paste('--read-calls-path', temp_out_read_calls, sep=" ") }

            system2(
                modkit_bin,
                args = c(
                    pass_ARGS,
                    paste('--region', region, sep = " "),
                    temp_pass_out_read_calls,
                    bamfile,
                    temp_out_extract_table
                ),
                wait = TRUE
            )
        }

        # Combine temporary outputs and cleanup:
        if (!is.null(out_extract_table)) {
            system(paste0("head -1 ", temp_out_extract_table, " > ",
                          pass_out_extract_table))
            system(paste0("cat ", tempdir,
                          "/*.etbl | grep -v '^read_id' | sort -u >> ",
                          pass_out_extract_table))
        }

        if (!is.null(out_read_calls)) {
            system(paste0("head -1 ", temp_out_read_calls, " > ",
                          out_read_calls))
            system(paste0("cat ", tempdir,
                          "/*.rdcl | grep -v '^read_id' | sort -u >> ",
                          out_read_calls))
        }

        unlink(tempdir, recursive=TRUE)
    }

    # return results
    # --------------------------------------------------------------------------
    return(c('extract-table' = ifelse(is.null(out_extract_table),
                                      NA, normalizePath(out_extract_table)),
              'read-calls' = ifelse(is.null(out_read_calls),
                                    NA, normalizePath(out_read_calls)),
              'run-log' = ifelse(is.null(out_log_file),
                                 NA, normalizePath(out_log_file))
    ))
}

