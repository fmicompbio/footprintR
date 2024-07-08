#' Wrapper around the `modkit extract` command to extract read-level base modification information from modBAM files into tab-separated values table(s)
#'
#' @description
#' This function is a wrapper around the `modkit extract` sub-command to extract read-level base modification information
#' from modBAM files into tab-separated values table(s).
#' For more information on available `modkit extract` arguments and output tables specifciation see https://nanoporetech.github.io/modkit/intro_extract.html
#'
#' @param modkit_bin Character vector specifying the path of the \code{modkit}
#'     binary. 
#' @param bamfile Character vector specifying the path to a \code{modBAM}
#'     file. An indexed BAM files can significanly speed up certain extract
#'     operations.
#' @param regions A \code{\link[GRanges]{GRanges}} object specifying which genomic regions to extract the reads from.
#'     Note that the reads are not trimmed to the boundaries of the specified ranges.
#'     As a result returned positions will typically extend out of the specified regions
#' @param num_reads Number of reads to  extract. Note that this is extracted reads per specified genomic range.
#'     When N genomic ranges are specified the number of extracted reads will be at most num_reads x N
#' @param out_extract_table Character vector specifying the path for the `extract table`. Can be NULL 
#'     if only a `read-calls` table is needed.
#' @param out_read_calls Character vector specifying the path for the `read-calls table`. Can be NULL 
#'     if only an `extract tabe` is needed.
#' @param out_log_file Character vector specifying the path to file to write the command run log. Can be NULL 
#'     if no log file is needed (not recommended).
#' @param modkit_args Character vector with additional `modkit extract` arguments. 
#'     Please refer to the `modkit extract` documentation for a complete list of possible arguments.
#' @param tempdir_base Character vector specifying the path to create the `modkit_temp` directory.
#'     This temporary directory is only created when multiple genomic regions are specified
#'
#' @return A named character vector specifying the paths to the generated table and log files.
#'     Has slots for 'extract-table', 'read-calls' and 'run-log'.
#'
#' @author Panagiotis Papapasaikas
#'
#' @examples
#' \dontrun{
#'
#' suppressWarnings(
#'    GR <- c( GRanges("chr1:12340678-12345678"), 
#'             GRanges("chr2:12340678-12345678")
#'    )
#' )
#'
#' modkit_args=c(
#'    '-t 8', '--force',
#'    '--mapped-only', '--edge-filter 100',
#'    '-p  0.2',
#'    '--mod-threshold m:0.8',
#'    '--mod-threshold a:0.5'
#' )
#'
#'
#' # Produce both  `extract table` and `read calls` files for multiple GRanges
#' modkitExtract(bamfile = BAMF,num_reads = 10,regions=GR[1:2], out_extract_table = "test.etbl",out_read_calls = "test.rdcl",modkit_args=modkit_args )
#' # Produce only  `extract table` file
#' modkitExtract(bamfile = BAMF,num_reads = 10,regions=GR[1], out_extract_table = "test.etbl",out_read_calls = NULL,modkit_args=modkit_args )
#' # Produce only  `read calls` file
#' modkitExtract(bamfile = BAMF,num_reads = 10,regions=GR[1], out_extract_table = NULL,out_read_calls = "test.rdcl",modkit_args=modkit_args )
#' }
#' 
#'
#'
#' @seealso [`modkit` software](https://nanoporetech.github.io/modkit),
#'     [`modkit extract` documentation and tabulated output formats specification ](https://nanoporetech.github.io/modkit/intro_extract.html),
#'     \code{\link[GRanges]{GRanges}} for the object used to specify genomic regions,
#'
#' @import GenomicRanges
#'
#' @export
#' 
modkitExtract <- function(modkit_bin="/tungstenfs/groups/gbioinfo/Appz/ONT_modkit/modkit_v0.3.0_centos7_x86_64/dist/modkit",
                          bamfile=NULL,
                          regions=NULL,
                          num_reads=NULL,
                          out_extract_table = NULL,
                          out_read_calls = NULL,
                          out_log_file=NULL,
                          modkit_args=NULL,
                          tempdir_base=tmpdir() ) {
    
    # digest arguments
    .assertScalar( x=modkit_bin, type = "character")
    if (!file.exists(modkit_bin)) {
        stop("modkit binary not present at: ", normalizePath(modkit_bin))
    }
    .assertScalar(x = bamfile, type = "character", allowNULL = TRUE)
    if (!file.exists(bamfile)) {
        stop("BAM file not found at: ", normalizePath(bamfile))
    }
    
    .assertScalar(x = num_reads, type = "numeric", allowNULL= TRUE)
    
    .assertScalar(x = out_read_calls, type = "character", allowNULL = TRUE)
    
    .assertScalar(x = out_extract_table, type = "character", allowNULL = TRUE)
    
    .assertScalar(x = out_log_file, type = "character", allowNULL = TRUE)
    
    .assertScalar(x = tempdir_base, type = "character")
    
    
    
    
    if (!is.null(regions) && !inherits(regions,"GRanges")){
        stop("`regions` should be a GRanges object")
    }
    
    #Convert GRanges to `chr:start-end` format
    if (!is.null(regions)){
        regions_char <- as.character(reduce(GRanges(regions, strand="*")))
    }
    
    
    
    
    
    #character vector of arguments to be passed to modkit.
    #first argument should always be the modkit command. In this case 'extract'
    pass_ARGS <- c('extract', '--suppress-progress')
    
    #Prepare --num-reads argument
    if(!is.null(num_reads)){
        pass_ARGS <- c(pass_ARGS, paste("--num-reads",num_reads, sep=" ") )
        if (is.null(regions)){
            pass_ARGS <- c(pass_ARGS, "--ignore-index") #in absence of a specified region --num-reads is only respected iff bam index is ignored
        }
    }
    
    #Prepare --log-fileparh argument
    if(!is.null(out_log_file)){
        pass_ARGS <- c(pass_ARGS, paste("--log-filepath",out_log_file, sep=" ") )
        print(c( "Specified path to run log:", normalizePath(out_log_file, mustWork=FALSE) ))
        if (file.exists(out_log_file)){
            print("Watning: Specified `out_log_file` already exists. The log will be appened to the existing file")
        }
    }
    
    #combine function exposed args with user-passed modkit args:
    pass_ARGS <- c(pass_ARGS,modkit_args)
    
    
    #Prepare separately --read-calls-path argument
    if(!is.null(out_read_calls) ){
        pass_out_read_calls <- paste('--read-calls-path',out_read_calls,sep=" ")
        print(c( "Specified path to read-calls table:", normalizePath(out_read_calls, mustWork=FALSE) ) )
    } else{
        pass_out_read_calls <- NULL
    }
    
    
    #Prepare <OUT_PATH> argument
    if(is.null(out_extract_table) ){
        pass_out_extract_table <-  'null'
    } else{
        pass_out_extract_table <- out_extract_table
        print(c( "Specified path to extract table:", normalizePath(out_extract_table, mustWork=FALSE) ))
    }
    

    
    #Execute system modkit command(s):
    
    if (is.null(regions)){
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
    
    
    
    if (!is.null(regions) & length(regions)==1 ){
        system2( 
            modkit_bin,
            args=c(
                pass_ARGS,
                paste('--region',regions_char,sep=" "),
                pass_out_read_calls,
                bamfile,
                pass_out_extract_table
            ),
            wait=TRUE
        )    
    }
    
    
    
    if (!is.null(regions) & length(regions) > 1 ){
        #Prepare tempdir:
        tempdir <- paste0(tempdir_base,"/modkit_temp" )
        dir.create(tempdir, showWarnings = TRUE)
        file.remove(list.files(tempdir, full.names = TRUE))
        
        for(region in regions_char ){
            
            temp_out_extract_table <- ifelse(is.null(out_extract_table),'null', tempfile(pattern=region, fileext= '.etbl', tmpdir = tempdir )   ) 
            temp_out_read_calls <- tempfile(pattern=region, fileext= '.rdcl', tmpdir = tempdir )
            temp_pass_out_read_calls <- if(is.null(out_read_calls)) NULL else {paste('--read-calls-path',temp_out_read_calls,sep=" " )}
            
            system2( 
                modkit_bin,
                args=c(
                    pass_ARGS,
                    paste('--region',region,sep=" "),
                    temp_pass_out_read_calls,
                    bamfile,
                    temp_out_extract_table
                ),
                wait=TRUE
            ) 
            
        }
        
        #Combine temporary outputs and cleanup:
        if (!is.null(out_extract_table)){
            system( paste0 ("head -1 ", temp_out_extract_table, " > ", pass_out_extract_table ) )
            system ( paste0( "cat ", tempdir, "/*.etbl | grep -v '^read_id' | sort -u >> ", pass_out_extract_table) )
        }
        
        if (!is.null(out_read_calls)){
            system( paste0 ("head -1 ", temp_out_read_calls, " > ", out_read_calls ) )
            system (paste0( "cat ", tempdir, "/*.rdcl | grep -v '^read_id' | sort -u >> ", out_read_calls ) )
        }
        
        
        
        unlink(tempdir, recursive=TRUE)
    }
    
    
    return( c('extract-table'=ifelse(is.null(out_extract_table),NA,normalizePath(out_extract_table)) ,
              'read-calls'=ifelse(is.null(out_read_calls),NA,normalizePath(out_read_calls)) ,
              'run-log'=ifelse(is.null(out_log_file),NA,normalizePath(out_log_file)) 
    )  
    )
    
}

