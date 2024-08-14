#include <htslib/sam.h>
#include <string>
#include <vector>
#include <Rcpp.h>
// using namespace Rcpp;

// Helper functions
// -----------------------------------------------------------------------------

std::vector<int> read_to_reference_pos(const bam1_t *alignment, const std::vector<int> &read_positions) {
    const uint32_t *cigar = bam_get_cigar(alignment);  // CIGAR array
    int ref_pos = alignment->core.pos;  // Reference position (0-based)
    int read_index = 0;  // Current index in the read sequence

    std::vector<int> ref_positions(read_positions.size(), -1);  // Initialize all positions to -1

    size_t read_pos_index = 0;

    // Iterate over the CIGAR operations
    for (int i = 0; i < alignment->core.n_cigar && read_pos_index < read_positions.size(); i++) {
        int op_len = bam_cigar_oplen(cigar[i]);  // Length of the operation
        int op = bam_cigar_op(cigar[i]);  // CIGAR operation type

        switch (op) {
        case BAM_CMATCH:  // Match or mismatch (M)
        case BAM_CEQUAL:  // Sequence match (=)
        case BAM_CDIFF:   // Sequence mismatch (X)
            // Check each read position within this segment
            while (read_pos_index < read_positions.size() &&
                   read_index + op_len > read_positions[read_pos_index]) {
                ref_positions[read_pos_index] = ref_pos + (read_positions[read_pos_index] - read_index);
                read_pos_index++;
            }
            ref_pos += op_len;
            read_index += op_len;
            break;

        case BAM_CINS:  // Insertion (I)
            // Handle positions within the insertion
            if (read_index + op_len > read_positions[read_pos_index]) {
                // The current read position is within an insertion
                while (read_pos_index < read_positions.size() &&
                       read_index + op_len > read_positions[read_pos_index]) {
                    ref_positions[read_pos_index] = -1;  // No corresponding reference position
                    read_pos_index++;
                }
            }
            read_index += op_len;
            break;

        case BAM_CDEL:  // Deletion (D)
        case BAM_CREF_SKIP:  // Reference skip (N)
            ref_pos += op_len;
            break;

        case BAM_CSOFT_CLIP:  // Soft clipping (S)
            // Handle positions within the soft clipping
            if (read_index + op_len > read_positions[read_pos_index]) {
                // The current read position is within a soft-clipped region
                while (read_pos_index < read_positions.size() &&
                       read_index + op_len > read_positions[read_pos_index]) {
                    ref_positions[read_pos_index] = -1;  // No corresponding reference position
                    read_pos_index++;
                }
            }
            read_index += op_len;
            break;

        case BAM_CHARD_CLIP:  // Hard clipping (H)
        case BAM_CPAD:        // Padding (P)
            // Hard clipping and padding do not consume any positions in the read or reference
            break;

        default:
            Rcpp::Rcerr << "Unknown CIGAR operation: " << op << std::endl;
        return ref_positions;
        }
    }

    return ref_positions;
}


//' Read base modifications from a bam file.
//'
//' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
//' section 1.7) and return a list of information on modified bases.
//'
//' @param inname_str Character scalar with name of the input bam file.
//' @param regions Character vector specifying the region(s) for which
//'     to extract overlapping reads, in the form \code{"chr:start-end"}
//' @param modbase Character scalar defining the modified base to extract.
//' @param verbose Logical scalar. If \code{TRUE}, report on progress.
//'
//' @examples
//' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
//'                           package = "footprintR")
//' res <- read_modbam(modbamfile, "chr1:6940000-6955000", TRUE)
//' str(res)
//'
//' @seealso https://samtools.github.io/hts-specs/SAMtags.pdf describing the
//'     SAM ML and MM tags for base modifications.
//'     Helpful examples are available in
//'      https://github.com/samtools/htslib/blob/develop/samples/modstate.c
//'     Documentation of the htslib C API is available in
//'      https://github.com/samtools/htslib/blob/develop/htslib/sam.h
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
 // [[Rcpp::export]]
Rcpp::List read_modbam(std::string inname_str,
                       std::vector<std::string> regions,
                       char modbase,
                       bool verbose = false) {
    // variable declarations
    int c = 0, i = 0, j = 0, r = 0, impl = 0;
    // int pos = 0;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    uint8_t *data = NULL;
    hts_base_mod_state *ms = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    unsigned int regcnt = 0, alncnt = 0;
    char canonical = '0', **regions_c = NULL;

    std::vector<std::string> read_name;
    std::vector<char> modified_base;
    std::vector<char> canonical_base;
    std::vector<char> strand;
    std::vector<std::string> ref_name;
    std::vector<int> read_position;
    std::vector<int> ref_position;
    std::vector<double> call_prob;

    const char* inname = inname_str.c_str();

    // // get unmodified base corresponding to modbase
    // switch (modbase) {
    // case 'm':
    // case 'h':
    // case 'f':
    // case 'c':
    // case 'C':
    //     canonical = 'C';
    //     break;
    // case 'g':
    // case 'e':
    // case 'b':
    // case 'T':
    //     canonical = 'T';
    //     break;
    // case 'U':
    //     canonical = 'U';
    //     break;
    // case 'a':
    // case 'A':
    //     canonical = 'A';
    //     break;
    // case 'o':
    // case 'G':
    //     canonical = 'G';
    //     break;
    // case 'n':
    // case 'N':
    // default:
    //     canonical = 'N';
    // }

    // initialize bam data storage
    if (!(bamdata = bam_init1())) {
        Rprintf("Failed to initialize bamdata\n");
        goto end;
    }
    if (!(ms = hts_base_mod_state_alloc())) {
        printf("Failed to allocate state memory\n");
        goto end;
    }

    // open input file
    if (verbose) {
        Rprintf("opening input file %s\n", inname);
    }
    if (!(infile = sam_open(inname, "r"))) {
        Rprintf("Could not open input file %s\n", inname);
        goto end;
    }

    // load index file
    if (!(idx = sam_index_load(infile, inname))) {
        Rprintf("Failed to load the index for %s\n", inname);
        goto end;
    }

    // read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        Rprintf("Failed to read header from file %s\n", inname);
        goto end;
    }

    // convert regions to C arrays
    regcnt = (unsigned int) regions.size();
    regions_c = (char**) calloc(regcnt, sizeof(char*));
    for (i = 0; i < regcnt; i++) {
        regions_c[i] = (char*) regions[i].c_str();
    }

    // create multi-region iterator
    if (!(iter = sam_itr_regarray(idx, in_samhdr, regions_c, regcnt))) {
        printf("Failed to get iterator\n");
        goto end;
    }

    // iterate over regions
    if (verbose) {
        Rprintf("reading alignments overlapping any of %u regions\n", regcnt);
    }
    // read overlapping alignments using iterator
    while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
        alncnt++;

        // process alignment
        // ... get sequence
        i = 0; // position in read sequence
        data = bam_get_seq(bamdata);

        // ... parse base modifications
        if (bam_parse_basemod(bamdata, ms)) {
            Rprintf("Failed to parse the base mods (read %s)\n",
                    bam_get_qname(bamdata));
            goto end;
        }
        hts_base_mod mod[5] = {{0}};  //for ATCGN

        // ... process read if modifications of the right type are present
        //     bam_mods_query_type:
        //     - returns 0 on success, -1 if not found
        //     - also fills out `canonical` (unmodified base) and `impl`
        //       (boolean for whether unlisted positions should be implicitly
        //        assumed to be unmodified, or require an explicit score and
        //        should be considered as unknown)
        if (bam_mods_query_type(ms, modbase, NULL, &impl, &canonical) == 0) {
            // ... loop over sequence positions i
            for (; i < bamdata->core.l_qseq; i++) {
                if ((r = bam_mods_at_next_pos(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]))) <= -1) {
                    // r: number of found modifications (>=1, 0 or -1 if failed)
                    Rprintf("Failed to get modifications (read %s)\n",
                            bam_get_qname(bamdata));
                    goto end;

                } else if (r > (sizeof(mod) / sizeof(mod[0]))) {
                    Rprintf("More modifications than this app can handle, update the app\n");
                    goto end;

                } else if (!r && impl) {
                    // no modification at position i
                    if (seq_nt16_str[bam_seqi(data, i)] == canonical) {
                        // found omitted, likely unmodified base -> add to results
                        read_name.push_back(bam_get_qname(bamdata));
                        read_position.push_back(i);
                        ref_name.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                        modified_base.push_back('-');
                        canonical_base.push_back(canonical);
                        strand.push_back(bam_is_rev(bamdata) ? '-' : '+');
                        call_prob.push_back(-1.0); // special value of -1.0 indicates inferred unmodified base
                    }
                }
                //modifications
                for (j = 0; j < r; j++) {
                    if (mod[j].modified_base == modbase) {
                        // found modified base of the right type -> add to results
                        read_name.push_back(bam_get_qname(bamdata));
                        read_position.push_back(i);
                        ref_name.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                        modified_base.push_back((char) mod[j].modified_base);
                        canonical_base.push_back((char) mod[j].canonical_base);
                        strand.push_back(mod[j].strand ? '-' : '+');
                        // qual of `N` corresponds to call probability
                        //     in [N/256, (N+1)/256] -> store midpoint
                        call_prob.push_back(((double) mod[j].qual + 0.5) / 256.0);
                    }
                }
            }
        }

        // while ((r = bam_next_basemod(bamdata, ms, mod,
        //                              sizeof(mod)/sizeof(mod[0]), &pos)) >= 0) {
        //     // r: number of found modifications
        //     //   (typically 1, or 0 at the end, or -1 if failed)
        //
        //     // non-modified bases omitted from BAM
        //     for (; i < bamdata->core.l_qseq && i < pos; i++) {
        //         if (seq_nt16_str[bam_seqi(data, i)] == canonical) {
        //             // found omitted base -> add to results
        //             read_name.push_back(bam_get_qname(bamdata));
        //             read_position.push_back(i);
        //             ref_name.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
        //             modified_base.push_back('-');
        //             canonical_base.push_back(canonical);
        //             strand.push_back(bam_is_rev(bamdata) ? '-' : '+');
        //             call_prob.push_back(-1.0); // special value of -1.0 indicates inferred unmodified base
        //         }
        //     }
        //
        //     // modifications
        //     for (j = 0; j < r; j++) {
        //         // Rprintf("%c %c %c %d\n", mod[j].canonical_base, mod[j].strand ? '-' : '+', mod[j].modified_base, pos);
        //         // append modification to vectors
        //         if (mod[j].modified_base == modbase) {
        //             read_name.push_back(bam_get_qname(bamdata));
        //             read_position.push_back(pos);
        //             ref_name.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
        //             modified_base.push_back((char) mod[j].modified_base);
        //             canonical_base.push_back((char) mod[j].canonical_base);
        //             strand.push_back(mod[j].strand ? '-' : '+');
        //             // qual of `N` corresponds to call probalitiy in [N/256, (N+1)/256] -> use midpoint
        //             call_prob.push_back(((double) mod[j].qual + 0.5) / 256.0);
        //         }
        //     }
        //
        //     // skip the modification already processed
        //     if (i == pos)
        //         i++;
        //
        //     // non-modified bases after the last modified base omitted from the BAM file
        //     if (!r) {
        //         for (; i < bamdata->core.l_qseq; ++i) {
        //             if (seq_nt16_str[bam_seqi(data, i)] == canonical) {
        //                 // found omitted base -> add to results
        //                 read_name.push_back(bam_get_qname(bamdata));
        //                 read_position.push_back(i);
        //                 ref_name.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
        //                 modified_base.push_back('-');
        //                 canonical_base.push_back(canonical);
        //                 strand.push_back(bam_is_rev(bamdata) ? '-' : '+');
        //                 call_prob.push_back(-1.0); // special value of -1.0 indicates inferred unmodified base
        //             }
        //         }
        //         break;
        //     }
        // }
        // if (r <= -1) {
        //     Rprintf("Failed to get modifications (read %s)\n", bam_get_qname(bamdata));
        //     goto end;
        // }

        // ... convert read positions to reference coordinates
        std::vector<int> read_position_converted = read_to_reference_pos(bamdata, read_position);
        ref_position.reserve(ref_position.size() + read_position_converted.size());
        ref_position.insert(ref_position.end(),
                            read_position_converted.begin(),
                            read_position_converted.end());
        read_position.clear();
    }
    if (c != -1) {
        Rprintf("Error while reading from %s - aborting\n", inname);
        goto end;
    }
    if (verbose) {
        Rprintf("read %u alignments\n", alncnt);
    }

    end:
        // create return list
        Rcpp::List res = Rcpp::List::create(
            Rcpp::_["read_name"] = read_name,
            Rcpp::_["ref_name"] = ref_name,
            Rcpp::_["ref_position"] = ref_position,
            Rcpp::_["modified_base"] = modified_base,
            Rcpp::_["canonical_base"] = canonical_base,
            Rcpp::_["strand"] = strand,
            Rcpp::_["call_prob"] = call_prob);

        //cleanup
        if (regions_c) {
            free((void*) regions_c);
            regions_c = NULL;
        }
        if (in_samhdr) {
            sam_hdr_destroy(in_samhdr);
        }
        if (infile) {
            sam_close(infile);
        }
        if (bamdata) {
            bam_destroy1(bamdata);
        }
        if (ms) {
            hts_base_mod_state_free(ms);
        }
        if (iter) {
            sam_itr_destroy(iter);
        }
        if (idx) {
            hts_idx_destroy(idx);
        }
        return res;
}
