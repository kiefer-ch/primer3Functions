################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
# Uses large parts of https://gist.github.com/al2na/8540391
#
################################################################################

library("readr")
library("dplyr")
library("tidyr")
library("purrr")

################################################################################

#' call primer3 for a given set of DNAstringSet object
#'
#' @param seq The DNA string for which primer pairs shall be found.
#' @param seq_id The sequence id of the DNA string. Just anything at the moment.
#' @param product_size_range default: c(90-120)
#' @param Tm melting temprature parameters default: c(60, 62.5, 65)
#' @param primer_size_range default: c(15, 18, 21)
#' @param seq_target A character string of format
#'   "target,length target,length".
#' @param ovl_target A character string of space separated overlap targets.
#' @author CHK modified from Altuna Akalin modified Arnaud Krebs' original
#'   function. (For Altuna Akalins original original function see:
#'   https://gist.github.com/al2na/8540391).
#'
#' @export
call.primer3 <- function(seq, seq_id, product_size_range = c(90, 120),
    Tm = c(60, 62.5, 65), primer_size_range = c(15, 18, 21),
    seq_target = NULL, ovl_target = NULL, n = 10) {

    p3.input <- tempfile()
    p3.output <- tempfile()

    seq <- gsub('\n', '', seq)

    product_size_range <- paste(product_size_range[1], product_size_range[2],
        sep = '-')

    if (is.null(seq_target) & is.null(ovl_target)) {
        write(
            paste(sprintf("SEQUENCE_ID=%s\n", seq_id),
                sprintf("SEQUENCE_TEMPLATE=%s\n", as.character(seq)),
                "PRIMER_TASK=generic\n",
                "PRIMER_PICK_LEFT_PRIMER=1\n",
                "PRIMER_PICK_INTERNAL_OLIGO=0\n",
                "PRIMER_PICK_RIGHT_PRIMER=1\n",
                "PRIMER_EXPLAIN_FLAG=0\n",
                "PRIMER_PAIR_MAX_DIFF_TM=3\n",
                sprintf("PRIMER_MIN_SIZE=%s\n", primer_size_range[1]),
                sprintf("PRIMER_OPT_SIZE=%s\n", primer_size_range[2]),
                sprintf("PRIMER_MAX_SIZE=%s\n", primer_size_range[3]),
                sprintf("PRIMER_MIN_TM=%s\n", Tm[1]),
                sprintf("PRIMER_OPT_TM=%s\n", Tm[2]),
                sprintf("PRIMER_MAX_TM=%s\n", Tm[3]),
                sprintf("PRIMER_NUM_RETURN=%i\n", n),
                sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n", product_size_range),
                "=",
            sep = ''),
        p3.input)
    } else if (!is.null(seq_target)) {
        seq_target <- paste(seq_target, collapse = ' ')
        write(
            paste(sprintf("SEQUENCE_ID=%s\n", seq_id),
                sprintf("SEQUENCE_TEMPLATE=%s\n", as.character(seq)),
                "PRIMER_TASK=generic\n",
                "PRIMER_PICK_LEFT_PRIMER=1\n",
                "PRIMER_PICK_INTERNAL_OLIGO=0\n",
                "PRIMER_PICK_RIGHT_PRIMER=1\n",
                "PRIMER_EXPLAIN_FLAG=0\n",
                "PRIMER_PAIR_MAX_DIFF_TM=3\n",
                sprintf("PRIMER_MIN_SIZE=%s\n", primer_size_range[1]),
                sprintf("PRIMER_OPT_SIZE=%s\n", primer_size_range[2]),
                sprintf("PRIMER_MAX_SIZE=%s\n", primer_size_range[3]),
                sprintf("PRIMER_MIN_TM=%s\n", Tm[1]),
                sprintf("PRIMER_OPT_TM=%s\n", Tm[2]),
                sprintf("PRIMER_MAX_TM=%s\n", Tm[3]),
                sprintf("PRIMER_NUM_RETURN=%i\n", n),
                sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n", product_size_range),
                sprintf("SEQUENCE_TARGET=%s\n", seq_target),
                "=",
            sep = ''),
        p3.input)
    } else if (!is.null(ovl_target)) {
        ovl_target <- paste(ovl_target, collapse = ' ')
        write(
            paste(sprintf("SEQUENCE_ID=%s\n", seq_id),
                sprintf("SEQUENCE_TEMPLATE=%s\n", as.character(seq)),
                "PRIMER_TASK=generic\n",
                "PRIMER_PICK_LEFT_PRIMER=1\n",
                "PRIMER_PICK_INTERNAL_OLIGO=0\n",
                "PRIMER_PICK_RIGHT_PRIMER=1\n",
                "PRIMER_EXPLAIN_FLAG=0\n",
                "PRIMER_PAIR_MAX_DIFF_TM=3\n",
                "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=7\n",
                sprintf("PRIMER_MIN_SIZE=%s\n", primer_size_range[1]),
                sprintf("PRIMER_OPT_SIZE=%s\n", primer_size_range[2]),
                sprintf("PRIMER_MAX_SIZE=%s\n", primer_size_range[3]),
                sprintf("PRIMER_MIN_TM=%s\n", Tm[1]),
                sprintf("PRIMER_OPT_TM=%s\n", Tm[2]),
                sprintf("PRIMER_MAX_TM=%s\n", Tm[3]),
                sprintf("PRIMER_NUM_RETURN=%i\n", n),
                sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s\n", product_size_range),
                sprintf("SEQUENCE_OVERLAP_JUNCTION_LIST=%s\n", ovl_target),
                "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=5\n",
                "=",
            sep = ''),
        p3.input)
    } else {
        cat("nur eins von beiden")
    }

    system(paste("primer3_core", p3.input, '>', p3.output))
    out <- read_delim(p3.output, delim = '=', col_names = FALSE)
    unlink(c(p3.input, p3.output))

    out
}


#' Format primer3 Output File
#'
#' @param primer3_output Output file from primer3.
#'
#' @export
format.primer3 <- function(primer3_output) {
    id <- primer3_output %>%
        filter(X1 == "SEQUENCE_ID") %>%
        pull(X2)

    skip <- primer3_output %>%
        mutate(n = row_number()) %>%
        filter(X1 == "PRIMER_PAIR_NUM_RETURNED") %>%
        pull(n)

    out <- primer3_output %>%
        filter(row_number() > skip) %>%
        filter(grepl("TM", X1) | grepl("SEQUENCE", X1) | grepl("GC", X1)) %>%
        mutate(primer_pair = substr(sub("^.*T_", '', X1), 1, 1)) %>%
        mutate(type = sub("^.*T_._", '', X1)) %>%
        mutate(side = sub("_.*$", '', sub('^PRIMER_', '', X1))) %>%
        mutate(names = paste(type, side, sep = '_')) %>%
        select(-X1, -type, -side) %>%
        spread(names, X2)

    if (nrow(out) == 0) {
        return(NA)
    } else {
    out %>%
        mutate(LENGTH_LEFT = nchar(SEQUENCE_LEFT),
            LENGTH_RIGTH = nchar(SEQUENCE_RIGHT)) %>%
            return()
    }
}


#' Get A Single Sequence From A Multi Sequence Fasta File
#'
#' @param fasta Fasta to look in for sequence.
#' @param seq_id Sequence to look for.
#'
#' @example
#' \dontrun{
#'  get.seq("/home/christoph/Dokumente/ngs/transcriptomes/QUICKNONCODE_ens92_ncd5.fa.gz", "QNCMUST00114190")
#' }
#'
#' @export
get.seq <- function(fasta, seq_id) {
    stopifnot(file.exists(fasta))

    if (substr(fasta, nchar(fasta[1]) - 2, nchar(fasta)) == ".gz") {
        con <- gzfile(fasta, "r")
    } else {
        con <- file(fasta, "r")
    }

    seq <- character(length = 10)
    match <- paste0("^>", seq_id)
    i <- 1
    z <- 10
    while (TRUE) {
        line <- scan(file = con, what = "character", nlines = 1, quiet = TRUE,
            sep = '\n')
        if (grepl(match, line)) {
            break
        }
    }
    while (TRUE) {
        line <- scan(file = con, what = "character", nlines = 1, quiet = TRUE,
            sep = '\n')
        if (substr(line, 1, 1) == '>') {
            break
            close(con)
        } else {
            seq[i] <- line
            i <- i + 1
            if (i == z) {
                length(seq) <- length(seq) + 10
                z <- z + 10
            }
        }
    }
    seq <- seq[!is.na(seq)]
    paste(seq, collapse = '')
}

get.seq2 <- function (fasta, seq_id) {
    stopifnot(file.exists(fasta))

    read_file_raw(fasta)


}

#' Get A Transcript Out Of A gtf File
#'
#' @param gtf Gtf to look in for sequence.
#' @param seq_id Sequence to look for.
#'
#' @example
#' \dontrun{
#' get.transcript("/home/christoph/Dokumente/ngs/transcriptomes/QUICKNONCODE_ens92_ncd5.gtf", "QNCMUST00114190")
#' }
#'
#' @export
get.transcript <- function(gtf, seq_id) {
    grep_output <- tempfile()

    if (substr(gtf, nchar(gtf[1]) - 2, nchar(gtf)) == ".gz") {
        system(paste("zgrep", seq_id, gtf, '>', grep_output))
    } else {
        system(paste("grep", seq_id, gtf, '>', grep_output))
    }

    df <- read_tsv(grep_output, col_names = FALSE) %>%
        select(type = X3, start = X4, end = X5, strand = X7)

    strand <- unique(df$strand)

    if (strand == '+') {
        df <- df %>%
            arrange(start)
    } else {
        df <- df %>%
            arrange(desc(start))
    }

    df <- df %>%
        filter(type == "exon") %>%
        mutate(start = as.numeric(.$start), end = as.numeric(.$end)) %>%
        mutate(length = end - start + 1) %>%
        mutate(junction = cumsum(length))

    unlink(grep_output)

    df
}

#' Get The Exon Exon Junctions For A Transcript From A GTF
#'
#' @param gtf Gtf to look in.
#' @param seq_id Transcript to look for.
#'
#' @example
#' \dontrun{
#' get.junctions("/home/christoph/Dokumente/ngs/transcriptomes/QUICKNONCODE_ens92_ncd5.gtf", "QNCMUST00114190")
#' }
#'
#' @export
get.junctions <- function(gtf, seq_id) {
    junctions <- get.transcript(gtf, seq_id)$junction
    junctions[1:(length(junctions) - 1)]
}

#' Get Primer Pairs Given A Gtf And A Transcript ID
#'
#' `get.primer` returns one list of primer pairs per exon-exon junction. One of
#'  the primers in each pair will be spanning the junction.
#'
#' @param seq_id Sequence to look for.
#' @param transcriptome Transcriptome fasta file.
#' @param gtf Gtf to look in for sequence.
#' @param type One of either "overlap" or "target"
#' @param product_size_range default: '90-120'
#' @param Tm melting temprature parameters default: c(60, 62.5, 65)
#' @param primer_size_range default: c(15, 18, 21)
#' @example
#' \dontrun{
#' get.primer("QNCMUST00000012", "/home/christoph/Dokumente/ngs/transcriptomes/QUICKNONCODE_ens92_ncd5.fa.gz", "/home/christoph/Dokumente/ngs/transcriptomes/QUICKNONCODE_ens92_ncd5.gtf")
#' }
#'
#' @export
get.primer <- function(seq_id, transcriptome, gtf, type = "overlap",
    product_size_range = c(90, 120), Tm = c(60, 62.5, 65),
    primer_size_range = c(15, 18, 21)) {

    stopifnot(file.exists(transcriptome), file.exists(gtf),
        type %in% c("taget", "overlap"))

    seq <- get.seq(transcriptome, seq_id)
    junctions <- get.junctions(gtf, seq_id)

    if (type == "overlap") {
        junctions %>%
            map(call.primer3, seq = seq, seq_id = seq_id, product_size_range = product_size_range,
                Tm = Tm, primer_size_range = primer_size_range,
                seq_target = NULL) %>%
            map(format.primer3) %>%
            return(.)
    } else if (type == "target") {
        paste0(junctions, ",0") %>%
            map(call.primer3, seq = seq, seq_id = seq_id, product_size_range = product_size_range,
                Tm = Tm, primer_size_range = primer_size_range) %>%
            map(format.primer3) %>%
            return(.)
    }
}
