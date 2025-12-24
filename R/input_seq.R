#' Import region / constructing regioned_dna object
#'
#' Constructing \code{regioned_dna} from DNAStringSet. Optionally input a
#' \code{region} data.frame to define restricted amino-acid region for mutation.
#'
#' @param object Filepath or DNAstringSet. The input sequences is suggested to
#'   be in open reading frame(ORF).
#' @param region \code{NA}. A data.frame specifying paticular regions (positions
#'   in amino acid sequence) that is allowed to be mutated in the sequences.
#'   Both \code{1 / 0} or \code{TRUE / FALSE} encoding is OK. Alternatively, a
#'   numeric vector of 1-based codon positions to be mutated can be supplied
#'   and will be applied to all sequences.
#'   Examples below for reference.
#' @param ... ...
#' @return A regioned_dna-class object
#' @seealso \code{\link{get_cu}}, \code{\link{get_du}},
#'   \code{\link{get_region}}, \code{\link{get_dna}}
#' @examples
#' # Creating a input_seq class directly from system file
#' filepath <- system.file("extdata", "example.fasta", package = "SynMut")
#' rgd.seq <- input_seq(filepath)
#'
#' # Optionally input with region dataframe
#' filepath.fasta <- system.file("extdata", "example.fasta", package = "SynMut")
#' fp.csv <- system.file("extdata", "target_regions.csv", package = "SynMut")
#' region <- read.csv(fp.csv)
#' rgd.seq <- input_seq(filepath.fasta, region)
#'
#' # Creating from exsisting DNAStringSet object
#' seq <- Biostrings::DNAStringSet("ATCGATCGA")
#' rgd.seq <- input_seq(seq)
#'
#' # Input with codon position vector (1-based)
#' rgd.seq <- input_seq(seq, region = c(1, 3))
#'
#' @name input_seq
#' @rdname input_seq-methods
#' @import methods
#' @import BiocGenerics
#' @exportMethod input_seq
setGeneric(
    name = "input_seq",
    def = function(object, region = NA, ...) standardGeneric("input_seq")
)

#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings oligonucleotideFrequency
#' @rdname input_seq-methods
setMethod(
    f = "input_seq",
    signature = signature(object = "character"),
    definition = function(object, region) {
        dnaseq <- readDNAStringSet(filepath = object)
        dnaseq <-
            append(dnaseq, DNAStringSet('atg')) #helper sequence
        gernerate_rgd(dnaseq, region)
    }
)

#' @rdname input_seq-methods
setMethod(
    f = "input_seq",
    signature = signature(object = "DNAStringSet"),
    definition = function(object, region) {
        dnaseq <-
            append(object, DNAStringSet('atg')) #helper sequence
        gernerate_rgd(dnaseq, region)
    }
)

#' @rdname input_seq-methods
setMethod(
    f = "input_seq",
    signature = signature(object = "DNAString"),
    definition = function(object, region) {
        dnaseq <- #helper sequence
            append(DNAStringSet(object), DNAStringSet('atg'))
        gernerate_rgd(dnaseq, region)
    }
)


# helper function ---------------------------------------------------------

gernerate_rgd <- function(dnaseq, region){
    if (all(is.na(region))) {
        return(new(
            "regioned_dna",
            dnaseq = dnaseq,
            region = list(NA)
        ))
    } else {
        if (is.numeric(region) && is.vector(region) && length(region) > 0) {
            positions <- unique(as.integer(region))
            if (any(is.na(positions)) || any(positions < 1)) {
                stop("codon positions must be positive integers")
            }
            seq.lengths <- vapply(dnaseq, length, numeric(1))
            seq.codons <- seq.lengths / 3
            seq.codons <- seq.codons[seq_len(length(seq.codons) - 1)]
            region <- lapply(seq.codons, function(n.codons) {
                if (any(positions > n.codons)) {
                    stop("codon positions exceed sequence length")
                }
                mask <- rep(FALSE, n.codons)
                mask[positions] <- TRUE
                mask
            })
            return(new(
                "regioned_dna",
                dnaseq = dnaseq,
                region = c(region, list(TRUE))
            ))
        }
        if (!is(region, "data.frame")) {
            stop("the region input must be in data.frame format")
        }
        region <- as.list(region)
        region <-
            lapply(region, function(x) {
                as.logical(x[!is.na(x)])
            })
        return(new(
            "regioned_dna",
            dnaseq = dnaseq,
            region = c(region, list(TRUE))
        ))
    }
}
