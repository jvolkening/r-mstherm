#' @title Create a new MSThermExperiment.
#'
#' @description \code{MSThermExperiment} creates a new experiment object from
#'   a set of filenames or data frames.
#'
#' @param control data frame or filename of tab-delimited table describing the
#'   experimental setup and locations of data on disk (see Details)
#' @param annotations data frame or filename to tab-delimited table containing
#'   protein names and annotations (usually functional descriptions but can be
#'   any text
#'
#' @details Both parameters can take either a data frame or a tab-delimited
#'   filename on disk (which will be read into a data frame). "control" should
#'   contain columns with the following headers (in any order):
#'     \describe{
#'        \item{"name"}{Unique identifier of a single replicate}
#'        \item{"sample"}{Sample name that a replicate belongs to}
#'        \item{"data_file"}{Path to file on disk containing the quantification
#'            data}
#'        \item{"meta_file"}{Path to file on disk containing the labeling
#'            metadata}
#'     }
#'   The "meta_file" should be tab-delimited text and contain two columns
#'   labeled "channel" and "temp".
#'   The "data_file" should be tab-delimited text and contain, at a minimum,
#'   the following columns:
#'   \describe{
#'      \item{"peptide"}{Sequence of the matched peptide in single-letter IUPAC}
#'      \item{"protein"}{Protein or protein group to which the peptide belongs}
#'      \item{"..."}{One column per isobaric channel, containing absolute
#'         quantification values. Column names must match those in the
#'         "channel" column of the meta file, with the exception that R will
#'         automatically convert any name not compatible with its syntax
#'         rules. To be safe, use only letters, digits, underscores, and
#'         periods in channel names and never start with a digit (e.g. use "TMT.126"
#'         rather than "126")}
#'   }
#'   The following columns can also be utilized for filtering if included (all
#'   others will simply be ignored):
#'   \describe{
#'      \item{"coelute_inf"}{Calculated precursor co-isolation interference
#'         (0.0-1.0)}
#'      \item{"score"}{Score assigned by the processing software to the PSM}
#'   }
#'   "annotations" should contain two columns with the headers "name" and
#'   "annotation". "name" should match the protein names in the data files,
#'   and "annotation" can contain any text (generally a functional description)
#'
#' @return An MSThermExperiment object
#'
#' @examples
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#'
#' @export

MSThermExperiment <- function(control, annotations) {

    self <- structure(
        list(
            samples = list(),
            annot   = list()
        ),
        class = "MSThermExperiment"
    )
       
    # read in and populate annotations if present
    if (!missing(annotations)) {
        self$annot <- .to_dataframe(annotations)
    }

    # return empty object if no control file/data frame specified
    if (missing(control)) { return(self) }

   
    wd <- getwd()
    if (is.character(control)) {
        wd <- dirname(normalizePath(control))
    }
    control <- .to_dataframe(control)

    sample_names <- unique(control$sample)
    self <- add_samples( expt=self, samples = lapply(
        sample_names, function(s) {
            sample <- MSThermSample(s)
            replicate_names <- control$name[control$sample == s]
            sample <- add_replicates(sample=sample, replicates= lapply(
                replicate_names, function(r) {
                    repl <- MSThermReplicate(
                        name      = r,
                        data = control$data_file[control$name == r],
                        meta = control$meta_file[control$name == r],
                        wd   = wd
                    )
                }
            ))
        }
    ))

    return(self)

}



#----------------------------------------------------------------------------#
# These constructors are not called directly and thus are not documented in
# detail
#----------------------------------------------------------------------------#

MSThermSample <- function(name) {

    self <- structure( list( name = name,replicates = list() ),
        class = "MSThermSample" )

    return(self)

}



MSThermReplicate <- function(name, data, meta, wd) {

    # data and meta paths are relative to the control filename, if given,
    # so temporarily switch to that directory
    wd_old <- getwd()
    setwd(wd)

    data <- .to_dataframe(data)
    meta <- .to_dataframe(meta)
    meta <- meta[order(meta$temp,decreasing=F),]

    # and switch back
    setwd(wd_old)

    self <- structure(
        list(
            meta = meta,
            data = data,
            name = name
        ),
        class = "MSThermReplicate"
    )


    return( self )

}



#' @title Summarize MSThermResult object.
#'
#' @description Print a summary of an MSThermResult, including samples and
#' parameters.
#'
#' @param object an MSThermResult object
#' @param ... additional arguments passed to or from other functions
#'
#' @return Nothing
#'
#' @examples
#' # see model_protein() for an example
#'
#' @export

summary.MSThermResult <- function(object, ...) {

    cat(paste("Name:",object$name),"\n",sep='')
    cat("Samples:\n");
    for (i in object$sample_names) {
        cat("    ",i,"\n",sep='')
    }
    cat("Parameters:\n");
    for (i in names(object$parameters)) {
        cat("    ",i,": ",object$parameters[[i]],"\n",sep='')
    }

    return(invisible(NULL))

}



#' @title Summarize MSThermResultSet object.
#'
#' @description Print a summary of an MSThermResultSet, including samples and
#' parameters.
#'
#' @param object an MSThermResultSet object
#' @param ... additional arguments passed to or from other functions
#'
#' @return Nothing
#'
#' @examples
#' # see model_experiment() for an example
#'
#' @export

summary.MSThermResultSet <- function(object, ...) {

    result <- object[[1]]
    cat("Samples:\n");
    for (i in result$sample_names) {
        cat("    ",i,"\n",sep='')
    }
    cat("Parameters:\n");
    for (i in names(result$parameters)) {
        if (i != 'protein') {
            cat("    ",i,": ",result$parameters[[i]],"\n",sep='')
        }
    }

    return(invisible(NULL))

}



# populate MSThermSample object with child replicates
add_replicates <- function(sample, replicates) {

    names(replicates) <- sapply(replicates, '[[', "name")
    sample$replicates = c(sample$replicates, replicates)
    return(sample)

}



# populate MSThermExperiment object with child samples
add_samples <- function(expt, samples) {

    names(samples) <- sapply(samples, '[[', "name")
    expt$samples = c(expt$samples, samples)
    return(expt)

}



# takes a value that could be either data frame or string. If string, assumes
# a file path and tries to read into a data frame. Returns a data frame or
# throws and exception.
.to_dataframe <- function(str) {
    if (is.character(str)) {
        str <- read.delim(str,stringsAsFactors=F,header=T,comment.char="#")
    }
    if (! is.data.frame(str)) {
        stop("input must be a valid filename or data.frame")
    }
    return(str)
}



# generates an annotation/description string for an input protein or protein
# group ID, splitting the input and combining the output with the given record
# separator if necessary
gen_description <- function(expt, protein, sep='|') {

    parts <- strsplit(protein,sep,fixed=T)[[1]]
    desc <- sapply(parts, function(x) ifelse( x %in% expt$annot$name,
        expt$annot$annotation[expt$annot$name == x],
        ''))
    return(paste(desc,collapse=sep))

}


