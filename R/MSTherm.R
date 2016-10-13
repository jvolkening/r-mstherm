#' @importFrom stats coefficients density dnorm loess mad median nls.control
#'  p.adjust pnorm quantile var
#' @importFrom graphics par plot points rect mtext lines curve legend arrows
#'  abline
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nls2 nls2
#' @importFrom plotrix addtable2plot

#----------------------------------------------------------------------------#
# Generate ratio values from absolute quantitation
#
# data   : vector of absolute quantitation values
# method : how to calculate reference value
#----------------------------------------------------------------------------#

abs_to_ratio <- function(data,method='first') {

    denom <- switch( method,
        first  = data[1],
        max    = max(data),
        top3   = mean( data[ order(data,decreasing=T)[1:3] ] ),
        near   = median( data[ data > data[1]*0.8 ] ),
        compat = {
                m <- mean( data[ order(data,decreasing=T)[1:3] ] )
                b <- mean( data[ data > m*0.8 & data < m*1.2 ] )
                ifelse(is.na(b),m,b)
            },
        stop("Invalid method type",call.=T)
    )

    return( data/denom )

}


#----------------------------------------------------------------------------#
# Generate protein ratio profile from matrix of spectrum quantification values
#
# data   : data frame or matrix of quantification values, one row per spectrum
# method : how to roll up spectrum data into protein-level data
# ...    : parameters passed to 'abs_to_ratio'
#----------------------------------------------------------------------------#

gen_profile <- function( data, method='sum', method.denom='first' ) {

    summarized <- switch( method,
        sum = apply(data, 2, sum),
        median = apply(data, 2, median),
        # for these, the inner apply() will transpose, so the outer is applied
        # to rows rather than columns
        ratio.median = 
            apply( apply(data,1,abs_to_ratio,method=method.denom),1,median ),
        ratio.mean = 
            apply( apply(data,1,abs_to_ratio,method=method.denom),1,mean ),
        stop("Invalid method type", call.=T)
    )
    return( abs_to_ratio(summarized, method=method.denom) )

}

#' Create a new MSThermExperiment
#'
#' Creates a new experiment object from a set of filenames or dataframes
#'
#' @param control dataframe or filename to tab-delimited table describing the
#'    experimental setup and locations of data on disk (see Details)
#' @param annotations dataframe or filename to tab-delimited table containing
#'   protein names and annotations (usually functional descriptions but can be
#'   any text
#'
#' @details Both parameters can take either a dataframe or a filename on disk
#'   (which will be read into a dataframe). "control" should contain
#'   columns with the following headers (in any order):
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
#'   the following columns (others will currently be ignored):
#'   \describe{
#'      \item{"peptide"}{Sequence of the matched peptide in single-letter IUPAC}
#'      \item{"protein"}{Protein or protein group to which the peptide belongs}
#'      \item{"coelute_inf"}{Calculated precursor co-isolation interference
#'         (0.0-1.0)}
#'      \item{"..."}{One column per isobaric channel, containing absolute
#'         quantification values. Column names must match those in the
#'         "channel" column of the meta file, with the exception that R will
#'         automatically convert any name not compatible with its syntax
#'         rules. To be safe, use only letters, digits, underscores, and
#'         periods in channel names and never start with a digit (e.g. use "TMT.126"
#'         rather than "126")}
#'   }
#'   "annotations" should contain two columns with the headers "name" and
#'   "annotation". "name" should match the protein names in the data files,
#'   and "annotation" can contain any text (generally a functional description)
#'
#' @return An MSThermExperiment object
#'
#' @examples
#'\dontrun{
#' expt    <- MSThermExperiment(control="control.tsv",annotations="annots.tsv")
#'}
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

    # return empty object if no control file/dataframe specified
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

add_replicates <- function(sample,replicates) {

    names(replicates) <- sapply(replicates, '[[', "name")
    sample$replicates = c(sample$replicates, replicates)
    return(sample)

}

add_samples <- function(expt,samples) {

    names(samples) <- sapply(samples, '[[', "name")
    expt$samples = c(expt$samples, samples)
    return(expt)

}

.to_dataframe <- function(str) {
    if (is.character(str)) {
        str <- read.delim(str,stringsAsFactors=F,header=T,comment.char="#")
    }
    if (! is.data.frame(str)) {
        stop("input must be a valid filename or data.frame")
    }
    return(str)
}
    
sigShift <- function( df, repl1, repl2, bin ) {

    deltas <- df[[paste0(repl2,'.tm')]] - df[[paste0(repl1,'.tm')]]
    m <- median(deltas,na.rm=T)
    d <- mad(deltas,na.rm=T)
    print(m)
    print(d)
    z <- (deltas - m)/d
    p <- 2*pnorm(-abs(z))
    q <- p.adjust(p,method="BH")
    plot(density(deltas,na.rm=T),xlim=c(-10,10))
    x <- NULL; rm(x) # silence R CMD check noise due to curve() call below
    curve(dnorm(x,mean=m,sd=d),col="red",add=T)
    return(q)

} 
    
#' MSResultSet to data frame
#'
#' Populates a dataframe with information from an MSResultSet, one row per
#' protein/group
#'
#' @param x an MSResultSet object
#' @param ... additional arguments passed to or from other functions
#'
#' @return A dataframe populated with relevant information per result
#'
#' @examples
#'\dontrun{
#' df <- as.data.frame(expt)
#' write.table(df, "results.tsv")
#'}
#'
#' @export

as.data.frame.MSThermResultSet <- function( x, ... ) {

    df <- data.frame(
        row.names = sapply(x, '[[', "name")
    )
    df$annotation <- sapply(x, '[[', "annotation")

    repl_lists <- lapply(x, function(d)
        unique(sapply( d$series, '[[', "name" )))
    repl_names <- unique(unlist(repl_lists))
    repl_names <- repl_names[order(repl_names)]

    for (r in repl_names) {
        x <- x[1]$series[[r]][['x']]

        for(col in c("tm","psm","inf","slope","k","plat","r2","rmsd")) {
            df[[paste0(r,'.',col)]]  <- sapply(x, function(v) v$series[[r]][[col]])
        }
    }

    return(df)

} 

#' Re-normalize based on Tm
#'
#' Normalizes each replicate of an experiment based on linear regression
#' of calculated Tm (corrects for remaining systematic error)
#'
#' @param expt An MSThermExperiment object
#' @param res An MSThermResultSet object
#'
#' @return An MsThermExperiment object with re-normalized data slots
#'
#' @examples
#'\dontrun{
#' expt <- normalize_to_tm(expt, res)
#'}
#'
#' @export

normalize_to_tm <- function( expt, res ) {

    repl_names <- sapply(res[[1]]$series, '[[', "name")

    # filter low-quality fits
    c <- 1
    m.r2 <- matrix(ncol=length(repl_names),nrow=length(res))
    for (r in repl_names) {
        m.r2[,c] <- sapply(res, function(v) v$series[[r]][['r2']])
        c <- c + 1
    }
    min.r2 <- apply(m.r2,1,min)
    res <- res[min.r2 > .98 & ! is.na(min.r2)]
    
    # determine baseline replicate
    if (length(repl_names) < 3) {
        baseline <- repl_names[1]
    }
    else {
        c <- 1
        m <- matrix(ncol=length(repl_names))
        for (r in repl_names) {
            m[,c] <- sapply(res, function(v) v$series[[r]][['tm']])
            c <- c + 1
        }
        meds <- apply(m,1,median)
        c <- 1
        diffsums = c()
        for (r in repl_names) {
            diffsums[c] <- sum(abs(m[,c]-meds))
        }
        baseline <- repl_names[which.min(diffsums)]
    }
    tm.b <- as.vector(sapply(res, function(v) v$series[[baseline]][['tm']]))
    sample_names <- sapply(expt$samples, '[[', "name")
    for (r in repl_names) {
        print(r)
        if (r == baseline) {
            next
        }
        tm.u <- as.vector(sapply(res, function(v) v$series[[r]][['tm']]))
        #l <- lm(tm.b ~ tm.u)
        l <- median(tm.b - tm.u)

        for (s in sample_names) {
            print(s)
            r_names <- sapply(expt$samples[[s]]$replicates, '[[', "name")
            print(r_names)
            for (r2 in r_names) {
                print(r2)
                if (r2 == r) {
                    print(l)
                    #expt$samples[[s]]$replicates[[r]]$meta$temp <- expt$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]
                    expt$samples[[s]]$replicates[[r]]$meta$temp <- expt$samples[[s]]$replicates[[r]]$meta$temp + l
                }
            }
        }

    }

    return( expt )

}

#' Normalize to a spike-in standard
#'
#' Normalizes each replicate of an experiment based on a given spike-in
#' protein standard (assumed to be present in equimolar amounts in each
#' channel)
#'
#' @param expt An MSThermExperiment object
#' @param protein ID of a protein to normalize against
#'
#' @return An MsThermExperiment object with normalized data slots
#'
#' @examples
#'\dontrun{
#' expt <- normalize_to_std(expt, "cRAP_ALBU_BOVIN")
#'}
#'
#' @export

normalize_to_std <- function( expt, protein ) {

    n_replicates <- length(unlist(lapply(expt$samples,
        "[[", "replicates"),recursive=F))
    n_rows <- 1;
    while (n_rows*(n_rows+1) < n_replicates) {
        n_rows <- n_rows + 1
    }
    par(mfrow=c(n_rows,n_rows))

    n_samples <- length(expt$samples)
        
    for (i_sample in 1:n_samples) {

        sample <- expt$samples[[i_sample]]
        n_replicates <- length(sample$replicates)

        for (i_replicate in 1:n_replicates) {
            
            replicate <- sample$replicates[[i_replicate]]
            expt$samples[[i_sample]]$replicates[[i_replicate]]$data <-
                norm_to_std(replicate,protein)

        }

    }

    return( expt )

}

#' Normalize to a profile
#'
#' Normalizes each replicate of an experiment based on a given spike-in
#' protein standard (assumed to be present in equimolar amounts in each
#' channel)
#'
#' @param expt An MSThermExperiment object
#' @param profile A vector of relative values
#' @param model Whether to fit scale factors to model
#'
#' @return An MsThermExperiment object with normalized data slots
#'
#' @examples
#'\dontrun{
#' expt <- normalize_to_profile(expt, conc, model=T)
#'}
#'
#' @export

normalize_to_profile <- function( expt, profile, model ) {

    n_replicates <- length(unlist(lapply(expt$samples,
        "[[", "replicates"),recursive=F))
    n_rows <- 1;
    while (n_rows*(n_rows+1) < n_replicates) {
        n_rows <- n_rows + 1
    }
    par(mfrow=c(n_rows,n_rows))

    n_samples <- length(expt$samples)
        
    for (i_sample in 1:n_samples) {

        sample <- expt$samples[[i_sample]]
        n_replicates <- length(sample$replicates)

        for (i_replicate in 1:n_replicates) {
            
            replicate <- sample$replicates[[i_replicate]]
            expt$samples[[i_sample]]$replicates[[i_replicate]]$data <-
                norm_to_profile(replicate,profile, model)

        }

    }

    return( expt )

}

gen_description <- function(expt, protein, sep='|') {

    parts <- strsplit(protein,sep,fixed=T)[[1]]
    desc <- sapply(parts, function(x) ifelse( x %in% expt$annot$name,
        expt$annot$annotation[expt$annot$name == x],
        ''))
    return(paste(desc,collapse=sep))

}

#' Model single protein
#'
#' Model a single protein from an MSThermExperiment object
#'
#' @param expt An MSThermExperiment object
#' @param protein ID of the protein to model
#' @param min_rep_psm Minimum number of spectral matches required for each
#'   replicate to model protein
#' @param min_smp_psm Minimum number of spectral matches required for each
#'   sample to model protein
#' @param min_tot_psm Minimum number of spectral matches required across all
#'   replicates to model protein
#' @param max_inf Maximum co-isolation interference level allowed to include a
#'   spectrum in protein-level quantification
#' @param min_score minimum score allowed to include a
#'   spectrum in protein-level quantification
#' @param max_score maximum score allowed to include a
#'   spectrum in protein-level quantification
#' @param smooth (t/F) Perform loess smoothing on the data prior to modeling
#' @param method Protein quantification method to use (see Details)
#' @param method.denom Method used to calculate denominator of abundance
#'   (see Details)
#' @param trim (t/F) Trim all lower data points less than the abundance maximum
#' @param bootstrap (T/F) Perform bootstrap analysis to determine confidence
#'   intervals (slow)
#' @param min_bs_psms Minimum number of spectral matches required to perform
#'   bootstrapping
#' @param merge_reps Treat replicates as overlapping temperature series and
#'   and merge into single series (EXPERIMENTAL!)
#' @param annot_sep Symbol used to separate protein group IDs (used for
#'   retrieval of annotations) (default: '|')
#'
#' @details Valid quantification methods include:
#'   \describe{
#'      \item{"sum"}{Absolute sums for all spectra in each channel are used
#'         to calculate ratios}
#'      \item{"median"}{Median absolute intensities for spectra in each channel are used
#'         to calculate ratios}
#'      \item{"ratio.median"}{Ratios are calculated for each spectrum and
#'         the median ratio for each channel is used}
#'      \item{"ratio.mean"}{Ratios are calculated for each spectrum and
#'         the mean ratio for each channel is used}
#'   }
#' Valid denominator methods include:
#'   \describe{
#'      \item{"first"}{Use the value of the lowest temperature}
#'      \item{"max"}{Use the largest value}
#'      \item{"top3"}{Use the mean of the three largest values}
#'      \item{"near"}{Add description}
#'   }
#'
#' @return MSThermResult object
#'
#' @examples
#'\dontrun{
#' p01  <- model_protein(expt,protein="protein01",smooth=T,bootstrap=F)
#'}
#'
#' @export

model_protein <- function( expt, protein,
  min_rep_psm  = 0,
  min_smp_psm  = 0,
  min_tot_psm  = 0,
  max_inf      = 1,
  min_score,
  max_score,
  smooth       = 0,
  method       = 'sum',
  method.denom = 'near',
  trim         = 0,
  bootstrap    = 0,
  min_bs_psms  = 8,
  merge_reps   = 0,
  annot_sep    = '|'
) {

    self <- structure(
        list(
            name         = protein,
            series       = list(),
            sample_names = c(),
            parameters   = as.list(environment())
        ),
        class = "MSThermResult"
    )
    self$parameters[['expt']] <- NULL

    self$annotation <- gen_description(expt, protein, sep=annot_sep)

    psm_tot   <- 0 # track total PSMs for protein
    n_samples <- length(expt$samples)
        
    for (i_sample in 1:n_samples) {

        psm_smp <- 0 # track total PSMs for sample

        sample <- expt$samples[[i_sample]]
        n_replicates <- length(sample$replicates)
        self$sample_names[i_sample] <- sample$name

        # Track combined data for using "merge_reps"
        merged_profiles <- c()
        merged_temps   <- c()
        merged_splits  <- c(0)
        merged_psms    <- 0

        # Process each replicate
        for (i_replicate in 1:n_replicates) {
            
            replicate <- sample$replicates[[i_replicate]]
           
            # Set score range if not yet defined
            if ("score" %in% colnames(replicate$data)) {
                if (missing(max_score)) {
                    max_score = max(replicate$data$score)
                }
                if (missing(min_score)) {
                    min_score = min(replicate$data$score)
                }
            }

            temps <- replicate$meta$temp

            # Track global min and max temperature points for the protein
            if (is.null(self$tmin) || self$tmin > min(temps)) {
                self$tmin <- min(temps)
            }
            if (is.null(self$tmax) || self$tmax < max(temps)) {
                self$tmax <- max(temps)
            }

            # Pull out matching data points passing coisolation and score
            # thresholds
            sub <- replicate$data[which(replicate$data$protein == protein
                & replicate$data$coelute_inf <= max_inf),];
            if ("score" %in% colnames(sub)) {
                sub <- sub[which( sub$score >= min_score & sub$score <= max_score ),]
            }

            # Reorder channels based on metadata
            quant_columns <- match(replicate$meta$channel,colnames(sub))
            quant <- sub[,quant_columns]

            # Filter out rows with NA or with all zero values
            ok    <- apply(quant,1,function(v) {all(!is.na(v)) & any(v>0)})
            sub   <- sub[ok,]
            quant <- quant[ok,]

            n_psms  <- nrow(sub)

            # Update PSM totals and check cutoffs
            psm_tot <- psm_tot + n_psms
            psm_smp <- psm_smp + n_psms
            if (n_psms < min_rep_psm) {
                return(NULL)
            }

            # Obviously don't try to model if no rows pass filtering
            if (n_psms < 1) {
                next
            }

            profile <- gen_profile(quant,method,method.denom=method.denom)

            # Updated merged variables for use with merge_reps
            merged_profiles <- c(merged_profiles,profile)
            merged_temps    <- c(merged_temps,temps)
            merged_splits   <- c(merged_splits, length(merged_temps))
            merged_psms     <- merged_psms + n_psms

            fit <- try_fit(profile, temps, trim=trim, smooth=smooth)
            fit$is.fitted <- !is.null(fit)

            # calculate weighted co-inf
            sums       <- apply(quant, 1, sum)
            fit$inf    <- sum(sub$coelute_inf * sums) / sum(sums)

            # keep track of other data for later use
            fit$psm    <- n_psms
            fit$name   <- replicate$name
            fit$sample <- sample$name
            fit$x      <- temps
            fit$y      <- profile

            if (fit$is.fitted) {

                #bootstrap if asked
                bs <- c()
                iterations <- 20
                bs.ratios <- matrix(nrow=iterations,ncol=length(profile))
                fit.count <- 0
                if (bootstrap & nrow(quant) >= min_bs_psms) { 
                    for (i in 1:iterations) {
                        
                        quant.bs <- quant[sample(nrow(quant),nrow(quant),replace=T),]
                        profile.bs <- gen_profile(quant.bs,method,method.denom=method.denom)
                        bs.ratios[i,] <- profile.bs
                        fit.bs <- try_fit(profile.bs,temps,trim,smooth)
                        is.fitted <- !is.null(fit.bs)
                        if (is.fitted) {
                            bs[i] <- fit.bs$tm
                            fit.count <- fit.count + 1
                        }
                    }
                    fit$bs.lowers <- apply(bs.ratios,2,function(x) quantile(x,0.025))
                    fit$bs.uppers <- apply(bs.ratios,2,function(x) quantile(x,0.975))
                    if (fit.count > iterations * 0.8) {
                        fit$tm_CI <- quantile(bs,c(0.025,0.975),na.rm=T)
                    }
                }

            }

            # If unfitted, these variables still need to be set but should be NA
            else {
                fit$tm    <- NA
                fit$slope <- NA
                fit$k     <- NA
                fit$plat  <- NA
                fit$r2    <- NA
            }

            self$series[[replicate$name]] <- fit
        }

        # Filter on total PSMs for sample
        if (psm_smp < min_smp_psm) {
            return( NULL )
        }
        
        # Handle merging of replicates (currently ALPHA and intentionally
        # undocumented!)

        # HERE BE DRAGONS!
        if (merge_reps & length(merged_profiles)>0) {
            
            sfs <- c()
            for (i in 1:(length(merged_temps)-1)) {
                for (j in (i+1):length(merged_temps)) {
                    if (merged_temps[i] == merged_temps[j]) {
                        sfs <- c(sfs,merged_profiles[i]/merged_profiles[j])
                    }
                }
            }

            # Here we renormalize the second sample to the first,
            # ASSUMING only two samples of equal length !!!!
            x.merged <- merged_temps
            y.merged <- merged_profiles
            if (length(sfs)>0) {
                #sf <- sum(sfs)/length(sfs)
                sf <- median(sfs)
                e <- length(merged_profiles)
                p <- e/2+1
                merged_profiles[p:e] <- merged_profiles[p:e] * sf
                trim.l <- floor( (length(sfs)+1)/2 )
                trim.r <- floor( (length(sfs)+0)/2 )
                if (trim.l > 0) {
                    merged_profiles <- merged_profiles[-((p-trim.l):(p+trim.r-1))]
                    merged_temps    <- merged_temps[-((p-trim.l):(p+trim.r-1))]
                }
                x.merged <- merged_temps
                y.merged <- merged_profiles

                merged_profiles <- merged_profiles[order(merged_temps)]
                merged_profiles <- abs_to_ratio(merged_profiles,method=method.denom)
                merged_temps    <- merged_temps[order(merged_temps)]
            } 

            fit <- try_fit(merged_profiles,merged_temps,trim=trim,smooth=smooth)
            fit$is.fitted <- !is.null(fit)
            fit$x.merged <- x.merged
            fit$y.merged <- y.merged
            fit$psm      <- merged_psms
            fit$splits   <- merged_splits
            fit$name     <- sample$name
            fit$sample   <- sample$name
            fit$x        <- merged_temps
            fit$y        <- merged_profiles

            if (! fit$is.fitted) {
                fit$tm    <- NA
                fit$slope <- NA
                fit$k     <- NA
                fit$plat  <- NA
                fit$r2    <- NA
            }
            self[['merged']][[sample$name]] <- fit

        }

    }

    # Do final filtering on total PSMs for protein
    if (psm_tot < min_tot_psm) {
        return( NULL )
    }

    # If asked to merge reps, replace individual replicates with merged
    # profiles
    if (merge_reps) {
        self$series <- self$merged
        self$merged <- NULL
    }

    return( self )
}

#' Plot MSThermResultSet object
#'
#' Generate a series of denaturation plots for all results in an MSThermResultSet
#'
#' @param x an MSThermResultSet object
#' @param ... other parameters are passed through to plot.MSThermResult
#' 
#' @details Since this function makes multiple sequential calls to
#'   plot.MSThermResult, it is usually used in conjunction with a multipage
#'   graphics device such as \code{"pdf()"}. Otherwise each subsequent call
#'   will only overwrite the previous output.
#'
#' @return Nothing
#'
#' @examples
#'\dontrun{
#' res  <- model_experiment(expt)
#' plot(res)
#'}
#'
#' @export

plot.MSThermResultSet <- function(x, ...) {

    sapply(x, plot, ...)
    return(invisible(NULL))

}

#' Summarize MSThermResult object
#'
#' Print a summary of an MSThermResult, including samples and parameters
#'
#' @param object an MSThermResult object
#' @param ... additional arguments passed to or from other functions
#'
#' @return Nothing
#'
#' @examples
#'\dontrun{
#' m  <- model_protein(expt,"P38707")
#' summmary(m)
#'}
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

#' Summarize MSThermResultSet object
#'
#' Print a summary of an MSThermResultSet, including samples and parameters
#'
#' @param object an MSThermResultSet object
#' @param ... additional arguments passed to or from other functions
#'
#' @return Nothing
#'
#' @examples
#'\dontrun{
#' res  <- model_experiment(expt)
#' summmary(res)
#'}
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
    
#' Plot MSThermResult object
#'
#' Generate a denaturation plot for an modeled protein/group
#'
#' @param x An MSThermResult object
#' @param table (T/f) include table of per-replicate parameters
#' @param col array of colors used to plot samples
#' @param CI.points (T/F) plot temperature point confidence intervals
#' @param CI.Tm (T/F) plot Tm confidence intervals
#' @param ... other parameters passed through to plot()
#'
#' @return Nothing
#'
#' @examples
#'\dontrun{
#' m  <- model_protein(expt,"P38707")
#' plot(m)
#'}
#'
#' @export

plot.MSThermResult <- function(
    x,
    table=T,
    col,
    CI.points=T,
    CI.Tm=T,
    ...
    ) {

    # object argument must be called 'x' to satisfy S3 rules, but we'll use a
    # more descriptive variable onwards
    result <- x

    if (length(result$series) < 1) {
        return(NULL)
    }

    #library(RColorBrewer, quietly=T)
    #library(plotrix, quietly=T)

    if (missing(col)) {
        col <- brewer.pal(max(3,length(result$sample_names)),"Set1") 
    }
    if (length(col) < length(result$sample_names)) {
        stop("Length of color array must >= sample count")
    }
    colors2 <- paste0(col,"99");
    colors3 <- paste0(col,"33");

    par(mgp = c(2.2, 0.7, 0))
    par(mar = c(4, 4, 3.5, 3.5))
    def_args <- list(
        x    = 0,
        y    = 0,
        xlim = c(result$tmin,result$tmax),
        ylim = c(0,1.3),
        main = result$name,
        xlab=expression(temperature~(degree*C)),
        ylab="relative soluble fraction"
    )
    passed <- list(...)
    do.call("plot",
        c(passed,def_args[!names(def_args) %in% names(passed)])
    )
    mtext(result$annotation)

    for (i_series in 1:length(result$series)) {
        
        series <- result$series[[i_series]]
        if (is.null(series)) { next () }

        i_sample <- which(result$sample_names == series$sample)

        # plot TM confidence intervals if requested
        if (CI.Tm) {
            if (!is.null(series$tm_CI)) {
                rect(series$tm_CI[1],-2,series$tm_CI[2],2,col=colors3[i_sample],border=F)
            }
        }
        if (series$is.fitted) {
            curve(sigmoid(series$plat, series$k, series$tm, x), col=col[i_sample], lwd=2, add=T)
            abline(v=series$tm,col=col[i_sample])
        }

        #merged_splits <- series$splits
        #for (i in 1:(length(merged_splits)-1)) {
            #x <- series$x.merged[(merged_splits[i]+1):merged_splits[i+1]]
            #y <- series$y.merged[(merged_splits[i]+1):merged_splits[i+1]]
            lines(series, lty=2, col=col[i_sample])
            points(series, pch=1, cex=0.8, col=col[i_sample])
        #}

        # plot point confidence intervals if requested
        if (CI.points) {
            if (!is.null(series$bs.lowers)) {
                for (i in 1:length(series$bs.lowers)) {
                    j <- (result$tmax - result$tmin)/100
                    q <- jitter(series$x[i],j)
                    if (series$bs.lowers[i] < series$bs.uppers[i]) {
                        arrows(q, series$bs.lowers[i], q, series$bs.uppers[i], code=3,
                            angle=90, length=0.02,col=colors2[i_sample])
                    }
                    else {
                        lines(c(q - 0.02,q + 0.02),rep(series$bs.lowers[i],2),col=colors2[i_sample])
                    }
                }
            }
        }

    }

    l.dims <- legend("topright",legend=result$sample_names,fill=col,inset=0.02,cex=0.9,bg="white")

    if (table) {

        tm.mean <- suppressWarnings( mean( sapply(result$series, '[[', "tm"), na.rm=T ) )
        lims <- par('usr')
        if (!is.na(tm.mean) & tm.mean < (result$tmax+result$tmin)/2) {
            t.x <- lims[2] - ( lims[2] - lims[1] )*0.02
            t.y <- l.dims$rect$top - l.dims$rect$h - ( lims[4] - lims[3] )*0.02
            just.x <- 1
            just.y <- 0
        }
        else {
            t.x <- lims[1] + ( lims[2] - lims[1] )*0.02
            t.y <- lims[3] + ( lims[4] - lims[3] )*0.02
            just.x <- 0
            just.y <- 1
        }

        l <- length(result$series)*3
        tbl <- matrix(rep(0,length(result$series)*3),ncol=3)
        tbl[,1] <- sapply(result$series, function(x) ifelse(is.null(x$psm),NA,x$psm))
        tbl[,2] <- round(sapply(result$series, function(x) ifelse(is.null(x$tm),NA,x$tm)),1)
        tbl[,3] <- round(sapply(result$series, function(x) ifelse(is.null(x$slope),NA,x$slope)),2)
        colnames(tbl) <- c("PSM",expression(T[m]),"Slp")
        rownames(tbl) <- sapply(result$series, '[[', "name")
        addtable2plot(t.x,t.y,table=tbl,bty="o",lwd=1,hlines=T,xjust=just.x,yjust=just.y,display.rownames=T,xpad=0.4,ypad=1.0,cex=0.7,bg="#FFFFFF77")
    }



}

#' Model MSThermExperiment
#'
#' Model multiple proteins from an MSThermExperiment object
#'
#' @param expt An MSThermExperiment object
#' @param proteins A vector of protein IDs to model (default is all
#'   proteins). 
#' @param np Number of parallel jobs to start (default = number of available
#'   processors)
#' @param ... Parameters passed to model_protein()
#'
#' @return MSThermResultSet object
#'
#' @examples
#'\dontrun{
#' res  <- model_proteins(expt,proteins=c("Q03262","P32582"))
#' res  <- model_proteins(expt,np=1)
#'}
#'
#' @export

model_experiment <- function(expt,proteins,np,...) {

    # parallel processing details
    #suppressMessages(library(doParallel))
    #suppressMessages(library(foreach))
    np <- ifelse( !missing(np), np, detectCores() )
    cl <- makeCluster(np,outfile="")
    registerDoParallel(cl,cores=np)

    if (missing(proteins)) {
        protein_lists <- lapply(expt$samples,function(d)
            unique(sapply(d$replicates, function(l) {l$data$protein})) )
        # proteins <- Reduce(union, protein_lists)
        proteins <- unique(unlist(protein_lists))
    }

    pb <- txtProgressBar(min=0,max=length(proteins),style=3)
    i <- NULL; rm(i) # silence R CMD check noise due to foreach() call below
    results <-
    foreach(i=1:length(proteins),.export=c(
        "model_protein",
        "abs_to_ratio",
        "gen_profile",
        "try_fit",
        "sigmoid",
        "sigmoid.d1",
        "gen_description"
    )) %dopar% {
        protein <- proteins[i]
        setTxtProgressBar(pb,i)
        tryCatch( {
            res <- model_protein(expt,protein,...)
            return(res)
        },error = function(e) {print(e)})
    }
    stopCluster(cl)
    close(pb)
    names(results) <- sapply(results, '[[', "name")
    self <- structure(
        results[!sapply(results,is.null)],
        class = "MSThermResultSet"
    )
    return( self )

}

#' Export MSThermResultSet to an SQLite database
#'
#' Exports and MSThermResultSet object to a new SQLite database file.
#' Each model (specific to a given replicate and protein) is exported as an
#' individual record. The schema used for the 'data' table can be seen in the
#' code below.
#'
#' @param res An MSThermResultSet object
#' @param file Path to the output sqlite database to be created
#'
#' @return Nothing
#'
#' @examples
#'\dontrun{
#' write.sqlite(res, "/path/to/db.sqlite")
#'}
#'
#' @export

write.sqlite <- function( res, file ) {

    #library("RSQLite")

    if (requireNamespace("RSQLite")) {

        con <- RSQLite::dbConnect(RSQLite::SQLite(), file)

        create <- "CREATE TABLE data (
            protein TEXT NOT NULL,
            replicate TEXT NOT NULL,
            sample TEXT NOT NULL,
            m REAL NOT NULL,
            k REAL NOT NULL,
            p REAL NOT NULL,
            x BLOB NOT NULL,
            y BLOB NOT NULL,
            psm INTEGER NOT NULL,
            inf REAL NOT NULL,
            slope REAL NOT NULL,
            r2 REAL NOT NULL,
            el BLOB,
            eu BLOB,
            ml REAL,
            mu REAL,
            PRIMARY KEY (protein, replicate)
        );"

        rs <- RSQLite::dbSendQuery(con, create)
        RSQLite::dbClearResult(rs)

        for (p in res) {

            for (r in p$series) {

                # skip un-modeled series
                if (! r$is.fitted) {
                    next
                }

                has_pci <- ! is.null(r$bs.lowers)
                has_tci <- ! is.null(r$tm_CI)

                # even though last two values are floats, treat as strings so that
                # 'NULL' can be substituted when necessary
                sql <- sprintf("INSERT INTO data VALUES ('%s', '%s', '%s', %f, %f,
                    %f, '%s', '%s', %d, %f, %f, %f, %s, %s, %s, %s )",
                    p$name,
                    r$name,
                    r$sample,
                    r$tm,
                    r$k,
                    r$plat,
                    x_to_str(r$x),
                    y_to_str(r$y),
                    r$psm,
                    r$inf,
                    r$slope,
                    r$r2,
                    if(has_pci) paste0("'",y_to_str( r$bs.lowers ),"'") else 'NULL',
                    if(has_pci) paste0("'",y_to_str( r$bs.uppers ),"'") else 'NULL',
                    if(has_tci) r$tm_CI[[1]] else 'NULL',
                    if(has_tci) r$tm_CI[[2]] else 'NULL'
                )

                rs <- RSQLite::dbSendQuery(con, sql)
                RSQLite::dbClearResult(rs)
            }
        }

        RSQLite::dbDisconnect(con)

        return(invisible(NULL))

    }
    else {
        warning("RSQLite not installed so this functionality is unavailable", call.=F) 
        return(FALSE)
    }

}

try_fit <- function(ratios,temps,trim,smooth) {

    #library(nls2, quietly=T)

    x <- temps
    y <- ratios

    if (!missing(trim) & trim) {
        x <- temps[which.max(temps):length(temps)]
        y <- ratios[which.max(temps):length(temps)]
    }

    if (smooth) {
        f <- loess(y ~ x, span=0.6)
        y <- f$fitted
    }

    fit <- list()

    st.coarse <- expand.grid(p=c(0,0.3),k=seq(0,4000,by=1000),m=seq(30,60,by=15))
    st.fine   <- expand.grid(p=c(0,0.3),k=seq(0,8000,by=200),m=seq(30,80,by=10))
    for (st in list(st.coarse,st.fine)) {
        tryCatch( {
            mod <- nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=st,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=5000))
            fit <-
            nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=F),algorithm="port",lower=c(0,1,10),upper=c(0.4,100000,100))
            obj <- list()
            obj$plat <- coefficients(fit)[1]
            obj$k <- coefficients(fit)[2]
            obj$tm <- coefficients(fit)[3]
            obj$slope <- sigmoid.d1(obj$plat,obj$k,obj$tm,obj$tm)
            y.fit <- sigmoid(obj$plat,obj$k,obj$tm,temps)
            obj$y.fit <- y.fit
            obj$resid <- ratios - y.fit
            obj$r2 <- 1-(sum(obj$resid^2)/(length(ratios)*var(ratios)))
            obj$rmsd <- sqrt( sum(obj$resid^2)/length(ratios) )
            return(obj)
        },error = function(e) {})
    }
    return(NULL)

}


sigmoid <- function(p,k,m,x) {

    (1-p)/(1+exp(-k*(1/x-1/m)))+p

}

sigmoid.d1 <- function(p,k,m,x) {

    -((1 - p) * (exp(-k * (1/x - 1/m)) * (k * (1/x^2)))/(1 + exp(-k * (1/x - 1/m)))^2)

}

norm_to_profile <- function(replicate,profile,model=T) {

    #library(RColorBrewer, quietly=T)
    #library(nls2, quietly=T)
    cols <- brewer.pal(3,"Set1")

    quant_columns <- match(replicate$meta$channel,colnames(replicate$data))

    temps <- replicate$meta$temp
    quant <- replicate$data[,quant_columns]

    std.ratios <- profile/profile[1]
    sums <- apply(as.matrix(quant),2,function(x){sum(as.numeric(x),na.rm=T)})
    ratios <- sums/sums[1]

    if(model) {
        st.coarse <- expand.grid(p=seq(0,0.4,by=0.1),k=seq(0,2000,by=100),m=seq(20,80,by=5))
        x <- temps[std.ratios < 1.2]
        y <- std.ratios[std.ratios < 1.2]
        mod <- nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=st.coarse,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=50000))
        fit <- nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=F,maxiter=50000),algorithm="port",lower=c(0,1,1),upper=c(0.5,2000,100))
        p <- coefficients(fit)[['p']]
        k <- coefficients(fit)[['k']]
        m <- coefficients(fit)[['m']]
        yfit <- sigmoid(p,k,m,temps)

        sf <- yfit/ratios
    }
    else {
        sf <- std.ratios/ratios
    }
    norm <- replicate$data
    norm[,quant_columns[order(replicate$meta$temp)] ] <- t(t(as.matrix(quant))*sf)
    corrected <- ratios*sf

    plot( temps,ratios,ylim=c(0,max(c(ratios,std.ratios))),xlab=expression(paste("temperature ", degree, "C")),ylab="relative fraction remaining",col=cols[3],main=replicate$name)
    points(temps,std.ratios,col=cols[1])
    points(temps,corrected,col=cols[2])
    if (model) {
        curve(sigmoid(p,k,m,x),col=cols[2],add=T)
    }

    return(norm)

}


norm_to_std <- function(replicate,protein) {

    #library(RColorBrewer, quietly=T)
    #library(nls2, quietly=T)
    cols <- brewer.pal(3,"Set1")

    quant_columns <- match(replicate$meta$channel,colnames(replicate$data))

    temps <- replicate$meta$temp
    quant <- replicate$data[,quant_columns]

    std <- matrix(nrow=0,ncol=ncol(quant)) 
    for (p in protein) {
        tmp <- quant[which(replicate$data$protein == p),]
        std <- rbind(std, tmp)
    }
    std.sums <- apply(std,2,sum)
    std.ratios <- std.sums/std.sums[1]
    sums <- apply(as.matrix(quant),2,function(x){sum(as.numeric(x),na.rm=T)})
    ratios <- sums/sums[1]

    corrected <- ratios/std.ratios

    st.coarse <- expand.grid(p=seq(0,0.4,by=0.1),k=seq(0,2000,by=100),m=seq(20,80,by=5))
    x <- temps[corrected < 1.2]
    y <- corrected[corrected < 1.2]
    mod <- nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=st.coarse,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=50000))
    fit <- nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=F,maxiter=50000),algorithm="port",lower=c(0,1,1),upper=c(0.5,2000,100))
    p <- coefficients(fit)[['p']]
    k <- coefficients(fit)[['k']]
    m <- coefficients(fit)[['m']]
    yfit <- sigmoid(p,k,m,temps)

    sf <- yfit/ratios
    norm <- replicate$data
    norm[,quant_columns[order(replicate$meta$temp)] ] <- t(t(as.matrix(quant))*sf)

    plot( temps,ratios,ylim=c(0,max(c(ratios,std.ratios))),xlab=expression(paste("temperature ", degree, "C")),ylab="relative fraction remaining",col=cols[3],main=replicate$name)
    points(temps,std.ratios,col=cols[1])
    for (i in 1:ncol(std)) {
        v <- c(1:nrow(std))
        for (j in 1:nrow(std)) {
            v[j] <- std[j,i]/std[j,1]
        }
        if (quantile(v,.25,na.rm=T) < quantile(v,.75,na.rm=T)) {
            arrows(temps[i],quantile(v,.25,na.rm=T),temps[i],quantile(v,.75,na.rm=T),code=3,angle=90,length=0.02,col=cols[1])
        }
    }
    points(temps,corrected,col=cols[2])
    curve(sigmoid(p,k,m,x),col=cols[2],add=T)

    return(norm)

}

str_to_y <- function (s) {

    y.conv <- as.integer(charToRaw(s))
    y.orig <- (y.conv-127)/64
    return(y.orig)

}

y_to_str <- function (v) {

    y.conv <- vapply(v, function(y) {
        min(255, round(y*64+127, 0))
    }, double(1))
    bs <- rawToChar(as.raw(y.conv))
    return(bs)

}

x_to_str <- function (v) {

    # check input range
    if (min(v) < 20 | max(v) > 83) {
        stop("temps must be between 20 and 83 oC")
    }

    x.conv <- vapply(v, function(x) {
        round( (x-19)/64*255, 0 )
    }, double(1))
    bs <- rawToChar(as.raw(x.conv))
    return(bs)

}
    
str_to_x <- function (s) {

    x.conv <- as.integer(charToRaw(s))
    x.orig <- x.conv/255*64+19;
    return(x.orig)

}
