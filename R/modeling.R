#' @title Model single protein.
#'
#' @description Model a single protein from an MSThermExperiment object.
#'
#' @param expt An MSThermExperiment object
#' @param protein ID of the protein to model
#' @param min_rep_psm Minimum number of spectral matches required for each
#'   replicate to model protein
#' @param min_smp_psm Minimum number of spectral matches required for each
#'   sample to model protein
#' @param min_tot_psm Minimum number of spectral matches required across all
#'   replicates to model protein
#' @param min_r2 Minimum R2 value to consider model (implies "only_modeled")
#' @param max_slope Maximum slope to consider model (implies "only_modeled")
#' @param min_reps Minimum number of modeled replicates for each sample to
#'   return protein
#' @param max_inf Maximum co-isolation interference level allowed to include a
#'   spectrum in protein-level quantification
#' @param min_score minimum score allowed to include a
#'   spectrum in protein-level quantification
#' @param max_score maximum score allowed to include a
#'   spectrum in protein-level quantification
#' @param only_modeled (t/F) Only consider modeled proteins
#' @param check_missing (t/F) Run simple test to filter out PSMs with missing
#'   quantification channels where values are expected
#' @param missing_cutoff Minimum fraction relative to surrounding data points
#'   used in the check for missing channels
#' @param smooth (t/F) Perform loess smoothing on the data prior to modeling
#' @param method Protein quantification method to use (see Details)
#' @param method.denom Method used to calculate denominator of abundance
#'   (see Details)
#' @param trim (t/F) Trim all lower data points less than the abundance maximum
#' @param bootstrap (T/F) Perform bootstrap analysis to determine confidence
#'   intervals (slow)
#' @param min_bs_psms Minimum number of spectral matches required to perform
#'   bootstrapping
#' @param annot_sep Symbol used to separate protein group IDs (used for
#'   retrieval of annotations) (default: '|')
#'
#' @details Valid quantification methods include:
#'     \describe{
#'        \item{"sum"}{use the sum of the spectrum values for each channel}
#'        \item{"median"}{use the median of the spectrum values for each
#'          channel}
#'        \item{"ratio.median"}{Like "median", but values for each spectrum
#'          are first converted to ratios according to "method.denom"
#'          channel}
#'        \item{"ratio.mean"}{Like "ratio.median" but using mean of ratios}
#'     }
#' Valid denominator methods include:
#'     \describe{
#'        \item{"first"}{Use the first value (lowest temperature point)
#'          (default)}
#'        \item{"max"}{Use the maximum value}
#'        \item{"top3"}{Use the mean of the three highest values}
#'        \item{"near"}{Use the median of all values greater than 80% of
#'          the first value}
#'     }
#'
#' @return MSThermResult object
#'
#' @examples
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#'
#' model   <- model_protein(expt, "P38707", smooth=TRUE, bootstrap=FALSE)
#' summary(model)
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
  annot_sep    = '|',
  max_slope    = 0,
  min_r2       = 0,
  min_reps     = 0,
  only_modeled = 0,
  check_missing = 0,
  missing_cutoff = 0.3
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

    if (min_r2 > 0 || max_slope < 0 ) {
        only_modeled = 1
    }
        
    for (i_sample in 1:n_samples) {

        psm_smp <- 0 # track total PSMs for sample
        worst_slope <- -1
        worst_r2    <- 1

        sample <- expt$samples[[i_sample]]
        n_replicates <- length(sample$replicates)
        self$sample_names[i_sample] <- sample$name

        n_reps <- 0

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

            if (nrow(sub) < 1) {
                next
            }

            # Reorder channels based on metadata
            quant_columns <- match(replicate$meta$channel,colnames(sub))
            quant <- sub[,quant_columns]

            # Filter out rows with NA or with all zero values
            ok <- apply(quant,1,function(v) {
                all(!is.na(v)) &
                  any(v>0)     &
                  ( (! check_missing) || is_consistent(v, missing_cutoff) )
            })
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
                if (fit$r2 < min_r2) {
                    next;
                }
                if (fit$slope > max_slope) {
                    next;
                }
                n_reps = n_reps + 1

            }

            # If unfitted, these variables still need to be set but should be NA
            else {
                if (only_modeled) {
                    next;
                }
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
        # Filter on minimum modeled replicates for sample
        if (n_reps < min_reps) {
            return( NULL )
        }

    }

    # Do final filtering on total PSMs for protein
    if (psm_tot < min_tot_psm) {
        return( NULL )
    }
    # Filter on worst slope
    if (worst_slope > max_slope) {
        return( NULL )
    }
    # Filter on worst R2
    if (worst_r2 < min_r2) {
        return( NULL )
    }

    return( self )
}



#' @title Model MSThermExperiment.
#'
#' @description Model multiple proteins from an MSThermExperiment object.
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
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#'
#' res     <- model_experiment(expt, bootstrap=FALSE, np=2)
#' summary(res)
#'
#' @export

model_experiment <- function(expt,proteins,np,...) {

    # parallel processing details
    np <- ifelse( !missing(np), np, detectCores() )
    registerDoParallel(cores=np)

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
        if (getTxtProgressBar(pb) < i) {
            setTxtProgressBar(pb,i)
        }
        tryCatch( {
            res <- model_protein(expt,protein,...)
            return(res)
        },error = function(e) {print(e)})
    }
    close(pb)
    names(results) <- sapply(results, '[[', "name")
    self <- structure(
        results[!sapply(results,is.null)],
        class = "MSThermResultSet"
    )

    stopImplicitCluster()

    return( self )

}
   


#' @title Convert absolute quantitation to relative ratios.
#'
#' @description \code{abs_to_ratio} takes a vector of absolute values and
#' returns a vector of ratios relative to some starting point.
#'
#' @param x vector of numeric absolute quantitation values
#' @param method method to use to determine starting value (denominator)
#'
#' @details The denominator used to calculate relative protein concentrations
#'   can affect the ability to model noisy data. In the theoretically ideal
#'   scenario, everything would be relative to the lowest temperature point.
#'   However, other methods can be used to help alleviate problems related to
#'   noise. Available methods include:
#'     \describe{
#'        \item{"first"}{Use the first value (lowest temperature point)
#'          (default)}
#'        \item{"max"}{Use the maximum value}
#'        \item{"top3"}{Use the mean of the three highest values}
#'        \item{"near"}{Use the median of all values greater than 80% of
#'          the first value}
#'     }
#'
#' @return A numeric vector of the same length as input
#' @keywords internal

abs_to_ratio <- function(x, method='first') {

    denom <- switch( method,
        first  = x[1],
        max    = max(x),
        top3   = mean( x[ order(x,decreasing=T)[1:3] ] ),
        near   = median( x[ x > x[1]*0.8 ] ),
        #compat = {
                #m <- mean( x[ order(x,decreasing=T)[1:3] ] )
                #b <- mean( x[ x > m*0.8 & x < m*1.2 ] )
                #ifelse(is.na(b),m,b)
            #},
        stop("Invalid method type",call.=T)
    )

    return( x/denom )

}



#' @title Generate protein ratio profile from spectrum quantification matrix.
#'
#' @description \code{gen_profile} takes a matrix of spectrum channel
#' quantification values belonging to a protein and "rolls them up" into a
#' vector of protein-level relative quantification values.
#'
#' @param x matrix of spectrum quantification values, one row per spectrum and
#' one column per channel
#' @param method method to use to "roll up" spectrum values to protein level
#' @param method.denom method used to determine ratio denominator, passed as
#' the "method" argument to \code{abs_to_ratio}
#'
#' @details The following methods for spectrum-to-protein conversion are
#' supported:
#'     \describe{
#'        \item{"sum"}{use the sum of the Spectrum values for each channel}
#'        \item{"median"}{use the median of the spectrum values for each
#'          channel}
#'        \item{"ratio.median"}{Like "median", but values for each spectrum
#'          are first converted to ratios according to "method.denom"
#'          channel}
#'        \item{"ratio.mean"}{Like "ratio.median" but using mean of ratios}
#'     }
#'
#' @return A numeric vector of the same length as the number of matrix columns
#' @keywords internal

gen_profile <- function( x, method='sum', method.denom='first' ) {

    summarized <- switch( method,
        sum    = apply(x, 2, sum),
        median = apply(x, 2, median),
        # for these, the inner apply() will transpose, so the outer is applied
        # to rows rather than columns
        ratio.median = 
            apply( apply(x,1,abs_to_ratio,method=method.denom),1,median ),
        ratio.mean = 
            apply( apply(x,1,abs_to_ratio,method=method.denom),1,mean ),
        stop("Invalid method type", call.=T)
    )
    return( abs_to_ratio(summarized, method=method.denom) )

}



# Perform the actual model fitting
try_fit <- function(ratios,temps,trim,smooth) {

    x <- temps
    y <- ratios

    if (!missing(trim) & trim) {
        x <- temps[which.max(temps):length(temps)]
        y <- ratios[which.max(temps):length(temps)]
    }

    if (smooth) {
        f <- loess(y ~ x, span=0.65)
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
            obj$plat  <- as.numeric(coefficients(fit)[1])
            obj$k     <- as.numeric(coefficients(fit)[2])
            obj$tm    <- as.numeric(coefficients(fit)[3])
            obj$slope <- as.numeric(sigmoid.d1(obj$plat,obj$k,obj$tm,obj$tm))
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



# The model
sigmoid <- function(p,k,m,x) {

    (1-p)/(1+exp(-k*(1/x-1/m)))+p

}



# first derivative
sigmoid.d1 <- function(p,k,m,x) {

    -((1 - p) * (exp(-k * (1/x - 1/m)) * (k * (1/x^2)))/(1 + exp(-k * (1/x - 1/m)))^2)

}

# look for probable missing values
is_consistent <- function(v,cutoff=0.3) {

    len <- length(v)

    if (len < 2) {
        return(1)
    }

    for (i in 1:(len-1)) {
        if (i == 1 & (v[i]/v[i+1]) < cutoff) {
            return(0)
        }
        if (i > 1) {
            if ((v[i]/v[i-1]) < cutoff & (v[i]/v[i+1]) < cutoff) {
                return(0)
            }
        }
    }

    return(1)

}
