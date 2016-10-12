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
#' expt  <- MSThermExperiment(control="expt_1.control",annotations="protein_annots.tsv")
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
        annotations <- .to_dataframe(annotations)
        if (! is.data.frame(annotations)) {
            stop("annotations must be a valid filename or data.frame")
        }
        self$annot <- annotations
    }

    # return empty object if no control file/dataframe specified
    if (missing(control)) { return(self) }

    control <- .to_dataframe(control)
    if (! is.data.frame(control)) {
        stop("control must be a valid filename or data.frame")
    }

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
                        meta = control$meta_file[control$name == r]
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

MSThermReplicate <- function(name, data, meta) {

    data <- .to_dataframe(data)
    meta <- .to_dataframe(meta)
    meta <- meta[order(meta$temp,decreasing=F),]

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
    curve(dnorm(x,mean=m,sd=d),col="red",add=T)
    return(q)

} 
    
#' MSResultSet to data frame
#'
#' Populates a dataframe with information from an MSResultSet, one row per
#' protein/group
#'
#' @param set An MSResultSet object
#'
#' @return A dataframe populated with relevant information per result
#'
#' @examples
#' df  <- as.data.frame(expt)
#' write.table(df, "results.tsv")
#'
#' @export

as.data.frame.MSThermResultSet <- function( set ) {

    df <- data.frame(
        row.names = sapply(set, '[[', "name")
    )
    df$annotation <- sapply(set, '[[', "annotation")

    repl_lists <- lapply(set, function(d)
        unique(sapply( d$series, '[[', "name" )))
    repl_names <- unique(unlist(repl_lists))
    repl_names <- repl_names[order(repl_names)]

    for (r in repl_names) {
        x <- set[1]$series[[r]][['x']]

        for(col in c("tm","psm","inf","slope","k","plat","r2","rmsd")) {
            df[[paste0(r,'.',col)]]  <- sapply(set, function(v) v$series[[r]][[col]])
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
#' expt <- normalize_to_tm(expt, res)
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
#' @param gname Name of a protein to normalize against
#'
#' @return An MsThermExperiment object with normalized data slots
#'
#' @examples
#' expt <- normalize_to_std(expt, "bovine_serum_albumin")
#'
#' @export

normalize_to_std <- function( expt, gname, ... ) {

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
                norm_to_std(replicate,gname)

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
#' expt <- normalize_to_profile(expt, conc, model=T)
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

gen_description <- function(expt, gname, sep='|') {

    parts <- strsplit(gname,sep,fixed=T)[[1]]
    desc <- sapply(parts, function(x) ifelse( x %in% expt$annot$name,
        expt$annot$annotation[expt$annot$name == x],
        ''))
    return(paste(desc,collapse=sep))

}

model_gene <- function( expt, gname,
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
  annot_sep    = '|'
) {

    self <- structure(
        list(
            name         = gname,
            series       = list(),
            sample_names = c(),
            parameters   = as.list(environment())
        ),
        class = "MSThermResult"
    )
    self$parameters[['expt']] <- NULL

    desc <- gen_description(expt, gname, sep=annot_sep)
    self$annotation <- desc

    res_cols_per_repl <- 28
    tbl_cols_per_repl <- 6

    replicate_total <- 0
    psm_tot <- 0
    n_samples <- length(expt$samples)
        
    for (i_sample in 1:n_samples) {

        psm_smp <- 0

        sample <- expt$samples[[i_sample]]
        n_replicates <- length(sample$replicates)
        self$sample_names[i_sample] <- sample$name

        for (i_replicate in 1:n_replicates) {

            
            replicate <- sample$replicates[[i_replicate]]
            replicate_total <- replicate_total + 1
            
            if ("score" %in% colnames(replicate$data)) {
                if (missing(max_score)) {
                    max_score = max(replicate$data$score)
                }
                if (missing(min_score)) {
                    min_score = min(replicate$data$score)
                }
            }

            temps <- replicate$meta$temp
            if (is.null(self$tmin) || self$tmin > min(temps)) {
                self$tmin <- min(temps)
            }
            if (is.null(self$tmax) || self$tmax < max(temps)) {
                self$tmax <- max(temps)
            }

            sub <- replicate$data[which(replicate$data$protein == gname
                & replicate$data$coelute_inf <= max_inf),];
            if ("score" %in% colnames(sub)) {
                sub <- sub[which( sub$score >= min_score & sub$score <= max_score ),]
            }

            quant_columns <- match(replicate$meta$channel,colnames(sub))
            quant <- sub[,quant_columns]

            ok <- apply(quant,1,function(v) {all(!is.na(v)) & any(v>0)})
            sub <- sub[ok,]
            quant <- quant[ok,]

            psm_tot <- psm_tot + nrow(sub)
            psm_smp <- psm_smp + nrow(sub)
            if (nrow(sub) < min_rep_psm) {
                return(NULL)
            }
            if (nrow(sub) < 1) {
                next
            }

            profile <- gen_profile(quant,method,method.denom=method.denom)
            fit <- try_fit(profile,temps,trim=trim,smooth=smooth)
            fit$is.fitted <- !is.null(fit)

            # calculate weighted co-inf
            sums       <- apply(quant,1,sum)
            fit$inf    <- sum(sub$coelute_inf * sums) / sum(sums)
            fit$psm    <- nrow(sub)
            fit$name   <- replicate$name
            fit$sample <- sample$name
            fit$x      <- temps
            fit$y      <- profile

            if (fit$is.fitted) {

                #bootstrap
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
            else {
                fit$tm <- NA
                fit$slope <- NA
                fit$k <- NA
                fit$plat <- NA
                fit$r2 <- NA
            }

            self$series[[replicate$name]] <- fit
        }

        if (psm_smp < min_smp_psm) {
            return( NULL )
        }

    }

    if (psm_tot < min_tot_psm) {
        return( NULL )
    }
    return( self )
}

#' Plot MSThermResultSet object
#'
#' Generate a series of denaturation plots for all results in an MSThermResultSet
#'
#' @param set An MSThermResultSet object
#' @param ... Other parameters are passed through to plot.MSThermResult
#' 
#' @details Since this function makes multiple sequential calls to
#'   plot.MSThermResult, it is usually used in conjunction with a multipage
#'   graphics device such as \code{"pdf()"}. Otherwise each subsequent call
#'   will only overwrite the previous output.
#'
#' @return Nothing
#'
#' @examples
#' res  <- model_genes(expt)
#' plot(res)
#'
#' @export

plot.MSThermResultSet <- function(set,...) {

    sapply(set, plot, ...)
    return(invisible(NULL))

}

#' Summarize MSThermResult object
#'
#' Print a summary of an MSThermResult, including samples and parameters
#'
#' @param result An MSThermResult object
#'
#' @return Nothing
#'
#' @examples
#' m  <- model_gene(expt,"P12345")
#' summmary(m)
#'
#' @export

summary.MSThermResult <- function(result) {

    cat(paste("Name:",result$name),"\n",sep='')
    cat("Samples:\n");
    for (i in result$sample_names) {
        cat("    ",i,"\n",sep='')
    }
    cat("Parameters:\n");
    for (i in names(result$parameters)) {
        cat("    ",i,": ",result$parameters[[i]],"\n",sep='')
    }

    return(invisible(NULL))

}

#' Summarize MSThermResultSet object
#'
#' Print a summary of an MSThermResultSet, including samples and parameters
#'
#' @param set An MSThermResultSet object
#'
#' @return Nothing
#'
#' @examples
#' res  <- model_genes(expt)
#' summmary(res)
#'
#' @export

summary.MSThermResultSet <- function(set) {

    result <- set[[1]]
    cat("Samples:\n");
    for (i in result$sample_names) {
        cat("    ",i,"\n",sep='')
    }
    cat("Parameters:\n");
    for (i in names(result$parameters)) {
        if (i != 'gname') {
            cat("    ",i,": ",result$parameters[[i]],"\n",sep='')
        }
    }

    return(invisible(NULL))

}
    
#' Plot MSThermResult object
#'
#' Generate a denaturation plot for an modeled protein/group
#'
#' @param result An MSThermResult object
#' @param table (T/F) Include table of per-replicate parameters
#' @param col Array of colors used to plot samples
#' @param CI.points (T/F) Plot temperature point confidence intervals
#' @param CI.Tm (T/F) Plot Tm confidence intervals
#'    intervals
#'
#' @return Nothing
#'
#' @examples
#' m  <- model_gene(expt,"P12345")
#' plot(m)
#'
#' @export

plot.MSThermResult <- function(result,
    table=T,
    col,
    CI.points=T,
    CI.Tm=T,
    ...) {

    if (length(result$series) < 1) {
        return(NULL)
    }

    library(RColorBrewer, quietly=T)
    library(plotrix, quietly=T)

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

        lines(series, lty=2, col=col[i_sample])

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
        tbl[,3] <- round(sapply(result$series, function(x) ifelse(is.null(x$inf),NA,x$inf)),2)
        colnames(tbl) <- c("PSM",expression(T[m]),"Inf")
        rownames(tbl) <- sapply(result$series, '[[', "name")
        addtable2plot(t.x,t.y,table=tbl,bty="o",lwd=1,hlines=T,xjust=just.x,yjust=just.y,display.rownames=T,xpad=0.4,ypad=1.0,cex=0.7,bg="#FFFFFF77")
    }


}

#' Model MSThermExperiment
#'
#' Model one or more proteins in an MSThermExperiment object
#'
#' @param expt An MSThermExperiment object
#' @param genes A vector of gene/protein IDs to model (default is all genes). 
#' @param np Number of parallel jobs to start (default = number of available
#'   processors)
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
#' @param smooth (T/F) Perform loess smoothing on the data prior to modeling
#' @param method Protein quantification method to use (see Details)
#' @param method.denom Method used to calculate denominator of abundance
#'   (see Details)
#' @param bootstrap (T/F) Perform bootstrap analysis to determine confidence
#'   intervals (slow)
#' @param min_bs_psms Minimum number of spectral matches required to perform
#'   bootstrapping
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
#' @return MSThermResultSet object
#'
#' @examples
#' res  <- model_genes(expt,genes=c("gene01","gene02"))
#' res  <- model_genes(expt,np=4)
#'
#' @export

model_genes <- function(expt,genes,np,...) {

    # parallel processing details
    suppressMessages(library(doParallel))
    suppressMessages(library(foreach))
    np <- ifelse( !missing(np), np, detectCores() )
    cl <- makeCluster(np,outfile="")
    registerDoParallel(cl,cores=np)

    if (!missing(genes)) {
        proteins <- genes
    }
    else {
        protein_lists <- lapply(expt$samples,function(d)
            unique(sapply(d$replicates, function(l) {l$data$protein})) )
        # proteins <- Reduce(union, protein_lists)
        proteins <- unique(unlist(protein_lists))
    }

    pb <- txtProgressBar(min=0,max=length(proteins),style=3)
        
    results <-
    foreach(i=1:length(proteins),.export=c(
        "model_gene",
        "abs_to_ratio",
        "gen_profile",
        "try_fit",
        "sigmoid",
        "sigmoid.d1",
        "gen_description"
    )) %dopar% {
        gene <- proteins[i]
        setTxtProgressBar(pb,i)
        tryCatch( {
            res <- model_gene(expt,gene,...)
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

try_fit <- function(ratios,temps,trim,smooth) {

    library(nls2, quietly=T)

    x <- temps
    y <- ratios

    if (!missing(trim) & trim > 1) {
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

    library(RColorBrewer, quietly=T)
    library(nls2, quietly=T)
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


norm_to_std <- function(replicate,gene) {

    library(RColorBrewer, quietly=T)
    library(nls2, quietly=T)
    cols <- brewer.pal(3,"Set1")

    quant_columns <- match(replicate$meta$channel,colnames(replicate$data))

    temps <- replicate$meta$temp
    quant <- replicate$data[,quant_columns]

    std <- matrix(nrow=0,ncol=ncol(quant)) 
    for (g in gene) {
        tmp <- quant[which(replicate$data$protein == g),]
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

