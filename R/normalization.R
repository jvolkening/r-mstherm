#' @title Re-normalize based on Tm.
#'
#' @description Normalizes each replicate of an experiment based on linear
#' regression of calculated Tm (corrects for remaining systematic error).
#'
#' @details An assumption can be made in most TPP experiments using a single
#' organism that the Tm of most proteins should not be changing. However,
#' global shifts have been observed between replicates, presumably due to
#' technical variance, which complicate data interpretation. This method
#' attempts to remove this source of error by doing a bootstrap
#' renormalization of the quantification values based on pairwise linear
#' regression between calculated Tms of each replicate. A reference set of Tms
#' is calculated based on all replicates, and each replicate is normalized to
#' this based on the calculated slope and intercept of the input data.
#'
#' @param expt An MSThermExperiment object
#' @param res An MSThermResultSet object
#' @param mode Regression mode
#' @param min_r2 Minimum R2 for protein to use Tm in regression
#' @param max_slope Maximum slope for protein to use Tm in regression
#'
#' @return An MsThermExperiment object with re-normalized data slots
#'
#' @examples
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#' res     <- model_experiment(expt, smooth=TRUE, bootstrap=FALSE, np=2)
#'
#' expt    <- normalize_to_tm(expt, res)
#'
#' @export

normalize_to_tm <- function( expt, res, mode=1, min_r2=0.95, max_slope=-0.06 ) {

    df.tm <- data.frame(row.names = sapply(res, "[[", "name"))
    df.r2 <- data.frame(row.names = sapply(res, "[[", "name"))
    df.sl <- data.frame(row.names = sapply(res, "[[", "name"))
    repl_lists <- lapply(res, function(d) unique(sapply(d$series, 
        "[[", "name")))
    repl_names <- unique(unlist(repl_lists))
    repl_names <- repl_names[order(repl_names)]
    for (r in repl_names) {
        s <- sapply(res, function(v) v$series[[r]][["tm"]])
        s[sapply(s, is.null)] <- NA
        df.tm[[r]] <- unlist(s, use.names = F)
        s <- sapply(res, function(v) v$series[[r]][["r2"]])
        s[sapply(s, is.null)] <- NA
        df.r2[[r]] <- unlist(s, use.names = F)
        s <- sapply(res, function(v) v$series[[r]][["slope"]])
        s[sapply(s, is.null)] <- NA
        df.sl[[r]] <- unlist(s, use.names = F)
    }
    min.r2  <- apply(df.r2, 1, min)
    max.sl  <- apply(df.sl, 1, max)
    df.tm <- df.tm[(min.r2 > min_r2) & !is.na(min.r2) & (max.sl < max_slope) & !is.na(max.sl), ]
    mean.tm <- apply(df.tm, 1, function(v) mean(v,na.rm=T))
    print(paste("rows used:",nrow(df.tm)))
    if (length(repl_names) < 3) {
        baseline <- repl_names[1]
    }
    else {
        meds <- apply(df.tm, 1, median)
        diffsums = c()
        c <- 1
        for (r in repl_names) {
            diffsums[c] <- sum(abs(df.tm[, c] - meds))
            c <- c + 1
        }
        baseline <- repl_names[which.min(diffsums)]
        print(paste("baseline:",baseline))
    }
    #tm.b <- sapply(res, function(v) v$series[[baseline]][["tm"]])
    #tm.b[sapply(tm.b, is.null)] <- NA
    #tm.b <- unlist(tm.b)
    sample_names <- sapply(expt$samples, "[[", "name")
    tm.b <- df.tm[[baseline]]

    if (mode == 3) {
        tm.b <- tm.b[ mean.tm < quantile(mean.tm, 0.2)
                   | mean.tm > quantile(mean.tm, 0.8) ]
        print(paste("final count:",length(tm.b)))
    }

    for (r in repl_names) {
        if (r == baseline) {
            next
        }
        #tm.u <- sapply(res, function(v) v$series[[r]][["tm"]])
        #tm.u[sapply(tm.u, is.null)] <- NA
        #tm.u <- unlist(tm.u)
        tm.u <- df.tm[[r]]
        if (mode == 3) {
            tm.u <- tm.u[ mean.tm < quantile(mean.tm, 0.2)
                    | mean.tm > quantile(mean.tm, 0.8) ]
        }
        if (mode == 1 | mode == 3) {
            print("Mode 1/3 reg")
            l <- lm(tm.b ~ tm.u)
            tm.u <- tm.u * l$coefficients[2] + l$coefficients[1]
            l2 <- lm(tm.b ~ tm.u)
        }
        else if (mode == 2) {
            print("Mode 2 reg")
            l <- lm(tm.u ~ tm.b)
            tm.u <- (tm.u - l$coefficients[1]) / l$coefficients[2]
            l2 <- lm(tm.u ~ tm.b)
        }
        else {
            print("bad mode")
            return(expt)
        }
        print(paste(l$coefficients[2],l$coefficients[1]))
        print(paste(l2$coefficients[2],l2$coefficients[1]))
        print("------")

        for (s in sample_names) {
            r_names <- sapply(expt$samples[[s]]$replicates, '[[', "name")
            for (r2 in r_names) {
                if (r2 == r) {
                    if (mode == 1 | mode == 3) {
                        print("Mode 1/3 norm")
                        expt$samples[[s]]$replicates[[r]]$meta$temp <- expt$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]
                    }
                    else if (mode == 2) {
                        print("Mode 2 norm")
                        expt$samples[[s]]$replicates[[r]]$meta$temp <- (expt$samples[[s]]$replicates[[r]]$meta$temp - l$coefficients[1]) / l$coefficients[2]
                    }
                }
            }
        }
    }
    return(expt)

}



#' @title Normalize to a spike-in standard.
#'
#' @description Normalizes each replicate of an experiment based on a given
#' spike-in protein standard (assumed to be present in equimolar amounts in
#' each channel).
#'
#' @param expt an MSThermExperiment object
#' @param protein ID of a protein to normalize against
#' @param model whether to fit scale factors to model
#' @param plot (T/f) whether to display a summary plot
#'
#' @return An MsThermExperiment object with normalized data slots
#'
#' @examples
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#'
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#'
#' @export

normalize_to_std <- function(
    expt,
    protein,
    model=T,
    plot=T
) {

    n_replicates <- length(unlist(lapply(expt$samples,
        "[[", "replicates"),recursive=F))
    n_rows <- 1;
    while (n_rows*(n_rows+1) < n_replicates) {
        n_rows <- n_rows + 1
    }
    par(mfrow=c(n_rows,n_rows+1))

    n_samples <- length(expt$samples)
        
    for (i_sample in 1:n_samples) {

        sample <- expt$samples[[i_sample]]
        n_replicates <- length(sample$replicates)

        for (i_replicate in 1:n_replicates) {
            
            replicate <- sample$replicates[[i_replicate]]
            expt$samples[[i_sample]]$replicates[[i_replicate]]$data <-
                norm_to_std(replicate,protein,model,plot)

        }

    }

    return( expt )

}



#' @title Normalize to a profile.
#'
#' @description Normalizes an MSThermReplicate based on a pre-determined
#' vector of relative abundances
#'
#' @param replicate an MSThermReplicate object
#' @param profile a vector of relative values
#' @param model whether to fit scale factors to model
#' @param plot (T/f) whether to display a summary plot
#'
#' @return An MsThermReplicate object with normalized data slots
#'
#' @examples
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#'
#' profile <- c(50.0, 50.5, 47.5, 42.0, 37.0, 25.0, 16.0, 11.5, 10.5, 10.0)
#' expt$samples$Control$replicates$C1 <- normalize_to_profile(
#'    expt$samples$Control$replicates$C1, profile, plot=FALSE
#' )
#'
#' @export

normalize_to_profile <- function(
    replicate,
    profile,
    model=T,
    plot=T
) {

    replicate$data <- norm_to_profile(replicate, profile, model, plot)
    return( replicate )

}



# performs the actual profile normalization
norm_to_profile <- function(
    replicate,
    profile,
    model=T,
    plot=T
) {

    cols <- brewer.pal(3,"Set1")

    quant_columns <- match(replicate$meta$channel,colnames(replicate$data))

    temps <- replicate$meta$temp
    quant <- replicate$data[,quant_columns]

    std.ratios <- profile/profile[1]
    sums <- apply(as.matrix(quant),2,function(x){sum(as.numeric(x),na.rm=T)})
    ratios <- sums/sums[1]

    if(model) {
        sigmoid.formula    <- as.formula(paste("y ~ ", sigmoid))
        sigmoid.d1.formula <- as.formula(paste("y ~ ", sigmoid.d1))
        st.coarse <- expand.grid(p=c(0,0.3),k=seq(0,4000,by=1000),m=seq(30,60,by=15))
        x <- temps[std.ratios < 1.2]
        y <- std.ratios[std.ratios < 1.2]
        mod <- nls2(sigmoid.formula,data=list(x=x,y=y),start=st.coarse,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=50000))
        fit <- nls2(sigmoid.formula,data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=F),algorithm="port",lower=c(0,1,10),upper=c(0.4,100000,100))
        fit <- try_fit(y,x,trim=F,smooth=F)
        if (is.null(fit)) {
            model <- 0
        }
        #p <- fit$plat
        #k <- fit$k
        #m <- fit$tm
        #p <- coefficients(fit)[['p']]
        #k <- coefficients(fit)[['k']]
        #m <- coefficients(fit)[['m']]
    }
    if (model) {
        yfit <- eval(sigmoid,list(p=fit$plat, k=fit$k, m=fit$tm,x=temps))
        sf <- yfit/ratios
    }
    else {
        sf <- std.ratios/ratios
    }
    norm <- replicate$data
    norm[,quant_columns[order(replicate$meta$temp)] ] <- t(t(as.matrix(quant))*sf)
    corrected <- ratios*sf

    if (plot) {

        plot( temps,ratios,ylim=c(0,max(c(ratios,std.ratios))),xlab=expression(paste("temperature ", degree, "C")),ylab="relative fraction remaining",col=cols[3],main=replicate$name)
        points(temps,std.ratios,col=cols[1])
        points(temps,corrected,col=cols[2])
        if (model) {
            x.curve <- seq(min(temps), max(temps), by=0.1)
            y.curve <- eval(sigmoid,list(p=fit$plat, k=fit$k, m=fit$tm,x=x.curve))
            lines(x.curve, y.curve, col=cols[2])
            #curve(sigmoid(p,k,m,x),col=cols[2],add=T)
        }

    }

    return(norm)

}



# Perfoms the actual normalization to a spike-in standard (see
# "normalize_to_std" for the public method)
norm_to_std <- function(
    replicate,
    protein,
    model=T,
    plot=T
) {

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

    if (model) {
        sigmoid.formula    <- as.formula(paste("y ~ ", sigmoid))
        sigmoid.d1.formula <- as.formula(paste("y ~ ", sigmoid.d1))
        st.coarse <- expand.grid(p=c(0,0.3),k=seq(0,4000,by=1000),m=seq(30,60,by=15))
        x <- temps[corrected < 1.2]
        y <- corrected[corrected < 1.2]
        mod <- nls2(sigmoid.formula,data=list(x=x,y=y),start=st.coarse,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=5000))
        fit <- nls2(sigmoid.formula,data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=F),algorithm="port",lower=c(0,1,10),upper=c(0.4,100000,100))
        fit <- try_fit(y,x,trim=F,smooth=F)
        if (is.null(fit)) {
            model <- 0
        }
        #p <- fit$plat
        #k <- fit$k
        #m <- fit$tm
        #p <- coefficients(fit)[['p']]
        #k <- coefficients(fit)[['k']]
        #m <- coefficients(fit)[['m']]
    }
    if (model) {
        yfit <- eval(sigmoid,list(p=fit$plat, k=fit$k, m=fit$tm,x=temps))
        sf <- yfit/ratios
    }
    else {
        sf <- corrected/ratios
    }

    norm <- replicate$data
    norm[,quant_columns[order(replicate$meta$temp)] ] <- t(t(as.matrix(quant))*sf)

    if (plot) {

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
        if (model) {
            x.curve <- seq(min(temps), max(temps), by=0.1)
            y.curve <- eval(sigmoid,list(p=fit$plat, k=fit$k, m=fit$tm,x=x.curve))
            lines(x.curve, y.curve, col=cols[2])
            #curve(sigmoid(p,k,m,x),col=cols[2],add=T)
        }

    }

    return(norm)

}
