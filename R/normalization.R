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



#' @title Normalize to a spike-in standard.
#'
#' @description Normalizes each replicate of an experiment based on a given
#' spike-in protein standard (assumed to be present in equimolar amounts in
#' each channel).
#'
#' @param expt an MSThermExperiment object
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



#' @title Normalize to a profile.
#'
#' @description Normalizes an MSThermReplicate based on a pre-determined
#' vector of relative abundances
#'
#' @param replicate An MSThermReplicate object
#' @param profile A vector of relative values
#' @param model Whether to fit scale factors to model
#'
#' @return An MsThermReplicate object with normalized data slots
#'
#' @examples
#'\dontrun{
#' expt$ <- normalize_to_profile(expt, concentrations, model=T)
#'}
#'
#' @export

normalize_to_profile <- function( replicate, profile, model ) {

    replicate$data <- norm_to_profile(replicate,profile, model)
    return( replicate )

}



# performs the actual profile normalization
norm_to_profile <- function(replicate,profile,model=T) {

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



# Perfoms the actual normalization to a spike-in standard (see
# "normalize_to_std" for the public method)
norm_to_std <- function(replicate,protein) {

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
