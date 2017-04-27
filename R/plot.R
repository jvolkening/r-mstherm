#' @title Plot MSThermResult object.
#'
#' @description Generate a denaturation plot for an modeled protein/group.
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
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#' res     <- model_experiment(expt, bootstrap=FALSE, np=2)
#'
#' # plot single MSThermResult
#' plot(res$P38707)
#'
#' # plot all proteins (e.g. to pdf device, one-per-page)
#' plot(res)
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

        #l <- length(result$series)*4
        tbl <- matrix(rep(0,length(result$series)*4),ncol=4)
        tbl[,1] <- sapply(result$series, function(x) ifelse(is.null(x$psm),NA,x$psm))
        tbl[,2] <- round(sapply(result$series, function(x) ifelse(is.null(x$tm),NA,x$tm)),1)
        tbl[,3] <- round(sapply(result$series, function(x) ifelse(is.null(x$slope),NA,x$slope)),2)
        tbl[,4] <- round(sapply(result$series, function(x) ifelse(is.null(x$r2),NA,x$r2)),2)
        colnames(tbl) <- c("PSM",expression(T[m]),"Slp","R2")
        rownames(tbl) <- sapply(result$series, '[[', "name")
        addtable2plot(t.x,t.y,table=tbl,bty="o",lwd=1,hlines=T,xjust=just.x,yjust=just.y,display.rownames=T,xpad=0.4,ypad=1.0,cex=0.7,bg="#FFFFFF77")
    }

}



#' @title Plot MSThermResultSet object.
#'
#' @description Generate a series of denaturation plots for all results in an
#' MSThermResultSet.
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
#' # see plot.MSThermResult for an example
#'
#' @export

plot.MSThermResultSet <- function(x, ...) {

    sapply(x, plot, ...)
    return(invisible(NULL))

}
