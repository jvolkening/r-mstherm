# currently unused
#sigShift <- function( df, repl1, repl2, bin ) {
#
#   deltas <- df[[paste0(repl2,'.tm')]] - df[[paste0(repl1,'.tm')]]
#   m <- median(deltas,na.rm=T)
#   d <- mad(deltas,na.rm=T)
#   print(m)
#   print(d)
#   z <- (deltas - m)/d
#   p <- 2*pnorm(-abs(z))
#   q <- p.adjust(p,method="BH")
#   plot(density(deltas,na.rm=T),xlim=c(-10,10))
#   x <- NULL; rm(x) # silence R CMD check noise due to curve() call below
#   curve(dnorm(x,mean=m,sd=d),col="red",add=T)
#   return(q)
