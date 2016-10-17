#' @title MSResultSet to data frame.
#'
#' @description Populates a data frame with information from an MSResultSet,
#' one row per protein/group
#'
#' @param x an MSResultSet object
#' @param ... additional arguments passed to or from other functions
#'
#' @return A data frame populated with relevant information per result
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



#' @title Export MSThermResultSet to an SQLite database.
#'
#' @description Exports and MSThermResultSet object to a new SQLite database file.
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



#----------------------------------------------------------------------------#
# These functions are used for converting given data series in low resolution
# to and from byte strings for storage in an SQL database
#----------------------------------------------------------------------------#

# stores abundance ratios as signed 8-bit values clamped between -4 and +4 and
# returned as a string
y_to_str <- function (v) {

    y.conv <- vapply(v, function(y) {
        min(255, round(y*64+127, 0))
    }, double(1))
    bs <- rawToChar(as.raw(y.conv))
    return(bs)

}

# converts the output of "y_to_str" back to a numeric vector of abundance
# ratios
str_to_y <- function (s) {

    y.conv <- as.integer(charToRaw(s))
    y.orig <- (y.conv-127)/64
    return(y.orig)

}

# stores temperatures as unsigned 8-bit values returned as a string 
# input values must be between 20 and 83 (within expected Celcius temperature
# range for TPP experiments). This could be modified to accept a wider range
# of temperatures with a corresponding loss in precision of stored data or
# else using more bytes per values
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
    
# converts the output of "x_to_str" back to a numeric vector of temperatures
str_to_x <- function (s) {

    x.conv <- as.integer(charToRaw(s))
    x.orig <- x.conv/255*64+19;
    return(x.orig)

}
