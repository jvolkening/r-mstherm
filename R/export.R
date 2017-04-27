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
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#' res     <- model_experiment(expt, bootstrap=FALSE, np=2)
#'
#' df <- as.data.frame(res)
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
        #x <- x[1]$series[[r]][['x']]

        for(col in c("tm","psm","inf","slope","k","plat","r2","rmsd")) {
            s <- sapply(x, function(v) v$series[[r]][[col]])
            s[sapply(s, is.null)] <- NA
            df[[paste0(r,'.',col)]]  <- unlist(s, use.names=F)
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
#' control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
#' annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")
#' expt    <- MSThermExperiment(control, annotations=annots)
#' expt    <- normalize_to_std(expt, "cRAP_ALBU_BOVIN", plot=FALSE)
#' res     <- model_experiment(expt, bootstrap=FALSE, np=2)
#'
#' fn <- tempfile(fileext = ".sqlite")
#' write.sqlite(res, fn)
#' unlink(fn) # for example only
#'
#' @export

write.sqlite <- function( res, file ) {

    if (requireNamespace("RSQLite")) {

        con <- RSQLite::dbConnect(RSQLite::SQLite(), file)

        create <- "CREATE TABLE proteins (
            id INTEGER PRIMARY KEY NOT NULL,
            name TEXT NOT NULL,
            description TEXT
        );"
        rs <- RSQLite::dbSendQuery(con, create)
        RSQLite::dbClearResult(rs)

        create <- "CREATE TABLE samples (
            id INTEGER PRIMARY KEY NOT NULL,
            name TEXT NOT NULL
        );"
        rs <- RSQLite::dbSendQuery(con, create)
        RSQLite::dbClearResult(rs)

        create <- "CREATE TABLE replicates (
            id INTEGER PRIMARY KEY NOT NULL,
            name TEXT NOT NULL
        );"
        rs <- RSQLite::dbSendQuery(con, create)
        RSQLite::dbClearResult(rs)

        create <- "CREATE TABLE temps (
            id INTEGER PRIMARY KEY NOT NULL,
            value BLOB NOT NULL
        );"
        rs <- RSQLite::dbSendQuery(con, create)
        RSQLite::dbClearResult(rs)

        create <- "CREATE TABLE data (
            protein INTEGER NOT NULL,
            replicate INTEGER NOT NULL,
            sample INTEGER NOT NULL,
            m REAL NOT NULL,
            k REAL NOT NULL,
            p REAL NOT NULL,
            x INTEGER NOT NULL,
            y BLOB NOT NULL,
            psm INTEGER NOT NULL,
            inf INTEGER NOT NULL,
            slope REAL NOT NULL,
            r2 INTEGER NOT NULL,
            el BLOB,
            eu BLOB,
            ml REAL,
            mu REAL,
            PRIMARY KEY (protein, replicate)
        );"
        rs <- RSQLite::dbSendQuery(con, create)
        RSQLite::dbClearResult(rs)

        i.p <- 0
        i.s <- 0
        i.r <- 0
        i.x <- 0

        used.p <- c()
        used.s <- c()
        used.r <- c()
        used.x <- c()

        for (p in res) {

            p.id <- match(p$name, used.p, nomatch=0) - 1
            if (p.id < 0) {
                p.id = i.p
                i.p <- i.p + 1
                used.p <- append(used.p, p$name)

                df <- data.frame(id=p.id,name=p$name,annot=p$annotation)
                sql <- "INSERT INTO proteins VALUES (:id, :name, :annot)"
                rs <-RSQLite::dbGetQuery(con, sql, df)

            }

            for (r in p$series) {

                # skip un-modeled series
                if (! r$is.fitted) {
                    next
                }

                x_blob <- x_to_str(r$x)
                x_tag  <- rawToChar(x_blob)
                x.id <- match(x_tag, used.x, nomatch=0) - 1
                if (x.id < 0) {
                    x.id = i.x
                    i.x <- i.x + 1
                    used.x <- append(used.x, x_tag)

                    df <- data.frame(id=x.id,value=I(list(x_blob)))
                    sql <- "INSERT INTO temps VALUES (:id, :value)"
                    rs <-RSQLite::dbGetQuery(con, sql, df)

                }

                r.id <- match(r$name, used.r, nomatch=0) - 1
                if (r.id < 0) {
                    r.id = i.r
                    i.r <- i.r + 1
                    used.r <- append(used.r, r$name)
                    
                    df <- data.frame(id=r.id,name=r$name)
                    sql <- "INSERT INTO replicates VALUES (:id, :name)"
                    rs <-RSQLite::dbGetQuery(con, sql, df)

                }

                s.id <- match(r$sample, used.s, nomatch=0) - 1
                if (s.id < 0) {
                    s.id = i.s
                    i.s <- i.s + 1
                    used.s <- append(used.s, r$sample)

                    df <- data.frame(id=s.id,name=r$sample)
                    sql <- "INSERT INTO samples VALUES (:id, :name)"
                    rs <-RSQLite::dbGetQuery(con, sql, df)
                }

                has_pci <- ! is.null(r$bs.lowers)
                has_tci <- ! is.null(r$tm_CI)
                has_inf <- ! is.null(r$inf)

                y_blob <- y_to_str(r$y)

                # even though last two values are floats, treat as strings so that
                # 'NULL' can be substituted when necessary
                df <- data.frame(
                    pid = p.id,
                    rid = r.id,
                    sid = s.id,
                    tm  = r$tm,
                    k   = r$k,
                    p   = r$plat,
                    xid = x.id,
                    y   = I(list(y_blob)),
                    psm = r$psm,
                    inf = if(has_inf) as.integer( r$inf * 255 ) else 'NULL',
                    slp = r$slope,
                    r2  = as.integer( r$r2 * 255 ),
                    v1  = if(has_pci) I(list(y_to_str( r$bs.lowers ))) else 'NULL',
                    v2  = if(has_pci) I(list(y_to_str( r$bs.uppers ))) else 'NULL',
                    v3  = if(has_tci) r$tm_CI[[1]] else 'NULL',
                    v4  = if(has_tci) r$tm_CI[[2]] else 'NULL'
                )
                sql <- "INSERT INTO data VALUES (:pid, :rid, :sid, :tm, :k, :p, :xid, :y, :psm, :inf, :slp, :r2, :v1, :v2, :v3, :v4)"
                rs <-RSQLite::dbGetQuery(con, sql, df)

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
    bs <- as.raw(y.conv)
    return(bs)

}

# converts the output of "y_to_str" back to a numeric vector of abundance
# ratios
str_to_y <- function (s) {

    y.conv <- as.integer(s)
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
    bs <- as.raw(x.conv)
    return(bs)

}
    
# converts the output of "x_to_str" back to a numeric vector of temperatures
str_to_x <- function (s) {

    x.conv <- as.integer(s)
    x.orig <- x.conv/255*64+19;
    return(x.orig)

}
