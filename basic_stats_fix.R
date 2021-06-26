basic.stats_fix <- function (df_test, diploid = TRUE, digits = 4) 
{
    if (is.genind(df_test)) 
        #####FIX#####
        df_test <- genind2hierfstat_fix(df_test)
    loc.names <- names(df_test)[-1]
    if (length(table(df_test[, 1])) < 2) 
        df_test[dim(df_test)[1] + 1, 1] <- df_test[dim(df_test)[1], 1] + 
            1
    if (dim(df_test)[2] == 2) 
        df_test <- df_test.frame(df_test, dummy.loc = df_test[, 2])
    p <- pop.freq(df_test, diploid)
    n <- t(ind.count(df_test))
    if (diploid) {
        dum <- getal.b(df_test[, -1])
        Ho <- dum[, , 1] == dum[, , 2]
        sHo <- (1 - t(apply(Ho, 2, fun <- function(x) tapply(x, 
            df_test[, 1], mean, na.rm = TRUE))))
        mHo <- apply(sHo, 1, mean, na.rm = TRUE)
    }
    else {
        sHo <- NA
        mHo <- NA
    }
    sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
    sp2 <- matrix(unlist(sp2), nrow = dim(df_test[, -1])[2], byrow = TRUE)
    if (diploid) {
        Hs <- (1 - sp2 - sHo/2/n)
        Hs <- n/(n - 1) * Hs
        Fis = 1 - sHo/Hs
    }
    else {
        Hs <- n/(n - 1) * (1 - sp2)
        Fis <- NA
    }
    np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
    mn <- apply(n, 1, fun <- function(x) {
        np <- sum(!is.na(x))
        np/sum(1/x[!is.na(x)])
    })
    msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
    mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
    mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
    if (diploid) {
        mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
        Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
        mFis = 1 - mHo/mHs
    }
    else {
        mHs <- mn/(mn - 1) * (1 - msp2)
        Ht <- 1 - mp2 + mHs/mn/np
        mFis <- NA
    }
    Dst <- Ht - mHs
    Dstp <- np/(np - 1) * Dst
    Htp = mHs + Dstp
    Fst = Dst/Ht
    Fstp = Dstp/Htp
    Dest <- Dstp/(1 - mHs)
            #####FIX#####
    res <- as.data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst, 
        Fstp, mFis, Dest))
    names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp", "Fst", 
        "Fstp", "Fis", "Dest")
    if (diploid) {
        rownames(sHo) <- loc.names
        rownames(Fis) <- loc.names
    }
    is.na(res) <- do.call(cbind, lapply(res, is.infinite))
    overall <- apply(res, 2, mean, na.rm = TRUE)
    overall[7] <- overall[4]/overall[3]
    overall[8] <- overall[6]/overall[5]
    overall[9] <- 1 - overall[1]/overall[2]
    overall[10] <- overall[6]/(1 - overall[2])
    names(overall) <- names(res)
    if (!diploid) {
        overall[-2] <- NA
    }
    all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round, 
        digits), Ho = round(sHo, digits), Hs = round(Hs, digits), 
        Fis = round(Fis, digits), perloc = round(res, digits), 
        overall = round(overall, digits))
    class(all.res) <- "basic.stats"
    all.res
}

        #####FIX#####
genind2hierfstat_fix <- function (dat, pop = dat@strata$pop1_2) 
{
    if (!is.genind(dat)) 
        stop("dat must be a genind object. Exiting")
    if (is.null(pop)) {
        if (is.null(adegenet::pop(dat))) {
            stop("population factor must be defined")
        }
        else {
            pop <- adegenet::pop(dat)
        }
    }
    if (dat@type != "codom") 
        stop("data type must be codominant. Exiting")
    ploid <- unique(dat@ploidy)
    if (length(ploid) != 1) 
        stop("data must contain only diploids or only haploids. Exiting")
    if (ploid > 2L) 
        stop("Data must come from diploids or haploids. Exiting")
    alleles.name <- toupper(unique(unlist(adegenet::alleles(dat))))
    ids <- adegenet::indNames(dat)
    x <- if (!all(alleles.name %in% c("A", "C", "G", "T"))) {
        if (length(grep("[[:alpha:]|[:punct:]]", alleles.name)) > 
            0) {
            max.length <- max(sapply(adegenet::alleles(dat), 
                length))
            digits <- floor(log10(max.length))
            dat <- as.matrix(adegenet::genind2df(dat, sep = "", 
                usepop = FALSE, oneColPerAll = TRUE))
            dat <- apply(dat, 2, function(a) ifelse(a == "NA", 
                NA, a))
            do.call(cbind, lapply(seq(ploid, ncol(dat), by = ploid), 
                function(end) {
                  start <- end - ploid + 1
                  allelesid <- sort(unique(as.vector(dat[, start:end])))
                  apply(dat[, start:end, drop = FALSE], 1, function(gntp) {
                    if (any(is.na(gntp))) 
                      return(NA)
                    gntp <- match(gntp, allelesid)
                    gntp <- formatC(gntp, digits = digits, flag = "0", 
                      mode = "integer")
                    as.integer(paste(gntp, collapse = ""))
                  })
                }))
        }
        else {
            tmp1 <- unique(nchar(unlist(dat@all.names)))
            if (length(tmp1) > 1) {
                dig <- max(tmp1)
                allnames <- lapply(dat@all.names, formatC, width = dig, 
                  flag = "0", mode = "integer")
                dat@all.names <- allnames
            }
            dat <- adegenet::genind2df(dat, sep = "", usepop = FALSE)
            do.call(cbind, lapply(dat, as.integer))
        }
    }
    else {
        dat <- adegenet::genind2df(dat, sep = "", usepop = FALSE)
        do.call(cbind, lapply(dat, function(a) {
            a <- gsub("[aA]", "1", a)
            a <- gsub("[cC]", "2", a)
            a <- gsub("[gG]", "3", a)
            a <- gsub("[tT]", "4", a)
            as.integer(a)
        }))
    }
    x <- data.frame(pop = pop, x)
    rownames(x) <- ids
    if (is.factor(pop) & nlevels(pop) == 1) {
        dum1 <- dim(dat)[1]
        x[dum1 + 1, ] <- NA
        x[, 1] <- factor(c(pop, "dumpop"))
        rownames(x) <- c(ids, "dumind")
    }
    return(x)
}