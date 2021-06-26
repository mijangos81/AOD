
write.bayescan <- function (dat = dat, diploid = TRUE,file_name) {
    nloc <- dim(dat)[2] - 1
    npop <- length(table(dat[, 1]))
    alc.dat <- allele.count(dat, diploid)
    nal <- unlist(lapply(alc.dat, function(x) dim(x)[1]))
    nindx <- sapply(alc.dat, function(x) apply(x, 2, sum))
    cat(paste("[loci]=", nloc,"\r\n",sep = ""),file = file_name)
    cat(paste("\r\n"), file =file_name, append = TRUE)
    cat(paste("[populations]=", npop,"\r\n",sep = ""),file = file_name, append = TRUE)
    for (ip in 1:npop) {
        cat(paste("\r\n"), file = file_name, append = TRUE)
        cat(paste("[pop]=", ip,"\r\n",sep = ""),file = file_name, append = TRUE)
        for (il in 1:nloc) {
            tow <- c(il,nindx[ip, il], nal[il], alc.dat[[il]][,ip])
            tow <- c(tow," ")
            cat(tow,"\r\n",file = file_name, append = TRUE)
        }
    }
}

# Function taken from Whitlock and Lotterhos 2014
ReturnCorrectedPVal.BS <- function(bayescan_ouput){
    ### BayeScan outputs a file with q-values, but not p-values, already corrected for
    ### this function returns corrected p values for a Bayescan outfile
    
    bs <- bayescan_ouput
    p.valC <- as.numeric(1-bs$POST_PROB)
    num.obs <- length(bs$POST_PROB)
    cent <- mean(as.numeric(bs$FST))
    
    Tail <- rep(NA, num.obs)
    Tail[bs$FST<cent] <- "L"
    Tail[bs$FST>cent] <- "R"
    Tail <- as.factor(Tail)
    
    Bonf.CO <- 0.05/num.obs
    Bonf <- rep(FALSE, num.obs)
    Bonf[p.valC<Bonf.CO] <- TRUE
    
    FDR.01 <- rep(FALSE, num.obs)
    FDR.01[bs$Q_VALUE<0.01] <-  TRUE
    
    FDR.05 <- rep(FALSE, num.obs)
    FDR.05[bs$Q_VALUE<0.05] <-  TRUE
    
    out <- cbind(bs, tail=Tail, p.val.tail=p.valC, Bonf, FDR.01, FDR.05)
    return(out)
}

run_bayescan <- function (data, n = 5000, thin = 10, nbp = 20, pilot = 5000,
burn = 50000, pr_odds, subsample = NULL, iteration.subsample = 1,
parallel.core = number_cores, bayescan.path = bayescan.path,time_report)
{
    # cat("#######################################################################\n")
    # cat("###################### radiator::run_bayescan #########################\n")
    # cat("#######################################################################\n")
    timing <- proc.time()
    res <- list()
    if (!file.exists(bayescan.path)) {
        rlang::abort("Path to BayeScan install is not valid")
    }
    if (missing(data))
    rlang::abort("Input file missing")
    if (missing(pr_odds))
    rlang::abort("Prior odds for the neutral model is missing.\n                             No shortcut with default here, sorry.\n                             Please read the BayeScan manual...")
    file.date <- format(Sys.time(), "%Y%m%d@%H%M%S")
    if ((!is.null(subsample))) {
        folder.message <- stringi::stri_join("radiator_bayescan_subsampling_",
        file.date, sep = "")
    }
    else {
        folder.message <- stringi::stri_join("radiator_bayescan_",
        file.date, sep = "")
    }
    # dir_data <- basename(dirname(data))
    # path.folder <- stringi::stri_join(getwd(),"/",dir_data ,"/", folder.message,sep = "")
    path.folder <- stringi::stri_join(getwd(), "/", folder.message,sep = "")
    dir.create(file.path(path.folder))
    # message("\nFolder created: \n", folder.message)
    if (!is.null(subsample)) {
        message("Subsampling: selected")
        data.type <- radiator::detect_genomic_format(data = data)
        if (is.vector(data)) {
            if (data.type != "fst.file") {
                rlang::abort("Using subsample argument requires a tidy data frame saved by\n             radiator::tidy_genomic_data function")
            }
            else {
                data <- radiator::tidy_genomic_data(data = data,
                monomorphic.out = FALSE, common.markers = FALSE,
                verbose = FALSE)
            }
        }
        else {
            columns.tidy <- colnames(data)
            want <- c("GT_VCF", "GT_VCF_NUC", "GT", "GT_BIN")
            want.check <- TRUE %in% (unique(want %in% columns.tidy))
            want.more <- c("MARKERS", "INDIVIDUALS", "POP_ID")
            want.more.check <- isTRUE(unique(want.more %in% columns.tidy))
            is.tidy <- isTRUE(unique(c(want.check, want.more.check)))
            if (!is.tidy)
            rlang::abort("A tidy data frame object required")
        }
        ind.pop.df <- dplyr::distinct(.data = data, POP_ID, INDIVIDUALS)
        strata.stats <- ind.pop.df %>% dplyr::group_by(POP_ID) %>%
        dplyr::tally(.) %>% dplyr::mutate(STRATA = stringi::stri_join(POP_ID,
        n, sep = " = "))
        n.pop <- dplyr::n_distinct(ind.pop.df$POP_ID)
        n.ind <- dplyr::n_distinct(ind.pop.df$INDIVIDUALS)
        message("Number of populations: ", n.pop)
        message("Number of individuals: ", n.ind)
        message("Number of ind/pop:\n", stringi::stri_join(strata.stats$STRATA,
        collapse = "\n"))
        message("Number of markers: ", dplyr::n_distinct(data$MARKERS))
        if (subsample == "min") {
            subsample <- ind.pop.df %>% dplyr::group_by(POP_ID) %>%
            dplyr::tally(.) %>% dplyr::filter(n == min(n)) %>%
            dplyr::ungroup(.) %>% dplyr::select(n) %>% purrr::flatten_int(.)
            message("\nSubsample used: ", subsample)
        }
        subsample.list <- purrr::map(.x = 1:iteration.subsample,
        .f = subsampling_data, ind.pop.df = ind.pop.df, subsample = subsample)
        subsampling.individuals <- dplyr::bind_rows(subsample.list)
        readr::write_tsv(x = subsampling.individuals, path = file.path(path.folder,
        "radiator_bayescan_subsampling_individuals.tsv"),
        col_names = TRUE, append = FALSE)
        res$subsampling.individuals <- subsampling.individuals
    }
    else {
        iteration.subsample <- 1
    }
    if (is.null(subsample)) {
        res <- bayescan_one(data = data, n = n, thin = thin,
        nbp = nbp, pilot = pilot, burn = burn, pr_odds = pr_odds,
        parallel.core = parallel.core, path.folder = path.folder,
        file.date = file.date, bayescan.path = bayescan.path)
    }
    else {
        subsample.bayescan <- purrr::map(.x = subsample.list,
        .f = bayescan_one, data = data, n = n, thin = thin,
        nbp = nbp, pilot = pilot, burn = burn, pr_odds = pr_odds,
        subsample = subsample, iteration.subsample = iteration.subsample,
        parallel.core = parallel.core, path.folder = path.folder,
        file.date = file.date, bayescan.path = bayescan.path)
        cat("\n\n#######################################################################\n")
        message("Summarizing subsampling results...")
        res$bayescan.all.subsamples <- purrr::map_df(subsample.bayescan,
        "bayescan") %>% dplyr::select(-BAYESCAN_MARKERS)
        readr::write_tsv(x = res$bayescan.all.subsamples, path = file.path(path.folder,
        "bayescan.all.subsamples.tsv"), col_names = TRUE,
        append = FALSE)
        iteration.number <- dplyr::n_distinct(res$bayescan.all.subsamples$ITERATIONS)
        markers.summary <- dplyr::ungroup(res$bayescan.all.subsamples) %>%
        dplyr::select(MARKERS) %>% dplyr::group_by(MARKERS) %>%
        dplyr::tally(.)
        markers.whitelist <- dplyr::filter(markers.summary, n ==
        iteration.number) %>% dplyr::distinct(MARKERS)
        markers.all.iterations <- nrow(markers.whitelist)
        total.unique.markers <- dplyr::n_distinct(markers.summary$MARKERS)
        proportion.keeper <- round(markers.all.iterations/total.unique.markers,
        2)
        message("BayeScan subsampling summary: ")
        message("    number of unique markers: ", total.unique.markers)
        message("    keeping markers common in all iterations: ",
        markers.all.iterations, " (= ", proportion.keeper,
        ")")
        bayescan.all.subsamples.filtered <- dplyr::left_join(markers.whitelist,
        res$bayescan.all.subsamples, by = "MARKERS")
        res$selection.accuracy <- bayescan.all.subsamples.filtered %>%
        dplyr::group_by(MARKERS, SELECTION) %>% dplyr::tally(.)
        readr::write_tsv(x = res$selection.accuracy, path = file.path(path.folder,
        "selection.accuracy.tsv"), col_names = TRUE, append = FALSE)
        res$accurate.markers <- dplyr::ungroup(res$selection.accuracy) %>%
        dplyr::filter(n == iteration.number) %>% dplyr::distinct(MARKERS,
        .keep_all = TRUE) %>% dplyr::select(-n)
        readr::write_tsv(x = res$accurate.markers, path = file.path(path.folder,
        "accurate.markers.tsv"), col_names = TRUE, append = FALSE)
        accurate.markers.summary <- res$accurate.markers %>%
        dplyr::group_by(SELECTION) %>% dplyr::tally(.)
        accurate.markers.number <- nrow(res$accurate.markers)
        res$accuracy.summary <- tibble::data_frame(total = total.unique.markers,
        `found in all iterations` = markers.all.iterations,
        `not accurate` = markers.all.iterations - accurate.markers.number,
        accurate = accurate.markers.number, `accurate + neutral` = accurate.markers.summary$n[accurate.markers.summary$SELECTION ==
        "neutral"], `accurate + balancing` = accurate.markers.summary$n[accurate.markers.summary$SELECTION ==
        "balancing"], `accurate + diversifying` = accurate.markers.summary$n[accurate.markers.summary$SELECTION ==
        "diversifying"]) %>% tidyr::gather(key = "ACCURACY_MARKERS",
        value = "N") %>% dplyr::mutate(PROP = N/total.unique.markers)
        readr::write_tsv(x = res$accuracy.summary, path = file.path(path.folder,
        "accuracy.summary.tsv"), col_names = TRUE, append = FALSE)
        res$bayescan.summary <- dplyr::left_join(dplyr::select(res$accurate.markers,
        MARKERS), res$bayescan.all.subsamples, by = "MARKERS") %>%
        dplyr::group_by(MARKERS) %>% dplyr::summarise_if(.tbl = .,
        .predicate = is.numeric, .funs = mean) %>% dplyr::select(-ITERATIONS) %>%
        dplyr::mutate(SELECTION = factor(dplyr::if_else(ALPHA >=
        0 & Q_VALUE <= 0.05, "diversifying", dplyr::if_else(ALPHA >=
        0 & Q_VALUE > 0.05, "neutral", "balancing"))),
        PO_GROUP = factor(dplyr::if_else(LOG10_PO > 2,
        "decisive", dplyr::if_else(LOG10_PO > 1.5,
        "very strong", dplyr::if_else(LOG10_PO >
        1, "strong", dplyr::if_else(LOG10_PO >
        0.5, "substantial", "no evidence")))),
        levels = c("no evidence", "substantial", "strong",
        "very strong", "decisive"), ordered = TRUE)) %>%
        dplyr::ungroup(.) %>% dplyr::mutate(FST_GROUP = dplyr::ntile(FST,
        5), FST_GROUP = dplyr::if_else(FST_GROUP == 1, "0-20%",
        dplyr::if_else(FST_GROUP == 2, "20-40%", dplyr::if_else(FST_GROUP ==
        3, "40-60%", dplyr::if_else(FST_GROUP == 4, "60-80%",
        "80-100%"))))) %>% dplyr::arrange(FST)
        readr::write_tsv(x = res$bayescan.summary, path = file.path(path.folder,
        "bayescan.summary.tsv"), col_names = TRUE, append = FALSE)
        # res$bayescan.summary.plot <- plot_bayescan(res$bayescan.summary)
        # ggplot2::ggsave(filename = file.path(path.folder, "bayescan.summary.plot.pdf"),
        #     plot = res$bayescan.summary.plot, width = 30, height = 15,
        #     dpi = 600, units = "cm", useDingbats = FALSE)
        res$selection.summary <- res$bayescan.summary %>% dplyr::group_by(SELECTION,
        PO_GROUP) %>% dplyr::tally(.) %>% dplyr::rename(MARKERS = n)
        readr::write_tsv(x = res$selection.summary, path = file.path(path.folder,
        "selection.summary.tsv"), col_names = TRUE, append = FALSE)
        message("Generating blacklist and whitelists for all iterations")
        all.markers <- dplyr::distinct(markers.summary, MARKERS)
        res$whitelist.markers.positive.selection <- res$bayescan.summary %>%
        dplyr::filter(SELECTION == "diversifying" & PO_GROUP !=
        "no evidence") %>% dplyr::distinct(MARKERS) %>%
        dplyr::arrange(MARKERS)
        if (nrow(res$whitelist.markers.positive.selection) >
        0) {
            readr::write_tsv(x = res$whitelist.markers.positive.selection,
            path = file.path(path.folder, "whitelist.markers.positive.selection.tsv"))
            positive <- TRUE
            message("    whitelist positive/directional selection: generated")
        }
        else {
            message("    whitelist positive/directional selection: not generated")
            positive <- FALSE
        }
        res$whitelist.markers.neutral.selection <- res$bayescan.summary %>%
        dplyr::filter(SELECTION == "neutral") %>% dplyr::distinct(MARKERS) %>%
        dplyr::arrange(MARKERS)
        if (nrow(res$whitelist.markers.neutral.selection) > 0) {
            readr::write_tsv(x = res$whitelist.markers.neutral.selection,
            path = file.path(path.folder, "whitelist.markers.neutral.selection.tsv"))
            neutral <- TRUE
            message("    whitelist neutral selection: generated")
        }
        else {
            message("    whitelist neutral selection: not generated")
            neutral <- FALSE
        }
        if (neutral && positive) {
            res$whitelist.markers.neutral.positive.selection <- res$bayescan.summary %>%
            dplyr::filter(SELECTION == "neutral" | (SELECTION ==
            "diversifying" & PO_GROUP != "no evidence")) %>%
            dplyr::distinct(MARKERS) %>% dplyr::arrange(MARKERS)
            readr::write_tsv(x = res$whitelist.markers.neutral.positive.selection,
            path = file.path(path.folder, "whitelist.markers.neutral.positive.selection.tsv"))
            message("    whitelist neutral and positive/directional selections: generated")
        }
        else {
            message("    whitelist neutral and positive/directional selections: not generated")
        }
        res$blacklist.markers.balancing.selection <- res$bayescan.summary %>%
        dplyr::filter(SELECTION == "balancing") %>% dplyr::distinct(MARKERS) %>%
        dplyr::arrange(MARKERS)
        if (nrow(res$blacklist.markers.balancing.selection) >
        0) {
            readr::write_tsv(x = res$blacklist.markers.balancing.selection,
            path = file.path(path.folder, "blacklist.markers.balancing.selection.tsv"))
            balancing <- TRUE
            message("    blacklist balancing selection: generated")
        }
        else {
            message("    blacklist balancing selection: not generated")
            balancing <- FALSE
        }
        if (neutral && positive && balancing) {
            res$whitelist.markers.without.balancing.positive <- dplyr::anti_join(all.markers,
            res$blacklist.markers.balancing.selection, by = "MARKERS") %>%
            dplyr::anti_join(res$whitelist.markers.positive.selection,
            by = "MARKERS") %>% readr::write_tsv(x = res$whitelist.markers.without.balancing.positive,
            path = file.path(path.folder, "whitelist.markers.without.balancing.positive.tsv"))
            message("    whitelist without balancing and positive selection: generated")
        }
        if (neutral && balancing && !positive) {
            res$whitelist.markers.without.balancing.positive <- dplyr::anti_join(all.markers,
            res$blacklist.markers.balancing.selection, by = "MARKERS")
            readr::write_tsv(x = res$whitelist.markers.without.balancing.positive,
            path = file.path(path.folder, "whitelist.markers.without.balancing.positive.tsv"))
            message("    whitelist without balancing and positive selection: generated")
        }
    }
    timing <- proc.time() - timing
     if(time_report==T){message("\nComputation time: ", round(timing[[3]]), " sec")}
   # cat("############################## completed ##############################\n")
    return(res)
}

detect_genomic_format <- function (data){
    if (!is.vector(data)) {
        if (tibble::has_name(data, "INDIVIDUALS")) {
            data.type <- "tbl_df"
        }
        else {
            data.type <- class(data)[1]
            if (!data.type %in% c("genind", "genlight", "gtypes",
            "SeqVarGDSClass"))
            rlang::abort("Input file not recognised")
        }
    }
    else {
        data.type <- suppressWarnings(readLines(con = data, n = 1L))
        file.ending <- stringi::stri_sub(str = data, from = -4,
        to = -1)
        if (identical(data.type, "##fileformat=VCF") || file.ending ==
        ".vcf") {
            data.type <- "vcf.file"
        }
        if (file.ending == ".tped") {
            data.type <- "plink.file"
            if (!file.exists(stringi::stri_replace_all_fixed(str = data,
            pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
                rlang::abort("Missing tfam file with the same prefix as your tped")
            }
            return(data.type)
        }
        if (stringi::stri_detect_fixed(str = data.type, pattern = "POP_ID") |
        stringi::stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") |
        stringi::stri_detect_fixed(str = data.type, pattern = "MARKERS") |
        stringi::stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
            data.type <- "tbl_df"
            return(data.type)
        }
        if (stringi::stri_detect_fixed(str = data.type, pattern = "Catalog")) {
            data.type <- "haplo.file"
        }
        if (file.ending == ".gen") {
            data.type <- "genepop.file"
        }
        if (file.ending == ".dat") {
            data.type <- "fstat.file"
        }
        if (TRUE %in% (c(".rad", ".gds") %in% file.ending)) {
            if (stringi::stri_detect_fixed(str = data.type, pattern = "COREARRAY")) {
                data.type <- "gds.file"
            }
            else {
                data.type <- "fst.file"
            }
        }
        dart.temp <- check_dart(data)
        if (dart.temp$data.type %in% c("dart", "silico.dart")) {
            data.type <- dart.temp$data.type
        }
    }
    return(data.type)
}

tidy_genomic_data <- function (data, strata = NULL, filename = NULL, parallel.core = 7, verbose = FALSE, ...){
    if (verbose) {
        cat("################################################################################\n")
        cat("######################### radiator::tidy_genomic_data ##########################\n")
        cat("################################################################################\n")
    }
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    if (verbose)
    message("Execution date/time: ", file.date)
    old.dir <- getwd()
    opt.change <- getOption("width")
    options(width = 70)
    timing <- proc.time()
    res <- list()
    on.exit(setwd(old.dir), add = TRUE)
    on.exit(options(width = opt.change), add = TRUE)
    on.exit(timing <- proc.time() - timing, add = TRUE)
    on.exit(if (verbose) message("\nComputation time, overall: ",
    round(timing[[3]]), " sec"), add = TRUE)
    on.exit(if (verbose) cat("######################### tidy_genomic_data completed ##########################\n"),
    add = TRUE)
    rad.dots <- radiator_dots(func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(), args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error",
    .check_assign = TRUE), keepers = c("path.folder",
    "parameters", "keep.allele.names", "blacklist.id",
    "whitelist.markers", "filter.common.markers", "filter.monomorphic",
    "vcf.metadata", "vcf.stats", "blacklist.genotypes",
    "internal"), deprecated = c("maf.thresholds", "common.markers",
    "max.marker", "monomorphic.out", "snp.ld", "filter.call.rate",
    "filter.markers.coverage", "filter.markers.missing",
    "number.snp.reads", "mixed.genomes.analysis", "duplicate.genomes.analysis",
    "maf.data", "hierarchical.levels", "imputation.method",
    "pred.mean.matching", "num.tree", "pop.levels", "pop.labels",
    "pop.select"), verbose = FALSE)
    if (missing(data))
    rlang::abort("data is missing")
    path.folder <- generate_folder(f = path.folder, rad.folder = "radiator_tidy_genomic",
    internal = internal, file.date = file.date, verbose = verbose)
    write_rad(data = rad.dots, path = path.folder, filename = stringi::stri_join("radiator_tidy_genomic_data_args_",
    file.date, ".tsv"), tsv = TRUE, internal = internal,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose)
    skip.tidy.wide <- FALSE
    data.type <- radiator::detect_genomic_format(data)
    whitelist.markers <- read_whitelist(whitelist.markers, verbose)
    blacklist.id <- read_blacklist_id(blacklist.id, verbose)
    strata.df <- read_strata(strata = strata, pop.id = TRUE,
    blacklist.id = blacklist.id, verbose = verbose) %$% strata
    if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
        if (!"SeqVarTools" %in% utils::installed.packages()[,
        "Package"]) {
            rlang::abort("Please install SeqVarTools for this option:\n\n                   install.packages(\"BiocManager\")\n                   BiocManager::install(\"SeqVarTools\")")
        }
        if (data.type == "gds.file") {
            data <- radiator::read_rad(data, verbose = verbose)
        }
        data <- gds2tidy(gds = data, parallel.core = parallel.core)
        data.type <- "tbl_df"
    }
    if (data.type == "vcf.file") {
        if (verbose)
        message("Importing and tidying the VCF...")
        input <- radiator::tidy_vcf(data = data, strata = strata.df,
        parallel.core = parallel.core, verbose = verbose,
        whitelist.markers = whitelist.markers, filter.monomorphic = filter.monomorphic,
        filter.common.markers = filter.common.markers, filename = NULL,
        vcf.metadata = vcf.metadata, vcf.stats = vcf.stats,
        gt.vcf.nuc = TRUE, gt.vcf = TRUE, gt = TRUE, gt.bin = TRUE,
        path.folder = path.folder, internal = TRUE, tidy.check = FALSE)
        biallelic <- radiator::detect_biallelic_markers(input)
    }
    if (data.type == "plink.file") {
        if (verbose)
        message("Importing the PLINK files...")
        input <- tidy_plink(data = data, strata = strata.df,
        verbose = verbose, whitelist.markers = whitelist.markers,
        blacklist.id = blacklist.id)
        biallelic <- input$biallelic
        input <- input$input
    }
    if (data.type == "haplo.file") {
        if (verbose)
        message("Importing STACKS haplotype file")
        strata.df <- strata_haplo(strata = strata.df, data = data,
        blacklist.id = blacklist.id)
        want <- tibble::tibble(INFO = "CATALOG", COL_TYPE = "c") %>%
        dplyr::bind_rows(dplyr::select(strata.df, INFO = INDIVIDUALS) %>%
        dplyr::mutate(COL_TYPE = rep("c", n()), INFO = clean_ind_names(INFO)))
        haplo.col.type <- readr::read_tsv(file = data, n_max = 1,
        na = "-", col_names = FALSE, col_types = readr::cols(.default = readr::col_character())) %>%
        tidyr::gather(data = ., key = DELETE, value = INFO) %>%
        dplyr::mutate(INFO = clean_ind_names(INFO)) %>% dplyr::select(-DELETE) %>%
        dplyr::mutate(INFO = clean_ind_names(INFO)) %>% dplyr::left_join(want,
        by = "INFO") %>% dplyr::mutate(COL_TYPE = stringi::stri_replace_na(str = COL_TYPE,
        replacement = "_")) %>% dplyr::select(COL_TYPE)
        haplo.col.type[1, 1] <- "c"
        haplo.col.type <- purrr::flatten_chr(haplo.col.type) %>%
        stringi::stri_join(collapse = "")
        input <- readr::read_tsv(file = data, col_names = TRUE,
        na = "-", col_types = haplo.col.type)
        colnames(input) <- stringi::stri_replace_all_fixed(str = colnames(input),
        pattern = c("# Catalog ID", "Catalog ID", "# Catalog Locus ID"),
        replacement = c("LOCUS", "LOCUS", "LOCUS"), vectorize_all = FALSE)
        if (rlang::has_name(input, "Seg Dist")) {
            input <- dplyr::select(.data = input, -`Seg Dist`)
        }
        n.catalog.locus <- dplyr::n_distinct(input$LOCUS)
        n.individuals <- ncol(input) - 1
        message("\nNumber of loci in catalog: ", n.catalog.locus)
        message("Number of individuals: ", n.individuals)
        input <- tidyr::gather(data = input, key = "INDIVIDUALS",
        value = "GT_VCF_NUC", -LOCUS)
        input$INDIVIDUALS <- radiator::clean_ind_names(input$INDIVIDUALS)
        if (!is.null(whitelist.markers)) {
            input <- filter_whitelist(data = input, whitelist.markers = whitelist.markers)
        }
        if (verbose)
        message("\nScanning for consensus markers...")
        consensus.markers <- dplyr::filter(input, GT_VCF_NUC ==
        "consensus") %>% dplyr::distinct(LOCUS)
        if (length(consensus.markers$LOCUS) > 0) {
            input <- suppressWarnings(dplyr::anti_join(input,
            consensus.markers, by = "LOCUS"))
            readr::write_tsv(consensus.markers, "radiator.tidy.genomic.data.consensus.markers.tsv")
        }
        if (verbose)
        message("    number of consensus markers removed: ",
        dplyr::n_distinct(consensus.markers$LOCUS))
        consensus.markers <- NULL
        if (!is.null(strata)) {
            input %<>% join_strata(strata = strata.df)
            check.ref <- TRUE
        }
        if (verbose)
        message("Scanning for artifactual genotypes...")
        input <- input %>% dplyr::mutate(POLYMORPHISM = stringi::stri_count_fixed(GT_VCF_NUC,
        "/"))
        blacklist.paralogs <- input %>% dplyr::filter(POLYMORPHISM >
        1) %>% dplyr::select(LOCUS, INDIVIDUALS)
        if (verbose)
        message("    number of genotypes with more than 2 alleles: ",
        length(blacklist.paralogs$LOCUS))
        if (length(blacklist.paralogs$LOCUS) > 0) {
            input <- input %>% dplyr::mutate(GT_VCF_NUC = replace(GT_VCF_NUC,
            which(POLYMORPHISM > 1), NA)) %>% dplyr::select(-POLYMORPHISM)
            readr::write_tsv(blacklist.paralogs, "blacklist.genotypes.paralogs.tsv")
        }
        blacklist.paralogs <- NULL
        if (verbose)
        message("Calculating REF/ALT alleles...")
        input <- input %>% dplyr::mutate(GT_VCF_NUC = dplyr::if_else(POLYMORPHISM ==
        0, stringi::stri_join(GT_VCF_NUC, "/", GT_VCF_NUC),
        GT_VCF_NUC, missing = "./."), GT_VCF_NUC = dplyr::if_else(stringi::stri_detect_fixed(GT_VCF_NUC,
        "N"), "./.", GT_VCF_NUC)) %>% dplyr::select(-POLYMORPHISM)
        input.temp <- radiator::calibrate_alleles(data = input,
        biallelic = FALSE, parallel.core = parallel.core,
        verbose = verbose)
        input <- input.temp$input
        input.temp <- NULL
        biallelic <- FALSE
        input <- dplyr::rename(input, LOCUS = MARKERS)
    }
    if (data.type == "genepop.file") {
        if (verbose)
        message("Tidying the genepop file ...")
        input <- radiator::tidy_genepop(data = data, tidy = TRUE)
        skip.tidy.wide <- TRUE
    }
    if (data.type == "dart") {
        if (verbose)
        message("Tidying DArT data...")
        input <- radiator::read_dart(data = data, strata = strata,
        verbose = FALSE, parallel.core = parallel.core, tidy.dart = TRUE)
        skip.tidy.wide <- TRUE
    }
    if (data.type == "genind") {
        if (verbose)
        message("Tidying the genind object ...")
        input <- radiator::tidy_genind(data = data, gds = FALSE,
        keep.allele.names = keep.allele.names)
        data <- NULL
        input <- input %>% dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
        .funs = clean_ind_names) %>% dplyr::mutate_at(.tbl = .,
        .vars = "POP_ID", .funs = clean_pop_names)
        skip.tidy.wide <- TRUE
    }
    if (data.type == "genlight") {
        if (verbose)
        message("Tidying the genlight object ...")
        input <- radiator::tidy_genlight(data = data, gds = FALSE)
        data <- NULL
        input <- input %>% dplyr::mutate_at(.tbl = ., .vars = "INDIVIDUALS",
        .funs = clean_ind_names) %>% dplyr::mutate_at(.tbl = .,
        .vars = "POP_ID", .funs = clean_pop_names)
        biallelic <- TRUE
        skip.tidy.wide <- TRUE
    }
    if (data.type == "gtypes") {
        if (verbose)
        message("Tidying the gtypes object ...")
        input <- tidy_gtypes(data) %>% dplyr::mutate_at(.tbl = .,
        .vars = "INDIVIDUALS", .funs = clean_ind_names) %>%
        dplyr::mutate_at(.tbl = ., .vars = "POP_ID", .funs = clean_pop_names)
        data <- NULL
        skip.tidy.wide <- TRUE
    }
    if (data.type == "fst.file") {
        if (verbose)
        message("Importing the fst.file...")
        input <- read_rad(data = data)
        data.type <- "tbl_df"
        skip.tidy.wide <- TRUE
    }
    if (data.type == "tbl_df" || skip.tidy.wide) {
        if (!skip.tidy.wide) {
            if (verbose)
            message("Importing the data frame ...")
            input <- radiator::tidy_wide(data = data, import.metadata = TRUE)
            data <- NULL
        }
        if (!is.null(whitelist.markers)) {
            input <- filter_whitelist(data = input, whitelist.markers = whitelist.markers)
        }
        if (!is.null(blacklist.id)) {
            if (verbose)
            message("Filtering with blacklist of individuals")
            input %<>% dplyr::filter(!INDIVIDUALS %in% blacklist.id$INDIVIDUALS)
            check.ref <- TRUE
        }
        if (!is.null(strata)) {
            input %<>% join_strata(strata = strata.df)
            check.ref <- TRUE
        }
        if (rlang::has_name(input, "REF")) {
            check.ref <- FALSE
        }
        else {
            check.ref <- TRUE
        }
        if (check.ref) {
            input.temp <- radiator::calibrate_alleles(data = input)
            input <- input.temp$input
            biallelic <- input.temp$biallelic
            input.temp <- NULL
        }
        else {
            biallelic <- radiator::detect_biallelic_markers(data = input)
        }
    }
    if (!is.null(strata)) {
        strata.df <- generate_strata(input, pop.id = TRUE)
    }
    else {
        filter.common.markers <- FALSE
    }
    if (is.null(blacklist.genotypes)) {
        if (verbose)
        message("Erasing genotype: no")
    }
    else {
        input <- filter_blacklist_genotypes(data = input, blacklist.genotypes = blacklist.genotypes,
        verbose = verbose)
    }
    blacklist.id <- whitelist.markers <- whitelist.markers.ind <- NULL
    want <- blacklist.genotypes <- NULL
    filters.parameters <- radiator_parameters(generate = TRUE,
    initiate = TRUE, update = FALSE, parameter.obj = parameters,
    data = input, path.folder = path.folder, file.date = file.date,
    internal = FALSE, verbose = verbose)
    input <- filter_common_markers(data = input, filter.common.markers = filter.common.markers,
    verbose = verbose, path.folder = path.folder, parameters = filters.parameters,
    internal = TRUE)
    input <- filter_monomorphic(data = input, filter.monomorphic = filter.monomorphic,
    verbose = verbose, path.folder = path.folder, parameters = filters.parameters,
    internal = TRUE)
    if (!is.null(strata)) {
        input %<>% dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
    }
    else {
        input %<>% dplyr::arrange(INDIVIDUALS, MARKERS)
    }
    if (!is.null(filename)) {
        tidy.name <- stringi::stri_join(filename, ".rad")
        message("\nWriting tidy data set:\n", tidy.name)
        write_rad(data = input, path = file.path(path.folder,
        tidy.name))
    }
    n.markers <- length(unique(input$MARKERS))
    if (rlang::has_name(input, "CHROM")) {
        n.chromosome <- length(unique(input$CHROM))
    }
    else {
        n.chromosome <- "no chromosome info"
    }
    n.individuals <- length(unique(input$INDIVIDUALS))
    if (!is.null(strata))
    n.pop <- length(unique(input$POP_ID))
    if (verbose) {
        cat("################################### RESULTS ####################################\n")
        if (!is.null(filename)) {
            message("Tidy data written in global environment and working directory")
        }
        else {
            message("Tidy data written in global environment")
        }
        message("Data format: ", data.type)
        if (biallelic) {
            message("Biallelic data")
        }
        else {
            message("Multiallelic data")
        }
        message("\nTidy genomic data:")
        message("    Number of markers: ", n.markers)
        message("    Number of chromosome/contig/scaffold: ",
        n.chromosome)
        if (!is.null(strata))
        message("    Number of strata: ", n.pop)
        message("    Number of individuals: ", n.individuals)
    }
    return(input)
}

bayescan_one <- function (x = NULL, data, n = 5000, thin = 10, nbp = 20, pilot = 5000,
burn = 50000, pr_odds, subsample = NULL, iteration.subsample = 1,
parallel.core = number_cores, path.folder,
file.date, bayescan.path = bayescan.path){
    res <- list()
    if (!is.null(subsample)) {
        subsample.id <- unique(x$SUBSAMPLE)
        message("\nBayeScan, subsample: ", subsample.id, "\n")
        path.folder.subsample <- stringi::stri_join(path.folder,
        "/bayescan_subsample_", subsample.id)
        dir.create(file.path(path.folder.subsample))
        folder.message <- stringi::stri_join("bayescan_subsample_",
        subsample.id)
        message("Subsampling folder created: ", folder.message)
    }
    else {
        path.folder.subsample <- path.folder
    }
    output.folder <- stringi::stri_join("-od ", path.folder.subsample)
    log.file <- stringi::stri_join(path.folder.subsample, "/radiator_bayescan_",
    file.date, ".log")
    # message("For progress, look in the log file: radiator_bayescan_",
    # file.date, ".log")
    all.trace <- "-all_trace "
    parallel.core.bk <- parallel.core
    parallel.core <- stringi::stri_join("-threads ", parallel.core)
    n <- stringi::stri_join("-n ", n)
    thin <- stringi::stri_join("-thin ", thin)
    nbp <- stringi::stri_join("-nbp ", nbp)
    pilot <- stringi::stri_join("-pilot ", pilot)
    burn <- stringi::stri_join("-burn ", burn)
    pr.odds <- stringi::stri_join("-pr_odds ", pr_odds)
    if (!is.null(subsample)) {
        bayescan.filename <- stringi::stri_join("radiator_bayescan_subsample_",
        subsample.id)
        bayescan.sub <- radiator::write_bayescan(data = dplyr::semi_join(data,
        x, by = c("POP_ID", "INDIVIDUALS")), parallel.core = parallel.core.bk,
        filename = bayescan.filename)
        x <- NULL
        data <- stringi::stri_join(bayescan.filename, ".txt")
    }
    # message("Copying input BayeScan file in folder")
    link.problem <- stringi::stri_detect_fixed(str = data, pattern = getwd())
    if (link.problem) {
        new.data <- stringi::stri_replace_all_fixed(str = data,
        pattern = getwd(), replacement = "", vectorize_all = FALSE)
        new.data <- stringi::stri_join(path.folder.subsample,
        new.data)
    }
    else {
        new.data <- stringi::stri_join(path.folder.subsample,
        "/", data)
    }
    file.copy(from = data, to = new.data)
    if (!is.null(subsample)) {
        pop.dictionary <- bayescan.sub$pop.dictionary
        markers.dictionary <- bayescan.sub$markers.dictionary
    }
    else {
        pop.dic.file <- stringi::stri_replace_all_fixed(str = data,
        pattern = ".txt", "_pop_dictionary")
        pop.dic.file <- list.files(path = getwd(), pattern = pop.dic.file)
        markers.dic.file <- stringi::stri_replace_all_fixed(str = data,
        pattern = ".txt", "_markers_dictionary")
        markers.dic.file <- list.files(path = getwd(), pattern = markers.dic.file)
        if (length(pop.dic.file) > 0) {
            pop.dictionary <- readr::read_tsv(file = pop.dic.file,
            col_types = "ci")
            markers.dictionary <- readr::read_tsv(file = markers.dic.file,
            col_types = "ci")
            file.copy(from = pop.dic.file, to = stringi::stri_join(path.folder.subsample,
            "/", pop.dic.file))
            file.copy(from = markers.dic.file, to = stringi::stri_join(path.folder.subsample,
            "/", markers.dic.file))
        }
        else {
            pop.dictionary <- markers.dictionary <- NULL
        }
    }
    pop.dic.file <- markers.dic.file <- bayescan.sub <- NULL
    command.arguments <- paste(new.data, output.folder, all.trace,
    parallel.core, n, thin, nbp, pilot, burn, pr.odds)
    system2(command = bayescan.path, args = command.arguments,
    stderr = log.file, stdout = log.file)
    # message("Importing BayeScan results")
    res$bayescan <- suppressWarnings(readr::read_table2(file = list.files(path = path.folder.subsample,
    pattern = "_fst.txt", full.names = TRUE), skip = 1, col_names = c("BAYESCAN_MARKERS",
    "POST_PROB", "LOG10_PO", "Q_VALUE", "ALPHA", "FST"),
    col_types = c("iddddd"))) %>% dplyr::mutate(Q_VALUE = dplyr::if_else(Q_VALUE <=
    1e-04, 1e-04, Q_VALUE), Q_VALUE = round(Q_VALUE, 4),
    POST_PROB = round(POST_PROB, 4), LOG10_PO = round(LOG10_PO,
    4), ALPHA = round(ALPHA, 4), FST = round(FST, 6),
    SELECTION = factor(dplyr::if_else(ALPHA >= 0 & Q_VALUE <=
    0.05, "diversifying", dplyr::if_else(ALPHA >= 0 &
    Q_VALUE > 0.05, "neutral", "balancing"))), LOG10_Q = log10(Q_VALUE))
    if (!is.null(markers.dictionary)) {
        res$bayescan <- dplyr::right_join(markers.dictionary,
        res$bayescan, by = "BAYESCAN_MARKERS")
    }
    else {
        res$bayescan <- dplyr::mutate(res$bayescan, MARKERS = BAYESCAN_MARKERS)
    }
    res$bayescan <- res$bayescan %>% dplyr::mutate(PO_GROUP = factor(dplyr::if_else(LOG10_PO >
    2, "decisive", dplyr::if_else(LOG10_PO > 1.5, "very strong",
    dplyr::if_else(LOG10_PO > 1, "strong", dplyr::if_else(LOG10_PO >
    0.5, "substantial", "no evidence")))), levels = c("no evidence",
    "substantial", "strong", "very strong", "decisive"),
    ordered = TRUE)) %>% dplyr::ungroup(.) %>% dplyr::mutate(FST_GROUP = dplyr::ntile(FST,
    5), FST_GROUP = dplyr::if_else(FST_GROUP == 1, "0-20%",
    dplyr::if_else(FST_GROUP == 2, "20-40%", dplyr::if_else(FST_GROUP ==
    3, "40-60%", dplyr::if_else(FST_GROUP == 4, "60-80%",
    "80-100%"))))) %>% dplyr::arrange(FST)
    if (!is.null(subsample)) {
        res$bayescan <- dplyr::mutate(res$bayescan, ITERATIONS = rep(subsample.id,
        n()))
    }
    radiator.markers <- dplyr::distinct(res$bayescan, MARKERS)
    radiator.markers <- dplyr::filter(radiator.markers, !is.na(MARKERS))
    radiator.markers <- unique(stringi::stri_detect_fixed(str = radiator.markers$MARKERS,
    pattern = "__"))
    if (radiator.markers) {
        message("Detected SNP and LOCUS information in markers")
        res$bayescan <- res$bayescan %>% tidyr::separate(data = .,
        col = MARKERS, into = c("CHROM", "LOCUS", "POS"),
        sep = "__", remove = FALSE, extra = "warn")
        whitelist.multiple.snp <- dplyr::distinct(res$bayescan,
        LOCUS, POS) %>% dplyr::group_by(LOCUS) %>% dplyr::tally(.) %>%
        dplyr::filter(n > 1) %>% dplyr::select(LOCUS)
        markers.more.snp <- nrow(whitelist.multiple.snp)
        n.markers <- dplyr::n_distinct(res$bayescan$MARKERS)
        if (markers.more.snp > 0) {
            message("Detected markers > 1 SNP per LOCUS...")
            message("    total number of markers: ", n.markers)
            message("    markers with >1 SNPs/LOCUS: ", markers.more.snp,
            " (", round(markers.more.snp/n.markers, 2), ")")
            message("\nCalculating accuracy within LOCUS...")
            locus.accuracy <- dplyr::left_join(whitelist.multiple.snp,
            res$bayescan, by = "LOCUS") %>% dplyr::select(-BAYESCAN_MARKERS,
            -MARKERS) %>% dplyr::group_by(LOCUS, SELECTION) %>%
            dplyr::tally(.) %>% dplyr::group_by(LOCUS) %>%
            dplyr::mutate(SNP_NUMBER = sum(n), ACCURACY = dplyr::if_else(n ==
            SNP_NUMBER, "accurate", "not accurate"))
            accurate.locus <- locus.accuracy %>% dplyr::filter(ACCURACY ==
            "accurate") %>% dplyr::select(LOCUS, SELECTION,
            SNP_NUMBER)
            n.accurate.locus <- nrow(accurate.locus)
            n.not.accurate.locus <- markers.more.snp - n.accurate.locus
            message("Number of locus accurate: ", n.accurate.locus,
            " (", round(n.accurate.locus/markers.more.snp,
            2), ")")
            message("Number of locus NOT accurate: ", n.not.accurate.locus,
            " (", round(n.not.accurate.locus/markers.more.snp,
            2), ")")
            res$accurate.locus.summary <- accurate.locus %>%
            dplyr::group_by(SELECTION) %>% dplyr::tally(.)
            res$whitelist.accurate.locus <- dplyr::distinct(accurate.locus,
            LOCUS) %>% dplyr::arrange(LOCUS)
            readr::write_tsv(x = res$whitelist.accurate.locus,
            path = file.path(path.folder.subsample, "whitelist.accurate.locus.tsv"))
            res$blacklist.not.accurate.locus <- locus.accuracy %>%
            dplyr::filter(ACCURACY == "not accurate") %>%
            dplyr::distinct(LOCUS) %>% dplyr::arrange(LOCUS)
            readr::write_tsv(x = res$blacklist.not.accurate.locus,
            path = file.path(path.folder.subsample, "blacklist.not.accurate.locus.tsv"))
            res$accuracy.snp.number <- locus.accuracy %>% dplyr::distinct(LOCUS,
            SNP_NUMBER, ACCURACY) %>% dplyr::group_by(SNP_NUMBER,
            ACCURACY) %>% dplyr::tally(.)
            readr::write_tsv(x = res$accuracy.snp.number, path = file.path(path.folder.subsample,
            "accuracy.snp.number.tsv"))
            # res$accuracy.snp.number.plot <- ggplot2::ggplot(res$accuracy.snp.number,
            #     ggplot2::aes(y = n, x = SNP_NUMBER, fill = ACCURACY)) +
            #     ggplot2::geom_bar(stat = "identity") + ggplot2::labs(y = "Number of locus") +
            #     ggplot2::labs(x = "Number of SNPs per locus") +
            #     ggplot2::theme(axis.title.x = ggplot2::element_text(size = 12,
            #       family = "Helvetica", face = "bold"), axis.title.y = ggplot2::element_text(size = 12,
            #       family = "Helvetica", face = "bold"), legend.title = ggplot2::element_text(size = 12,
            #       family = "Helvetica", face = "bold"), legend.text = ggplot2::element_text(size = 12,
            #       family = "Helvetica", face = "bold"), strip.text.x = ggplot2::element_text(size = 12,
            #       family = "Helvetica", face = "bold"))
            # ggplot2::ggsave(filename = file.path(path.folder.subsample,
            #     "accuracy.snp.number.plot.pdf"), plot = res$accuracy.snp.number.plot,
            #     width = 20, height = 15, dpi = 600, units = "cm",
            #     useDingbats = FALSE)
            res$not.accurate.summary <- locus.accuracy %>% dplyr::filter(ACCURACY ==
            "not accurate") %>% dplyr::group_by(LOCUS) %>%
            dplyr::summarise(SELECTION_TYPE_ON_LOCUS = stringi::stri_join(SELECTION,
            collapse = " <-> ")) %>% dplyr::group_by(SELECTION_TYPE_ON_LOCUS) %>%
            dplyr::tally(.) %>% dplyr::rename(LOCUS_NUMBER = n) %>%
            dplyr::mutate(PROP = round(LOCUS_NUMBER/sum(LOCUS_NUMBER),
            4), PROP_TOTAL_MARKERS = round(LOCUS_NUMBER/markers.more.snp,
            4))
            readr::write_tsv(x = res$not.accurate.summary, path = file.path(path.folder.subsample,
            "not.accurate.summary.tsv"))
        }
    }
    radiator.markers <- NULL
    res$whitelist.markers.positive.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "diversifying" & PO_GROUP !=
    "no evidence") %>% dplyr::distinct(MARKERS) %>% dplyr::arrange(MARKERS)
    if (!is.null(subsample)) {
        res$whitelist.markers.positive.selection <- dplyr::mutate(res$whitelist.markers.positive.selection,
        ITERATIONS = rep(subsample.id, n()))
    }
    readr::write_tsv(x = res$whitelist.markers.positive.selection,
    path = file.path(path.folder.subsample, "whitelist.markers.positive.selection.tsv"))
    res$whitelist.markers.neutral.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "neutral") %>% dplyr::distinct(MARKERS) %>%
    dplyr::arrange(MARKERS)
    if (!is.null(subsample)) {
        res$whitelist.markers.neutral.selection <- dplyr::mutate(res$whitelist.markers.neutral.selection,
        ITERATIONS = rep(subsample.id, n()))
    }
    readr::write_tsv(x = res$whitelist.markers.neutral.selection,
    path = file.path(path.folder.subsample, "whitelist.markers.neutral.selection.tsv"))
    res$whitelist.markers.neutral.positive.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "neutral" | (SELECTION ==
    "diversifying" & PO_GROUP != "no evidence")) %>%
    dplyr::distinct(MARKERS) %>% dplyr::arrange(MARKERS)
    if (!is.null(subsample)) {
        res$whitelist.markers.neutral.positive.selection <- dplyr::mutate(res$whitelist.markers.neutral.positive.selection,
        ITERATIONS = rep(subsample.id, n()))
    }
    readr::write_tsv(x = res$whitelist.markers.neutral.positive.selection,
    path = file.path(path.folder.subsample, "whitelist.markers.neutral.positive.selection.tsv"))
    res$blacklist.markers.balancing.selection <- res$bayescan %>%
    dplyr::filter(SELECTION == "balancing") %>% dplyr::distinct(MARKERS) %>%
    dplyr::arrange(MARKERS)
    if (!is.null(subsample)) {
        res$blacklist.markers.balancing.selection <- dplyr::mutate(res$blacklist.markers.balancing.selection,
        ITERATIONS = rep(subsample.id, n()))
    }
    readr::write_tsv(x = res$blacklist.markers.balancing.selection,
    path = file.path(path.folder.subsample, "blacklist.markers.balancing.selection.tsv"))
    selection <- dplyr::group_by(res$bayescan, SELECTION, PO_GROUP) %>%
    dplyr::tally(.) %>% dplyr::rename(MARKERS = n)
    if (!is.null(subsample)) {
        selection <- dplyr::mutate(selection, ITERATIONS = rep(subsample.id,
        n()))
    }
    readr::write_tsv(x = selection, path = file.path(path.folder.subsample,
    "selection.summary.tsv"))
    # message("Generating plot")
    # res$bayescan.plot <- plot_bayescan(res$bayescan)
    # if (!is.null(subsample)) {
    #     temp.name <- stringi::stri_join("bayescan_plot_", subsample.id,
    #         ".pdf")
    #     ggplot2::ggsave(filename = file.path(path.folder.subsample,
    #         temp.name), plot = res$bayescan.plot, width = 30,
    #         height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
    # }
    # else {
    #     ggplot2::ggsave(filename = file.path(path.folder.subsample,
    #         "bayescan_plot.pdf"), plot = res$bayescan.plot, width = 30,
    #         height = 15, dpi = 600, units = "cm", useDingbats = FALSE)
    # }
    if (!is.null(subsample)) {
        temp.name <- stringi::stri_join("bayescan_", subsample.id,
        ".tsv")
        readr::write_tsv(x = res$bayescan, path = file.path(path.folder.subsample,
        temp.name))
    }
    else {
        readr::write_tsv(x = res$bayescan, path = file.path(path.folder.subsample,
        "bayescan.tsv"))
    }
    res$selection.summary <- selection
    if (!is.null(markers.dictionary))
    res$markers.dictionary <- markers.dictionary
    if (!is.null(pop.dictionary))
    res$pop.dictionary <- pop.dictionary
    return(res)
}

clean_ind_names <- function (x){
    x <- stringi::stri_replace_all_fixed(str = as.character(x),
    pattern = c("_", ":", " ", ","), replacement = c("-",
    "-", "", ""), vectorize_all = FALSE)
    x <- stringi::stri_replace_all_regex(str = x, pattern = "\\s+",
    replacement = "", vectorize_all = FALSE)
}

plot_bayescan <- function (data) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = LOG10_Q, y = FST)) +
    ggplot2::geom_point(ggplot2::aes(colour = PO_GROUP, shape = FST_GROUP)) +
    ggplot2::scale_shape_manual(name = "FST quantile group",
    values = c(5, 2, 3, 4, 1)) + ggplot2::scale_colour_manual(name = "Model choice",
    values = c("darkred", "yellow", "orange", "green", "forestgreen")) +
    ggplot2::labs(x = "Log10(Q_VALUE)") + ggplot2::labs(y = "Fst") +
    ggplot2::geom_vline(xintercept = c(log10(0.05)), color = "black") +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 16,
    family = "Helvetica", face = "bold"), legend.title = ggplot2::element_text(size = 16,
    family = "Helvetica", face = "bold"), legend.text = ggplot2::element_text(size = 16,
    family = "Helvetica", face = "bold"), legend.position = "right")
    return(plot)
}
