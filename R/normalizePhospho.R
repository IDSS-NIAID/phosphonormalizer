#' @title Pairwise Normalization of MS-based phosphoproteomic data
#' @description This function compensates for the bias introduced in global phosphorylation in the sample after using median normalization.
#'
#' @param enriched The enriched data with the type data.frame or MSnSet, which should contain the sequence, modification of the sequence with their phosphorylation site and their abundances across samples
#' @param non.enriched The non-enriched data with the type data.frame or MSnSet, which should contain the sequence, modification of the sequence with their phosphorylation site and their abundances across samples
#' @param phospho a string that shows the term that represents phosphorylation in the modification column of the data. If it is not assigned, "Phospho" will be used as the default value
#' @param samplesCols A data.frame with two columns, with the column names enriched and non.enriched, of type numeric or integer, which must contain the column number of samples that hold the abundances
#' @param modseqCols A data.frame with two columns, with the names enriched and non.enriched, of type numeric or integer, which must contain the column number of samples that hold the sequence and modifications of the peptides
#' @param techRep a factor that holds information about columns order and the technical replicates of the samples
#' @param plot.fc This parameter if set plots the fold change distribution before and after pairwise normalization. controls and samples should be set as named vectors in a list (look at the example)
#' @param prop_good Numeric proportion of good values to keep (see Details).
#' @param llod Numeric lower limit of detection. Any value at or below this level will be considered missing.
#'
#' @details
#' It is shown that global median normalization can introduce bias in the fold change of global phosphorylation between samples. It is suggested that by taking the non-enriched data into consideration, this bias could be compensated (Kauko et al. 2015).
#'
#' By default, only rows with complete data will be analyzed, but this can be relaxed by setting the `prop_good` parameter to a lower value.
#' Setting the value at 0.5, for example, will allow samples with up to half of the replicates missing or below the lower limit of detection.
#'
#' @return A data.frame with the normalized values and their sequence and modification.
#'
#' @references https://doi.org/10.1093/bioinformatics/btx573, https://www.nature.com/articles/srep13099
#' @seealso [MSnbase](https://bioconductor.org/packages/release/bioc/html/MSnbase.html)
#' @examples
#' #Specifying the column numbers of abundances in the original data.frame,
#' #from both enriched and non-enriched runs
#' samplesCols <- data.frame(enriched=3:17, non.enriched=3:17)
#'
#' #Specifying the column numbers of sequence and modification in the original data.frame,
#' #from both enriched and non-enriched runs
#' modseqCols <- data.frame(enriched = 1:2, non.enriched = 1:2)
#'
#' #The samples and their technical replicates
#' techRep <- factor(x = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5))
#'
#' #Call the function
#' norm <- normalizePhospho(enriched = enriched.rd, non.enriched = non.enriched.rd,
#'                          samplesCols = samplesCols, modseqCols = modseqCols, techRep = techRep,
#'                          plot.fc = list(control = c(1,2), samples = c(3,4,5)))
#' head(norm)
#'
#' @export
#' @importFrom plyr ddply
#' @importFrom graphics abline boxplot
#' @importFrom matrixStats rowMaxs colMedians
#' @importFrom methods is
#' @importFrom stats median
normalizePhospho <- function(enriched, non.enriched, phospho = NULL, samplesCols, modseqCols, techRep, plot.fc=NULL,
                             prop_good = 1, llod = 0)
{
    # take care of annoying no visible binding note
    if(FALSE)
        mod <- modSeq <- NULL

    ##### Argument checks #####

    #Check if all the necessary arguments are present
    if(missing(enriched)) stop("The function parameter (enriched) is missing!")
    if(missing(non.enriched)) stop("The function parameter (non.enriched) is missing!")
    if(missing(samplesCols)) stop("The function parameter (samplesCols) is missing!")
    if(missing(modseqCols)) stop("The function parameter (modseqCols) is missing!")
    if(missing(techRep)) stop("The function parameter (techRep) is missing!")

    #Check the type of the arguments
    class.dfCols <- function(df, types)
    {
        all(apply(df, MARGIN = 2, function(x) class(x)[1]) %in% types)
    }
    if(!is(enriched, "data.frame") & !is(enriched, "MSnSet"))
        stop("The enriched parameter must be of either types of data.frame or MSnSet")

    if(!is(non.enriched, "data.frame") & !is(non.enriched, "MSnSet"))
        stop("The non.enriched parameter must be of either types of data.frame or MSnSet")

    if(is(enriched, "MSnSet"))
    {
        enriched.eset <- enriched
        enriched <- cbind(MSnbase::fData(enriched)[, modseqCols$enriched], MSnbase::exprs(enriched)[, samplesCols$enriched])
        non.enriched <- cbind(MSnbase::fData(non.enriched)[, modseqCols$non.enriched], MSnbase::exprs(non.enriched)[, samplesCols$non.enriched])
        modseqCols <- data.frame(enriched = c(1,2), non.enriched = c(1,2))
        samplesCols <- data.frame(enriched = 3:ncol(enriched), non.enriched = 3:ncol(non.enriched))
    }

    if(!is.null(phospho) & !methods::is(phospho, "character"))
        stop("The phospho parameter must be of type of character")

    if(!methods::is(samplesCols, "data.frame") |
       ncol(samplesCols) != 2 | !class.dfCols(samplesCols, c("integer","numeric")) |
       all(!colnames(samplesCols) %in% c("enriched", "non.enriched")))
        stop("The samplesCols must be data.frame with two columns, with the column names enriched and non.enriched
             , of type numeric or integer, which must contain the column number of samples that hold the abundances")

    if(!methods::is(modseqCols, "data.frame") |
       ncol(modseqCols) != 2 | !class.dfCols(modseqCols, c("integer","numeric")) |
       all(!colnames(modseqCols) %in% c("enriched", "non.enriched")))
        stop("The modseqCols must be data.frame with two columns, with the names enriched and non.enriched
             , of type numeric or integer, which must contain the column number of samples that hold the sequence
             and modifications of the peptides")

    if(any(!class.dfCols(enriched[, modseqCols$enriched], "character") ,
           !class.dfCols(non.enriched[, modseqCols$non.enriched], "character")))
        stop("The sequence and modification columns that is specified are not of type charachter!")

    if(any(!class.dfCols(enriched[, samplesCols$enriched], c("numeric", "integer")) ,
           !class.dfCols(non.enriched[, samplesCols$non.enriched], c("numeric", "integer"))))
        stop("The samples specified are not of type of numeric")


    ##### Missingness checks and Filters #####

    # convert values at or below `llod` to NA
    enriched[,samplesCols$enriched] <- replace(enriched[,samplesCols$enriched],
                                               enriched[,samplesCols$enriched] <= llod,
                                               NA)

    non.enriched[,samplesCols$non.enriched] <- replace(non.enriched[,samplesCols$non.enriched],
                                                       non.enriched[,samplesCols$non.enriched] <= llod,
                                                       NA)


    # drop values within each sample with too many missing or low values
    for(tr in levels(techRep))
    {
        tr_num <- techRep[techRep == tr] |> unique() |> as.numeric()

        drop <- check_missingness(enriched,
                                  samplesCols$enriched[techRep == tr_num],
                                  prop_good,
                                  llod)

        enriched[drop, samplesCols$enriched[techRep == tr_num]] <- NA

        drop <- check_missingness(non.enriched,
                                  samplesCols$non.enriched[techRep == tr_num],
                                  prop_good,
                                  llod)

        non.enriched[drop, samplesCols$non.enriched[techRep == tr_num]] <- NA
    }


    # keep these for later
    enriched.original.mat <- as.matrix(enriched[, samplesCols$enriched])
    seqMod <- enriched[, modseqCols$enriched]


    # drop rows with too many missing or low values
    too_many_missing <- apply(enriched[,samplesCols$enriched], 1, function(.x) all(is.na(.x) | .x == 0))
    enriched <- enriched[!too_many_missing,]

    too_many_missing <- apply(non.enriched[,samplesCols$non.enriched], 1, function(.x) all(is.na(.x) | .x == 0))
    non.enriched <- non.enriched[!too_many_missing,]


    # housekeeping of seq and mod columns
    colnames(enriched)[modseqCols$enriched] <- c("seq", "mod")
    colnames(non.enriched)[modseqCols$non.enriched] <- c("seq", "mod")

    enriched[,"modSeq"] <- paste(enriched$seq, enriched$mod,sep = ", ")
    non.enriched[,"modSeq"] <- paste(non.enriched$seq, non.enriched$mod,sep = ", ")


    # filter out phospho modifications
    if(is.null(phospho))
    {
        enriched <- enriched[grepl(pattern = "Phospho", x = enriched$mod),]
        non.enriched <- non.enriched[grepl(pattern = "Phospho", x = non.enriched$mod),]
    } else {
        enriched <- enriched[grepl(pattern = phospho, x = enriched$mod),]
        non.enriched <- non.enriched[grepl(pattern = phospho, x = non.enriched$mod),]
    }


    # column sums for each combination of seq and mod
    enriched     <- plyr::ddply(    enriched, plyr::.(seq, mod, modSeq),
                                function(df) colSums(df[, samplesCols$enriched    ], na.rm = TRUE))
    non.enriched <- plyr::ddply(non.enriched, plyr::.(seq, mod, modSeq),
                                function(df) colSums(df[, samplesCols$non.enriched], na.rm = TRUE))

    # convert any 0s back to NAs
    enriched[enriched == 0] <- NA
    non.enriched[non.enriched == 0] <- NA


    # find intersection between enriched and non-enriched
    inter <- intersect(non.enriched$modSeq, enriched$modSeq)
    stopifnot(length(inter) > 0)
    enriched.olp.idx <- which(enriched$modSeq %in% inter)
    non.enriched.olp.idx <- which(non.enriched$modSeq %in% inter)

    enriched.mat <- enriched[enriched.olp.idx, ]
    non.enriched.mat <- non.enriched[non.enriched.olp.idx, ]


    # sort rows by modSeq and convert to numeric matrix
    enriched.mat     <-     enriched.mat[order(    enriched.mat$modSeq),]
    non.enriched.mat <- non.enriched.mat[order(non.enriched.mat$modSeq),]

                                                     # drop seq, mod and modSeq columns
    enriched.mat     <- as.matrix(    enriched.mat[, -c(1:3)])
    non.enriched.mat <- as.matrix(non.enriched.mat[, -c(1:3)])


    ##### Ratios and Averages #####

    if(length(inter) > 1) {
        # calculate ratios
        ratios <- non.enriched.mat/enriched.mat
        colnames(ratios) <- as.character(techRep)

        # double check for non-finite values
        ratios[!is.finite(ratios)] <- NA

        # this is where we will collect average ratios for each condition in techRep
        ratios.avg <- matrix(nrow = nrow(ratios), ncol = length(levels(techRep)),
                             dimnames = list(NULL, levels(techRep)))

        # calculate average ratios for each condition in techRep
        for (tr in levels(techRep)) {
            if(nrow(ratios.avg) == 1) {
                ratios.avg[,tr] <- mean(ratios[,colnames(ratios) == tr], na.rm = TRUE)
            } else {
                ratios.avg[,tr] <- rowMeans(ratios[,colnames(ratios) == tr], na.rm = TRUE)
            }
        }

        # calculate max fold change
        # i.e. find the phosphopeptides with the largest variation between total and enriched samples
        # we will use these to filter outlier rows from `ratios`
        max.fc <- log2(ratios.avg[,1]) - log2(ratios.avg[,2])

        if(ncol(ratios.avg) > 2)
        {
            for (i in 2:(ncol(ratios.avg)-1)) {
                for (j in (i+1):(ncol(ratios.avg))) {
                    max.fc <- matrixStats::rowMaxs(cbind(max.fc, log2(ratios.avg[,i]) - log2(ratios.avg[,j])), na.rm = TRUE)
                }
            }
        }

        # check for non-finite values in max.fc, stemming from NAs
        max.fc[!is.finite(max.fc)] <- NA

        # remove outlier rows from `ratios`
        boxp <- boxplot(max.fc, plot = FALSE)
        ratios <- ratios[!(is.na(max.fc) | max.fc > max(boxp$stats)),]

        # calculate rowMeans on log10 scale
        # this will be used to center the ratios
        lratios <- log10(ratios)
        if(methods::is(lratios, "matrix") | methods::is(lratios, "data.frame")) {
            col.sub <- rowMeans(lratios, na.rm = TRUE)
        } else {
            col.sub <- mean(lratios, na.rm = TRUE)
        }

        # center ratios
        lratios.norm <- lratios - col.sub

        # calculate column medians on log10 scale
        # these will be used as the normalization factors for each sample
        if(methods::is(lratios, "matrix") | methods::is(lratios, "data.frame")) {
            factors <- 10^(matrixStats::colMedians(lratios.norm, na.rm = TRUE))
        } else {
            factors <- 10^lratios.norm
        }
    } else {
        factors <- as.numeric(non.enriched.mat/enriched.mat)
    }

    enriched.normalized.mat <- t(t(enriched.original.mat) * factors)
    if(!is.null(plot.fc)) {
        for(i in plot.fc$control) {
            tr_i <- techRep == i
            for(j in plot.fc$samples) {
                tr_j <- techRep == j
                a.original <- rowMeans(log2(enriched.original.mat[,tr_i]+1),na.rm=TRUE)
                b.original <- rowMeans(log2(enriched.original.mat[,tr_j]+1),na.rm=TRUE)
                fc.original <- a.original - b.original
                a.normnalized <- rowMeans(log2(enriched.normalized.mat[,tr_i]+1),na.rm=TRUE)
                b.normnalized <- rowMeans(log2(enriched.normalized.mat[,tr_j]+1),na.rm=TRUE)
                fc.normnalized <- a.normnalized - b.normnalized
                boxplot(cbind(fc.original, fc.normnalized), range=1.5, outline=FALSE, main=paste0("Peptide log fold changes", " (sample ", j, " vs sample ", i, ")"), names=c("Median normalized","Pairwise normalized"))
                abline(h=0, lty=2)}
        }

    }
    message(paste0("The number of peptides in the intersect is: ", length(inter)))
    message(paste0(length(plot.fc$control) * length(plot.fc$samples), " plots generated. Browse through them."))
    data.frame(seqMod, t(t(enriched.original.mat) * factors))
}
