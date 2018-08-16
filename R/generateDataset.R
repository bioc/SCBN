#' Generate simulation data for different species
#' @description To generate RNA-seq genes between different species.
#' @usage generateDataset(commonTags=15000, uniqueTags=c(1000, 3000),
#'                        unmapped=c(4000, 2000),group=c(1, 2),
#'                        libLimits=c(.9, 1.1)*1e6, empiricalDist=NULL,
#'                        genelength, randomRate=1/100,
#'                        pDifferential=.05, pUp=.5, foldDifference=2)
#' @param commonTags The number of genes have the same expression level.
#' @param uniqueTags The number of genes only expressed in one species.
#' @param unmapped The number of genes only in one species.
#' @param group The number of species.
#' @param libLimits The limits for two species.
#' @param empiricalDist Define where to take random sample from (empirical
#' distribution OR random exponential), if NULL, the reads take from random
#' exponential.
#' @param genelength A vector of gene length for each gene of two species.
#' @param randomRate The parameter for exponential distribution.
#' @param pDifferential The propotion of differential expression genes.
#' @param pUp The probably for the reads in first species fold than the second
#' species.
#' @param foldDifference The fold for fold expression genes.
#' @return list(.) A list of output, "DATAN" represents the read counts for the
#' first species, "DATAM" represents the read counts for the second species,
#' "trueFactors" represents the true scaling factor for data, "group" represents
#' the number of species, "libSizes" represents the library size for data,
#' "differentialInd" represents the ID for differential expression genes,
#' "commonInd" represents the ID for common expression genes.
#' @examples
#' data(orthgenes)
#' orthgenes[, 6:9] <- round(orthgenes[, 6:9])
#' orthgenes1 <- orthgenes[!(is.na(orthgenes[,6])|is.na(orthgenes[,7])|
#'                        is.na(orthgenes[,8])|is.na(orthgenes[,9])), ]
#' sim_data <- generateDataset(commonTags=5000, uniqueTags=c(100, 300),
#'                             unmapped=c(400, 200),group=c(1, 2),
#'                             libLimits=c(.9, 1.1)*1e6,
#'                             empiricalDist=orthgenes1[, 6],
#'                             genelength=orthgenes1[, 2], randomRate=1/100,
#'                             pDifferential=.05, pUp=.5, foldDifference=2)
#' @export


generateDataset <- function(commonTags=15000, uniqueTags=c(1000, 3000),
                            unmapped=c(4000, 2000), group=c(1, 2),
                            libLimits=c(.9, 1.1)*1e6, empiricalDist=NULL,
                            genelength, randomRate=1/100,
                            pDifferential=.05, pUp=.5, foldDifference=2) {

  group <- as.factor(group)
  stopifnot(length(group) == length(uniqueTags))
  stopifnot(nlevels(group) == 2 )

  if (is.null(empiricalDist))
      exampleCounts <- ceiling(rexp(commonTags, rate=randomRate))
  else
      exampleCounts <- empiricalDist
  exampleLambda <- exampleCounts/sum(exampleCounts)/genelength

  nLibraries <- length(uniqueTags)
  libSizes <- runif(nLibraries, min=libLimits[1], max=libLimits[2] )

  en <- commonTags+cumsum(uniqueTags)
  st <- c(commonTags+1, en[-nLibraries]+1)

  LAMBDA <- matrix(0, nrow=max(en), ncol=nLibraries)
  lengthsampled <- matrix(0, nrow=max(en), ncol=nLibraries)
  idd <- seq_len(length(exampleLambda))
  idd1 <- sample(idd, commonTags, replace=TRUE)

  lengthsampled[seq_len(commonTags), ] <- genelength[idd1]
  LAMBDA[seq_len(commonTags), ] <- exampleLambda[idd1]

  for(i in seq_len(nLibraries))
      if(uniqueTags[i] > 0) {
         uniqueind<-sample(idd, uniqueTags[i])
         LAMBDA[st[i]:en[i], i] <- exampleLambda[uniqueind]
         lengthsampled[st[i]:en[i], ] <- genelength[uniqueind]
      }
  g <- group == levels(group)[1]
  LAMBDA <- cbind(LAMBDA[, g]*lengthsampled[, 1],
                  LAMBDA[, !g]*lengthsampled[, 2])

  ind <- seq_len(floor(pDifferential*commonTags))
  if (length(ind)>0) {
      fcDir <- sample(c(-1, 1), length(ind), prob=c(1-pUp, pUp), replace=TRUE)
      LAMBDA[ind, g] <- LAMBDA[ind, !g]*exp(log(foldDifference)/2*fcDir)
      LAMBDA[ind, !g] <- LAMBDA[ind, !g]*exp(log(foldDifference)/2*(-fcDir))
  }
  lamd <- cbind(lengthsampled,LAMBDA)

  uniqind1 <- sample(seq_len(length(lamd[, 1])), unmapped[1])
  uniqind2 <- sample(seq_len(length(lamd[, 1])), unmapped[2])
  lamd1 <- lamd[uniqind1, ]
  lamd2 <- lamd[uniqind2, ]
  humanlamd <- rbind(lamd[, c(1, 3)], lamd1[, c(1, 3)])
  mouselamd <- rbind(lamd[, c(2, 4)], lamd2[, c(2, 4)])

  sampFactors <- c(sum(humanlamd[, 2]), sum(mouselamd[, 2]))
  MEANh <- humanlamd[, 2]/sampFactors[1]*libSizes[1]
  MEANm <- mouselamd[, 2]/sampFactors[2]*libSizes[2]

  obsh <- rpois(length(MEANh), lambda=MEANh)
  obsm <- rpois(length(MEANm), lambda=MEANm)
  DATAH <- cbind(humanlamd[, 1], obsh)
  DATAM <- cbind(mouselamd[, 1], obsm)
  colnames(DATAH) <- c("genelength", "count")
  colnames(DATAM) <- c("genelength", "count")

  rsum <- DATAH[seq_len(length(LAMBDA[, 1])), 2]+
          DATAM[seq_len(length(LAMBDA[, 1])), 2]
  zero <- sum(rsum[ind] == 0)
  allzero <- sum(rsum[seq_len(commonTags)] == 0)
  uniquezero <- sum(rsum == 0)
  DATAH <- rbind(DATAH[which(rsum!=0), ],
                 DATAH[(length(LAMBDA[, 1])+1):length(DATAH[, 1]), ])
  DATAM <- rbind(DATAM[which(rsum!=0), ],
                 DATAM[(length(LAMBDA[, 1])+1):length(DATAM[, 1]), ])

  trueFactors <- sampFactors/sampFactors[1]
  list(DATAH=DATAH,DATAM=DATAM, trueFactors=trueFactors, group=group,
       libSizes=libSizes, commonInd=seq_len((commonTags-allzero)),
       differentialInd=c(seq_len((length(ind)-zero)),
       (commonTags-allzero+1):(nrow(lamd)-uniquezero)))
}
