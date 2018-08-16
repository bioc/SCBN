#' Iteration to find the optimal value
#' A iteration process to compute the normalization factor to identify
#' difference expression(DE) of genes between different species
#' @describeIn To obtain the optimal normalization value.
#' @usage Iter_optimal(scale, orth_gene, hkind, a)
#' @param a P-value cutoff in iteration process to find the optimal
#' normalization factor.
#' @return factor Computed normalization factor.
#' @examples
#' data(sim_data)
#' scale <- MediancalcNorm(orth_gene=sim_data, hkind=1:1000)
#' Iter_optimal(scale=scale, orth_gene=sim_data, hkind=1:1000, a=0.05)
#' @export

Iter_optimal <- function(scale, orth_gene, hkind, a) {
  if (scale > 0.5) {
      fW1 <- seq(scale-0.5, scale + 0.5, 0.1)
  } else {
      fW1 <- seq(0.02, scale +0.5, 0.1)
  }
  n1 <- length(fW1)
  fdr1 <- rep(0, n1)
  for(j in seq_len(n1)){
      sMm1 <- sageTestNew(orth_gene[, 2], orth_gene[, 4], orth_gene[, 1],
                          orth_gene[, 3], sum(orth_gene[, 2]),
                          sum(orth_gene[, 4]), fW1[j])
      fdr1[j] <- sum(sMm1[hkind] < a)/length(hkind)
  }
  fdr11 <- abs(fdr1-a)
  counts1 <- length(which(fdr11 == min(fdr11)))
  if (counts1 %% 2 == 1) {
      fw2 <- median(fW1[which(fdr11 == min(fdr11))])
  } else {
      fw2 <- fW1[which(fdr11 == min(fdr11))][counts1/2+1]
  }
  if (fw2 > 0.5) {
      fW3 <- seq(fw2-0.5, fw2+0.5, 0.05)
  } else {
      fW3 <- seq(0.02, fw2+0.5, 0.05)
  }
  n2 <- length(fW3)
  fdr2 <- rep(0,n2)
  for(j in seq_len(n2)){
      sMm2 <- sageTestNew(orth_gene[, 2], orth_gene[, 4], orth_gene[,1],
                          orth_gene[, 3], sum(orth_gene[, 2]),
                          sum(orth_gene[, 4]), fW3[j])
      fdr2[j] <- sum(sMm2[hkind] < a)/length(hkind)
  }
  fdr22 <- abs(fdr2-a)
  counts2 <- length(which(fdr22 == min(fdr22)))
  if (counts2%%2 == 1) {
      fw4 <- median(fW3[which(fdr22 == min(fdr22))])
  } else {
      fw4 <- fW3[which(fdr22 == min(fdr22))][counts2/2+1]
  }
  if (fw4 > 0.25) {
      fW5 <- seq(fw4-0.25,fw4+0.25,0.005)
  } else {
     fW5 <- seq(0.02, fw4+0.25, 0.005)
  }
  n3 <- length(fW5)
  fdr3 <- rep(0,n3)
  for (j in seq_len(n3)) {
       sMm3 <- sageTestNew(orth_gene[, 2], orth_gene[, 4], orth_gene[, 1],
                           orth_gene[, 3], sum(orth_gene[, 2]),
                           sum(orth_gene[, 4]), fW5[j])
       fdr3[j] <- sum(sMm3[hkind] < a)/length(hkind)
  }
  fdr33 <- abs(fdr3-a)
  counts3 <- length(which(fdr33 == min(fdr33)))
  if (counts3%%2 == 1) {
      fw6 <- median(fW5[which(fdr33 == min(fdr33))])
  } else {
      fw6 <- fW5[which(fdr33 == min(fdr33))][counts3/2+1]
  }
  factor <- fw6
  return(factor)
}
