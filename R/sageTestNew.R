#' Compute the false discovery rate
#' Compute the p-value for each orthologous genes between different species
#' @describeIn To obtain the p-value for each orthologous genes between
#' different species.
#' @usage sageTestNew(x, y, lengthx, lengthy, n1, n2, scale)
#' @param x The read counts for the first species.
#' @param y The read counts for the second species.
#' @param lengthx The gene length for the first species.
#' @param lengthy The gene length for the second species.
#' @param n1 The total read counts for the first species.
#' @param n2 The total read counts for the second species.
#' @param scale A value for normalization factor.
#' @return p_value P-values for each orthologous genes between different
#' species.
#' @examples
#' data(sim_data)
#' orth_gene <- sim_data
#' hkind <- 1:1000
#' scale <- MediancalcNorm(orth_gene=orth_gene, hkind=hkind)
#' x <- orth_gene[, 2]
#' y <- orth_gene[, 4]
#' lengthx <- orth_gene[, 1]
#' lengthy <- orth_gene[, 3]
#' n1 <- sum(x)
#' n2 <- sum(y)
#' p_value <- sageTestNew(x, y, lengthx, lengthy, n1, n2, scale)
#' @export


sageTestNew<-function (x, y, lengthx, lengthy, n1, n2, scale)
{
  if (any(is.na(x)) || any(is.na(y)))
      stop("missing values not allowed")
  if (any(x < 0) || any(y < 0))
      stop("x and y must be non-negative")
  if (length(x) != length(y))
      stop("x and y must have same length")
  x <- round(x)
  y <- round(y)
  nn <- round(n1)
  mm <- round(n2)

  size <- x + y
  pValue <- rep(1, length(size))
  rate <- scale*lengthx*nn/lengthy/mm
  prob <- rate/(1+rate)

  if (any(big <- size > 10000)) {
      ibig <- (seq_along(x))[big]
      for (i in ibig)
           pValue[i] <- chisq.test(matrix(c(x[i], y[i],
                                           nn-x[i], mm-y[i]), 2, 2))$p.value
  }

  size0 <- size[size > 0 & !big]
  x1 <- x[size > 0 & !big]
  prob1 <- prob[size > 0 & !big]
  pValue1 <- rep(1, length(size0))
  if (length(size0))
      for(i in seq_len(length(size0))){
          p <- dbinom(0:size0[i], prob=prob1[i], size=size0[i])
          o <- order(p)
          cumsump <- cumsum(p[o])[order(o)]
          pValue1[i] <- cumsump[x1[i] + 1]
    }

  pValue[size > 0 & !big] <- pValue1
  return(pValue)
}
