#' Compute the normalization factor to identify difference expression
#' (DE) of genes between different species
#' @description To normalize read counts and identify difference expression(DE)
#' of orthologous genes between different species.
#' @usage SCBN(orth_gene, hkind, a=0.05)
#' @param orth_gene Matrix or data.frame containing read counts and gene length
#' for each orthologous gene between different species. The first and third
#' column containing gene length, the second and the fourth column containing
#' read counts.
#' @param hkind A vector shows conserved genes position in orthologous genes.
#' @param a P-value cutoff in iteration process to find the optimal
#' normalization factor.
#' @return list(.) A list of computed normalization factors, "median_val"
#' represents foctors computed by median methods, "scbn_val" represents factors
#' computed by SCBN methods.
#' @examples
#' data(sim_data)
#' SCBN(orth_gene=sim_data, hkind=1:1000, a=0.05)
#' @export

SCBN <- function(orth_gene, hkind, a=0.05)
{
  if (all(!is.na(orth_gene)) == FALSE) {
      stop("The dataset of orthologous genes has NA values.", call.=TRUE)
  } else if (all(hkind %in% (seq_len(nrow(orth_gene)))) == FALSE){
             stop("The conserved genes are not included in dataset of
                   orthologous genes.", call.=TRUE)
  }

  scale <- MediancalcNorm(orth_gene, hkind)
  scbn_val <- Iter_optimal(scale=scale, orth_gene=orth_gene, hkind=hkind, a=a)
  list(median_val=scale, scbn_val=scbn_val)
}
