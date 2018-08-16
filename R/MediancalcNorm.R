#' Set the initial value
#' Using median method to compute the normalization factor to identify
#' difference expression (DE) of genes between different species
#' @describeIn To get scaling factor for different species.
#' @usage MediancalcNorm(orth_gene, hkind)
#' @param orth_gene Matrix or data.frame containing read counts and gene
#' length for each orthologous gene between different species. The first and
#' third column containing gene length, the second and the fourth column
#' containing read counts.
#' @param hkind A vector shows conserved genes position in orthologous genes.
#' @return scale Computed Normalization factor.
#' @references Brawand D, Soumillon M, Necsulea A, Julien P, Csardi G, Harrigan
#' P, et al. The evolution of gene expression levels in mammalian organs.
#' Nature. 2011;478:343-348.
#' @examples
#' data(sim_data)
#' MediancalcNorm(orth_gene=sim_data, hkind=1:1000)
#' @export

MediancalcNorm <- function(orth_gene, hkind)
{
  RPKMh <- orth_gene[, 2]/orth_gene[, 1]*(10^9/sum(orth_gene[, 2]))
  RPKMm <- orth_gene[, 4]/orth_gene[, 3]*(10^9/sum(orth_gene[, 4]))
  cRPKMh <- RPKMh[hkind]
  cRPKMm <- RPKMm[hkind]
  scale <- median(cRPKMh)/median(cRPKMm)
  return(scale)
}
