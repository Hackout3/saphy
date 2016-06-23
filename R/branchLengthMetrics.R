#' Produces a dataframe of tree branch length statistics
#'
#' \code{branchLengthMetrics} returns a list of imabalance statistics
#'
#' @param tree a phylogenetic tree (as a \code{phylo} object)
#'
#' @author John Lees (\email{lees.john6@@gmail.com})
#'
#' @return A list with 5 elements relating to tip lengths
#'
#' @export

branchLengthMetrics <- function(tree) {

  ls <- list()
  ls$meanTipLength <- mean(saphy::tiplength(tree, tree$tip.label))
  ls$meanInternalLength <- mean(tree$edge.length[tree$edge[-match(tree$tip.label, tree$tip.label),2]], na.rm = T)
  ls$tipRatio <- ls$meanTipLength/ls$meanInternalLength

  # Clustering - expect fewer/less well defined for transmission/Ne growth?
  # An alternative would be to use cmdscale + DBSCAN:
  # tree_mds <- cmdscale(cophenetic.phylo(tree))
  # kNNdistplot(tree_mds)
  # plot(tree_mds, col=factor(dbscan(tree_mds, 15)$cluster))
  # or hclust, but I didn't find much difference between transmission and
  # no transmission for a range of agglomeration methods
  # clusts <- hclust(dist(cophenetic.phylo(sim_tree$tree)), method="ward.D2")

  # Dists with long branches are normally distributed
  # with clusters there is a peak below the mean which looks like a cutoff
  distances <- cophenetic.phylo(tree)[upper.tri(cophenetic.phylo(tree))]
  ls$clustsShapiro <- shapiro.test(dists)$statistic

  # Larger for positive skew -> close tips; 1 for symmetric -> equally distant tips
  quantiles <- quantile(distances, probs = c(0.25,0.5,0.75))
  ls$clustsQuantileRatio <- as.numeric((quantiles[2]-quantiles[1])/(quantiles[3]-quantiles[2]))

  return(ls)

}
