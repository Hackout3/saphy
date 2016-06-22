#' Produces a list of permuted trees with the same tip and internal node times for bootstrapping purposes
#'
#' \code{permuteTrees} returns the length of a tip of a tree, given either (a) the tip name or
#' (b) the index of the tip (in the order of the tip labels in the tree).
#'
#' @param tree a phylogenetic tree (as a \code{phylo} object)
#' @param reps the number of permuted trees
#'
#' @author OJ Watson (\email{o.watson15@@imperial.ac.uk})
#'
#' @return A list of permuted trees
#'
#' @export

permuteTrees <- function(tree, reps = 100) {

  bs.trees <-  rep(tree,reps)

  bs.trees <- lapply(bs.trees,FUN = treeImbalance::getSimTree)

  return(bs.trees)
}
