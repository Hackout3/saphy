#' Produces a list of a list of sequentially growing trees given a simulated tree.
#'
#' \code{prunePermute} returns a list of length equal to the number of nodes in the given tree
#'
#' @param tree a phylogenetic tree (as a \code{phylo} object)
#' @param permutes the number of permutations to execute for the tree
#'
#' @author OJ Watson (\email{o.watson15@@imperial.ac.uk})
#'
#' @return A list of length equal to the number of nodes in the given tree
#'
#' @export

prunePermute <- function(tree,permutes=10) {

  nms <- names(timeprune(tree)$tipdates)
  p.tree <- permuteTrees(tree,permutes)
  n.tips <- tree$Nnode + 1

  lapply(p.tree,timeprune)

  prunes <- lapply(lapply(p.tree,timeprune),with,trees)


  same.prunes <- list()
  length(same.prunes) <- n.tips-1

  for (i in 1:(n.tips-1)){
    same.prunes[[i]] <- lapply(prunes, `[[`, i)
  }

  return(list("seqPrunes" = same.prunes,
              "seqTipNames" = nms))

}
