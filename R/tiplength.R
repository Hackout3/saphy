#' Extracts the length of a tip from a tree
#' 
#' \code{tiplength} returns the length of a tip of a tree, given either (a) the tip name or
#' (b) the index of the tip (in the order of the tip labels in the tree).
#' 
#' @param tree a phylogenetic tree (as a \code{phylo} object)
#' @param tipname the tip name, as a character string, or a numeric index
#' @return The tip length (as a \code{double}).
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#' @export
tiplength <- function(tree,tipname){
  if(is.character(tipname)){
    idx <- match(tipname,tree$tip.label)
  }else{
    idx <- tipname
  }
  idx2 <- match(idx,tree$edge[,2])
  return(tree$edge.length[idx2])
}