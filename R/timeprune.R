#' Prune tree in time order
#'
#' \code{timeprune} prunes a tree by removing tips one at a time, in reverse time order.
#' If the time is omitted and the tree is rooted, then the root-to-tip distance is used as
#' a proxy of sampling time.
#'
#' @param tr a tree, as a \code{phylo} object
#' @param tipdates a vector of sampling times or dates
#' @param reverse a boolean (default: \code{FALSE}) determining whether the list of trees should be returned
#' in reverse time
#'
#' @return A list with tipdates and trees.
#'
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#'
#' @export
timeprune <- function(tr,tipdates=NULL,reverse=FALSE){
  tl <- tr$tip.label
  if(!is.null(tipdates)){
    td <- tipdates
  }else{
    if(is.null(tipdates)&ape::is.rooted(tr)){
      td <- adephylo::distRoot(tr)
    }else{
      stop("Tipdates not specified and tree is unrooted.")
    }
  }
  ntips <- length(tr$tip.label)
  o <- order(td,decreasing=TRUE)
  tl <- tl[o]
  trs <- list()
  for(i in 1:(ntips-1)){
    trs[[tl[i]]] <- tr
    if(i<(ntips-1)){
      mintip <- tl[i]
      tr <- ape::drop.tip(tr,mintip,trim.internal=TRUE)
    }
  }
  if(!reverse){
    trs <- rev(trs)
    td <- rev(td[o])
  }
  else{
    td <- td[o]
  }
  return(list(tipdates=td,trees=trs))
}
