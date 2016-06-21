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
timeprune <- function(tr,tipdates=NULL,reverse=FALSE){
  tl <- tr$tip.label
  if(!is.null(tipdates)){
    td <- tipdates
  }else{
    if(is.null(tipdates)&is.rooted(tr)){
      td <- distRoot(tr)
    }else{
      stop("Tipdates not specified and tree is unrooted.")
    }
  }
  ntips <- length(tr$tip.label)
  o <- order(td,decreasing=TRUE)
  trs <- list()
  for(i in 1:(ntips-2)){
    mintip <- tl[o[i]]
    tr <- drop.tip(tr,mintip,trim.internal=TRUE)
    trs[[i]] <- tr
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
