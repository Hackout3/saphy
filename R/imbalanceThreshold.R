#' Produces a list of a confidence intervals given a tree through tip node topology permutation
#'
#' \code{imbalanceThreshold} returns a list of thresholds
#'
#' @param tree a phylogenetic tree (as a \code{phylo} object)
#'
#' @param reps the number of permutations to execute for the tree
#'
#' @param ci confidence interval. Default = 0.25, i.e. 95 percent
#'
#' @author OJ Watson (\email{o.watson15@@imperial.ac.uk}) Yu Luo (\email{yu.luo3@@mail.mcgill.ca})
#'
#' @return A list of length 12 for th first 12 imbalance metrics
#'
#' @export
imbalanceThreshold <- function(tree,reps=5,ci=0.025) {

  thre<-NULL
  bstrees<-permuteTrees(tree,reps=reps)
  bslist<-lapply(bstrees,imbalanceMetrics)
  treemetrics<-imbalanceMetrics(tree)
  nms <- names(treemetrics)[1:12]
  for (j in 1:12){
    ql<-NULL
    qu<-NULL
    bs<-NULL
    for (k in 1:reps){
      metr<-bslist[[k]]
      bs<-unlist(c(bs,metr[j]))
      }
    ql<-c(ql,quantile(bs,prob=c(ci)));qu<-c(qu,quantile(bs,prob=c(1-ci)))
    thre<-cbind(thre,rbind(qu,ql))

  }
  colnames(thre) <- nms
  return(thre)
  }
