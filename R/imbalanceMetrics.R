#' Produces a dataframe of tree imbalance statistics
#'
#' \code{imbalanceMetrics} returns a list of imabalance statistics
#'
#' @param tree a phylogenetic tree (as a \code{phylo} object)
#'
#' @author OJ Watson (\email{o.watson15@@imperial.ac.uk})
#'
#' @return A list of 14 imbalance metrics. TO DO. Improve local branching index utility
#'
#' @export

imbalanceMetrics <- function(tree) {

    ls <- list()
    xx <- as.treeshape(tree)
    ntips <- tree$Nnode + 1
    ntot <- n + tree$Nnode
    che <- treeImbalance::Ncherries(tree)

    ls$Colless.yule  <- apTreeshape::colless(xx,norm="yule")
    ls$Colless.pda  <- apTreeshape::colless(xx,"pda")
    ls$Sackin.yule  <- apTreeshape::sackin(xx,"yule")
    ls$Sackin.pds  <- apTreeshape::sackin(xx,"pda")
    ls$McKenzie <- abs(che - ntips/3)/sqrt(2 * ntips/45)

    ls$I1 <- treeImbalance::I1(tree)
    ls$I2 <- treeImbalance::I2(tree)
    ls$Ic <- treeImbalance::Ic(tree)
    ls$meanIprime <- treeImbalance::meanIprime(tree)
    ls$M <- treeImbalance::M(tree)
    ls$B1 <- treeImbalance::B1(tree)
    ls$B2 <- treeImbalance::B2(tree)
    ls$lbi.tips <- treeImbalance::lbi(tree)[1:ntips]
    ls$lbi.nodes <- treeImbalance::lbi(tree)[(n+1):ntot]

    return(ls)

}
