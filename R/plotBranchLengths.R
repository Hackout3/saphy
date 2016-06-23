#' Plots the average branch lengths over time
#'
#' \code{plot.branch.lengths} returns a list of imabalance statistics
#'
#' @param tree a phylogenetic tree
#' @param permutes number of permutation (passed to prunePermute)
#'
#' @author John Lees (\email{lees.john6@@gmail.com})
#'
#' @export

plot.branch.lengths <- function(tree, permutes=10) {
  pruned_trees <- saphy::prunePermute(tree, permutes)

  lengths <- lapply(pruned_trees$seqPrunes, lapply, saphy::tiplength, pruned_trees$seqTipNames)
  mean_lengths <- lapply(lengths, lapply, mean, na.rm=T)

  x <- seq(1,length(mean_lengths),1)
  y <- rep(NA, length(mean_lengths))
  for (permutation in 1:length(mean_lengths[[1]]))
  {
    for (i in 1:length(mean_lengths))
    {
      y[i] <- mean_lengths[[i]][[permutation]][1]
    }
    plot(x, y, xlab="Node", ylab="Mean tip length", main=paste0("Permutation ", permutation))
    readline(prompt="Press [enter] to view next plot")
  }

}
