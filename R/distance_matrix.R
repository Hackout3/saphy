#' Computes a pairwise distance matrix
#'
#' \code{compute.distance.matrix} calculates a pairwise distance matrix based on a multiple sequence
#' alignment. It is basically a thin wrapper to \code{dist.dna} in \code{ape}, with the default model
#' being TN93 rather than K80 and pairwise deletions considered, which returns a matrix
#' (rather than either a \code{dist} object or a matrix).
#'
#' @param aln a multiple sequence alignment
#' @param model a model of nucleotide substitution (default=TN93)
#' @param gamma the gamma parameter of rate variation
#'
#' @return A distance matrix (as a standard R matrix)
#'
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#'
#' @export
compute.distance.matrix <- function(aln, model = "TN93", gamma = FALSE){
  m <- dist.dna(aln,model=model,gamma=gamma,pairwise.deletion = TRUE,base.freq=NULL,as.matrix=TRUE)
  return(m)
}

#' Calculate pairwise distance between two sequences
#'
#' \code{pairwise.distance} calculates the pairwise distance between two sequences.
#'
#' @param s1 a nucleotide sequence
#' @param s2 a second nucleotide sequence. The sequences should be of equal length
#' @param model  model of nucleotide substitution (default=TN93)
#' @param gamma The gamma parameter of rate variation
#'
#' @return The distance (as a double).
#'
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#'
#' @export
pairwise.distance <- function(s1, s2, model = "TN93", gamma = FALSE){
  aln <- rbind(s1,s2)
  d <- dist.dna(aln,model=model,gamma=gamma,pairwise.deletion=TRUE,base.freq=NULL)
  return(as.double(d))
}

#' Updates a pre-existing distance matrix
#'
#' \code{expand.distance.matrix} updates a pairwise distance matrix with distances calculated between
#' the existing alignment and an additional sequence (assumed to be aligned)
#'
#' @param aln a multiple sequence alignment
#' @param s a nucleotide sequence
#' @param distmat a distance matrix (as a standard matrix). If \code{NULL}, then a distance matrix is computed.
#' @param model a model of nucleotide substitution (default=TN93)
#' @param gamma the gamma parameter of rate variation
#'
#' @return A list with two elements: distmat, a distance matrix (as a standard R matrix), and aln, an updated
#' alignment.
#'
#' @author Simon Frost (\email{sdwfrost@@gmail.com})
#'
#' @export
expand.distance.matrix <- function(aln, s, distmat = NULL, model = "TN93", gamma = FALSE){
  if(is.null(distmat)){
    distmat <- dist.dna(aln,model=model,gamma=gamma,pairwise.deletion=TRUE,base.freq=NULL,as.matrix=TRUE)
  }
  naln <- dim(aln)[1]
  rw <- rep(0.0,naln+1)
  for(i in 1:naln){
    rw[i] <- pairwise.distance(aln[i,],s,model=model,gamma=gamma)
  }
  distmat <- rbind(distmat,rw[1:naln])
  distmat <- cbind(distmat,rw)
  dimnames(distmat)[[1]][naln+1] <-row.names(s)
  dimnames(distmat)[[2]][naln+1] <-row.names(s)
  return(list(aln=rbind(aln,s),distmat=distmat))
}
