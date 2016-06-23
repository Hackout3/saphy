#' Simulate tree and sequences in two-population SIR-model
#'
#' \code{simulatetree} does a deterministic simulation of a two-population SIR-model,
#' takes samples and simulates a phylogenetic tree. For all this, it makes use of phydynR
#' (https://github.com/emvolz-phylodynamics/phydynR). Then it simulates DNA-sequences for these samples,
#' with simSeq (\pkg{phangorn}).
#'
#' @param N the population sizes of the two populations
#' @param beta an array of 2x2 transmission rate matrices. The third dimension is the number of matrices, which is the
#' number of time intervals with different transmission rates. Transmission is modelled with frequency-dependent transmission
#' terms (division by N).
#' @param gamma the recovery rate
#' @param mu the death rate. Births balance deaths.
#' @param I0 initial numbers of infecteds. The rest is assumed susceptible.
#' @param t1 a vector with times at which the transmission rates change to the next matrix in \code{beta}.
#' The last time is the duration of the simulation
#' @param res resolution of the epidemic simulation per transmission rate interval (ODE solver)
#' @param sampletimes indicates if sampling occurs at random times proportional to incidence (default \code{"incidence"}) or
#' at random uniform times (\code{"regular"})
#' @param samplesize the number of samples to be taken
#' @param sampleperiod the time interval from which the sampling times should be taken
#' @param sampledist indicates how samples should be divided over the two populations, apart from bias. The default \code{"incidence"} indicates
#' proportionality to incidence; \code{"popsize"} indicates proportionality to population size; \code{"fixed"} indicates
#' no proportionality, i.e. fifty-fifty
#' @param samplebias indicates the extra weight (multiplication factor) given for sampling the first population
#' relative to the second. If equal to 1, the weight is completely determined by \code{sampledist}.
#' @param ... arguments to be passed to simSeq (\pkg{phangorn})
#'
#' @return A list with a tree (class \code{phylo}), and a set of sequences (class \code{phyDat}), and the
#' time course of the prevalence (and susceptible populations).
#'
#' @author Don Klinkenberg (\email{don@@xs4all.nl})
#'
#' @export
simulatetree <- function(N = c(100, 10000),
                          beta = array(c(1.5,0,0,2, 1.5,0,0.1,2), dim = c(2,2,2)), gamma = 1, mu = .05,
                          I0 = c(1, 1),
                          t1 = c(50, 100), res = 1000,
                          sampletimes = c("incidence", "regular"),
                          samplesize = 200,
                          sampleperiod = c(0, tail(t1, 1)),
                          sampledist = c("incidence", "popsize", "fixed"),
                          samplebias = 1, ...) {
  if(length(dim(beta)) != 3) {
    dim(beta) <- c(2, 2, length(beta)/4)
  }
  if(length(t1) != dim(beta)[3]) {
    stop("number of transmission matrices in beta should match length of t1")
  }

  ### Build the demographic process for constant transmission rates, i.e. the transmission model
  dp.piece <- phydynR::build.demographic.process(
    births = matrix(c('parms$beta[1,1] * S1 * I1 / parms$N[1]', 'parms$beta[1,2] * S1 * I2 / parms$N[1]',
                      'parms$beta[2,1] * S2 * I1 / parms$N[2]', 'parms$beta[2,2] * S2 * I2 / parms$N[2]'), nrow = 2,
                    dimnames = list(c("I1","I2"),c("I1","I2"))),
    nonDemeDynamics = c(S1 = 'parms$mu * parms$N[1] - parms$beta[1,1] * S1 * I1 / parms$N[1] - parms$beta[1,2] * S1 * I2 / parms$N[1] - parms$mu * S1',
                        S2 = 'parms$mu * parms$N[2] - parms$beta[2,1] * S2 * I1 / parms$N[2] - parms$beta[2,2] * S2 * I2 / parms$N[2] - parms$mu * S2'),
    deaths = c(I1 = 'parms$gamma * I1 + parms$mu * I1',
               I2 = 'parms$gamma * I2 + parms$mu * I2'),
    parameterNames = c('beta', 'N', 'gamma', 'mu'),
    rcpp = FALSE, sde = FALSE
  )

  ### Make epi-curve, incidence curve, and initial states for each piecewise constant-rate transmission model
  epicurve <- dp.piece(theta = list(beta = beta[,,1], gamma = gamma, mu = mu, N = N),
                       x0 = c(I1 = I0[1], I2 = I0[2], S1 = N[1] - I0[1], S2 = N[2] - I0[2]),
                       t0 = 0, t1 = t1[1], res = res)[[5]]
  initialstates <- c()
  if(length(t1) > 1) {
    for(i in 2:length(t1)) {
      initialstates <- rbind(initialstates, tail(epicurve, 1))
      epicurve <- head(epicurve,-1)
      epicurve <- rbind(epicurve,
                        dp.piece(theta = list(beta = beta[,,i], gamma = gamma, mu = mu, N = N),
                                 x0 = c(I1 = initialstates[i-1,2], I2 = initialstates[i-1,3],
                                        S1 = initialstates[i-1,4], S2 = initialstates[i-1,5]),
                                 t0 = t1[i-1], t1 = t1[i], res = res)[[5]])
    }
  }
  epicurveoutput <- epicurve
  class(epicurveoutput) <- c("deSolve", "matrix")
  epicurve <- head(epicurve[epicurve[,1] >= sampleperiod[1] & epicurve[,1] <= sampleperiod[2],], -1)
  incfunc <- function(epirow, b_arr, t) {
    b <- b_arr[,,sum(epirow[1]>=t)]
    c(sum(epirow[4]*b[,1]*epirow[2:3]), sum(epirow[5]*b[,2]*epirow[2:3]))/N
  }
  inccurve <- apply(epicurve, 1, incfunc, b_arr = beta, t = c(0, t1))

  ### Determine sample times at times of curve
  if(sampletimes[1] == "incidence") {
    sampleTimes <- sample(epicurve[,1], samplesize, replace = TRUE, prob = colSums(inccurve))
  } else {
    sampleTimes <- sample(epicurve[,1], samplesize, replace = TRUE)
  }

  ### Determine sample population
  if(sampledist[1] == "incidence") {
    sampleprobs <- inccurve[, match(sampleTimes,epicurve[,1]) ] * c(samplebias, 1)
  } else if (sampledist[1] == "popsize") {
    sampleprobs <- matrix(rep(N,samplesize),nrow = 2) * c(samplebias, 1)
  } else {
    sampleprobs <- matrix(.5, nrow = 2, ncol = samplesize) * c(samplebias, 1)
  }
  sampledpops <- apply(sampleprobs, 2, sample, x = 2, size = 1, replace = FALSE)
  sampleStates <- matrix(0, nrow = samplesize, ncol = 2, dimnames = list(NULL, c("I1", "I2")))
  sampleStates[cbind(1:samplesize,sampledpops)] <- 1

  ### Add uniform increment to sample times, and sort times and populations
  add_incr <- function(sampleTime) {
    runif(1, sampleTime, epicurve[epicurve[,1] > sampleTime, 1][1])
  }
  sampleTimes <- sapply(sampleTimes, add_incr)
  sampleStates <- sampleStates[order(sampleTimes),]
  sampleTimes <- sort(sampleTimes)

  ### Build the complete demographic process, i.e. the transmission model
  if(dim(beta)[3] == 1) {
    beta <- beta[,,1]
    dp <- dp.piece
  } else {
    dp <- function(theta, x0, t0, t1, res = 1000, integrationMethod = "lsoda") {
      # time interval by time interval
      dpsections <- list()

      # first time interval start with original initial conditions
      dpsections[[1]] <- dp.piece(theta = list(beta = theta$beta[,,1], gamma = theta$gamma, mu = theta$mu,
                                               N = theta$N),
                                  x0 = x0, t0 = 0, t1 = theta$t[1], res = res, integrationMethod = integrationMethod)

      # other time intervals start with earlier obtained 'initial conditions'
      for(i in 2:length(theta$t)) {
        dpsections[[i]] <- dp.piece(theta = list(beta = theta$beta[,,i], gamma = theta$gamma, mu = theta$mu,
                                                 N = theta$N),
                                    x0 = c(I1 = initialstates[i-1,2], I2 = initialstates[i-1,3],
                                           S1 = initialstates[i-1,4], S2 = initialstates[i-1,5]),
                                    t0 = theta$t[i-1], t1 = theta$t[i],
                                    res = res, integrationMethod = integrationMethod)
      }

      # build a single object from all sections
      result <- list(
        times = tail(dpsections,1)[[1]]$times,
        births = tail(dpsections,1)[[1]]$births,
        migrations = tail(dpsections,1)[[1]]$migrations,
        sizes = tail(dpsections,1)[[1]]$sizes,
        dpsections[[1]][[5]])
      for(i in 2:length(dpsections)) {
        result$times <- c(result$times, tail(dpsections,i)[[1]]$times[-1])
        result$births <- c(result$births, tail(dpsections,i)[[1]]$births[-1])
        result$migrations <- c(result$migrations, tail(dpsections,i)[[1]]$migrations[-1])
        result$sizes <- c(result$sizes, tail(dpsections,i)[[1]]$sizes[-1])
        result[[5]] <- rbind(result[[5]], dpsections[[i]][[5]][-1,])
      }

      return(result)

    }

  }


  ### Simulate the phylogenetic tree
  cotree <- phydynR::sim.co.tree(theta = list(beta = beta, gamma = gamma, mu = mu, N = N,
                                              initialstates = initialstates, t = t1),
                                 dp,
                                 x0 = c(I1 = I0[1], I2 = I0[2], S1 = N[1] - I0[1], S2 = N[2] - I0[2]),
                                 t0 = 0, res = res, sampleTimes = sampleTimes,
                                 sampleStates = sampleStates)
  class(cotree) <- "phylo"

  ### Simulate the sequences
  seqs <- phangorn::simSeq(cotree, ...)
  return(list(tree = cotree,
              sequences = seqs,
              epicurve = epicurveoutput))

}


