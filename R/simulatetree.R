simulate <- function(N = c(100, 10000),
                     params = list(beta = matrix(c(2,0,0,2),nrow=2), gamma = 1, mu = .05),
                     I0 = c(1,1),
                     t1 = 100,
                     sampletimes = c("incidence", "regular"),
                     samplesize = 200,
                     sampleperiod = c(0, t1),
                     sampledist = c("incidence", "popsize", "fixed"),
                     samplebias = 1) {

  ### Build the demographic process, i.e. the transmission model
  dp <- build.demographic.process(
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

  ### Make incidence curves rather than prevalence curves
  epicurve <- dp(theta = c(params, list(N = N)),
                x0 = c(I1 = I0[1], I2 = I0[2], S1 = N[1] - I0[1], S2 = N[2] - I0[2]),
                t0 = 0, t1 = t1)[[5]]
  epicurve <- epicurve[epicurve[,1] > sampleperiod[1] & epicurve[,1] < sampleperiod[2] - t1/1000,]
  incfunc <- function(epirow, b) {
    c(sum(epirow[4]*b[,1]*epirow[2:3]), sum(epirow[5]*b[,2]*epirow[2:3]))/N
  }

  inccurve <- apply(epicurve, 1, incfunc, b = params$beta)

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
  sampleTimes <- sampleTimes + runif(samplesize, 0, t1/1000)
  sampleStates <- sampleStates[order(sampleTimes),]
  sampleTimes <- sort(sampleTimes)

  cotree <- sim.co.tree(theta = c(params, list(N = N)),
                        dp,
                        x0 = c(I1 = I0[1], I2 = I0[2], S1 = N[1] - I0[1], S2 = N[2] - I0[2]),
                        t0 = 0, sampleTimes = sampleTimes,
                        sampleStates = sampleStates)

  seqs <- simSeq(cotree)
  return(list(tree = cotree,
              sequences = seqs,
              sampledpopulations = sampleStates))

}


