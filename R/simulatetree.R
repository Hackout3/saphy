simulate <- function(N = c(100, 10000),
                     params = list(beta = matrix(c(2,0,0,2),nrow=2), gamma = 1, mu = .05),
                     t1 = 100,
                     sampletimes = c("incidence", "regular"),
                     samplesize = 200,
                     sampleperiod = c(0, t1),
                     sampledist = c("incidence", "popsize", "fixed"),
                     samplebias = ) {

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


  ###
  show.demographic.process(dp, theta = list(beta = params$beta,
                                                 N = N, gamma = params$gamma, mu = params$mu),
                           x0 = c(S1 = N[1] - 1, I1 = 1, S2 = N[2] - 1, I2 = 1),
                           t0 = 0, t1 = 100)

}

#TO PRODUCE
# annotated tree
#



