# JAGS model from CCTpack
mccrmid <- function() {
  for (l in 1:nobs) {
    Y[l,3] ~ dnorm((a[Y[l,1]]*T[Y[l,2],Om[Y[l,1]]])+b[Y[l,1]], pow(a[Y[l,1]]*E[Y[l,1]]*exp(lam[Y[l,2],Om[Y[l,1]]]),-2))
  }

  #Parameters
  for (i in 1:nresp){
    Om[i] ~ dcat(pi)
    Elog[i] ~ dnorm(Emu[Om[i]],Etau[Om[i]])
    E[i] <- exp(Elog[i])
    alog[i] ~ dnorm(amu[Om[i]],atau[Om[i]])T(-2.3,2.3)
    a[i] <- exp(alog[i])
    b[i] ~ dnorm(bmu[Om[i]],btau[Om[i]])
  }

  for (k in 1:nitem){
    for (v in 1:V){
      T[k,v] ~ dnorm(Tmu[v],Ttau[v])
      lam[k,v] ~ dnorm(lammu[v],lamtau[v])T(-2.3,2.3)
    }
  }

  pi[1:V] ~ ddirch(L)

  #Hyperparameters
  for (v in 1:V) {
    L[v] <- 1
    Tmu[v] ~ dnorm(0,.25)
    Tsig[v] ~ dunif(.25,3)
    Ttau[v] <- pow(Tsig[v],-2)
    Emu[v] ~ dnorm(0,.01)
    Etau[v] ~ dgamma(.01, .01)
    amu[v] <- 0
    atau[v] ~ dgamma(.01, .01)T(.01,)
    bmu[v] <- 0
    btau[v] ~ dgamma(.01, .01)
    lammu[v] <- 0
    lamsig[v] ~ dunif(.25, 2)
    lamtau[v] <- pow(lamsig[v],-2)
  }
}

# assumes V is 1
mccrmid_one <- function() {
  for (l in 1:nobs) {
    Y[l,3] ~ dnorm(
      (a[Y[l,1]] * T[Y[l,2], Om[Y[l,1]]]) + b[Y[l,1]],
      pow(a[Y[l,1]] * E[Y[l,1]] * exp(lam[Y[l,2], Om[Y[l,1]]]), -2)
    )
  }

  #Parameters
  for (i in 1:nresp){
    Om[i] ~ dcat(pi)
    Elog[i] ~ dnorm(Emu[Om[i]],Etau[Om[i]])
    E[i] <- exp(Elog[i])
    alog[i] ~ dnorm(amu[Om[i]],atau[Om[i]])T(-2.3,2.3)
    a[i] <- exp(alog[i])
    b[i] ~ dnorm(bmu[Om[i]],btau[Om[i]])
  }

  for (k in 1:nitem) {
    T[k] ~ dnorm(Tmu, Ttau)
    lam[k] ~ dnorm(lammu, lamtau)T(-2.3,2.3)
  }

  L <- 1
  Tmu ~ dnorm(0,.25)
  Tsig ~ dunif(.25,3)
  Ttau <- pow(Tsig[v],-2)
  Emu ~ dnorm(0,.01)
  Etau ~ dgamma(.01, .01)
  amu <- 0
  atau ~ dgamma(.01, .01)T(.01, )
  bmu <- 0
  btau ~ dgamma(.01, .01)
  lammu <- 0
  lamsig ~ dunif(.25, 2)
  lamtau <- pow(lamsig[v], -2)

}
