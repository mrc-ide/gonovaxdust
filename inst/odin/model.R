# Key to indices
# i = Activity group
## 1: Low
## 2: High
# j = Vaccination status
## 1: Unvaccinated
## 2: Bexsero
## 3: Waned

n_group <- 2
n_vax   <- user(1)

## calibrate time-varying parameters
# tt runs from t0 = 2009, to t10 = 2019
dt <- 1 / 365
initial(time) <- 0
update(time) <- (step + 1) * dt

## What we really want is min(step + 1, length(beta_step)) but that's not
## supported by odin (it could be made to support this).
update(beta) <- if (as.integer(step) >= length(beta_step))
  beta_step[length(beta_step)] else beta_step[step + 1]
update(eta[1]) <- if (as.integer(step) >= length(eta_l_step))
  eta_l_step[length(eta_l_step)] else eta_l_step[step + 1]
update(eta[2]) <- if (as.integer(step) >= length(eta_h_step))
  eta_h_step[length(eta_h_step)] else eta_h_step[step + 1]

## Core equations for transitions between compartments:

update(U[, ]) <- U[i, j] + n_xU[i, j] - n_UI[i, j] - n_Ux[i, j] + n_AU[i, j] +
  n_TU[i, j] + sum(wU[i, j, ]) -
  sum(n_vbe[i, j, ]) - sum(n_vod[i, j, ]) - sum(n_vos[i, j, ])

update(I[, ]) <- I[i, j] + n_UI[i, j] - n_I[i, j] + sum(wI[i, j, ])

update(A[, ]) <- A[i, j] + n_IA[i, j] - n_A[i, j] + sum(wA[i, j, ])

update(S[, ]) <- S[i, j] + n_IS[i, j] - n_S[i, j] + sum(wS[i, j, ])

update(T[, ]) <- T[i, j] + n_ST[i, j] + n_AT[i, j] - n_T[i, j] + sum(wT[i, j, ])

## Update population size
update(N[, ]) <- U[i, j] + I[i, j] + A[i, j] + S[i, j] + T[i, j]

n_xU[, 1] <- rpois(enr * q[i] * dt)

# calculate mixing matrix, probability of infection and force of infection
C[, ] <- (1 - vei[j]) * (I[i, j] + A[i, j] + S[i, j])
prop_C[] <- sum(C[i, ]) / sum(N[i, ])
Np[]    <- sum(N[i, ]) * p[i]

foi_LH[] <- prop_C[i] * Np[i] / sum(Np[])
update(lambda[]) <- p[i] * beta * (epsilon * prop_C[i] + (1 - epsilon) * sum(foi_LH[]))

## Individual probabilities of transition:
r_UI[, ] <- lambda[i] * (1 - vea[j])
r_UU[, ] <- eta[i]
r_AT[, ] <- eta[i]
r_AU[, ] <- nu / (1 - ved[j])
r_ST[, ] <- mu
r_TU[, ] <- rho

r_U[, ] <- r_UI[i, j] + r_UU[i, j] + exr
r_I[, ] <- sigma + exr
r_A[, ] <- r_AT[i, j] + r_AU[i, j] + exr
r_S[, ] <- r_ST[i, j] + exr
r_T[, ] <- r_TU[i, j] + exr

## Draws from binomial distributions for numbers leaving each compartments
n_U[, ] <- rbinom(U[i, j], 1 - exp(-r_U[i, j] * dt))
n_I[, ] <- rbinom(I[i, j], 1 - exp(-r_I[i, j] * dt))
n_A[, ] <- rbinom(A[i, j], 1 - exp(-r_A[i, j] * dt))
n_S[, ] <- rbinom(S[i, j], 1 - exp(-r_S[i, j] * dt))
n_T[, ] <- rbinom(T[i, j], 1 - exp(-r_T[i, j] * dt))

# Draw the number of leavers from each compartment
n_Ux[, ] <- rbinom(n_U[i, j], exr / r_U[i, j])
n_Ix[, ] <- rbinom(n_I[i, j], exr / r_I[i, j])
n_Ax[, ] <- rbinom(n_A[i, j], exr / r_A[i, j])
n_Sx[, ] <- rbinom(n_S[i, j], exr / r_S[i, j])
n_Tx[, ] <- rbinom(n_T[i, j], exr / r_T[i, j])

# Draw the numbers of transitions between compartments
n_UI[, ] <- rbinom(n_U[i, j] - n_Ux[i, j],
                   r_UI[i, j] / (r_UI[i, j] + r_UU[i, j]))
n_UU[, ] <- n_U[i, j] - n_Ux[i, j] - n_UI[i, j]
n_IS[, ] <- rbinom(n_I[i, j] - n_Ix[i, j], (1 - ves[j]) * psi)
n_IA[, ] <- n_I[i, j] - n_Ix[i, j] - n_IS[i, j]
n_AT[, ] <- rbinom(n_A[i, j] - n_Ax[i, j],
                   r_AT[i, j] / (r_AT[i, j] + r_AU[i, j]))
n_AU[, ] <- n_A[i, j] - n_Ax[i, j] - n_AT[i, j]
n_ST[, ] <- n_S[i, j] - n_Sx[i, j]
n_TU[, ] <- n_T[i, j] - n_Tx[i, j]

# Vaccination
## time-varying switch
vax_switch <- if (as.integer(step) >= length(vax_step))
  vax_step[length(vax_step)] else vax_step[step + 1]
## at screening
n_vos[, , ] <- rbinom(n_UU[i, k], vos[i, j, k] * vax_switch)
## on diagnosis
n_vod[, , ] <- rbinom(n_TU[i, k], vod[i, j, k] * vax_switch)
## on entry - no switch as background rate
n_vbe[, , ] <- rbinom(n_xU[i, k], vbe[i, j, k])

# Waning (inter-stratum transition) occurs after inter-compartment transition
wU[, , ] <- rbinom(U[i, k] - n_UI[i, j] - n_Ux[i, j], w[j, k])
wI[, , ] <- rbinom(I[i, k] - n_I[i, k], w[j, k])
wA[, , ] <- rbinom(A[i, k] - n_A[i, k], w[j, k])
wS[, , ] <- rbinom(S[i, k] - n_S[i, k], w[j, k])
wT[, , ] <- rbinom(T[i, k] - n_T[i, k], w[j, k])

## outputs
update(cum_incid[, ])      <- n_UI[i, j]
update(cum_diag_a[, ])     <- n_AT[i, j]
update(cum_diag_s[, ])     <- n_ST[i, j]
update(cum_treated[, ])    <- n_TU[i, j]
update(cum_screened[, ])   <- n_UU[i, j]
update(cum_vaccinated[, ]) <- n_vos[i, j, j] + n_vod[i, j, j] + n_vbe[i, j, j]

# aggregated time series for fitting mcmc
update(tot_treated) <- sum(cum_treated)
update(tot_attended) <- sum(cum_treated) + sum(cum_screened)
update(entrants) <- sum(n_xU)
update(leavers) <- sum(n_Ux) + sum(n_Ix) + sum(n_Ax) + sum(n_Sx) + sum(n_Tx)


## Set up compartments
## Initial states are all 0 as we will provide a state vbector
initial(U[, ]) <- U0[i, j]
initial(I[, ]) <- I0[i, j]
initial(A[, ]) <- A0[i, j]
initial(S[, ]) <- S0[i, j]
initial(T[, ]) <- T0[i, j]

U0[, ] <- user()
I0[, ] <- user()
A0[, ] <- user()
S0[, ] <- user()
T0[, ] <- user()

initial(cum_incid[, ])      <- 0
initial(cum_diag_a[, ])     <- 0
initial(cum_diag_s[, ])     <- 0
initial(cum_treated[, ])    <- 0
initial(cum_screened[, ])   <- 0
initial(cum_vaccinated[, ]) <- 0
initial(tot_treated) <- 0
initial(tot_attended) <- 0
initial(entrants) <- 0
initial(leavers) <- 0
initial(beta) <- 0
initial(eta[]) <- 0
initial(N[, ]) <- U0[i, j] + I0[i, j] + A0[i, j] + S0[i, j] + T0[i, j]
initial(lambda[]) <- 0

# set up dimensions of compartments
dim(U) <- c(n_group, n_vax)
dim(I) <- c(n_group, n_vax)
dim(A) <- c(n_group, n_vax)
dim(S) <- c(n_group, n_vax)
dim(T) <- c(n_group, n_vax)

dim(U0) <- c(n_group, n_vax)
dim(I0) <- c(n_group, n_vax)
dim(A0) <- c(n_group, n_vax)
dim(S0) <- c(n_group, n_vax)
dim(T0) <- c(n_group, n_vax)

dim(C)    <- c(n_group, n_vax)
dim(N)    <- c(n_group, n_vax)
dim(n_xU) <- c(n_group, n_vax)
dim(Np)     <- n_group
dim(prop_C) <- n_group
dim(foi_LH) <- n_group
dim(lambda) <- n_group


dim(r_UI) <- c(n_group, n_vax)
dim(r_AT) <- c(n_group, n_vax)
dim(r_AU) <- c(n_group, n_vax)
dim(r_ST) <- c(n_group, n_vax)
dim(r_TU) <- c(n_group, n_vax)
dim(r_UU) <- c(n_group, n_vax)

dim(r_U) <- c(n_group, n_vax)
dim(r_I) <- c(n_group, n_vax)
dim(r_A) <- c(n_group, n_vax)
dim(r_S) <- c(n_group, n_vax)
dim(r_T) <- c(n_group, n_vax)

dim(n_UI) <- c(n_group, n_vax)
dim(n_IA) <- c(n_group, n_vax)
dim(n_IS) <- c(n_group, n_vax)
dim(n_AT) <- c(n_group, n_vax)
dim(n_AU) <- c(n_group, n_vax)
dim(n_ST) <- c(n_group, n_vax)
dim(n_TU) <- c(n_group, n_vax)
dim(n_UU) <- c(n_group, n_vax)

dim(n_Ux) <- c(n_group, n_vax)
dim(n_Ix) <- c(n_group, n_vax)
dim(n_Ax) <- c(n_group, n_vax)
dim(n_Sx) <- c(n_group, n_vax)
dim(n_Tx) <- c(n_group, n_vax)

dim(n_U) <- c(n_group, n_vax)
dim(n_I) <- c(n_group, n_vax)
dim(n_A) <- c(n_group, n_vax)
dim(n_S) <- c(n_group, n_vax)
dim(n_T) <- c(n_group, n_vax)

dim(cum_incid)      <- c(n_group, n_vax)
dim(cum_diag_a)     <- c(n_group, n_vax)
dim(cum_diag_s)     <- c(n_group, n_vax)
dim(cum_treated)    <- c(n_group, n_vax)
dim(cum_screened)   <- c(n_group, n_vax)
dim(cum_vaccinated) <- c(n_group, n_vax)

## Parameters
p[] <- user()
q[] <- user()

enr <- user()
exr <- user()
beta_step[]  <- user()
eta_l_step[] <- user()
eta_h_step[] <- user()
epsilon <- user()
sigma   <- user()
psi     <- user()
nu      <- user()
mu      <- user()
rho     <- user()

## vaccination parameters
# vaccination routes
vbe[, , ] <- user()
vos[, , ] <- user()
vod[, , ] <- user()

# vaccine effects
vea[] <- user() # efficacy against acquisition
ved[] <- user() # efficacy against duration of infection
ves[] <- user() # efficacy against symptoms
vei[] <- user() # efficacy against infectiousness

w[, ]    <- user()
vax_step[]  <- user()

## par dimensions
dim(beta_step)  <- user()
dim(eta_l_step) <- user()
dim(eta_h_step) <- user()

dim(p)    <- n_group
dim(q)    <- n_group
dim(eta)  <- n_group
dim(vea)  <- n_vax
dim(ved)  <- n_vax
dim(ves)  <- n_vax
dim(vei)  <- n_vax
dim(vbe)   <- c(n_group, n_vax, n_vax)
dim(vod)   <- c(n_group, n_vax, n_vax)
dim(vos)   <- c(n_group, n_vax, n_vax)
dim(w)    <- c(n_vax, n_vax)
dim(vax_step) <- user()

dim(n_vbe) <- c(n_group, n_vax, n_vax)
dim(n_vos) <- c(n_group, n_vax, n_vax)
dim(n_vod) <- c(n_group, n_vax, n_vax)
dim(wU)   <- c(n_group, n_vax, n_vax)
dim(wI)   <- c(n_group, n_vax, n_vax)
dim(wA)   <- c(n_group, n_vax, n_vax)
dim(wS)   <- c(n_group, n_vax, n_vax)
dim(wT)   <- c(n_group, n_vax, n_vax)
