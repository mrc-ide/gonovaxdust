context("model (check)")


test_that("there are no infections when beta is 0", {
  pars <- model_params(gono_params = gono_params(), fix_pop = TRUE)
  pars$beta_step[] <- 0

  mod <- model$new(pars, step = 0, n_particles = 10, seed = 1L)
  y <- mod$run(step_end = 365)
  y <- mod$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_incid == 0))
  expect_true(all(unlist(y) >= 0))
  expect_true(all(apply(y$N, 3, sum) == pars$N0))

})

test_that("there are no symptomatic infections when psi = 0", {
  pars <- model_params(gono_params = gono_params(), fix_pop = TRUE)
  pars$psi <- 0

  mod <- model$new(pars, step = 0, n_particles = 10, seed = 1L)
  y <- mod$run(step_end = 365)
  y <- mod$transform_variables(y)

  expect_true(any(y$I != 0))
  expect_true(all(y$S == 0))
  expect_true(all(y$cum_diag_s == 0))
  expect_true(all(y$cum_diag_a[-1, , ] > 0))
  expect_true(all(unlist(y) >= 0))
  expect_true(all(apply(y$N, 3, sum) == pars$N0))
})

test_that("there are no asymptomatic infections when psi = 1", {
  pars <- model_params(gono_params = gono_params(), fix_pop = TRUE)
  pars$psi <- 1
  pars$S0[, ] <- pars$A0[, ]
  pars$A0[, ] <- 0

  mod <- model$new(pars, step = 0, n_particles = 10, seed = 1L)
  y <- mod$run(step_end = 10)
  y <- mod$transform_variables(y)
  expect_true(any(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$cum_diag_a == 0))
  expect_true(all(y$cum_diag_s[-1, , ] > 0))
  expect_true(all(unlist(y) >= 0))
  expect_true(all(apply(y$N, 3, sum) == pars$N0))
})

test_that("there are no infections when A0 = 0", {
  pars <- model_params(gono_params = gono_params(), fix_pop = TRUE)
  pars$A0[, ] <- 0

  mod <- model$new(pars, step = 0, n_particles = 10, seed = 1L)
  y <- mod$run(step_end = 10)
  y <- mod$transform_variables(y)

  expect_true(all(y$I == 0))
  expect_true(all(y$A == 0))
  expect_true(all(y$S == 0))
  expect_true(all(unlist(y) >= 0))
})

test_that("no-one is treated when mu and eta = 0", {
  pars <- model_params(gono_params = gono_params(), fix_pop = TRUE)
  pars$mu <- pars$eta_h_step[] <- pars$eta_l_step[] <-  0

  mod <- model$new(pars, step = 0, n_particles = 10, seed = 1L)
  y <- mod$run(step_end = 10)
  y <- mod$transform_variables(y)

  expect_true(all(y$T == 0))
  expect_true(all(y$cum_treated == 0))
  expect_true(all(unlist(y) >= 0))
})

test_that("the foi is calculated correctly", {

  vei <- 0.123
  vax_params <- gonovax:::vax_params_xvwv(uptake = 0.5, dur = 1,
                                strategy = "VoA", vei = vei)
  pars <- model_params(gono_params = gono_params(),
                       vax_params = vax_params, fix_pop = TRUE)
  expect_true(length(pars$beta_step) > 0)


  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 10)
  y <- mod$transform_variables(y)

  # run forward one more step
  y1 <- mod$run(step_end = 11)
  y1 <- mod$transform_variables(y1)

  # unpack parameters
  pL <- pars$p[1]
  pH <- pars$p[2]
  NL <- colSums(y$N[1, , ])
  NH <- colSums(y$N[2, , ])
  C <- y$I + y$A + y$S
  CL <- c((1 - vax_params$vei) %*% C[1, , ])
  CH <- c((1 - vax_params$vei) %*% C[2, , ])
  eps <- pars$epsilon
  beta <- pars$beta_step[10]

  np <- pL * NL + pH * NH
  npL <- pL * NL / np
  npH <- pH * NH / np

  # calculate FOI
  foi_cross <- (1 - eps) * (npL * CL / NL + npH * CH / NH)
  foi_L <- pL * beta * (eps * CL / NL + foi_cross)
  foi_H <- pH * beta * (eps * CH / NH + foi_cross)

  expect_equal(y1$lambda[1, ], foi_L)
  expect_equal(y1$lambda[2, ], foi_H)

})

test_that("time-varying eta works as expected", {
  gono_pars <- gono_params()
  inc <- seq(1, 2, length.out = length(gono_pars$eta_h_step))
  gono_pars$eta_h_step <- gono_pars$eta_h_step * inc
  gono_pars$eta_l_step <- gono_pars$eta_l_step * inc / 2

  pars <- model_params(gono_params = gono_pars, fix_pop = TRUE)
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 10)
  y <- mod$transform_variables(y)

  expect_true(all(y$eta[1, ] == gono_pars$eta_l_step[10]))
  expect_true(all(y$eta[2, ] == gono_pars$eta_h_step[10]))

  # check can switch off screening in a group
  pars$eta_l_step[] <- 0
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y1 <- mod$run(step_end = 10)
  y1 <- mod$transform_variables(y1)
  expect_equal(sum(y1$cum_screened[1, , ]), 0)
  expect_true(all(y1$cum_screened[2, , ] > 0))
})


test_that("Bex model runs with no vaccination", {

  params0 <- model_params(gono_params = gono_params())
  mod0 <- model$new(params0, step = 0, n_particles = 5, seed = 1L)
  y0 <- mod0$run(step_end = 10)
  y0 <- mod0$transform_variables(y0)

  params1 <- model_params(gono_params = gono_params(),
                          vax_params = gonovax:::vax_params_xvwv(vbe = 0))
  mod1 <- model$new(params1, step = 0, n_particles = 5, seed = 1L)
  y1 <- mod1$run(step_end = 10)
  y1 <- mod1$transform_variables(y1)

  # check that nil vaccination gives same results as before
  expect_true(all(y1$U[, 1, , drop = FALSE] == y0$U))
  expect_true(all(y1$I[, 1, , drop = FALSE] == y0$I))
  expect_true(all(y1$A[, 1, , drop = FALSE] == y0$A))
  expect_true(all(y1$S[, 1, , drop = FALSE] == y0$S))
  expect_true(all(y1$T[, 1, , drop = FALSE] == y0$T))

  expect_true(all(y1$N[, 2, ] == 0))
  expect_equal(apply(y1$N, 3, sum), apply(y0$N, 3, sum))

})

test_that("Bex model runs with vbe", {

  # with perfect efficacy
  pars <- model_params(gono_params = gono_params(),
                       vax_params = gonovax:::vax_params_xvwv(vbe = 1, vea = 1),
                       fix_pop = TRUE)
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 10)
  y <- mod$transform_variables(y)

  expect_equal(y$U,
               array(c(505325, 89205, 293, 52, 0, 0,
                       505298, 89218, 303, 58, 0, 0,
                       505336, 89208, 270, 48, 0, 0,
                       505315, 89209, 290, 48, 0, 0,
                       505329, 89207, 282, 47, 0, 0), dim = c(2L, 3L, 5L)))
  # check all entrants are being vaccinated
  expect_equivalent(colSums(y$cum_vaccinated[, 1, ]),  y$entrants)
  # check no-one else is
  expect_true(all(y$cum_vaccinated[, -1, ] == 0))
  # check no compartments are leaking
  expect_true(all(apply(y$N, 3, sum) == 6e5))
  # check there are infections in unvaccinated group
  expect_true(all(y$incid[, 1, ] > 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$incid[, 2, ] == 0))
})

test_that("Check vaccination on screening in Bex model", {

  # with perfect efficacy
  vp <- gonovax:::vax_params_xvwv(vbe = 0, uptake = 1, strategy = "VoS",
                                  vea = 1, dur = 1)
  pars <- model_params(gono_params = gono_params(), vax_params = vp,
                       fix_pop = TRUE)
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 365)
  y <- mod$transform_variables(y)

  expect_equal(y$U,
               array(c(330701, 57987, 117622, 20574, 61320, 10598,
                       330850, 57762, 117559, 20642, 61216, 10723,
                       330621, 57769, 117687, 20624, 61279, 10681,
                       330648, 57604, 117929, 20757, 61029, 10738,
                       330072, 57912, 117909, 20601, 61626, 10591),
                     dim = c(2L, 3L, 5L)))
  # check all screened U, W are being vaccinated, but not V
  expect_equivalent(y$cum_vaccinated[, -2, ],  y$cum_screened[, -2, ])
  expect_true(any(y$cum_screened[, 2, ] > 0))
  expect_true(all(y$cum_vaccinated[, 2, ] == 0))

  # check no compartments are leaking
  expect_true(all(apply(y$N, 3, sum) == 6e5))
  # check there are infections in unvaccinated group
  expect_true(all(y$incid[, 1, ] > 0))
  # check there are no infections in vaccinated group
  expect_true(all(y$incid[, 2, ] == 0))
})

test_that("Check vaccination on diagnosis in Bex model", {
  # with imperfect efficacy
  vp <- gonovax:::vax_params_xvwv(vbe = 0, uptake = 1, strategy = "VoD",
                                  vea = 0.5, dur = 1)
  pars <- model_params(gono_params = gono_params(), vax_params = vp,
                       fix_pop = TRUE)
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 365)
  y <- mod$transform_variables(y)

  expect_equal(y$U,
               array(c(508542, 86919, 516, 1136, 411, 508,
                       508460, 86713, 550, 1185, 443, 527,
                       508577, 86995, 522, 1059, 400, 496,
                       508570, 86904, 543, 1142, 412, 571,
                       508577, 86933, 521, 1100, 428, 512),
                     dim = c(2L, 3L, 5L)))
  # check all treated U, W people are being vaccinated (but not V)
  expect_equivalent(y$cum_vaccinated[, -2, ],  y$cum_treated[, -2, ])
  expect_true(any(y$cum_treated[, 2, ] > 0))
  expect_true(all(y$cum_vaccinated[, 2, ] == 0))

  # check no compartments are leaking
  expect_true(all(apply(y$N, 3, sum) == 6e5))
  # check there are infections in unvaccinated group
  expect_true(all(y$incid[, 1, ] > 0))
})

test_that("Check vaccination on attendance in Bex model", {
  # with imperfect efficacy
  vp <- gonovax:::vax_params_xvwv(vbe = 0, uptake = 1, strategy = "VoA",
                                  vea = 0.5, dur = 1)
  pars <- model_params(gono_params = gono_params(), vax_params = vp,
                       fix_pop = TRUE)
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 365)
  y <- mod$transform_variables(y)

  expect_equal(y$U,
               array(c(329335, 56568, 118180, 21147, 62073, 11102,
                       329848, 56600, 117973, 21141, 61702, 11032,
                       330278, 56546, 118171, 21172, 61136, 11012,
                       329395, 56936, 118459, 21234, 61743, 10839,
                       330027, 56283, 118013, 21592, 61541, 10944),
                     dim = c(2L, 3L, 5L)))
  # check all treated U, W people are being vaccinated (but not V)
  expect_equivalent(y$cum_vaccinated[, -2, ],
                    y$cum_treated[, -2, ] + y$cum_screened[, -2, ])
  expect_true(any(y$cum_treated[, 2, ] > 0))
  expect_true(any(y$cum_screened[, 2, ] > 0))
  expect_true(all(y$cum_vaccinated[, 2, ] == 0))

  # check no compartments are leaking
  expect_true(all(apply(y$N, 3, sum) == 6e5))
  # check there are infections in unvaccinated group
  expect_true(all(y$incid[, 1, ] > 0))
})

test_that("incidence time series output correctly", {
  ## check with single parameter set
  vp <- gonovax:::vax_params_xvwv(vbe = 0, uptake = 1, strategy = "VoA",
                                  vea = 0.5, dur = 1)
  pars <- model_params(gono_params = gono_params(), vax_params = vp,
                       fix_pop = TRUE)
  mod <- model$new(pars, step = 0, n_particles = 5, seed = 1L)
  y <- mod$run(step_end = 365)
  y <- mod$transform_variables(y)
  expect_equal(y$treated, y$cum_treated)
  expect_equal(y$screened, y$cum_screened)
  expect_equal(y$diag_a, y$cum_diag_a)
  expect_equal(y$diag_s, y$cum_diag_s)
  expect_equal(y$incid, y$cum_incid)
  expect_equal(y$vaccinated, y$cum_vaccinated)

  y1 <- mod$run(step_end = 365 + 10L)
  y1 <- mod$transform_variables(y1)

  expect_equal(y$treated + y1$treated, y1$cum_treated)
  expect_equal(y$screened + y1$screened, y1$cum_screened)
  expect_equal(y$diag_a + y1$diag_a, y1$cum_diag_a)
  expect_equal(y$diag_s + y1$diag_s, y1$cum_diag_s)
  expect_equal(y$incid + y1$incid, y1$cum_incid)
  expect_equal(y$vaccinated + y1$vaccinated, y1$cum_vaccinated)
})
