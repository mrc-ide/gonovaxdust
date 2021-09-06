##' @importFrom stats approx
##' @importFrom gonovax gono_params
gono_params <- function() {
  pars <- gonovax::gono_params(1)[[1]]
  tt <- pars$tt
  step <- seq(tt[1], tt[2], by = 1 / 365)
  pars$eta_h_step <- approx(x = tt, y = pars$eta_h_t, xout = step)$y
  pars$eta_l_step <- approx(x = tt, y = pars$eta_l_t, xout = step)$y
  pars$beta_step  <- approx(x = tt, y = pars$beta_t, xout = step)$y
  pars$tt <- pars$beta_t <- pars$eta_l_t <- pars$eta_h_t <- NULL
  pars
}

##' @importFrom gonovax model_params
model_params <- function(gono_params = NULL, demographic_params = NULL,
                         init_params = NULL, vax_params = NULL,
                         fix_pop = FALSE) {
  pars <- gonovax::model_params(gono_params, demographic_params,
                                init_params, vax_params)
  pars <- convert_vax_params(pars)
  pars$fix_pop <- fix_pop * 1

  pars
}

convert_vax_params <- function(pars) {

  tt <- pars$vax_t
  step <- seq(tt[1], tt[2], by = 1 / 365)
  pars$vax_step <- approx(x = tt, y = pars$vax_y, xout = step,
                          method = "constant")$y
  pars$vax_t <- pars$vax_y <- NULL

  pars$D <- diag(pars$w)
  pars$w <- sign(pars$w)

  if (pars$n_vax == 1) {
    f <- function(x) t(t(apply(x, 1, diag)))
  }  else {
    f <- function(x) t(apply(x, 1, diag))
  }
    pars$u_vbe <- f(pars$vbe)
    pars$u_vod <- f(pars$vod)
    pars$u_vos <- f(pars$vos)
    pars$vbe <- sign(pars$vbe)
    pars$vod <- sign(pars$vod)
    pars$vos <- sign(pars$vos)

  pars
}
