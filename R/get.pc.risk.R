
# returns risk estimated by a partly conditional model.
get.pc.risk <- function(dc, cohort, data.s, marker.names, use.timevarying){

    sp <- coxt.mix.wrapper(dc, cohort, data.s, marker.names, use.timevarying)

  if(is.null(sp)){
    return(NULL)
  }
  # return risk
  return(1-sp)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# coxt.mix
# wrapper to coxt.mix
coxt.mix.wrapper <- function(dc, cohort, data.s, marker.names, use.timevarying){

  meas.time.spline.basis <- ns(cohort$measurement.time,
                               df = dc$spline.s.df)
  z.covariate.matrix <- cbind(meas.time.spline.basis,
                              select(cohort, marker.names))

  # use the same knots and boundary knots to generate a basis of the measurement time for the test set
  meas.time.spline.basis.test <- ns(data.s$measurement.time,
                                    knots = attr(meas.time.spline.basis, "knots"),
                                    Boundary.knots = attr(meas.time.spline.basis, "Boundary.knots"),
                                    df = dc$spline.s.df)

  s.basis.eval.at.s <- predict(meas.time.spline.basis.test, dc$si)          # 1 x spline.s.df
  s.basis.eval.at.s.matrix <- VTM(s.basis.eval.at.s, dim(data.s)[1])        # n x spline.s.df
  z.covariate.matrix.test <- cbind(s.basis.eval.at.s.matrix,
                                   select(data.s, marker.names)) # n x (spline.s.df + 1)


  marker.index = which(is.element(names(z.covariate.matrix), marker.names ))[use.timevarying]


  if(any(use.timevarying)){
  out <- coxt.mix(x = cohort$t.star,
                      delta = cohort$status,
                      zc = z.covariate.matrix,
                      I.index = marker.index,
                      zc0 = z.covariate.matrix.test,
                      t0 = dc$pred.time,
                      mydf = dc$spline.t.df)


    out <- out$predp

  }else{

    #no time varying coefficients

    fit <- coxph(Surv(cohort$t.star, cohort$status) ~ . + cluster(cohort$id), data = z.covariate.matrix)

    names(z.covariate.matrix.test) <- names(z.covariate.matrix)
    out <- get.pc.risk_no_time_varying(dc, fit, z.covariate.matrix.test)

  }

}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to fit a cox proportional model with time-dependent and time-independent covariates

# INPUT
# (x, delta):       observed survival time subjected to right censoring
# z:                covariate matrix n by p
# constant.index:   vector indicate which covariates with constant regression coefficients
# control:          convergence criteria
# z0                covariate values at which the predictive probability is calculated
# t0                prediction time

# OUTPUT
# beta:       estimated constant regression coefficients
# predp:       survival probability
coxt.mix <- function(x, delta, zc, I.index, zc0, t0, mydf, control = 1e-10)
{
  # all covariates
  zc <- as.matrix(zc)

  # number of subjects
  n <- length(x)

  # number of parameters
  nc <- dim(zc)[2]

  # number of covariates in the prediction matrix
  nc0 <- dim(zc0)[1]

  # get fixed beta est
  beta.constant <- coxph(Surv(x, delta) ~ zc)$coef

  # survival times for subjects with events
  tt = x[delta == 1]

  # number of subjects with events
  ntt = length(tt)

  # covariates with time-varying effects
  zc.I = zc[, I.index, drop = FALSE]
  zc0.I = zc0[,I.index, drop = FALSE]

  # number of degrees of parameters for all the time varying covariates
  np.I <- length(I.index) * mydf

  # ***** temporary solution
  npp = np.I + nc

  # starting values for the time varying parameters
  beta.timevary <- rep(0.5, np.I)

  # spline basis for the event times
  ft = ns(tt, df = mydf)

  # placeholder for the beta's
  beta.est = c(beta.constant, beta.timevary)

  error <- 1

  fitting.iteration.counter <- 0 # how many iterations per fitting?

  # the fitting loop

  while(error > control)
  {
    # score function placeholder
    score <- 0

    # variance matrix placeholder
    v <- matrix(0, np.I + nc, np.I + nc)

    # predicted residual survival probability placeholder
    predp <- rep(0, nc0)

    # for each event time...
    for (i in 1:ntt) {
      the.t <- tt[i]          # current event time
      ft.t  <- ft[i, ]        # spline basis for current event time
     # zt    <- VTM(ft.t, n.rows =n)   # repeat the spline basis for current event time for all subjects (matrix n by mydf)

      # the interaction of the fixed covariate with time
      # does this work? n x mydf * n x length(I.index) - non-conformable arrays, no?
      # unless zc.I gets repped automatically to become mydf column matrix, in this case I.index has to equal 1.



      #we want to multiply ft.t by each row and column in zc.I
      # ft.t is a vector of length mydf
      # we repeat it once for each of our time-vary covariates

      zt <- rep(ft.t, ncol(zc.I))
      #expanding each column of zc.I to match the dimension of ft.t
      zc.I.repeat <- zc.I[, rep(1:ncol(zc.I), rep(mydf, ncol(zc.I)))]

      #multiply each row in zc.I.repeat by zt
      z.interaction <- sweep(zc.I.repeat, 2, zt, "*")


      #z.interaction <-  zc.I %>% mutate( zt*x)
      # all covariates, fixed and time-varying
      # columns: ns(meas.time), marker, marker:ns(current event time)
      zz <- cbind(zc, z.interaction)
      rm(z.interaction)

      # "new data", ie. covariates at the prediction time, fixed and time varying
      zc0.I.repeat <- zc0.I[, rep(1:ncol(zc0.I), rep(mydf, ncol(zc0.I)))]

      #multiply each row in zc.I.repeat by zt
      z0.interaction <- sweep(zc0.I.repeat, 2, zt, "*")
      zz0.t <- cbind(zc0, z0.interaction)
      rm(z0.interaction)

      #zz0.t <- cbind(zc0, VTM(ft.t, nc0) * zc0[, I.index])

      # n x num params %*% num params x 1 -> as vector -> length n
      effect <- as.vector(zz %*% beta.est)

      s0.temp = exp(effect) # length n, S^0(\beta, t)

      # sum exp(effect) among those at risk at the current event time
      # constant
      s0.t = sum(s0.temp[x >= the.t])

      # s0.temp gets repped to equal the same number of columns as zz
      # n x num params times length n -> transposed -> num params x n
      s1.temp = t(zz*s0.temp)  # params x n, S^1(\beta, t)

      # sum up zz*exp(effect) for those still at risk at current event time
      # vector of length p (num parameters)

      # this line causes an error sometimes... "dim(X) must have a positive length"
      # s1.t = apply(s1.temp[ , x >= the.t], 1, sum)
      s1.t = apply(as.matrix(s1.temp[ , x >= the.t], nrow = npp), 1, sum)

      # length p
      expectation.t = s1.t/s0.t

      # p x n %*% n x p with 0's for those not at risk
      # p x p
      s2.t = s1.temp %*% (zz * (x>= the.t)) # S^2(\beta, t)

      # information matrix (p x p)
      v = v + (-s2.t/s0.t + expectation.t %*% t(expectation.t))

      # score
      score = score + (zz[x == the.t,] - expectation.t)

      # if the current event time is less than prediction time
      if (the.t <= t0){ # shouldn't this also be > s?

        # predicted cumulative hazard
        predp = predp + exp(as.matrix(zz0.t) %*% matrix(beta.est, ncol = 1))/s0.t
      }

    }

    # variance (p x p)
    v <- solve(v)

    # betas (length p)
    beta.est <- as.vector(beta.est - v %*% score)
   cat(beta.est); cat("\n")
    error <- max(abs(score), na.rm = TRUE)
    fitting.iteration.counter <- fitting.iteration.counter + 1
    if(fitting.iteration.counter > 1e6) stop("Time varying coefficient estimates failing to converge after 1e6 iterations.")
  }

  # print(fitting.iteration.counter)
  return(list(beta.est = beta.est,
              predp = exp(-predp)
              ))
}

get.pc.risk_no_time_varying <- function(pred.time, fit, data.s){

  # dc  = list(pred.time = 2,
  #            si = 1,
  #            ti = 3,
  #            thresh = .5)
  # fit=checkmod
  # data.s=test.data.s

  sfit <- summary(survfit(formula = fit, newdata = data.s, se.fit = F),
                  times = pred.time)
  sp <- sfit$surv
  if(is.null(sp)){
    return(NULL)
  }
  # return risk
  sp
}


