


get.logistic.risk <- function(dc, cohort, data.s){
  data.tmp <- data.frame(ID = cohort$sub.id,
                         si = cohort$meas.time,
                         zi = cohort$marker,
                         xi = cohort$time,
                         di = cohort$status)

  glm.data <- preGLM.FUN(data.tmp, covariate.ix = 3, tau = dc$pred.time)

  glm.out  <- GLMt.FUN(glm.data,
                       s.vec = dc$si,
                       tau   = dc$pred.time,
                       y.vec = data.s$marker,
                       sim.control = dc)

  risk <- glm.out$predp.tau.s0.y0
  if(is.null(risk)){
    return(NULL)
  }
  return(risk)
}







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# example call:
# data.tmp <- data.frame(ID = cohort$sub.id, si = cohort$meas.time, zi = cohort$marker, xi = cohort$time, di = cohort$status)
# glm.data <- preGLM.FUN(data.tmp, covariate.index = 3, tau = dc$ti)
preGLM.FUN <- function(data, covariate.ix, tau, n.digits = 6){

  ID <- data$ID

  ## unpack data
  xi <- round(data$xi, digits = n.digits)    # surv time
  si <- round(data$si, digits = n.digits)    # meas.time
  zi <- data[ , covariate.ix, drop = FALSE]  # marker
  yi <- as.numeric(xi < si + tau)            # 1 if event before s + t

  ## extract censoring information for Ghat
  # index of first visit
  index   <- !duplicated(ID)
  xi.Ghat <- data$xi[index]
  di.Ghat <- data$di[index]
  ID.Ghat <- ID[index]

  # remove those who are censored and have been lost to followup before s + t
  rm.yes <- (1-data$di)*(data$xi < si + tau)

  working.dataset <- data.frame(ID = ID, yi = yi, Si = data$si, si = si, xi = xi, di = data$di, zi = zi)[rm.yes != 1, ]

  return(list(working.dataset = working.dataset, xi.Ghat = xi.Ghat, di.Ghat = di.Ghat, ID.Ghat = ID.Ghat, n.digits = n.digits))
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# example call:
# GLMt.FUN(glm.data, s.vec = dc$si, tau = dc$ti, y.vec = data.s$marker, sim.control = dc)
GLMt.FUN <- function(data, s.vec, tau, y.vec, vb = NULL, sim.control){
  yi <- data$working.dataset$yi # status of event before s + t
  xi <- data$working.dataset$xi # time
  di <- data$working.dataset$di # status
  si <- data$working.dataset$si # meas.time
  zi <- data$working.dataset$zi # marker
  bs.vec = ns(si, df = sim.control$spline.s.df)

  colnames(bs.vec) = paste(rep('s', sim.control$spline.s.df), 1:sim.control$spline.s.df, sep = '')

  # now zi becomes a vector of covariates (meas.time, marker)
  zi <- as.data.frame(cbind(bs.vec[1:length(xi), ], zi))

  # design matrix with intercept, spline function of measurement time, and marker
  ui <- as.matrix(cbind(rep(1, dim(zi)[1]), zi))

  xi.Ghat <- data$xi.Ghat
  di.Ghat <- data$di.Ghat

  wgt.IPW <- IPW.FUN(xi.Ghat, di.Ghat, xi, round(si+tau, digits = data$n.digits), di, si)


  fmla <- as.formula(paste("yi ~ ", paste(colnames(zi), collapse= "+")))

  fit <- glm(fmla, weights = wgt.IPW, data = zi, family = "binomial")

  beta  <- fit$coef
  score <- ui%*%beta
  predp.tau.i     <- g.logit(score)
  predp.tau.s0.y0 <- NULL

  if(!is.null(s.vec)){
    bs0 <- predict(bs.vec, s.vec)
    ns0 <- length(s.vec)
    ny0 <- length(y.vec)

    ui0 <- NULL
    for(j in 1:ny0){
      ui0 <- rbind(ui0, cbind(bs0, rep(y.vec[j], ns0)))
    }

    ui0 <- cbind(1, ui0)
    predp.tau.s0.y0 <- g.logit(ui0 %*% beta)
  }

  return(list(beta = beta, predp.tau.s0.y0 = predp.tau.s0.y0, predp.tau.i = predp.tau.i))
}


