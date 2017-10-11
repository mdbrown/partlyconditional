


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

  ID <- data[[1]]

  ## unpack data
  marker.ind <- 6:ncol(data)

  xi <- round(data[[2]], digits = n.digits)    # surv time
  si <- round(data[[4]], digits = n.digits)    # meas.time
  zi <- data[ , 6:ncol(data), drop = FALSE]  # marker
  yi <- as.numeric(xi < si + tau)            # 1 if event before s + t
  di <- data[[3]] #censoring status

  ## extract censoring information for Ghat
  # index of first visit
  index   <- !duplicated(ID)
  xi.Ghat <- data[[2]][index]
  di.Ghat <- di[index] #censoring status
  ID.Ghat <- ID[index]

  # remove those who are censored and have been lost to followup before s + t
  rm.yes <- (1-di)*(data[[2]] < si + tau)

  working.dataset <- data.frame(ID = ID, yi = yi, Si = data[[4]], si = si, xi = xi, di = di, zi = zi)[rm.yes != 1, ]

  return(list(working.dataset = working.dataset,
              xi.Ghat = xi.Ghat,
              di.Ghat = di.Ghat,
              ID.Ghat = ID.Ghat,
              n.digits = n.digits))
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

  colnames(bs.vec) = paste(rep('s', sim.control$spline.s.df),
                           1:sim.control$spline.s.df, sep = '')

  # now zi becomes a vector of covariates (meas.time, marker)
  zi <- as.data.frame(cbind(bs.vec[1:length(xi), ], zi))

  # design matrix with intercept, spline function of measurement time, and marker
  ui <- as.matrix(cbind(rep(1, dim(zi)[1]), zi))

  xi.Ghat <- data$xi.Ghat
  di.Ghat <- data$di.Ghat

  wgt.IPW <- IPW.FUN(xi.Ghat, di.Ghat, xi, round(si+tau, digits = data$n.digits), di, si)


  fmla <- as.formula(paste("yi ~ ", paste(colnames(zi), collapse= "+")))

  fit <- glm(fmla, weights = wgt.IPW, data = zi, family = "binomial")




  return(fit)
}



