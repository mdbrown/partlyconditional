

#' Fit a partly conditional GLM (logistic) model
#'
#' Fit a partly conditional (PC) logistic model. PC models are helpful predictive tools for (medical) contexts where long-term follow-up is available and interest lies in predicting patients’ risks for a future adverse outcome using repeatedly measured predictors over time. These methods model the risk of an adverse event conditional on survival up to a landmark time and information accrued by that time.
#'
#' @param id name of numeric subject id in data
#' @param stime name of survival time, may be repeated across subj. id.
#' @param status name of survival status indicator, generally 1 = bad outcome/death, 0 alive/censored.
#' @param measurement.time name of time of measurement from baseline.
#' @param markers character vector consisting of marker names to include.
#' @param data data.frame with id, stime, status, measurement time  and marker variables. Observations with missing data will be removed.
#' @param prediction.time numeric value for the prediction time of interest to fit the PC logistic model.
#'
#' @return
#'
#' An object of class "PC_GLM" which is a list containing:
#'
#' \item{model.fit }{ A 'glm' object . Please note that the estimates of standard error associated with the model coefficients DO NOT incorporate the variation due to marker smoothing using BLUPs. }
#' \item{variable.names }{vector of variable names used to fit the model. }
#' \item{call}{Function call. }
#'#'#' \item{ prediction.time}{Inputs from function call. }
#'
#' @examples
#' data(pc_data)
#'#'#log transform measurement time for use in models.
#'pc_data$log.meas.time <- log10(pc_data$meas.time + 1)
#'
#' pc.glm.1 <-  PC.GLM(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  predictors = c("log.meas.time", "marker_1", "marker_2"),
#'  prediction.time = 24,
#'  data = pc_data)
#'
#' pc.glm.1
#'
#'pc.glm.1$model.fit #direct access to the glm model object
#'
#'#see function BLUP to fit mixed effects models and obtain BLUP-smoothed predictors.
#'
#' @seealso \code{\link[partlyconditional]{BLUP}} \code{\link[partlyconditional]{PC.Cox}}
#' @import dplyr
#' @import survival
#' @importFrom nlme lme
#' @importFrom nlme VarCorr
#' @importFrom splines ns
#' @export


PC.GLM <- function( id,
                    stime,
                    status,
                    measurement.time,
                    predictors,
                    data,
                    prediction.time){

  call <- match.call()

  stopifnot(is.character(stime))
  stopifnot(is.character(status))
  #check data frame
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(prediction.time))
  stopifnot(length(prediction.time) ==1)


  #checkc to see if all variables are present in data
  tmpnames <- c(id, stime, status, measurement.time, predictors)
  if(!all(is.element(tmpnames, names(data)))) stop(paste("'", tmpnames[which(!is.element(tmpnames, names(data)))], "' cannot be found in data.frame provided", sep = ""))

  #only keep complete cases
  mycomplete <- complete.cases(data[,tmpnames]);

  #check for missing data and throw it out, print a warning
  if(nrow(data)!=sum(mycomplete)){
    warning(paste(nrow(data)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    data <- data[mycomplete,]
  }

  #remove any observations with survival time less than meas.time
  wonky.times <- data[[tmpnames[2]]] < data[[tmpnames[4]]]
  if(any(wonky.times)) {
    cat(paste0("... removing ",  sum(wonky.times) , " observations where ", tmpnames[2], " is greater than outcome variable ", tmpnames[4], ".\n"))
  }

  data <- data[!wonky.times,]


  #combine/filter data
  my.data <- data.frame(data[[id]],
                        data[[stime]],
                        data[[status]],
                        data[[measurement.time]],
                        data[[stime]] - data[[measurement.time]]
  )
  names(my.data) <- c(id, stime, status, measurement.time, "t.star")

  X.fit <- data[, predictors, drop = FALSE]

  my.data <- cbind(my.data, X.fit)

  #need to sort by meas.time

  #### preGLM fun ####
  glm.data <- preGLM.FUN(my.data,  tau = prediction.time)

  #### GLM fun ###

  # i paste the code from GLMfun here because
  # i want to allow for modeling of measurement time
  #without knots.
  ###########
  yi <- glm.data$working.dataset$yi # status of event before s + t
  xi <- glm.data$working.dataset$xi # time
  di <- glm.data$working.dataset$di # status
  si <- glm.data$working.dataset$si # meas.time


  xi.Ghat <- glm.data$xi.Ghat
  di.Ghat <- glm.data$di.Ghat

  #some of these are rounded by preGLM function
  #not sure of the purpose but I'm keeping the rounded versions
  #this is why it looks like i needlessly copy the data here (and above)
  #mb oct 11, 2017
  wgt.IPW <- IPW.FUN(xi.Ghat, di.Ghat, xi,
                     round(si+prediction.time, digits = 6),
                     di,
                     si)

  fmla <- as.formula(paste0("yi ~",  paste(predictors, collapse = " + ")))
  #fit model using blupss
  fit <- glm(fmla, weights = wgt.IPW, data = glm.data$working.dataset, family = "binomial")


  ################


  out <- list( model.fit = fit,
               call = call,
               variable.names = c(id, stime, status, measurement.time, predictors),
               prediction.time = prediction.time)
  class(out) <- "PC_GLM"

  out

}


#' print function for PC_GLM
#' @export

print.PC_GLM <- function(x, ...){

  cat("### Call:\n")
  print(x$call)
  cat("\n")


  cat("### Partly conditional Logistic model\n")
  cat("###  for prediction time:", x$prediction.time, "\n")

  print(summary(x$model.fit)$coef)
}



#' Predict method for partly conditional cox models.
#'
#' Estimate predicted absolute n-year risk from a PC.Cox model object.
#'
#' @param object object of class 'PC_cox' fit using the 'PC.Cox' function
#' @param newdata data.frame with new data for which to calculate predictions. All variables used to fit the PC.Cox model must be present. Observations with missing data will be removed.
#' @param prediction.time vector of prediction times (from marker measurement time) to estimate future risk. Prediction time is defined from time of measurement for each individual trajectory.
#'
#' @return A data.frame consisting of one row for each individual (grouped by id), with predicted risk estimates at the final observed measurement time for each individuals' marker trajectory. Risks are estimated for each 'prediction.time' from the last measurement time observed for each individual. Marker measurements, BLUP point estimates and spline bases for measurement times (if applicable) are provided on each observation as well.
#'
#'
#' @note Marker measurements repeated over time are only needed, otherwise, just providing the most recent marker measurements is adquate.
#'
#' @examples
#'data(pc_data)
#'
#'pc.model.1 <-  PC.Cox(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "log.meas.time",
#'  markers = c("marker", "marker_2"),
#'  data = pc_data,
#'  use.BLUP = c(FALSE, FALSE),
#'  knots.measurement.time = NA)
#'
#'pc.model.1
#'
#'newdata.subj.6 <- pc_data[pc_data$sub.id ==6,]
#'#last marker measured at time '54', so predictions
#'#will be conditional on surviving to time '54'
#'newdata.subj.6
#'
#'#estimate 12 and 24 month risk from month 54
#'#this pc model doesn't include BLUPs, so only the final marker
#'#measurements are considered.
#'predict(pc.model.1, newdata = newdata.subj.6, prediction.time = c(12, 24))
#'
#'
#'#fit a model using natural cubic splines to model measurement time
#'# and BLUPs to smooth marker measurements.
#'pc.model.2 <-  PC.Cox(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  markers = c("marker", "marker_2"),
#'  data = pc_data,
#'  use.BLUP = c(TRUE, TRUE),
#'  knots.measurement.time = 3)
#'
#'#BLUPs are used in this pc model so marker trajectory is needed for
#'#BLUP estimation and risk prediction calculation
#'predict(pc.model.2, newdata = newdata.subj.6, prediction.time = c(12, 24))
#'
#'
#' @export

predict.PC_GLM <- function(object, newdata, ...){


  stopifnot(is.data.frame(newdata))
  #check if newdata has the right variables
  if(!all(is.element(object$variable.names, names(newdata)) )){
    stop("variable(s): {", paste(object$names[which(!is.element(object$variable.names, names(newdata)))], collapse = ", "), "}, not found in 'newdata' ")
  }

  #only keep complete cases
  mycomplete <- complete.cases(newdata[,object$variable.names]);

  #check for missing data and throw it out, print a warning
  if(nrow(newdata)!=sum(mycomplete)){
    warning(paste(nrow(newdata)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    newdata <- newdata[mycomplete,]

  }

  #remove any observations with survival time less than meas.time
  wonky.times <-  newdata[[object$variable.names[2]]] < newdata[[object$variable.names[4]]]
  if(any(wonky.times)) {
    warning(paste0("... removing ",  sum(wonky.times) , " observations where ", object$variable.names[4], " is greater than outcome variable ", object$variable.names[2], ".\n"))
  }

  newdata <- newdata[!wonky.times,]


  my.data <- data.frame(id = newdata[[object$variable.names[1]]],
                        stime = newdata[[object$variable.names[2]]],
                        status = newdata[[object$variable.names[3]]],
                        measurement.time = newdata[[object$variable.names[4]]],
                        t.star = newdata[[object$variable.names[2]]] - newdata[[object$variable.names[4]]])

  names(my.data) <- c(object$variable.names[1:4], "t.star")
  my.data <- cbind(my.data, newdata[,object$variable.names[-c(1:4)], drop = FALSE])


  #return risk
  ##########
  aa <- data.frame( r = predict(object$model.fit,
                                  newdata = newdata,
                                  type = "response"))


  names(aa) <- paste0("risk_", object$prediction.time)

  out <- bind_cols(newdata, aa)

  return(out)





}
