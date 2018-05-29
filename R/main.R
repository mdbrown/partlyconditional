#' Fit a partly conditional Cox model
#'
#' Fit a partly conditional (PC) Cox model. PC models are helpful predictive tools for (medical) contexts where long-term follow-up is available and interest lies in predicting patients’ risks for a future adverse outcome using repeatedly measured predictors over time. These methods model the risk of an adverse event conditional on survival up to a landmark time and information accrued by that time. Methods to smooth markers through time using mixed effect models and BLUP estimates are also implemented.
#'
#' @param id name of numeric subject id in data
#' @param stime name of survival time, may be repeated across subj. id.
#' @param status name of survival status indicator, generally 1 = bad outcome/death, 0 alive/censored.
#' @param measurement.time name of time of measurement from baseline.
#' @param markers character vector consisting of marker names to include.
#' @param data data.frame with id, stime, status, measurement time  and marker variables. Observations with missing data will be removed.
#' @param use.BLUP a vector of logical variables, indicating for each marker whether best linear unbiased predictors (BLUPs) should be calculated. If true, a mixed effect model of the form `lme(marker ~ 1 + measurement.time, random = ~ 1 + measurement.time | id)`, will be fit univariately for each marker.
#' @param type.BLUP a vector of length equal to use.BLUP, indicating which type of BLUP information from the mixed effects model to use in the PC model. Options are  'fitted' (default) which uses the fitted BLUP estimates, 'intercept' for the BLUP intercept, 'slope' for the  or 'intercept_slope', which adds the BLUP intercept and slope to the PC model. 'intercept_slope' adds two terms to the PC model for each marker with corresponding use.BLUP slot equal to TRUE.
#' @param knots.measurement.time number of knots to use when modeling measurement.time using natural cubic splines (using function `ns`) in the PC Cox model. Set to 'NA' if no splines are to be used, and measurement.time will be included as a linear predictor in the PC Cox model.
#'
#'
#' @return
#'
#' An object of class "PC_cox" which is a list containing:
#'
#' \item{model.fit }{ A 'glm' object . Please note that the estimates of standard error associated with the model coefficients DO NOT incorporate the variation due to marker smoothing using BLUPs. }
#' \item{marker.blup.fit }{ A list of length equal to the number of markers. For each index where use.BLUP is TRUE, this list contains the 'nmle' object fit for the corresponding marker. }
#' \item{meas.time.spline}{ If knots.measurement.time is set, the measurement time spline basis matrix output from function 'ns'.  }
#' \item{call, variable.names, prediction.time, use.BLUP, knots.measurement.time}{Inputs from function call. }
#'#'
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
#'pc.model.1$model.fit #direct access to the coxph model object
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
#'pc.model.2
#'
#' @import dplyr
#' @import survival
#' @export

PC.Cox <- function( id,
                    stime,
                    status,
                    measurement.time,
                    predictors,
                    data,
                    additional.formula.pars = NULL){

  call <- match.call()

  stopifnot(is.character(stime))
  stopifnot(is.character(status))
  #check data frame
  stopifnot(is.data.frame(data))


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
  t.star <-  data[[stime]] - data[[measurement.time]]


  #sort by measurement time
  my.formula <- as.formula(paste0("Surv( t.star, ", status,") ~",
                                  paste(predictors, collapse = " + "),
                                    "+", "cluster(", id, ")",
                                    ifelse(is.null(additional.formula.pars), "", paste("+ ", additional.formula.pars)) ))
    #fit model using blups
    fit <- coxph(my.formula, data = data, x = TRUE)

    predictors.unique <- predictors[!is.element(predictors, c(id, stime, status, measurement.time))]
     out <- list( model.fit = fit,
               call = call,
               variable.names = c(id, stime, status, measurement.time, predictors.unique))
  class(out) <- "PC_cox"

  out

}

#' print function for PC_cox
#' @export

print.PC_cox <- function(x, ...){

  cat("### Call:\n")
  print(x$call)
  cat("\n")


  cat("### Partly conditional Cox model:\n")

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
#' @note Marker measurements repeated over time are only needed for PC models fit using smoothed BLUP marker values, otherwise, just providing the most recent marker measurements is adequate.
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
#'  urement.time = NA)
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
#'  urement.time = 3)
#'
#'#BLUPs are used in this pc model so marker trajectory is needed for
#'#BLUP estimation and risk prediction calculation
#'predict(pc.model.2, newdata = newdata.subj.6, prediction.time = c(12, 24))
#'
#'
#' @export

predict.PC_cox <- function(object, newdata,  prediction.time, ...){

  stopifnot(is.numeric(prediction.time))
  prediction.time = sort(prediction.time)

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

  #sort by measurement time
  id <- object$variable.names[[1]]
  measurement.time <-  object$variable.names[[4]]
  my.data <- my.data %>% arrange_(id, measurement.time)



  #calculate risk
  # we calculate risk for the test set, last observation before conditioning time

  newdata[["t.star"]] <-  newdata[[object$variable.names[2]]] - newdata[[object$variable.names[4]]]
  #return risk
  aa <- summary(survfit( object$model.fit, newdata = newdata, se.fit = F),
          times = prediction.time)
  aa <- as.data.frame(t(1-aa$surv))
  names(aa) <- paste0("risk_", prediction.time)

   out <- bind_cols(newdata, aa)

  return(out)





}



