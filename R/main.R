#' Fit a partly conditional Cox model
#'
#' Fit a partly conditional (PC) Cox model. PC models are helpful predictive tools for (medical) contexts where long-term follow-up is available and interest lies in predicting patients’ risks for a future adverse outcome using repeatedly measured predictors over time. These methods model the risk of an adverse event conditional on survival up to a landmark time and information accrued by that time. Methods to smooth markers through time using mixed effect models and BLUP estimates are also implemented.
#'
#' @param id name of numeric subject id in data
#' @param stime name of survival time, may be repeated across subj. id.
#' @param status name of survival status indicator, generally 1 = bad outcome/death, 0 alive/censored.
#' @param measurement.time name of time of measurement from baseline.
#' @param predictors character vector of names for predictors and covariates to include in the model.
#' @param data data.frame with id, stime, status, measurement time  and predictor variables. Observations with missing data will be removed.
#'
#'
#' @return
#'
#' An object of class "PC_cox" which is a list containing:
#'
#' \item{model.fit }{ A 'coxph' fit object . Please note that the estimates of standard error associated with the model coefficients DO NOT incorporate the variation due to marker smoothing using BLUPs. }
#' \item{variable.names }{vector of variable names used to fit the model. }
#' \item{call}{Function call. }
#'#'
#' @examples
#'
#'data(pc_data)
#'#log transform measurement time for use in models.
#'pc_data$log.meas.time <- log10(pc_data$meas.time + 1)
#'
#'pc.model.1 <-  PC.Cox(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  predictors = c("log.meas.time", "marker_1", "marker_2"),
#'  data = pc_data)
#'
#'pc.model.1
#'
#'pc.model.1$model.fit #direct access to the coxph model object
#'
#'#fit a model using
#'# BLUPs to smooth marker measurements.
#'
#'#fit mixed effects model to use for blups
#'myblup.marker1 <- BLUP(marker = "marker_1",
#'                       measurement.time = "meas.time",
#'                       fixed = c("log.meas.time"),
#'                       random = c("log.meas.time"),
#'                       id = "sub.id" ,
#'                       data = pc_data )
#'
#'#adding blup estimates to pc_data
#'fitted.blup.values.m1 <- predict(myblup.marker1, newdata = pc_data)
#'
#'#fitted.blup.values.m1 includes the fitted blup value
#'pc_data$marker_1_blup <- fitted.blup.values.m1$fitted.blup
#'
#'
#'pc.model.2 <-  PC.Cox(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  predictors = c("log.meas.time", "marker_1_blup", "marker_2"),
#'  data = pc_data)
#'
#'pc.model.2
#'
#' @seealso \code{\link[partlyconditional]{BLUP}} \code{\link[partlyconditional]{PC.GLM}}
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
    cat(paste0("... removing ",  sum(wonky.times) , " observations where outcome time", tmpnames[2], " is less than measurement time ", tmpnames[4], ".\n"))
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
#' @return A data.frame matching that of newdata but with extra columns containing predicted risk estimates (one for each prediction.time). Risks are estimated for each 'prediction.time' from the measurement time observed on each row.
#'
#'
#' @examples
#'data(pc_data)
#'
#'pc_data$log.meas.time <- log10(pc_data$meas.time + 1)
#'
#'pc.model.1 <-  PC.Cox(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  predictors = c("log.meas.time", "marker_1", "marker_2"),
#'  data = pc_data)
#'
#'pc.model.1
#'
#'newdata.subj.6 <- pc_data[pc_data$sub.id ==6,]
#'
#'#estimate 12 and 24 month risk for each measurement time
#'predict(pc.model.1,
#'        newdata = newdata.subj.6,
#'        prediction.time = c(12, 24))
#'
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



