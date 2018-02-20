

#' Fit a partly conditional GLM (logistic) model
#'
#' Fit a partly conditional (PC) logistic model. PC models are helpful predictive tools for (medical) contexts where long-term follow-up is available and interest lies in predicting patients’ risks for a future adverse outcome using repeatedly measured predictors over time. These methods model the risk of an adverse event conditional on survival up to a landmark time and information accrued by that time. Methods to smooth markers through time using mixed effect models and BLUP estimates are also implemented.
#'
#' @param id name of numeric subject id in data
#' @param stime name of survival time, may be repeated across subj. id.
#' @param status name of survival status indicator, generally 1 = bad outcome/death, 0 alive/censored.
#' @param measurement.time name of time of measurement from baseline.
#' @param markers character vector consisting of marker names to include.
#' @param data data.frame with id, stime, status, measurement time  and marker variables. Observations with missing data will be removed.
#' @param prediction.time numeric value for the prediction time of interest to fit the PC logistic model.
#' @param use.BLUP a vector of logical variables, indicating for each marker whether best linear unbiased predictors (BLUPs) should be calculated. If true, a mixed effect model of the form `lme(marker ~ 1 + measurement.time, random = ~ 1 + measurement.time | id)`, will be fit univariately for each marker.
#' @param type.BLUP a vector of length equal to use.BLUP, indicating which type of BLUP information from the mixed effects model to use in the PC model. Options are  'fitted' (default) which uses the fitted BLUP estimates, 'intercept' for the BLUP intercept, 'slope' for the  or 'intercept_slope', which adds the BLUP intercept and slope to the PC model. 'intercept_slope' adds two terms to the PC model for each marker with corresponding use.BLUP slot equal to TRUE.
#' @param knots.measurement.time number of knots to use when modeling measurement.time using natural cubic splines (using function `ns`) in the PC Cox model. Set to 'NA' if no splines are to be used, which measurement.time will be included as a linear predictor in the PC Cox model.
#'
#'
#' @return
#'
#' a 'PC_GLM' model fit object, with elements
#'
#'  #' @return
#'
#' An object of class "PC_GLM" which is a list containing:
#'
#' \item{model.fit }{ A 'glm' object . Please note that the estimates of standard error associated with the model coefficients DO NOT incorporate the variation due to marker smoothing using BLUPs. }
#' \item{marker.blup.fit }{ A list of length equal to the number of markers. For each index where use.BLUP is TRUE, this list contains the 'nmle' object fit for the corresponding marker. }
#' \item{meas.time.spline}{ If knots.measurement.time is set, the measurement time spline basis matrix output from function 'ns'.  }
#' \item{call, variable.names, prediction.time, use.BLUP, knots.measurement.time}{Inputs from function call. }
#'
#' @examples
#' data(pc_data)
#'
#' pc.glm.1 <-  PC.GLM(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "log.meas.time",
#'  markers = c("marker_1", "marker_2"),
#'  prediction.time = 24,
#'  data = pc_data,
#'  knots.measurement.time = NA) #no spline used
#'
#' pc.glm.1
#'
#'pc.glm.1$model.fit #direct access to the coxph model object
#'
#'#fit a model using natural cubic splines to model measurement time
#'# and BLUPs to smooth marker measurements.
#' pc.glm.2 <-  PC.GLM(
#'  id = "sub.id",
#'  stime = "time",
#'  status = "status",
#'  measurement.time = "meas.time",
#'  markers = c("marker_1", "marker_2"),
#'  use.BLUP = c(TRUE, TRUE),
#'  prediction.time = 24,
#'  data = pc_data,
#'  knots.measurement.time = 3) #spline with three knots
#'
#' pc.glm.2
#'
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
                    markers,
                    data,
                    prediction.time,
                    use.BLUP = rep(FALSE, length(markers)),
                    type.BLUP = rep("fitted", length(markers)),
                    knots.measurement.time = NULL){

  call <- match.call()

  stopifnot(is.character(stime))
  stopifnot(is.character(status))
  #check data frame
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(prediction.time))
  stopifnot(length(prediction.time) ==1)

  stopifnot(length(use.BLUP)== length(markers))
  stopifnot(length(type.BLUP) == length(markers))
  if(!all(is.element(type.BLUP,
                     c("fitted",  "intercept", "slope", "intercept_slope"))) ){
    stop("The elements of the vector 'type.BLUP' must match 'fitted', 'intercept', 'slope', or 'intercept_slope'. See documentation for more details.")
  }
  #checkc to see if all variables are present in data
  tmpnames <- c(id, stime, status, measurement.time, markers)
  if(!all(is.element(tmpnames, names(data)))) stop(paste("'", tmpnames[which(!is.element(tmpnames, names(data)))], "' cannot be found in data.frame provided", sep = ""))

  #only keep complete cases
  mycomplete <- complete.cases(data[,tmpnames]);

  #check for missing data and throw it out, print a warning
  if(nrow(data)!=sum(mycomplete)){
    warning(paste(nrow(data)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    data <- data[mycomplete,]

  }

  #combine/filter data
  my.data <- data.frame(data[[id]],
                        data[[stime]],
                        data[[status]],
                        data[[measurement.time]],
                        data[[stime]] - data[[measurement.time]]
  )
  names(my.data) <- c(id, stime, status, measurement.time, "t.star")

  X.fit <- data[, markers, drop = FALSE]
  marker.names <- NULL
  ### blups
  mle.fit = list()


  for(xxind in 1:length(markers)){

    #calculate blups if needed
    if(use.BLUP[xxind]){
      cat("...Calculating Best Linear Unbiased Predictors (BLUP's) for marker: ", markers[[xxind]]) ; cat("\n")
      #calculate blup for each marker

      if(!is.numeric(X.fit[[xxind]])){
        stop("Marker must be numeric for BLUP calculation.")
      }

      marker.name <- names(X.fit)[xxind]
      my.data[[marker.name]] <- X.fit[,xxind]
      # get fitted values for the marker in the training set
      m1 <- try(lme(as.formula(paste0(marker.name, "~ 1 +",  measurement.time)),
                    random = as.formula(paste("~ 1 +",measurement.time, "|", id)),
                    data = my.data))
      mle.fit[[xxind]] <- m1


      tmp <- get.lme.blup.fitted(my.data, m1,
                                 id = id,
                                 marker = marker.name,
                                 measurement.time = measurement.time)


      names(tmp) <- paste0( marker.name,
                            c("_BLUP_fitted","_BLUP_intercept", "_BLUP_slope"))
      X.fit <- cbind(X.fit, tmp)

      #choose which blup elements to include in model
      #using type.BLUP
      if(type.BLUP[xxind] == "fitted"){
        tmp.ind <- c(1)
      }else if(type.BLUP[xxind] == "intercept"){
        tmp.ind <- c(2)
      }else if(type.BLUP[xxind] == "slope"){
        tmp.ind <- c(3)
      }else if(type.BLUP[xxind] == "intercept_slope"){
        tmp.ind <- c(2,3)
      }

      marker.names <- c(marker.names, names(tmp)[tmp.ind] )

    }else{
      mle.fit[[xxind]] <- NA
      marker.names <-  c(marker.names, names(X.fit)[xxind])

    }

  }

  my.data <- cbind(my.data, X.fit)
  ## end blups


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


  #fit cox model using raw marker data
  if(is.numeric(knots.measurement.time)){
    ## spline for measurement time
    meas.time.spline.basis <- ns(si,
                                 df = knots.measurement.time)
    xx <- as.data.frame(meas.time.spline.basis)
    names(xx) = paste0("meas.time.spline.basis", 1:ncol(xx))

    glm.data$working.dataset <- cbind(glm.data$working.dataset, xx)



    fmla <- as.formula(paste0("yi ~", paste(names(xx), collapse = " + "), "+",
                                    paste(marker.names, collapse = " + ")))



    fit <- glm(fmla, weights = wgt.IPW, data = glm.data$working.dataset, family = "binomial")


  }else{


    names(glm.data$working.dataset)[4] <- measurement.time
    fmla <- as.formula(paste0("yi ~", measurement.time, "+", paste(marker.names, collapse = " + ")))
    #fit model using blupss
    fit <- glm(fmla, weights = wgt.IPW, data = glm.data$working.dataset, family = "binomial")
    meas.time.spline.basis <- NULL #nothing to put here
  }



  ################


  out <- list( model.fit = fit,
               marker.blup.fit = mle.fit,
               meas.time.spline = meas.time.spline.basis,
               call = call,
               variable.names = c(id, stime, status, measurement.time, markers),
               prediction.time = prediction.time,
               use.BLUP = use.BLUP,
               type.BLUP = type.BLUP,
               knots.measurement.time = knots.measurement.time)
  class(out) <- "PC_GLM"

  out

}


#' print function for PC_GLM
#' @export

print.PC_GLM <- function(x, ...){

  cat("### Call:\n")
  print(x$call)
  cat("\n")


  if(any(x$use.BLUP)){
    cat("### BLUPs fit for marker(s): ")
    for(i in 1:length(x$use.BLUP)){
      if(x$use.BLUP[i]) {
        cat("", x$variable.names[i+4],"")

        #print(x$marker.blup.fit[[i]])

      }
    }
    cat("\n   See x$marker.blup.fit for details on mixed effect model fits. \n\n")
  }


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
    warning(paste(nrow(data)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    data <- data[mycomplete,]

  }
  my.data <- data.frame(id = newdata[[object$variable.names[1]]],
                        stime = newdata[[object$variable.names[2]]],
                        status = newdata[[object$variable.names[3]]],
                        measurement.time = newdata[[object$variable.names[4]]],
                        t.star = newdata[[object$variable.names[2]]] - newdata[[object$variable.names[4]]])

  names(my.data) <- c(object$variable.names[1:4], "t.star")
  my.data <- cbind(my.data, newdata[,object$variable.names[-c(1:4)], drop = FALSE])

  #fit blups if needed
  if(any(object$use.BLUP)){

    for(i in 1:length(object$use.BLUP)){
      if(object$use.BLUP[i]){

        marker.name <- object$variable.names[i + 4]
        # get fitted values for the marker in the training set
        m1 <- object$marker.blup.fit[[i]]

        tmp <-  get.lme.blup.fitted(my.data,
                                    m1,
                                    id = object$variable.names[1],
                                    marker = marker.name,
                                    measurement.time = object$variable.names[4])
        names(tmp) <- paste0( marker.name,
                              c("_BLUP_fitted","_BLUP_intercept", "_BLUP_slope"))
        newdata <- cbind(newdata, tmp)

      }
    }

  }


  #measurement times:


  if(is.numeric(object$knots.measurement.time)){
    ## spline for measurement time
    meas.time.spline.basis <- object$meas.time.spline

    meas.time.spline.basis <- ns(my.data[[object$variable.names[4]]],
                                 knots = attr(meas.time.spline.basis, "knots"),
                                 Boundary.knots = attr(meas.time.spline.basis, "Boundary.knots"),
                                 df = object$knots.measurement.time)

    #fit model using blups
    xx <- as.data.frame(meas.time.spline.basis)
    names(xx) = paste0("meas.time.spline.basis", 1:ncol(xx))
    newdata <- cbind( newdata, xx)


  }



  #calculate risk
  # we calculate risk for the test set, last observation before conditioning time
  id_name <-  object$variable.names[1]
  measurement_name <- object$variable.names[4]

  my.data$ind <- 1:nrow(my.data)

  ind <- my.data %>%
    group_by_(id_name) %>%
    select( "ind",
            !!id_name,
            !!measurement_name) %>%
    top_n(n  = 1) %>%
    ungroup()

  d.test.s.last.obs <- newdata[ind$ind, ]

  #return risk
  ##########
    aa <- data.frame( r = predict(object$model.fit,
                                  newdata = d.test.s.last.obs,
                                  type = "response"))




  names(aa) <- paste0("risk_", object$prediction.time)

  out <- cbind(d.test.s.last.obs, aa)

  return(out)





}
