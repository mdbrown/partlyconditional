#' Fit a partly conditional Cox model
#'
#' later
#'
#' @param id name of numeric subject id in data
#' @param stime name of survival time, may be repeated across subj. id.
#' @param status name of survival status indicator, generally 1 = bad outcome/death, 0 alive/censored.
#' @param measurement.time name of time of measurement from baseline.
#' @param markers data.frame consisting of n x p markers with names.
#' @param data data.frame
#' @param use.BLUP a vector of logical variables, indicating for each marker whether best linear unbiased predictors (BLUPs) should be calculated. If true, a mixed effect model of the form `lme(marker ~ 1 + measurement.time, random = ~ 1 + measurement.time | id)`, will be fit univariately for each marker.
#' @param knots.measurement.time number of knots to use when modeling measurement.time using natural cubic splines (using function `ns`) in the partly conditional Cox model. Set to 'NA' if no splines are to be used.
#'
#'
#' @return a 'PC_cox' model fit object, consisting of
#'
#' @examples
#'
#'
#' @export
#' @import dplyr
#' @import nmle

PC.cox <- function( id,
                    stime,
                    status,
                    measurement.time,
                    markers,
                    data,
                    use.BLUP = rep(FALSE, length(markers)),
                    knots.measurement.time = NULL){

  call <- match.call()

  stopifnot(is.character(stime))
  stopifnot(is.character(status))
  #check data frame
  stopifnot(is.data.frame(data))

  stopifnot(length(use.BLUP)== length(markers))
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

    X.fit <- data[, markers]


    mle.fit = list()

    names(X.fit)[use.BLUP] <- paste0(names(X.fit)[use.BLUP], "_BLUP")

    for(xxind in 1:length(markers)){

      #calculate blups if needed
    if(use.BLUP[xxind]){
    cat("...Calculating Best Linear Unbiased Predictors (BLUP's) for marker: ", markers[[xxind]]) ; cat("\n")
    #calculate blup for each marker


      marker.name <- names(X.fit)[xxind]
      my.data[[marker.name]] <- X.fit[,xxind]
      # get fitted values for the marker in the training set
      m1 <- try(lme(as.formula(paste0(marker.name, "~ 1 +",  measurement.time)),
                    random = as.formula(paste("~ 1 +",measurement.time, "|", id)),
                    data = my.data))
      mle.fit[[xxind]] <- m1
      X.fit[[xxind]] <- get.lme.blup.fitted.1.covariate(my.data, m1,
                                                        id = id,
                                                        marker = marker.name,
                                                        measurement.time = measurement.time)

    }else{
      mle.fit[[xxind]] <- NA
    }

    }
    marker.names <- names(X.fit)
    my.data <- cbind(my.data, X.fit)


  #fit cox model using raw marker data
  if(is.numeric(knots.measurement.time)){
    ## spline for measurement time
    meas.time.spline.basis <- ns(my.data[[measurement.time]],
                                 df = knots.measurement.time)
    xx <- as.data.frame(meas.time.spline.basis)
    names(xx) = paste0("meas.time.spline.basis", 1:ncol(xx))

    my.data <- cbind(my.data, xx)

    my.formula <- as.formula(paste0("Surv( t.star, ", status,") ~", names(xx), "+",
                                    paste(marker.names, collapse = " + "),
                                    "+", "cluster(", id, ")"))
    #fit model using blups
    fit <- coxph(my.formula, data = my.data, x = TRUE)


  }else{

    my.formula <- as.formula(paste0("Surv( t.star, ", status,") ~", measurement.time, "+", paste(marker.names, collapse = " + ") ,
                              "+", "cluster(", id, ")"))
    #fit model using blups
    fit <- coxph(my.formula, data = my.data)
    meas.time.spline.basis <- NULL #nothing to put here
  }

  out <- list( model.fit = fit,
               marker.blup.fit = mle.fit,
               meas.time.spline = meas.time.spline.basis,
               call = call,
               variable.names = c(id, stime, status, measurement.time, markers),
               use.BLUP = use.BLUP,
               knots.measurement.time = knots.measurement.time)
  class(out) <- "PC_cox"

  out

}

#' print function for PC_cox

print.PC_cox <- function(x, ...){

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
  cat("### Partly conditional Cox model:\n")

  print(summary(x$model.fit)$coef)
}

#' Predict method for partly conditional cox models.
#'
#' Estimate predicted absolute n-year risk from a PC.cox model object.
#'
#' @param object object of class 'PC.cox' fit using the 'PC.cox' function
#' @param newdata data.frame with variables
#' @param prediction.time time from baseline of prediction. This time should exceed all measurement times found in newdata.
#'
#' @return
#'
#' @export

predict.PC_cox <- function(object, newdata, prediction.time, ...){

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
  if(nrow(data)!=sum(mycomplete)){
    warning(paste(nrow(data)-sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", sum(mycomplete)))
    data <- data[mycomplete,]

  }
  my.data <- data.frame(id = newdata[[object$variable.names[1]]],
                        stime = newdata[[object$variable.names[2]]],
                        status = newdata[[object$variable.names[3]]],
                        measurement.time = newdata[[object$variable.names[4]]],
                        t.star = newdata[[object$variable.names[2]]] - newdata[[object$variable.names[4]]])

  names(my.data) <- c(object$variable.names[1:4], "t.star")
  my.data <- cbind(my.data, newdata[,object$variable.names[-c(1:4)]])

  #fit blups if needed
  if(any(object$use.BLUP)){

    for(i in 1:length(object$use.BLUP)){
      if(object$use.BLUP[i]){

        marker.name <- object$variable.names[i + 4]
        # get fitted values for the marker in the training set
        m1 <- object$marker.blup.fit[[i]]

        newdata[[paste0(marker.name, "_BLUP")]] <- get.lme.blup.fitted.1.covariate(my.data,
                                                                                   m1,
                                                                                   id = object$variable.names[1],
                                                                                   marker = marker.name,
                                                                                   measurement.time = object$variable.names[4])


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
  aa <- summary(survfit( object$model.fit, newdata = d.test.s.last.obs, se.fit = F),
          times = prediction.time)
  aa <- as.data.frame(t(1-aa$surv))
  names(aa) <- paste0("risk_", prediction.time)

   out <- cbind(d.test.s.last.obs, aa)

  return(out)





}

