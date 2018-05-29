#' Calculate best linear unbiased predictors (BLUPs) for a mixed effects model.
#'
#' Smooth predictors/markers through time using mixed effect models and estimate BLUPs.
#'
#' @param marker name of continuous marker to be modeled using mixed effects model. This will be used as the outcome for the mixed effects model.
#' @param measurement.time name of time of measurement from baseline.
#' @param fixed name of variables to be used as fixed effects.
#' @param random name of variables to be used as random effects.
#' @param id character name of subject id in data
#' @param data data.frame consisting of id, marker, measurement.time, fixed and random effects used for modeling. Observations with missing data in any of these variables will be removed.
#'
#'
#' @return
#'
#' An object of class "PC_BLUP" which is a list containing:
#'
#' \item{model}{ A 'lme' fit object from the nlme package.}
#' \item{fit}{  A data.frame consisting of id, measurement time, marker, 'fitted.blup' giving the fitted blup estimates  for the marker at each measurement time for the complete data. There are also columns giving the fixed and random components of the blup for each variable in the mixed effect model. These could be used to calculate, for example, an blup estimate for the rate of change in marker values over time.
#'}
#' \item{id, marker, measurement.time, fixed, random}{Function inputs.}
#' @examples
#'
#'
#' @importFrom nlme lme
#' @importFrom nlme VarCorr
#' @export

BLUP <- function( marker, measurement.time, fixed, random = NULL, id, data){
  #checks
  stopifnot(is.data.frame(data))
  #make sure marker.name, fixed.effects and random.effects are in data
  #stopifnot(is.element(marker.name, names(data)))
  #stopifnot(is.element(random.effects, names(data)))
  #stopifnot(is.element(fixed.effects, names(data)))

  #stopifnot(is.numeric(data[[marker.name]]))

  #throw out missing data....print warning

  #end checks


  fixed.formula <- as.formula(paste0(marker, "~ 1 + ",
                                     paste(fixed, collapse = " + " )))


  if(is.null(random) ){
    message("fitting mixed effects model with random intercept (~ 1| id)")
      random.formula <- as.formula(paste0("~ 1 ",
                                        "|" ,id))
      random.formula.noid <- as.formula("~ 1")
      random <- NULL
    }else{
      random.formula <- as.formula(paste0("~ 1 + ",
                                     paste(random, collapse = " + "),
                                     "|" ,id))
      random.formula.noid <- as.formula(paste0("~ 1 + ",
                                               paste(random, collapse = " + ")))
  }
  #fit nmle model
  blup.model <- try(lme(fixed = fixed.formula, data = data, random = random.formula))

  #need to sort by measurement.time before calling get.lme.blup.fitted
  data$orig_order <- 1:dim(data)[1]

  data <- data %>% arrange_(id, measurement.time)

  #fit blups on provided data
  blup.fit <- get.lme.blup.fitted(data,
                                  lf = blup.model,
                                  id = id,
                                  marker = marker,
                                  random = random.formula.noid,
                                  fixed = fixed.formula)

  names(blup.fit$fixed) <-  c("intercept", all.vars(fixed.formula)[-1])
  names(blup.fit$random) <-  c("intercept", all.vars(random.formula.noid))

  blup.fit2 <-  cbind(data[,c(id, measurement.time, marker, "orig_order")],
                      fitted.blup = blup.fit$fitted,
                      fixed.coef = blup.fit$fixed,
                      random.coef = blup.fit$random)

  blup.fit2 <- blup.fit2 %>% arrange(orig_order) %>% select(-orig_order)


  out <- list(model = blup.model,
              fit = blup.fit2,
              id = id,
              marker = marker,
              measurement.time = measurement.time,
              fixed = fixed.formula,
              random = random.formula.noid)
  class(out) <- "PC_BLUP"
  return(out)
}


#' Predict method for BLUP object.
#'
#' Calculate best linear unbiased predictors (BLUPs) using a previously fit mixed effects model for new observations.
#'
#' @param x object of class `PC_BLUP` output from the BLUP function
#' @param newdata data.frame of new observations for which to obtain fitted BLUP estimates.
#' @param \dots ignored.
#'
#' @return A data.frame consisting of id, measurement time, marker, 'fitted.blup' giving the fitted blup estimates for the marker at each measurement time for newdata. There are also columns giving the fixed and random components of the blup for each variable in the mixed effect model. These could be used to calculate, for example, an blup estimate for the rate of change in marker values over time.
#'
#' @importFrom nlme lme
#' @importFrom nlme VarCorr
#' @export

predict.PC_BLUP <- function(x, newdata,  ...){



  measurement.time = x$measurement.time
  id = x$id
  #change id

  newdata$orig_order <- 1:dim(newdata)[1]
  newdata <- newdata %>% arrange_(id, measurement.time)

  #fit blups on provided data
  blup.fit <- get.lme.blup.fitted(data = newdata,
                                  lf = x$model,
                                  id = x$id,
                                  marker = x$marker,
                                  random = x$random,
                                  fixed = x$fixed)


  names(blup.fit$fixed) <- c("intercept", all.vars(x$fixed)[-1])
  names(blup.fit$random) <- c("intercept", all.vars(x$random))

  blup.fit2 <-  cbind(newdata[,c(x$id, x$measurement.time, x$marker, "orig_order")],
                         fitted.blup = blup.fit$fitted,
                         fixed.coef =  blup.fit$fixed,
                         random.coef = blup.fit$random)

  blup.fit2 <- blup.fit2 %>% arrange(orig_order) %>% select(-orig_order)

  blup.fit2

}
