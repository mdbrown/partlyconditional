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

  data <- data %>% arrange_(id, measurement.time)

  #fit blups on provided data
  blup.fit <- get.lme.blup.fitted(data,
                                  lf = blup.model,
                                  id = id,
                                  marker = marker,
                                  random = random.formula.noid,
                                  fixed = fixed.formula)
  out <- list(model = blup.model,
              fit = blup.fit,
              id = id,
              marker = marker,
              measurement.time = measurement.time,
              fixed = fixed.formula,
              random = random.formula.noid)
  class(out) <- "PC_BLUP"
  return(out)
}

#' @importFrom nlme lme
#' @importFrom nlme VarCorr
#' @export
#'
predict.PC_BLUP <- function(x, newdata,  ...){

  measurement.time = x$measurement.time
  id = x$id
  #change id

  orig_order <- 1:dim(newdata)[1]
  newdata <- newdata %>% arrange_(id, measurement.time)

  #fit blups on provided data
  blup.fit <- get.lme.blup.fitted(data = newdata,
                                  lf = x$model,
                                  id = x$id,
                                  marker = x$marker,
                                  random = x$random,
                                  fixed = x$fixed)

  blup.fit2 <-  cbind(newdata[,c(x$measurement.time, x$marker)],
                         fitted.blup = blup.fit$fitted,
                         fixed.coef =  blup.fit$fixed,
                         random.coef = blup.fit$random)

  blup.fit2 %>% arrange(orig_order)
}
