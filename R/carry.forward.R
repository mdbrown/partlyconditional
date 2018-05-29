
#' @export
#' @importFrom tidyr fill_
carry_forward <- function(data, id, measurement.time, new.times){
   #add new.times

   data <- data %>%
     group_by_(id) %>%
     do(pc_add_row(., id = id, measurement.time, new.times = new.times))

   #sort by measurement time
   data <- data %>% arrange_(id, measurement.time)

   #carry values forward
   data <- tidyr::fill_(data, names(data) , .direction = "down")
   data
}



pc_add_row <- function (.data, id, measurement.time, new.times )
{
  newdata <- tail(.data, 1)[rep(1, length(new.times)),]
  other_vars <- setdiff(names(newdata), c(id, measurement.time))
  newdata[other_vars] <- NA
  newdata[[measurement.time]] <- new.times

  #indicator for if the measurement time is already in .data
  already.present <- is.element(new.times, .data[[measurement.time]])

  newdata <- newdata[!already.present, ]

  rbind(.data, newdata)
}
