# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ghat.FUN <- function(time, status, new.time, type = "kaplan"){

  ranked.new.time <- rank(new.time)
  summary(survfit(Surv(time, status) ~ 1, se.fit = FALSE, type = type), sort(new.time))$surv[ranked.new.time]

}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# used as: IPW.FUN(xi.Ghat, di.Ghat, xi, round(si+tau, digits = data$n.digits), di, si)
IPW.FUN <- function(est.time, est.status, new.xi, new.s.tau, new.status, si, type = "kaplan"){

  numerator   <- new.status * (new.xi <= new.s.tau) + (new.xi >= new.s.tau)
  min.time    <- pmin(new.xi, new.s.tau)
  denominator <- Ghat.FUN(est.time, 1 - est.status, min.time, type = type)
  result      <- numerator/denominator
  result[is.na(result)] <- 0
  return(result)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in:
# vector.to.rep = vector to be replicated by row
# n.row = number of rows in the returned matrix made up of the vector vector.to.rep
# out:
# matrix n.rows x length(vector.to.rep)
VTM <- function(vector.to.rep, n.rows){
  matrix(vector.to.rep, ncol = length(vector.to.rep), nrow = n.rows, byrow = T)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g.logit <- function(x){
  return(1/(1 + exp(-x)))
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# keep only the subjects who survived up to time s
condition.on.time.s <- function(dc, data){

  sub.ids <- data$id
  sub.ids.event.by.s <- unique(sub.ids[data$stime <= dc$si]) # sub.ids are already sorted
  ix.event.by.s <- matrix(FALSE, nrow = length(sub.ids), ncol = 1)
  for(i in 1:length(sub.ids.event.by.s)){
    ix.event.by.s[sub.ids == sub.ids.event.by.s[i]] <- TRUE
  }
  return(data[!ix.event.by.s, ])
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# given a full dataset, return only data up to time s (as defined in the dc)
data.conditional.and.up.to.s <- function(dc, data){
  data.cond.s <- condition.on.time.s(dc, data)
  ix <- data.cond.s$measurement.time <= dc$si

  return(data.cond.s[ix, ])
}
