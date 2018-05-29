#########################################################################################################
# Marlena Maziarz & Yingye Zheng
# University of Washington
# May 15, 2016
#
# Code used in:
# On Longitudinal Prediction with Time-to-Event Outcome: Comparison of Modeling Options, Biometrics, 2016
#########################################################################################################


# returns BLUP fitted values, slopes and intercepts for all marker measurements.
get.lme.blup.fitted <- function(data, lf, id, marker,
                                random, fixed){


  n.obs <- dim(data)[1]
  #sub.ids.full <- data$sub.id.consec
  sub.ids <- unique(sort(data[[id]]))
  n.subjects <- length(sub.ids)

  # outcome
  yi <- matrix(data[[marker]], ncol = 1)
  # fixed effects
  xi <- model.matrix(fixed, data)
  # random effects
  zi <- model.matrix(random, data)

  # get the beta hat and variance estimates from the lme fit (lf)
  beta.hat <- as.matrix(lf$coeff$fixed, ncol = 1)

  # random effects variance (D below)
  # between individual variance in intercepts and slopes
  # VarCorr is in the nlme package, calculates the estimated variances and correlations
  # between random effect terms in LME model

  v <- VarCorr(lf)
  v.diag <- v[1:(nrow(v)-1), 1] #first col (variance estimates), makes up diag

  if(ncol(v) >= 3){
    #if there are random effects besides slope
    v.off <- v[2:(nrow(v) - 1), 3:(ncol(v)), drop = FALSE]

    #Create covariance matrix D
    D <- diag(v.diag)
    for(i in 1:nrow(D)){
      for(j in 1:ncol(D)){

        if(i > j){

          D[i,j] <- sqrt(D[i,i]*D[j,j])*as.numeric(v.off[i-1, j])
          D[j,i] <- D[i,j]
        }
      }
    }
  }else{
    v.off <- NA
    D <- as.matrix(as.numeric(v.diag), nrow = 1, ncol = 1)
  }

  lf.sigma.sq <- lf$sigma^2

  yi.hat <- rep(NA, n.obs)
  random <- matrix(NA,nrow= n.obs,  ncol = ncol(zi))
  fixed <- matrix(NA,nrow= n.obs,  ncol = ncol(xi))

  for(i in 1:n.subjects){

    ix.curr <- data[[id]] == sub.ids[i]
    n.curr <- sum(ix.curr)
    yi.curr <- matrix(yi[ix.curr], ncol = 1)
    xi.curr <- xi[ix.curr, ,drop = FALSE ]
    zi.curr <- zi[ix.curr, , drop = FALSE ]

    yi.hat.curr <- rep(NA, n.curr)
    random.effects.curr <-  matrix(NA, nrow = n.curr, ncol = ncol(zi.curr) )
    fixed.effects.curr <- matrix(NA, nrow = n.curr, ncol =ncol(xi.curr))
    # browser()


    for(j in 1:n.curr){

      yi.curr.s <- matrix(yi.curr[1:j], ncol = 1)
      yi.cc <- complete.cases(yi.curr.s)

      xi.curr.s <- xi.curr[1:j, , drop = FALSE]
      zi.curr.s <- zi.curr[1:j, , drop = FALSE]

      # within subject variability
      # cov(ei) = sigma^2*Identity matrix (dim ni x ni)
      R <- lf.sigma.sq * diag(sum(yi.cc))

      # DZi'(Ri + ZiDZi')^(-1)(Yi - XiBhat)
      # only use the known y's to get bi.hat.s

      xi.curr.s.cc <- xi.curr.s[yi.cc,, drop = FALSE]
      yi.curr.s.cc <- yi.curr.s[yi.cc,,drop = FALSE]
      zi.curr.s.cc <-zi.curr.s[yi.cc,, drop = FALSE]

      bi.hat.s <- D%*%t(zi.curr.s.cc)%*%solve(R + zi.curr.s.cc%*%D%*%t(zi.curr.s.cc))%*%(yi.curr.s.cc - xi.curr.s.cc%*%beta.hat)

      out <- xi.curr.s[j,]%*%beta.hat + zi.curr.s[j,]%*%bi.hat.s

      random.effects.curr[j, ] <- c(bi.hat.s) #c(0,rep(1, length(beta.hat)-1))%*%beta.hat + sum(bi.hat.s[-1])
      fixed.effects.curr[j,] <- c(beta.hat)
      ## need to figure out how to implement here
      #myint.curr[j] <- beta.hat[1] + bi.hat.s[1]
      yi.hat.curr[j] <- out

    }


    # fixed (xi%*%beta.hat) and random Y hats (xi%*%beta.hat + zi%*%bi.hat)
    # just need the random fitted Y
    yi.hat[ix.curr] <- yi.hat.curr
    fixed[ix.curr, ] <-  fixed.effects.curr
    random[ix.curr, ] <- random.effects.curr

    rm(ix.curr, n.curr, yi.curr, xi.curr, zi.curr, yi.hat.curr, R, bi.hat.s)
  }




  return(list(fitted = yi.hat,
              fixed = as.data.frame(fixed),
              random = as.data.frame(random)
                    #, blup.slope = myslopes
                    ))
}

