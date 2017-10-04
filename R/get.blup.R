#########################################################################################################
# Marlena Maziarz & Yingye Zheng
# University of Washington
# May 15, 2016
#
# Code used in:
# On Longitudinal Prediction with Time-to-Event Outcome: Comparison of Modeling Options, Biometrics, 2016
#########################################################################################################



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THIS IS NOT A GENERIC FUNCTION - it only works with the following LME model passed to it:
# data: dataframe with at least the following columns: marker, x1, x2 (fixed effect covariates), z1 (random effect covariate), and sub.id (subject ID)
# returns: BLUP fitted marker values
#
# assumes the LME model is:
# fixed: marker ~ 1 + x1 + x2 ---> yi ~ 1 + x1 + x2
# random: 1 + z1
# lf = linear mixed effects model fit (from lme() from nlme package) ie. the object returned by the following call:
# lf <- try(lme(marker ~ 1 + x1 + x2, random = ~ 1 + z1 | sub.id, data = data.training.set))
get.lme.blup.fitted.2.covariates <- function(data, lf){

    n.obs <- dim(data)[1]
    sub.ids.full <- data$sub.id
    sub.ids <- unique(sort(sub.ids.full))
    n.subjects <- length(sub.ids)

    # outcome
    yi <- matrix(data$marker, ncol = 1)
    # fixed effects
    xi <- matrix(cbind(1, data$x1, data$x2), ncol = 3, byrow = F)
    # random effects
    zi <- matrix(cbind(1, data$z1), ncol = 2, byrow = F)

    # get the beta hat and variance estimates from the lme fit (lf)
    beta.hat <- as.matrix(lf$coeff$fixed, ncol = 1)

    # random effects variance (D below)
    # between individual variance in intercepts and slopes
    v <- VarCorr(lf)
    v11 <- as.numeric(v[1,1])
    v22 <- as.numeric(v[2,1])
    v12 <- sqrt(v11 * v22) * as.numeric(v[2,3])
    D <- matrix(c(v11, v12, v12, v22), nrow = 2, byrow = T)

    lf.sigma.sq <- lf$sigma^2

    yi.hat <- rep(NA, n.obs)
    for(i in 1:n.subjects){
        ix.curr <- sub.ids.full == sub.ids[i]
        n.curr <- sum(ix.curr)
        yi.curr <- matrix(yi[ix.curr], ncol = 1)
        xi.curr <- matrix(xi[ix.curr, ], ncol = 3, byrow = F)
        zi.curr <- matrix(zi[ix.curr, ], ncol = 2, byrow = F)
        yi.hat.curr <- rep(NA, n.curr)
        for(j in 1:n.curr){
            yi.curr.s <- matrix(yi.curr[1:j], ncol = 1)
            xi.curr.s <- matrix(xi.curr[1:j, ], ncol = 3)
            zi.curr.s <- matrix(zi.curr[1:j, ], ncol = 2)
            # within subject variability
            # cov(ei) = sigma^2*Identity matrix (dim nixni)
            R <- lf.sigma.sq * diag(j)

            # DZi'(Ri + ZiDZi')^(-1)(Yi - XiBhat)
            bi.hat.s <- D%*%t(zi.curr.s)%*%solve(R + zi.curr.s%*%D%*%t(zi.curr.s))%*%(yi.curr.s - xi.curr.s%*%beta.hat)
            out <- xi.curr.s%*%beta.hat + zi.curr.s%*%bi.hat.s
            yi.hat.curr[j] <- out[j]
        }
        # fixed (xi%*%beta.hat) and random Y hats (xi%*%beta.hat + zi%*%bi.hat)
        # just need the random fitted Y
        yi.hat[ix.curr] <- yi.hat.curr
        rm(ix.curr, n.curr, yi.curr, xi.curr, zi.curr, yi.hat.curr, R, bi.hat.s)
    }
    return(yi.hat)
}






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BLUP fit with LME with 1 fixed and 1 random covariate
# THIS IS NOT A GENERIC FUNCTION - it only works with the following LME model passed to it:
# data: dataframe with at least the following columns: marker, x1 (fixed effect covariate), z1 (random effect covariate), and sub.id (subject ID)
# returns: BLUP fitted marker values
#
# assumes the LME model is:
# fixed: marker ~ 1 + x1 ---> yi ~ 1 + x1
# random: 1 + z1
# lf = linear mixed effects model fit (from lme() from nlme package), ie. the object returned by the following call:
# lf <- try(lme(marker ~ 1 + x1, random = ~ 1 + z1 | sub.id, data = data.training.set))

get.lme.blup.fitted.1.covariate <- function(data, lf, id, marker, measurement.time){


    n.obs <- dim(data)[1]
    sub.ids <- unique(sort(data[[id]]))
    n.subjects <- length(sub.ids)

    # outcome
    yi <- matrix(data[[marker]], ncol = 1)
    # fixed effects
    xi <- matrix(cbind(1, data[[measurement.time]]), ncol = 2, byrow = F)
    # random effects
    zi <- matrix(cbind(1, data[[measurement.time]]), ncol = 2, byrow = F)

    # get the beta hat and variance estimates from the lme fit (lf)
    beta.hat <- as.matrix(lf$coeff$fixed, ncol = 1)

    # random effects variance (D below)
    # between individual variance in intercepts and slopes
    v <- VarCorr(lf)
    v11 <- as.numeric(v[1,1])
    v22 <- as.numeric(v[2,1])
    v12 <- sqrt(v11 * v22) * as.numeric(v[2,3])
    D <- matrix(c(v11, v12, v12, v22), nrow = 2, byrow = T)

    lf.sigma.sq <- lf$sigma^2

    yi.hat <- rep(NA, n.obs)
    for(i in 1:n.subjects){
        ix.curr <- data[[id]] == sub.ids[i]
        n.curr <- sum(ix.curr)
        yi.curr <- matrix(yi[ix.curr], ncol = 1)
        xi.curr <- matrix(xi[ix.curr, ], ncol = 2, byrow = F)
        zi.curr <- matrix(zi[ix.curr, ], ncol = 2, byrow = F)
        yi.hat.curr <- rep(NA, n.curr)
        for(j in 1:n.curr){
            yi.curr.s <- matrix(yi.curr[1:j], ncol = 1)
            xi.curr.s <- matrix(xi.curr[1:j, ], ncol = 2)
            zi.curr.s <- matrix(zi.curr[1:j, ], ncol = 2)

            # within subject variability
            # cov(ei) = sigma^2*Identity matrix (dim ni x ni)
            R <- lf.sigma.sq * diag(j)

            # DZi'(Ri + ZiDZi')^(-1)(Yi - XiBhat)
            bi.hat.s <- D%*%t(zi.curr.s)%*%solve(R + zi.curr.s%*%D%*%t(zi.curr.s))%*%(yi.curr.s - xi.curr.s%*%beta.hat)
            out <- xi.curr.s%*%beta.hat + zi.curr.s%*%bi.hat.s
            yi.hat.curr[j] <- out[j]
        }
        # fixed (xi%*%beta.hat) and random Y hats (xi%*%beta.hat + zi%*%bi.hat)
        # just need the random fitted Y
        yi.hat[ix.curr] <- yi.hat.curr
        rm(ix.curr, n.curr, yi.curr, xi.curr, zi.curr, yi.hat.curr, R, bi.hat.s)
    }
    return(yi.hat)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# END OF GET THE BLUP FITTED MARKER VALUES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


