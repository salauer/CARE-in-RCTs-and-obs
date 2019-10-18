## Utilities for the covariate adjusted residuals estimator
## and its use in both randomized and observational settings

#' Find recent file
#' Find the most recent file with a given starting string in a folder.
#'
#' @param name_start character string, first letters in file name
#' @param path character string, path to folder of interest, end with "/"
#' @param exclude character string, patterns to exclude from the file names of interest
#'
#' @return character string, path to most recent file
#' @export
#'
#' @examples
find_recent_file <- function(name_start, path, exclude=NULL){
    if(substring(path, nchar(path))!="/")
        warning('Path does not end with a "/", problems may ensue.')
    ## view all files of that name at that path
    file_list <- list.files(path=path,
                            pattern=paste0(name_start, "*"))
    ## remove files with unwanted patterns
    if(!is.null(exclude)){
        for(i in 1:length(exclude))
            file_list <- file_list[!grepl(pattern = exclude[i], file_list)]
    }
    if(length(file_list)==0){
        warning('File not found')
        return(NA)
    }
    ## view file info
    file_info <- file.info(paste0(path, file_list))
    ## find most recent file
    most_recent_file <- paste0(path,
                               file_list[which.max(file_info$mtime)])
    return(most_recent_file)
}

#' Get inference from estimator
#'
#' @param name character, name of the estimator
#' @param IC numeric vector, influence curve values for estimator
#' @param pt numeric value, the mean estimate from estimator (can be used instead of `IC``)
#' @param se numeric value, the standard error from estimator (can be used instead of `IC`)
#' @param alpha numeric value, the alpha level used for calculating confidence intervals and coverage
#' @param Psi.F numeric value, the true value that the estimator is trying to estimate
#'
#' @return data.frame with point estimate, standard error, confidence interval upper and lower bounds, whether the true value was covered, the p-value of the estimate, and whether or not the null hypothesis that the estimate is equal to zero is rejected at the specified alpha level.
#' @export
#'
#' @examples
get.inference<- function(name, pt, IC, se, alpha=0.05, ATE){
    if(missing(pt)){
        ## point estimate from IC curve
        pt <- mean(IC)
    }
    if(missing(se)){
        ## variance of asy lin est is var(IC)/n
        IC_var <- var(IC)/length(IC)
        se <- sqrt(IC_var)
    } else{
        IC_var <- se^2
    }
    ## testing and CI
    cutoff <- qnorm(alpha/2, lower.tail=F)  #*** change to t-dst for small sampes
    CI.lo <- pt - cutoff*se
    CI.hi <- pt + cutoff*se
    ## test statistic and pvalue
    pval <- 2*pnorm(abs(pt/se), lower.tail=F)  #*** change to t-dist for small sampes
    ## confidence interval coverage
    cover<- ( CI.lo <= ATE & ATE <= CI.hi )
    ## reject the null
    reject <- pval < alpha
    return(data.frame(name=name,
                      pt,
                      var=IC_var,
                      se,
                      CI.lo,
                      CI.hi,
                      cover,
                      pval,
                      reject))
}

#' Draws from an inverse logit function
#'
#' @param x numeric vector to be converted into probabilities by the inverse logit function
#'
#' @return
#' @export
#'
#' @examples
rexpit <- function(x,...){
    rbinom(n=length(x), size=1, prob=plogis(x,...))
}

#' Generate data for the simulation study of the CARE paper
#'
#' @param n numeric, number of clusters
#' @param effect logical, did the bednets reduce mortality?
#' @param trial logical, is this a randomized trial (T) or an observational study (F)?
#'
#' @return
#' @export
#'
#' @examples
gen_sim_data <- function(n, trial=F, effect=T){
    UY <- runif(n)
    W1 <- rnorm(n)
    W2 <- rnorm(n)
    W3 <- runif(n)
    W4 <- rbinom(n, 1, plogis(.5))
    if(trial){
        A <- rexpit(rep(0.5,n))
    } else{
        A <- rexpit(1 -.75*W1 - 2*W4 + 0.5*W2)
    }
    # summary(pscore)
    Y0 <- gen_sim_Y(W1, W2, W3, W4, 0, UY)
    if(effect){
        Y1 <- gen_sim_Y(W1, W2, W3, W4, 1, UY)
    } else{
        Y1 <- Y0
    }
    Y <- ifelse(A, Y1, Y0)
    return(data.frame(W1, W2, W3, W4, A, Y1, Y0, Y))
}

#' Generate Y values for the simulation study
#'
#' @param W1,W2,W3,W4 numeric vectors, covariates that influence the exposure (A) and the outcome (Y)
#' @param A numeric vector indicating whether a unit was expose (1) or unexposed (0)
#' @param UY numeric vector that indicates the level at which a unit receives the outcome or not
#'
#' @return
#' @export
#'
#' @examples
gen_sim_Y <- function(W1,W2,W3,W4,A,UY){
    p.out <- plogis(-.25 + .5*W1  -1*W3 + 2*W4 - 1.25*A - .5*A*W3)
    UY < p.out
}

#' Effect estimates for the simulation study
#'
#' @param sample_dat data frame, sample of clusters used for the experiment
#' @param y_formula character string, the formula for the covariate adjusted-residuals estimator to estimate the outcome of each cluster
#' @param y_family character string, the probability distribution for the covariate adjusted-residuals estimator to estimate the outcome of each cluster (default is 'binomial')
#' @param a_formula character string, the formula for the inverse probability of treatment weighting estimator to estimate the likelihood of each cluster being exposed
#' @param ATE numeric, the true average treatment effect across all clusters in the population
#'
#' @return data frame with all of the information related to estimating the effect of the exposure
#'
#' @export
#'
#' @examples
do_estimation_sim <- function(sample_dat,
                              y_formula="W1+W3+W4",
                              y_family="binomial",
                              a_formula="W1+W4",
                              ATE){
    txt<- sample_dat$A==1
    con <- sample_dat$A==0
    A1bar<- mean(txt)
    A0bar <- 1- A1bar
    Y <- sample_dat$Y
    Y1bar <- mean(Y[txt])
    Y0bar <- mean(Y[con])

    # unadj
    unadj <- Y1bar-Y0bar
    IC_unadj <- (txt/A1bar)*(Y-Y1bar) - (con/A0bar)*(Y-Y0bar)+Y1bar-Y0bar-unadj
    inf_unadj <- get.inference(name="unadj", pt=unadj,
                               IC=IC_unadj, ATE=ATE)

    # CARE
    glm_formula <- paste0("Y~",y_formula)
    sample_dat$Yhat <- predict(glm(as.formula(glm_formula),
                                   family=y_family,
                                   data=sample_dat),
                               type='response')
    Yhat <- sample_dat$Yhat
    CARE <- mean((txt/A1bar - con/A0bar)*(Y - Yhat))
    IC_CARE <- (txt/A1bar - con/A0bar)*(Y - Yhat)-CARE
    inf_CARE <- get.inference(name="CARE", pt=CARE,
                              IC=IC_CARE, ATE=ATE)

    ## IPW + CARE-IPW
    ## get propensity score
    txt_formula <- paste0("A~", a_formula)

    ## fit propensity score model
    A_glm <- glm(as.formula(txt_formula),
                 family='binomial', data=sample_dat)
    Ahat1 <- A_glm$fitted.values
    Ahat1 <- pmax(Ahat1, .01)
    Ahat1 <- pmin(Ahat1, .99)
    Ahat0 <- (1-Ahat1)
    IPW <- mean((txt/Ahat1 - con/Ahat0)*Y)
    IC_IPW <- (txt/Ahat1 - con/Ahat0)*Y-IPW
    inf_IPW <- get.inference(name="IPW", pt=IPW,
                              IC=IC_IPW, ATE=ATE)

    CARE_IPW <- mean((txt/Ahat1-con/Ahat0)*(Y - Yhat))
    IC_CARE_IPW <- (txt/Ahat1-con/Ahat0)*(Y-Yhat)-CARE_IPW
    inf_CARE_IPW <- get.inference(name="CARE_IPW",
                                  pt=CARE_IPW,
                                  IC=IC_CARE_IPW,
                                  ATE=ATE)

    ## return a data frame with all relevant info
    return(cbind(ATE=ATE, Y1bar, Y0bar, y_formula, a_formula,
                 A=sum(txt),
                 rbind(inf_unadj,
                       inf_CARE,
                       inf_IPW,
                       inf_CARE_IPW)))
}
