## CARE Simulation
## Stephen A Lauer, Laura B Balzer, Nicholas G Reich
library(MASS)
library(dplyr)
library(tidyr)
library(doMC)
library(here)
registerDoMC()
setwd(here::here())
source("code/CARE-utilities.R")

set.seed(1)

population_n <- 100000
trial_n <- 96
n_sims <- 5000

## simulation with binary exposure and outcome
###############
# giant sample size
set.seed(1)
est_sim_FT <- gen_sim_data(population_n, trial=F, effect=T)
# true value causal parameter
ATE_FT <- mean(est_sim_FT$Y1 - est_sim_FT$Y0)

est_false_true <- foreach(j=seq(n_sims), .combine = rbind) %dopar% {
    if(j %% 1000 == 0)
        print(paste("Sim", j, "of", n_sims, Sys.time()))
    set.seed(j)
    sample_dat <- gen_sim_data(trial_n, trial=FALSE,
                                effect=TRUE)
    cbind(sim=j, trial=FALSE,
          suppressWarnings(
              do_estimation_sim(sample_dat = sample_dat,
                                 ATE=ATE_FT)))
}
est_false_true %>%
    group_by(trial, name) %>%
    summarise(mean(pt-ATE), var(pt), mean(var),
              mean(cover), mean(reject))

###############
set.seed(1)
est_sim_FF <- gen_sim_data(population_n, trial=F, effect=F)
# true value causal parameter
ATE_FF <- mean(est_sim_FF$Y1 - est_sim_FF$Y0)

est_false_false <- foreach(j=seq(n_sims), .combine = rbind) %dopar% {
    if(j %% 1000 == 0)
        print(paste("Sim", j, "of", n_sims, Sys.time()))
    set.seed(j)
    sample_dat <- gen_sim_data(trial_n, trial=FALSE,
                                effect=FALSE)
    cbind(sim=j, trial=FALSE,
          suppressWarnings(
              do_estimation_sim(sample_dat = sample_dat,
                                 ATE=ATE_FF)))
}
est_false_false %>%
    group_by(trial, name) %>%
    summarise(mean(pt-ATE), var(pt), mean(var),
              mean(cover), mean(reject))

###############
# giant sample size
set.seed(1)
est_sim_TT <- gen_sim_data(population_n, trial=T, effect=T)
# true value causal parameter
ATE_TT <- mean(est_sim_TT$Y1 - est_sim_TT$Y0)


est_true_true <- foreach(j=seq(n_sims), .combine = rbind) %dopar% {
    if(j %% 1000 == 0)
        print(paste("Sim", j, "of", n_sims, Sys.time()))
    set.seed(j)
    sample_dat <- gen_sim_data(trial_n, trial=TRUE,
                                effect=TRUE)
    cbind(sim=j, trial=TRUE,
          suppressWarnings(
              do_estimation_sim(sample_dat = sample_dat,
                                 ATE=ATE_TT)))
}
est_true_true %>%
    group_by(trial, name) %>%
    summarise(mean(pt-ATE), var(pt), mean(var),
              mean(cover), mean(reject))

###############
# giant sample size
set.seed(1)
est_sim_TF <- gen_sim_data(population_n, trial=T, effect=F)
# true value causal parameter
ATE_TF <- mean(est_sim_TF$Y1 - est_sim_TF$Y0)


est_true_false <- foreach(j=seq(n_sims), .combine = rbind) %dopar% {
    if(j %% 1000 == 0)
        print(paste("Sim", j, "of", n_sims, Sys.time()))
    set.seed(j)
    sample_dat <- gen_sim_data(trial_n, trial=TRUE,
                                effect=FALSE)
    cbind(sim=j, trial=TRUE,
          suppressWarnings(
              do_estimation_sim(sample_dat = sample_dat,
                                 ATE=ATE_TF)))
}
est_true_false %>%
    group_by(trial, name) %>%
    summarise(mean(pt-ATE), var(pt), mean(var),
              mean(cover), mean(reject))

sim_est <- bind_rows(est_false_true, est_false_false, est_true_false, est_true_true)

saveRDS(sim_est, file=paste0("data/simulation-output/sim-output-", Sys.Date(), ".rds"))
