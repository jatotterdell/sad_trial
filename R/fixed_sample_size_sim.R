# run_basic_sequential_sim.R

source("R/basic_sequential_sim.R")

# Set up a NULL and ALTERNATIVE scenario
pt <- ptrunc(0:56, spec = "logis", a = -0.5, b = 56, 17.5, 3)
p0 <- diff(c(0, pt))
m0 <- sum(p0*0:56)
s0 <- sqrt(sum(p0*(0:56 - m0)^2))
or <- c(1/16, 1/8, 1/4, 1/2)
p1 <- pomodm(p = p0, odds.ratio = or[1])
m1 <- sum(p1*0:56)
s1 <- sqrt(sum(p1*(0:56 - m1)^2))
p2 <- pomodm(p = p0, odds.ratio = or[2])
m2 <- sum(p2*0:56)
p3 <- pomodm(p = p0, odds.ratio = or[3])
m3 <- sum(p3*0:56)
p4 <- pomodm(p = p0, odds.ratio = or[4])
m4 <- sum(p4*0:56)


p_alt00 <- rbind(p0, p0, p0)
p_alt11 <- rbind(p0, p1, p1)
p_alt22 <- rbind(p0, p2, p2)
p_alt33 <- rbind(p0, p3, p3)
p_alt44 <- rbind(p0, p4, p4)
p_alt01 <- rbind(p0, p0, p1)
p_alt21 <- rbind(p0, p2, p1)
p_alt31 <- rbind(p0, p3, p1)
p_alt41 <- rbind(p0, p4, p1)
p_alt02 <- rbind(p0, p0, p2)
p_alt12 <- rbind(p0, p1, p2)
p_alt32 <- rbind(p0, p3, p2)
p_alt42 <- rbind(p0, p4, p2)
p_alt03 <- rbind(p0, p0, p3)
p_alt13 <- rbind(p0, p1, p3)
p_alt23 <- rbind(p0, p2, p3)
p_alt43 <- rbind(p0, p4, p3)
p_alt04 <- rbind(p0, p0, p4)
p_alt14 <- rbind(p0, p1, p4)
p_alt24 <- rbind(p0, p2, p4)
p_alt34 <- rbind(p0, p3, p4)
p_alt10 <- rbind(p0, p1, p0)
p_alt20 <- rbind(p0, p2, p0)
p_alt30 <- rbind(p0, p3, p0)
p_alt40 <- rbind(p0, p4, p0)

looks <- list(150, 200, 250)
p <- list(
  p_alt01, p_alt11, p_alt21, p_alt31, p_alt41,
  p_alt02, p_alt12, p_alt22, p_alt32, p_alt42,
  p_alt03, p_alt13, p_alt23, p_alt33, p_alt43,
  p_alt04, p_alt14, p_alt24, p_alt34, p_alt44,
  p_alt00, p_alt10, p_alt20, p_alt30, p_alt40)
epsilon <- list(0.85, 0.9, 0.95)
mu0 <- list(c(18, 0, 0))
# Sigma0 <- list(diag(10^2, 3))
Sigma0 <- list(diag(c(10^2, 3.9^2, 3.9^2)))
drop_rule <- list("one")
MID <- list(5)
NonInfMarg <- list(3.3)
scenarios <- scenario_list(n_list(p, epsilon, mu0, Sigma0, drop_rule, MID, NonInfMarg, looks))
scenarios_dt <- rbindlist(lapply(scenarios, function(x)
  as.data.table(lapply(x, function(a) ifelse(length(a)>1, list(a), I(a))))))[, scenario := 1:.N]
scenarios_dt[, effect1 := rep(
  c("null", "very large", "large", "moderate", "small",
    "null", "very large", "large", "moderate", "small",
    "null", "very large", "large", "moderate", "small",
    "null", "very large", "large", "moderate", "small",
    "null", "very large", "large", "moderate", "small"), 9)]
scenarios_dt[, effect2 := rep(c(
  "very large", "very large", "very large", "very large", "very large",
  "large", "large", "large", "large", "large",
  "moderate", "moderate", "moderate", "moderate", "moderate",
  "small", "small", "small", "small", "small",
  "null", "null", "null", "null", "null"), 9)]

lvls <-  c("null", "small", "moderate", "large", "very large")
scenarios_dt[, `:=`(effect1 = factor(effect1, levels = lvls),
                    effect2 = factor(effect2, levels = lvls))]

res <- rbindlist(lapply(scenarios, function(x) sim_replicate(1e3, trial_par = x, cores = 15)), idcol = "scenario")
res_sum <- res[, lapply(.SD, mean), keyby = .(scenario, analysis)]
dat_trials <- res[scenarios_dt, on = .(scenario)]
dat <- res_sum[scenarios_dt, on = .(scenario)]

saveRDS(scenarios_dt, file = "~/out_files/mental_health_sims/mid_fixed_ss_scenarios.rds")
saveRDS(res, file = "~/out_files/mental_health_sims/mid_fixed_ss_sims.rds")
saveRDS(res_sum, file = "~/out_files/mental_health_sims/mid_fixed_ss_summary.rds")


# Example trial 1
set.seed(123)
sc1 <- scenarios[[43]]
sc1 <- scenarios[[1]]
sc1$looks <- seq(30, 150, 30)
ex1 <- sim_trial(sc1)
saveRDS(ex1, file = "~/out_files/mental_health_sims/example_trial1.rds")

# Example trial 2
set.seed(763753)
sc2 <- scenarios[[45]]
sc2 <- scenarios[[3]]
sc2$looks <- seq(30, 150, 30)
ex2 <- sim_trial(sc2)
saveRDS(ex2, file = "~/out_files/mental_health_sims/example_trial2.rds")

# Example trial 3
set.seed(672356)
sc3 <- scenarios[[58]]
sc3 <- scenarios[[16]]
sc3$looks <- seq(30, 150, 30)
ex3 <- sim_trial(sc3)
saveRDS(ex3, file = "~/out_files/mental_health_sims/example_trial3.rds")

# Example trial 4
set.seed(124572)
sc <- scenarios[[75]]
sc$MID <- 0
sc$looks <- seq(30, 150, 30)
ex <- sim_trial(sc)
saveRDS(ex, file = "~/out_files/mental_health_sims/example_trial4.rds")

# Example trial 5
set.seed(13546)
sc <- scenarios[[64]]
sc$MID <- 0
sc$looks <- seq(30, 150, 30)
ex <- sim_trial(sc)
saveRDS(ex, file = "~/out_files/mental_health_sims/example_trial5.rds")
