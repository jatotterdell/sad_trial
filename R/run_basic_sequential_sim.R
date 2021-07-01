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


p <- rbind(p0, p0, p0)
p_alt11 <- rbind(p0, p1, p1)
p_alt22 <- rbind(p0, p2, p2)
p_alt33 <- rbind(p0, p3, p3)
p_alt44 <- rbind(p0, p4, p4)
p_alt21 <- rbind(p0, p2, p1)
p_alt31 <- rbind(p0, p3, p1)
p_alt41 <- rbind(p0, p4, p1)

looks <- list(c(30, 60, 90, 120, 150))
p <- list(p, p_alt11, p_alt22, p_alt33, p_alt44, p_alt21, p_alt31, p_alt41)
epsilon <- list(0.85, 0.9, 0.95)
mu0 <- list(c(20, 0, 0))
Sigma0 <- list(diag(10^2, 3))
# drop_rule <- list("one", "two", "joint")
drop_rule <- list("two")
scenarios <- scenario_list(n_list(looks, p, epsilon, mu0, Sigma0, drop_rule))
scenarios_dt <- rbindlist(lapply(scenarios, function(x)
  as.data.table(lapply(x, function(a) ifelse(length(a)>1, list(a), I(a))))))[, scenario := 1:.N]
scenarios_dt[, effect := factor(rep(c("Null", "Very Large", "Large", "Moderate", "Small"), 3),
                                levels = c("Null", "Small", "Moderate", "Large", "Very Large"))]

res <- rbindlist(lapply(scenarios, function(x) sim_replicate(2000, trial_par = x)), idcol = "scenario")
res_sum <- res[, lapply(.SD, mean), keyby = .(scenario, analysis)]

saveRDS(scenarios_dt, file = "~/out_files/mental_health_sims/basic_sequential_scenarios.rds")
saveRDS(res, file = "~/out_files/mental_health_sims/basic_sequential_sims.rds")
saveRDS(res_sum, file = "~/out_files/mental_health_sims/basic_sequential_summary.rds")
