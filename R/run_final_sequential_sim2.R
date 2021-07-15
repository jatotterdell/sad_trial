# run_basic_sequential_sim.R

source("R/basic_sequential_sim_final.R")

# Data Generation setup 1
pt <- ptrunc(0:56, spec = "logis", a = -0.5, b = 56, 17.49, 2.5)
p0 <- diff(c(0, pt))
m0 <- sum(p0*0:56)
s0 <- sqrt(sum(p0*(0:56 - m0)^2))
or <- c(1/63, 1/26.5, 1/(11.5), 1/5, 1/2.25)
p1 <- pomodm(p = p0, odds.ratio = or[1])
m1 <- sum(p1*0:56)
s1 <- sqrt(sum(p1*(0:56 - m1)^2))
p2 <- pomodm(p = p0, odds.ratio = or[2])
m2 <- sum(p2*0:56)
s2 <- sqrt(sum(p2*(0:56 - m2)^2))
p3 <- pomodm(p = p0, odds.ratio = or[3])
m3 <- sum(p3*0:56)
s3 <- sqrt(sum(p3*(0:56 - m3)^2))
p4 <- pomodm(p = p0, odds.ratio = or[4])
m4 <- sum(p4*0:56)
s4 <- sqrt(sum(p4*(0:56 - m4)^2))
p5 <- pomodm(p = p0, odds.ratio = or[5])
m5 <- sum(p5*0:56)
s5 <- sqrt(sum(p5*(0:56 - m5)^2))

round(c(m1, m2, m3, m4, m5) - m0, 2)

plot(0:56, p0, type = 'b', ylim = c(0, 0.15))
lines(0:56, p1, type = 'b', col = 2)
points(0:56, p2, type = 'b', col = 3)
points(0:56, p3, type = 'b', col = 4)
points(0:56, p4, type = 'b', col = 5)
points(0:56, p5, type = 'b', col = 5)

# # Data Generation setup 2 (skewed)
# pt <- ptrunc(0:56, spec = "sn", a = -0.5, b = 56, 9.825, 10, 3.5)
# p0 <- diff(c(0, pt))
# m0 <- sum(p0*0:56)
# s0 <- sqrt(sum(p0*(0:56 - m0)^2))
# or <- c(1/55, 1/18, 1/7.3, 1/3.45, 1/1.8)
# p1 <- pomodm(p = p0, odds.ratio = or[1])
# m1 <- sum(p1*0:56)
# s1 <- sqrt(sum(p1*(0:56 - m1)^2))
# p2 <- pomodm(p = p0, odds.ratio = or[2])
# m2 <- sum(p2*0:56)
# s2 <- sqrt(sum(p2*(0:56 - m2)^2))
# p3 <- pomodm(p = p0, odds.ratio = or[3])
# m3 <- sum(p3*0:56)
# s3 <- sqrt(sum(p3*(0:56 - m3)^2))
# p4 <- pomodm(p = p0, odds.ratio = or[4])
# m4 <- sum(p4*0:56)
# s4 <- sqrt(sum(p4*(0:56 - m4)^2))
# p5 <- pomodm(p = p0, odds.ratio = or[5])
# m5 <- sum(p5*0:56)
# s5 <- sqrt(sum(p5*(0:56 - m5)^2))
#
# round(c(m1, m2, m3, m4, m5) - m0, 2)
#
# plot(0:56, p0, type = 'b', ylim = c(0, 0.15))
# lines(0:56, p1, type = 'b', col = 2)
# points(0:56, p2, type = 'b', col = 3)
# points(0:56, p3, type = 'b', col = 4)
# points(0:56, p4, type = 'b', col = 5)
# points(0:56, p5, type = 'b', col = 6)

# Outcome distribution combinations
p_alt00 <- rbind(p0, p0, p0)
p_alt01 <- rbind(p0, p0, p1)
p_alt02 <- rbind(p0, p0, p2)
p_alt03 <- rbind(p0, p0, p3)
p_alt04 <- rbind(p0, p0, p4)
p_alt05 <- rbind(p0, p0, p5)
p_alt10 <- rbind(p0, p1, p0)
p_alt11 <- rbind(p0, p1, p1)
p_alt12 <- rbind(p0, p1, p2)
p_alt13 <- rbind(p0, p1, p3)
p_alt14 <- rbind(p0, p1, p4)
p_alt15 <- rbind(p0, p1, p5)
p_alt20 <- rbind(p0, p2, p0)
p_alt21 <- rbind(p0, p2, p1)
p_alt22 <- rbind(p0, p2, p2)
p_alt23 <- rbind(p0, p2, p3)
p_alt24 <- rbind(p0, p2, p4)
p_alt25 <- rbind(p0, p2, p5)
p_alt30 <- rbind(p0, p3, p0)
p_alt31 <- rbind(p0, p3, p1)
p_alt32 <- rbind(p0, p3, p2)
p_alt33 <- rbind(p0, p3, p3)
p_alt34 <- rbind(p0, p3, p4)
p_alt35 <- rbind(p0, p3, p5)
p_alt40 <- rbind(p0, p4, p0)
p_alt41 <- rbind(p0, p4, p1)
p_alt42 <- rbind(p0, p4, p2)
p_alt43 <- rbind(p0, p4, p3)
p_alt44 <- rbind(p0, p4, p4)
p_alt45 <- rbind(p0, p4, p5)
p_alt50 <- rbind(p0, p5, p0)
p_alt51 <- rbind(p0, p5, p1)
p_alt52 <- rbind(p0, p5, p2)
p_alt53 <- rbind(p0, p5, p3)
p_alt54 <- rbind(p0, p5, p4)
p_alt55 <- rbind(p0, p5, p5)

# post_sd <- uniroot(\(x) pnorm(0, -5, x) - 0.95, interval = c(1, 10))$root
# par(mar = c(4,4,1,1), cex = 0.95)
# curve(dnorm(x, -5, post_sd), -15, 5, xlab = "Difference mean HAM-A", ylab = "Density", n = 1e3,
#       main = "Treat 1 - Usual care")
# abline(v = -5, lty = 2)
# abline(v = 0, lty = 2)
# xseq <- seq(0, 10, 0.005)
# yseq <- dnorm(xseq, -5, post_sd)
# polygon(list(x = c(0, xseq, 10), y = c(0, yseq, 0)), col = "grey80", border = NA)
# text(2.5, 0.05, "Pr(d>0|y)=0.05", cex = 0.9)
# text(-7, 0.05, "E[d|y]=-5", cex = 0.9)
#
# post_sd <- uniroot(\(x) pnorm(0, -5, x) - 0.99, interval = c(1, 10))$root
# par(mar = c(4,4,1,1), cex = 0.95)
# curve(dnorm(x, -5, post_sd), -15, 5, xlab = "Difference mean HAM-A", ylab = "Density", n = 1e3,
#       main = "Treat 1 - Usual care")
# abline(v = -5, lty = 2)
# abline(v = 0, lty = 2)
# xseq <- seq(0, 10, 0.005)
# yseq <- dnorm(xseq, -5, post_sd)
# polygon(list(x = c(0, xseq, 10), y = c(0, yseq, 0)), col = "grey80", border = NA)
# text(2.5, 0.05, "Pr(d>0|y)=0.01", cex = 0.9)
# text(-7, 0.05, "E[d|y]=-5", cex = 0.9)


# Scenario setup
# Prior - suppose we want Pr(delta < -10) approx 0.05
pr_sd <- uniroot(\(x) pnorm(-10, 0, x) - 0.05, interval = c(1, 10))$root
looks <- list(seq(30, 150, 30))
p <- list(
  p_alt00, p_alt01, p_alt02, p_alt03, p_alt04, p_alt05,
  p_alt10, p_alt11, p_alt12, p_alt13, p_alt14, p_alt15,
  p_alt20, p_alt21, p_alt22, p_alt23, p_alt24, p_alt25,
  p_alt30, p_alt31, p_alt32, p_alt33, p_alt34, p_alt35,
  p_alt40, p_alt41, p_alt42, p_alt43, p_alt44, p_alt45,
  p_alt50, p_alt51, p_alt52, p_alt53, p_alt54, p_alt55)
mu0 <- list(c(18, 0, 0))
Sigma0 <- list(diag(c(10^2, pr_sd^2, pr_sd^2)))
epsilon <- list(c(0.99, 0.05, 0.025, 0.85))
Delta <- list(c(2.5, 5, 5))
stop_rule <- list(2)
drop_rule <- list(1)

scenarios <- scenario_list(n_list(p, mu0, Sigma0, epsilon, Delta, stop_rule, drop_rule, looks))
scenarios_dt <- rbindlist(lapply(scenarios, function(x)
  as.data.table(lapply(x, function(a) ifelse(length(a)>1, list(a), I(a))))))[, scenario := 1:.N]
scenarios_dt[, effect1 := rep(c(rep(0, 6), rep(10, 6), rep(8, 6), rep(6, 6), rep(4, 6), rep(2, 6)))]
scenarios_dt[, effect2 := rep(rep(c(0, 10, 8, 6, 4, 2), 6))]

RNGkind("L'Ecuyer-CMRG")
set.seed(26246)

res_trial <- rbindlist(pbmcapply::pbmclapply(scenarios, function(tt) sim_replicate(1e4, trial_par = tt), mc.cores = 12), idcol = "scenario")
res_fin   <- res_trial[, .SD[.N], by = .(scenario, trial)]
res_fin[, analysis := as.integer(analysis)]
res_sum   <- res_fin[, !c("status", "trial")][, lapply(.SD, mean), keyby = .(scenario)]
dat_trial <- res_trial[scenarios_dt, on = .(scenario)]
dat_fin   <- dat_trial[, .SD[.N], by = .(scenario, trial)]
dat_sum   <- res_sum[scenarios_dt, on = .(scenario)]

saveRDS(scenarios_dt, file = "~/out_files/mental_health_sims/basic_sequential_scenarios_final2.rds")
saveRDS(dat_trial, file = "~/out_files/mental_health_sims/basic_sequential_sims_final2.rds")
saveRDS(dat_sum, file = "~/out_files/mental_health_sims/basic_sequential_summary_final2.rds")
