

scenarios_dt <- readRDS(file = "~/out_files/mental_health_sims/basic_sequential_scenarios_final2.rds")
dat_trial <- readRDS(file = "~/out_files/mental_health_sims/basic_sequential_sims_final2.rds")
dat_sum <- readRDS(file = "~/out_files/mental_health_sims/basic_sequential_summary_final2.rds")


dat_trial[, analysis := as.integer(analysis)]
dat_trial[,
          `:=`(
            cycle = ifelse(shift(stp, fill = 0, type = "lag") == 1, shift(analysis, type = "lag"), analysis),
            final = ifelse(shift(stp, fill = 0, type = "lag") == 1 | analysis == .N, 1, 0),
            result = ifelse(shift(stp, fill = 0, type = "lag") == 1, shift(status, type = "lag"),
                            ifelse(analysis == .N, "max sample size", "continue"))
          ),
          by = .(scenario, trial)
]

dat_trial_complete <- dat_trial[
  , .(analysis, trial, effect1, effect2, eff1, eff2, fut1, fut2, inf1, inf2, heff1, heff2, drp0, drp1, drp2)][
    CJ(analysis = unique(dat_trial$analysis),
       trial = unique(dat_trial$trial),
       effect1 = seq(0,10,2),
       effect2 = seq(0,10,2)),
    on = c("analysis", "trial", "effect1", "effect2")
  ]
dat_trial_complete[, `:=`(
  eff1 = nafill(eff1, type = "locf"),
  eff2 = nafill(eff2, type = "locf"),
  inf1 = nafill(inf1, type = "locf"),
  inf2 = nafill(inf2, type = "locf"),
  fut1 = nafill(fut1, type = "locf"),
  fut2 = nafill(fut2, type = "locf"),
  heff1 = nafill(heff1, type = "locf"),
  heff2 = nafill(heff2, type = "locf"),
  drp0 = nafill(drp0, type = "locf"),
  drp1 = nafill(drp1, type = "locf"),
  drp2 = nafill(drp2, type = "locf")),
  by = .(trial, effect1, effect2)]


