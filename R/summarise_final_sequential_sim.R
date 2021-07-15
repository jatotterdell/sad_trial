library(data.table)
library(ggplot2)
library(colorspace)
library(gtable)
library(grid)
library(gridExtra)

add_facet_labs <- function(p, labelT = "", labelR = "") {
  g <- ggplotGrob(p)
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(g$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(g$layout, grepl("strip-t", name), select = t:r)
  # Add a new column to the right of current right strips,
  # and a new row on top of current top strips
  if(nrow(posR) > 0)
    width <- g$widths[max(posR$r)]    # width of current right strips
  if(nrow(posT) > 0)
    height <- g$heights[min(posT$t)]  # height of current top strips
  if(nrow(posR) > 0)
    g <- gtable_add_cols(g, height, max(posR$r))
  if(nrow(posT) > 0)
    g <- gtable_add_rows(g, height, min(posT$t)-1)

  # Construct the new strip grobs
  if(nrow(posR) > 0) {
    stripR <- gTree(name = "Strip_right", children = gList(
      rectGrob(gp = gpar(col = "black", fill = "grey90")),
      textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, fontface = 'bold', col = "grey10"))))
  }
  if(nrow(posT) > 0) {
    stripT <- gTree(name = "Strip_top", children = gList(
      rectGrob(gp = gpar(col = "black", fill = "grey90")),
      textGrob(labelT, gp = gpar(fontsize = 8.8, fontface = 'bold', col = "grey10"))))
  }

  # Position the grobs in the gtable
  if(nrow(posR) > 0) {
    g <- gtable_add_grob(g, stripR, t = min(posR$t)+1,
                         l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  }
  if(nrow(posT) > 0) {
    g <- gtable_add_grob(g, stripT, t = min(posT$t),
                         l = min(posT$l), r = max(posT$r), name = "strip-top")
  }

  # Add small gaps between strips
  if(nrow(posR) > 0)
    g <- gtable_add_cols(g, unit(1/5, "line"), max(posR$r))
  if(nrow(posT) > 0)
    g <- gtable_add_rows(g, unit(1/5, "line"), min(posT$t))
  return(g)
}

scenarios_dt <- readRDS(file = "~/out_files/mental_health_sims/basic_sequential_scenarios_final2.rds")
dat_trial <- readRDS(file = "~/out_files/mental_health_sims/basic_sequential_sims_final2.rds")
dat_sum <- readRDS(file = "~/out_files/mental_health_sims/basic_sequential_summary_final2.rds")

dat_trial[, analysis := as.integer(analysis)]
# Note that, if we stop at analysis 3, say, then analysis 4 is final, but
# is not the same as an interim analysis 4.
# Therefore, define cycle to be the same as analysis, but if it is final
# following a stopping rule, then cycle = interim at which stopping occurred.
# E.g. analysis = 3 and stop = 1 => at analysis 4, cycle = 3 and final = 1.
dat_trial[,
  `:=`(
    cycle = ifelse(shift(stp, fill = 0, type = "lag") == 1, shift(analysis, type = "lag"), analysis),
    final = ifelse(shift(stp, fill = 0, type = "lag") == 1 | analysis == .N, 1, 0),
    result = ifelse(shift(stp, fill = 0, type = "lag") == 1, shift(status, type = "lag"),
                    ifelse(analysis == .N, "max sample size", "continue"))
  ),
  by = .(scenario, trial)
]

dat_fin <- dat_trial[, .SD[.N], by = .(scenario, trial)]
dat_cycle <- dat_trial[, .SD[.N], by = .(scenario, trial, cycle)]

stopping_time <- dat_trial[, .N, keyby = .(cycle, effect1, effect2, result)]
stopping_time[, P := N / sum(N), by = .(cycle, effect1, effect2)]

ggplot(stopping_time,
       aes(cycle, P, fill = result)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(effect1 ~ effect2)

ggplot(stopping_time,
       aes(cycle, P, fill = result)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(effect1 ~ effect2)

# Note that some trial finish early, meaning that a straight summary by
# analysis is implicitly conditional on trials which make it to that analysis
# The alternative is to just carry the results forward,
# e.g. if a trials stops following analysis 3, then the results at analysis 4 and 5 are the same
# that is, if an assertion was made at analysis 3 and the trial stopped, then that assertion had also been
# made by analysis 4 and 5.
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

final_cycle <- dat_trial[final == 1, .N, keyby = .(cycle, effect1, effect2)]
final_cycle <- final_cycle[CJ(cycle = 1:5, effect1 = seq(0, 10, 2), effect2 = seq(0, 10, 2))]
final_cycle[is.na(N), N := 0]
final_cycle[, P := N/sum(N), keyby = .(effect1, effect2)]

prob_any_sup <- dat_trial[, .(mean(eff1 | eff2)), by = .(analysis, effect1, effect2)]
dcast(prob_any_sup, analysis + effect1 ~ effect2, value.var = "V1")

prob_any_sup <- dat_trial[, .(mean(eff1 | eff2)), by = .(cycle, effect1, effect2)]
dcast(prob_any_sup, cycle + effect1 ~ effect2, value.var = "V1")

# Stopping

dat_stp     <- dat_trial[stp == 1 & analysis < 5, .SD[1], by = .(scenario, trial)]
stp_pattern <- dat_stp[, .N, keyby = .(effect1, effect2, eff1, eff2, inf1, inf2, fut1, fut2, heff1, heff2)]
stp_pattern[, `:=`(
  both_fut = fut1 & fut2,
  both_heff = heff1 & heff2,
  one_eff_one_inf = (eff1 & inf2) | (eff2 & inf1)
)]
stp_pattern[, P := N / max(dat_trial$trial), keyby = .(effect1, effect2)]
stp_pattern_grp <- copy(stp_pattern)
stp_pattern_grp <- stp_pattern_grp[, .(N = sum(N)), keyby = .(effect1, effect2, both_fut, both_heff, one_eff_one_inf)]
stp_pattern_grp[, P := N / max(dat_trial$trial), keyby = .(effect1, effect2)]

stp_pattern[effect1 == 0 & effect2 == 0][rev(order(P))]
stp_pattern[effect1 == 2 & effect2 == 2][rev(order(P))]
stp_pattern[effect1 == 10 & effect2 == 10][rev(order(P))]

stp_pattern_grp[effect1 == 0 & effect2 == 0][rev(order(P))]
stp_pattern_grp[effect1 == 2 & effect2 == 2][rev(order(P))]
stp_pattern_grp[effect1 == 10 & effect2 == 10][rev(order(P))]

dat_nostp <- dat_trial[stp == 0 & analysis == 5, .SD[1], by = .(scenario, trial)]
nostp_pattern <- dat_nostp[, .N, keyby = .(effect1, effect2, eff1, eff2, inf1, inf2, fut1, fut2, heff1, heff2)]
nostp_pattern[, P := N / max(dat_trial$trial), keyby = .(effect1, effect2)]

nostp_pattern[effect1 == 0 & effect2 == 0][rev(order(P))]
nostp_pattern[effect1 == 2 & effect2 == 2][rev(order(P))]
nostp_pattern[effect1 == 10 & effect2 == 10][rev(order(P))]

ggplot(final_cycle,
       aes(cycle, P)) +
  geom_bar(stat = "identity") +
  facet_grid(effect1 ~ effect2)



#
# SUMMARISE THE PROBABILITY OF CRITERION
# BEING SATISFIED AT SOME POINT DURING THE TRIAL
# ================================================

# Probability of having ASSERTED for treatment 1
prob_trt1 <- dat_trial_complete[, .(
  `Effective\nP(d1 < 0) > 0.975` = mean(eff1),
  `Futile\nP(d1 < -2.5) < 0.05` = mean(fut1),
  `Inferior\nP(d1-d2 < 0) < 0.025` = mean(inf1),
  `Highly effective\nP(d1 < -5) > 0.85` = mean(heff1)),
  by = .(Analysis = analysis, effect1, effect2)]
prob_trt1 <- melt(prob_trt1, id.vars = c("Analysis", "effect1", "effect2"))
p_trt1 <- ggplot(prob_trt1,
       aes(effect1, effect2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 2) +
  facet_grid(variable ~ Analysis) +
  scale_fill_continuous_sequential("Green-Yellow", n_interp = 50) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(x = "Effect size treatment 1", y = "Effect size treatment 2",
       title = "Probability that criterion met for treatment 1 by analysis.") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
p_trt1 <- add_facet_labs(p_trt1, "Analysis", "Assertion")
grid.arrange(p_trt1)

# Probability of having ASSERTED for treatment 2
prob_trt2 <- dat_trial_complete[, .(
  `Effective\nP(d1 < 0) > 0.975` = mean(eff2),
  `Futile\nP(d1 < -2.5) < 0.05` = mean(fut2),
  `Inferior\nP(d1-d2 < 0) < 0.025` = mean(inf2),
  `Highly effective\nP(d1 < -5) > 0.85` = mean(heff2)),
  by = .(Analysis = analysis, effect1, effect2)]
prob_trt2 <- melt(prob_trt2, id.vars = c("Analysis", "effect1", "effect2"))
p_trt2 <- ggplot(prob_trt2,
                 aes(effect1, effect2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 2) +
  facet_grid(variable ~ Analysis) +
  scale_fill_continuous_sequential("Green-Yellow", n_interp = 50) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(x = "Effect size treatment 1", y = "Effect size treatment 2",
       title = "Probability that criterion met for treatment 2 by analysis.") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
p_trt2 <- add_facet_labs(p_trt2, "Analysis", "Assertion")
grid.arrange(p_trt2)

# Probability of assertion for both treatments
prob_both <- dat_trial_complete[, .(
  `Effective\nP(d1 < 0) > 0.975` = mean(eff1 & eff2),
  `Futile\nP(d1 < -2.5) < 0.05` = mean(fut1 & fut2),
  `Inferior\nP(d1-d2 < 0) < 0.025` = mean(inf1 & inf2),
  `Highly effective\nP(d1 < -5) > 0.85` = mean(heff1 & heff2)),
  by = .(Analysis = analysis, effect1, effect2)]
prob_both <- melt(prob_both, id.vars = c("Analysis", "effect1", "effect2"))
p_both <- ggplot(
  # prob_both[variable != "Inferior"],
  prob_both,
  aes(effect1, effect2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 2) +
  facet_grid(variable ~ Analysis) +
  scale_fill_continuous_sequential("Green-Yellow", n_interp = 50) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(x = "Effect size treatment 1", y = "Effect size treatment 2",
       title = "Probability that criterion met for both treatments by analysis.") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
p_both <- add_facet_labs(p_both, "Analysis", "Assertion")
grid.arrange(p_both)

# Probability of assertion for either treatments
prob_any <- dat_trial_complete[, .(
  `Effective\nP(d1 < 0) > 0.975` = mean(eff1 | eff2),
  `Futile\nP(d1 < -2.5) < 0.05` = mean(fut1 | fut2),
  `Inferior\nP(d1-d2 < 0) < 0.025` = mean(inf1 | inf2),
  `Highly effective\nP(d1 < -5) > 0.85` = mean(heff1 | heff2)),
  by = .(Analysis = analysis, effect1, effect2)]
prob_any <- melt(prob_any, id.vars = c("Analysis", "effect1", "effect2"))
p_any <- ggplot(prob_any,
                 aes(effect1, effect2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 2) +
  facet_grid(variable ~ Analysis) +
  scale_fill_continuous_sequential("Green-Yellow", n_interp = 50) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  coord_fixed() +
  theme(legend.position = "none") +
  labs(x = "Effect size treatment 1", y = "Effect size treatment 2",
       title = "Probability that criterion met for either treatment by analysis.") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"))
p_any <- add_facet_labs(p_any, "Analysis", "Assertion")
grid.arrange(p_any)

gg <- arrangeGrob(p_trt1, p_trt2, p_any, p_both, ncol = 2)
grid.arrange(gg)
ggsave("summary_plot2.pdf", grid.arrange(gg), dpi = 240, width = 420, height = 350, units = "mm")
