# ============================================================
# Refined resynthesis ratio template
# ------------------------------------------------------------
# Accounts for two things revealed by the real trajectory data:
#
# 1. PEAK WINDOW CONTAMINATION
#    The 90-180 min window is only a clean, uncontaminated
#    single-feed peak for the group whose refeed happens AFTER
#    it closes. In this design that's the Late group (refeed at
#    210 min). Early (refeed @ 90) and Mid (refeed @ 150) have
#    their 90-180 min reading contaminated by their own refeed,
#    so their "peak" values aren't trustworthy on their own -
#    the reference peak increment is borrowed from Late instead.
#
# 2. POST-REFEED LAG
#    The post-refeed response doesn't appear immediately either
#    (mirrors the same lag seen after the first feed) - so the
#    LATER post-refeed window (90-180 min after each subject's
#    own refeed) is used as the resynthesis value, not the
#    immediate one (0-90 min after refeed), which may understate
#    the fully-developed response.
#
# Time structure per group (mirrors the real data):
#   Early (refeed @ 90):  measured at 0, 90, 180, 270
#   Mid   (refeed @ 150): measured at 0, 90, 180, 240, 330
#   Late  (refeed @ 210): measured at 0, 90, 180, 300, 390
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# ============================================================
# 1. SIMULATE EXAMPLE DATA
# ============================================================

groups <- tibble::tibble(
  group      = c("Early", "Mid", "Late"),
  IFI        = c(90, 150, 210),
  post1_time = IFI + 90,   # first post-refeed window (0-90 min post-refeed)
  post2_time = IFI + 180   # later post-refeed window (90-180 min post-refeed)
)

n_per_group <- 7

simulate_subject <- function(group_name) {

  baseline <- rnorm(1, 0.045, 0.006)
  dip <- baseline - abs(rnorm(1, 0.010, 0.004))

  # "true" undisturbed single-feed peak increment - same underlying
  # biology for everyone, but only Late's measurement of it is clean
  true_peak_increment <- rnorm(1, 0.025, 0.006)

  # how much of the natural single-feed rise had already happened
  # by the time each group's own refeed occurred
  fraction_realised <- switch(group_name,
                               Early = runif(1, 0.55, 0.75),
                               Mid   = runif(1, 0.80, 0.95),
                               Late  = 1)

  refeed_bump <- rnorm(1, 0.010, 0.004)

  peak_window <- if (group_name == "Late") {
    baseline + true_peak_increment + rnorm(1, 0, 0.003)      # clean
  } else {
    baseline + true_peak_increment * fraction_realised +
      refeed_bump + rnorm(1, 0, 0.004)                       # contaminated
  }

  # post-refeed lag: post1 captures only part of the eventual rise
  resynthesis_increment <- rnorm(1, 0.030, 0.007)

  if (group_name == "Early") {
    # Early's refeed occurs exactly at the start of the 90-180 min
    # window, so peak_window and post1 are literally the same
    # single measurement, not two independent ones - matches the
    # real design where Early only has 4 distinct windows, not 5
    post1 <- peak_window
  } else {
    post1 <- baseline + resynthesis_increment * runif(1, 0.35, 0.55) +
      rnorm(1, 0, 0.004)
  }
  post2 <- baseline + resynthesis_increment + rnorm(1, 0, 0.004)

  tibble::tibble(baseline, dip, peak_window, post1, post2,
                 true_peak_increment)  # kept for checking the sim only -
                                       # you won't have this in real data
}

sim_list <- list()
for (g in seq_len(nrow(groups))) {
  grp <- groups$group[g]
  ifi <- groups$IFI[g]
  for (i in 1:n_per_group) {
    row <- simulate_subject(grp)
    row$subject <- paste0(grp, i)
    row$group <- grp
    row$IFI <- ifi
    sim_list[[length(sim_list) + 1]] <- row
  }
}

sim <- bind_rows(sim_list) %>%
  left_join(groups %>% select(group, post1_time, post2_time), by = "group") %>%
  mutate(group = factor(group, levels = c("Early", "Mid", "Late")))

# ============================================================
# 2. SANITY CHECK: are baselines comparable across groups?
# worth checking before borrowing a peak reference from one
# group and applying it to the others
# ============================================================

summary(aov(baseline ~ group, data = sim))
kruskal.test(baseline ~ group, data = sim)

ggplot(sim, aes(x = group, y = baseline)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs(title = "Baseline MPS by group - check for imbalance") +
  theme_minimal()

# ============================================================
# 3. REFERENCE PEAK INCREMENT
# borrowed from the one group with an uncontaminated peak window
# change reference_group if your own data means a different
# group is the clean one
# ============================================================

reference_group <- "Late"

ref_peak_increment <- sim %>%
  filter(group == reference_group) %>%
  summarise(ref = mean(peak_window - baseline)) %>%
  pull(ref)

ref_peak_increment

# ============================================================
# 4. RESYNTHESIS RATIO
# uses the LATER post-refeed window (post2) and each subject's
# own baseline - NOT each group's own (possibly contaminated)
# peak_window value
# ============================================================

sim <- sim %>%
  mutate(
    own_increment     = post2 - baseline,
    resynthesis_ratio = own_increment / ref_peak_increment
  )

sim %>%
  select(subject, group, IFI, baseline, peak_window, post1, post2,
         own_increment, resynthesis_ratio)

# ============================================================
# 5. VISUAL CHECKS
# ============================================================

ggplot(sim, aes(x = IFI, y = resynthesis_ratio)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Inter-feed interval (min)", y = "Resynthesis ratio",
       title = "Resynthesis ratio vs inter-feed interval",
       subtitle = "Dashed line = fully matched the reference peak increment") +
  theme_minimal()

traj_long <- sim %>%
  select(subject, group, IFI, baseline, dip, peak_window, post1, post2,
         post1_time, post2_time) %>%
  pivot_longer(cols = c(baseline, dip, peak_window, post1, post2),
               names_to = "window", values_to = "MPS") %>%
  mutate(time = case_when(
    window == "baseline"    ~ 0,
    window == "dip"         ~ 90,
    window == "peak_window" ~ 180,
    window == "post1"       ~ post1_time,
    window == "post2"       ~ post2_time
  ))

ref_peak_level <- ref_peak_increment + mean(sim$baseline[sim$group == reference_group])

ggplot(traj_long, aes(x = time, y = MPS, color = group, group = subject)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(aes(group = group), fun = mean, geom = "point", size = 2.5) +
  geom_hline(yintercept = ref_peak_level, linetype = "dashed", color = "grey40") +
  labs(x = "Time (min)", y = "MPS", color = "Group",
       title = "Simulated MPS trajectories with reference peak level") +
  theme_minimal()

# ============================================================
# 6. PRIMARY ANALYSIS
# ============================================================

m_ifi <- lm(resynthesis_ratio ~ IFI, data = sim)
summary(m_ifi)
confint(m_ifi)

kruskal.test(resynthesis_ratio ~ group, data = sim)
summary(aov(resynthesis_ratio ~ group, data = sim))

# sensitivity check: does controlling for baseline directly
# (rather than the ratio) change the picture?
summary(lm(post2 ~ baseline + IFI, data = sim))

# ============================================================
# To bring plasma AA back in (from the earlier AUC/rmcorr work):
# merge AA_at_refeed / post_AA_iAUC onto this `sim` data frame by
# subject, then correlate them against resynthesis_ratio exactly
# as in resynthesis_analysis_template.R
# ============================================================

