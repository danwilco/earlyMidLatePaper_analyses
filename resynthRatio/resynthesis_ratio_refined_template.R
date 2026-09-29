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
library(emmeans)   # pairwise post-hoc contrasts
library(FSA)        # Dunn's test (non-parametric post-hoc)
library(ggsignif)   # significance brackets on plots

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

ref_peak_n <- sim %>%
  filter(group == reference_group) %>%
  summarise(n_complete = sum(!is.na(peak_window) & !is.na(baseline)),
            n_total = n())
ref_peak_n   # check how many Late subjects this reference is actually based on

ref_peak_increment <- sim %>%
  filter(group == reference_group) %>%
  summarise(ref = mean(peak_window - baseline, na.rm = TRUE)) %>%
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
# 7. PLASMA AA - AUC HELPER FUNCTIONS
# same trapezoidal logic as calculate_auc.R, included here so
# this script runs standalone
# ============================================================

trap_auc <- function(t1, t2, c1, c2) (c1 + c2) / 2 * (t2 - t1)

window_auc <- function(time, conc, t_start, t_end, baseline = NULL) {
  keep <- time >= t_start & time <= t_end
  t <- time[keep]; c <- conc[keep]
  ord <- order(t); t <- t[ord]; c <- c[ord]

  if (length(t) == 0 || min(t) > t_start) {
    c_start <- approx(time, conc, xout = t_start)$y
    t <- c(t_start, t); c <- c(c_start, c)
  }
  if (max(t) < t_end) {
    c_end <- approx(time, conc, xout = t_end)$y
    t <- c(t, t_end); c <- c(c, c_end)
  }

  total_auc <- sum(sapply(seq_len(length(t) - 1), function(i)
    trap_auc(t[i], t[i + 1], c[i], c[i + 1])))

  if (!is.null(baseline)) {
    iauc <- total_auc - baseline * (t_end - t_start)
    return(c(total_auc = total_auc, iauc = iauc))
  }
  c(total_auc = total_auc)
}

# ============================================================
# 8. SIMULATE PLASMA AA TIME SERIES
# one curve per subject (every 30 min), with a two-bolus shape:
# rise/decay after the initial feed, and a second rise/decay
# after that subject's own refeed
# ============================================================

simulate_aa_series <- function(ifi, end_time) {
  t <- seq(0, end_time, by = 30)
  base <- 1800
  peak1 <- 3200
  peak2 <- 3000

  conc <- base +
    (peak1 - base) * exp(-((t - 60)^2) / (2 * 45^2)) +
    (peak2 - base) * exp(-((t - (ifi + 60))^2) / (2 * 45^2)) +
    rnorm(length(t), 0, 40)

  tibble::tibble(time = t, conc = round(conc, 1))
}

aa_list <- lapply(seq_len(nrow(sim)), function(i) {
  end_t <- sim$post2_time[i]
  aa <- simulate_aa_series(sim$IFI[i], end_t)
  aa$subject <- sim$subject[i]
  aa
})
aa_df <- bind_rows(aa_list)

# ============================================================
# 9. AA WINDOW METRICS PER SUBJECT
# mirrors the same windows used for MPS: dip, peak_window,
# AA-at-refeed, post1, post2
# ============================================================

aa_summary <- aa_df %>%
  left_join(sim %>% select(subject, IFI, post1_time, post2_time), by = "subject") %>%
  group_by(subject) %>%
  group_modify(~ {
    ifi <- unique(.x$IFI)
    p1t <- unique(.x$post1_time)
    p2t <- unique(.x$post2_time)
    baseline_val <- .x$conc[.x$time == 0][1]

    dip  <- window_auc(.x$time, .x$conc, 0, 90, baseline = baseline_val)
    peak <- window_auc(.x$time, .x$conc, 90, 180, baseline = baseline_val)
    p1   <- window_auc(.x$time, .x$conc, ifi, p1t, baseline = baseline_val)
    p2   <- window_auc(.x$time, .x$conc, p1t, p2t, baseline = baseline_val)

    tibble::tibble(
      AA_baseline   = baseline_val,
      AA_at_refeed  = approx(.x$time, .x$conc, xout = ifi)$y,
      AA_dip_iAUC   = dip["iauc"],
      AA_peak_iAUC  = peak["iauc"],
      AA_post1_iAUC = p1["iauc"],
      AA_post2_iAUC = p2["iauc"]
    )
  }) %>%
  ungroup()

# note: for the Early group, AA_peak_iAUC and AA_post1_iAUC use the
# exact same window bounds (90-180 min), so they'll come out
# identical - that's expected, mirroring the same overlap you
# confirmed for the MPS peak_window/post1 columns. Unlike the MPS
# simulation, no manual fix was needed here, since both are just
# the same integral computed over the same bounds from one
# underlying AA curve.

# ============================================================
# 10. MERGE WITH MPS DATA AND LINK TO THE RESYNTHESIS RATIO
# ============================================================

sim <- sim %>%
  left_join(aa_summary, by = "subject")

# (i) did AA clear back toward baseline before the refeed?
cor.test(sim$AA_at_refeed, sim$resynthesis_ratio, method = "pearson")
cor.test(sim$AA_at_refeed, sim$resynthesis_ratio, method = "spearman")

# (ii) does post-refeed AA exposure track the resynthesis ratio?
# (post2 window matches the MPS post2 window used in the ratio)
cor.test(sim$AA_post2_iAUC, sim$resynthesis_ratio, method = "pearson")
cor.test(sim$AA_post2_iAUC, sim$resynthesis_ratio, method = "spearman")

# (optional) does AA explain resynthesis beyond IFI alone?
# exploratory only, given n=15-24 total across groups
summary(lm(resynthesis_ratio ~ IFI + AA_at_refeed, data = sim))

ggplot(sim, aes(x = AA_at_refeed, y = resynthesis_ratio, color = group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40") +
  labs(x = "Plasma AA at refeed", y = "Resynthesis ratio",
       title = "Does AA clearance before refeed predict resynthesis?") +
  theme_minimal()

# ============================================================
# 11. BASELINE-ADJUSTED SENSITIVITY ANALYSIS
# ------------------------------------------------------------
# resynthesis_ratio already corrects for each subject's OWN
# baseline, but if baseline differs SYSTEMATICALLY by group
# (e.g. a significant kruskal.test(baseline ~ group)), that's
# still worth checking directly - both because baseline may
# affect how much room there is to increase (a ceiling/floor
# effect), and because it raises a broader question about
# whether the groups are comparable on other unmeasured factors.
# ============================================================

# does baseline itself predict the resynthesis ratio, regardless
# of group? if not, the baseline imbalance is less likely to be
# quietly driving the group comparison
cor.test(sim$baseline, sim$resynthesis_ratio, method = "spearman")

# ANCOVA-style check: does IFI still predict resynthesis_ratio
# once baseline is explicitly controlled for?
m_adjusted <- lm(resynthesis_ratio ~ IFI + baseline, data = sim)
summary(m_adjusted)
confint(m_adjusted)

# compare to the unadjusted model (m_ifi, from Section 6) - if the
# IFI estimate barely moves once baseline is added, baseline isn't
# doing much confounding; if it shifts a lot, baseline was masking
# or inflating the IFI effect
summary(m_ifi)

# visual check: does baseline relate to the resynthesis ratio in a
# way that could plausibly be driving the group differences?
ggplot(sim, aes(x = baseline, y = resynthesis_ratio, color = group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "grey40") +
  labs(x = "Baseline MPS", y = "Resynthesis ratio",
       title = "Does baseline relate to the resynthesis ratio?") +
  theme_minimal()

# ------------------------------------------------------------
# optional: a proportional (fold-change) version of the ratio,
# in case the relationship between baseline and increment is
# multiplicative rather than a fixed additive amount - useful as
# a robustness check alongside the main (additive) ratio
# ------------------------------------------------------------

ref_fold_increase <- sim %>%
  filter(group == reference_group) %>%
  summarise(ref = mean(peak_window / baseline, na.rm = TRUE)) %>%
  pull(ref)

sim <- sim %>%
  mutate(resynthesis_ratio_fold = (post2 / baseline) / ref_fold_increase)

cor.test(sim$resynthesis_ratio, sim$resynthesis_ratio_fold, method = "spearman")
summary(lm(resynthesis_ratio_fold ~ IFI, data = sim))

# ============================================================
# 12. MISSING DATA CHECK
# ------------------------------------------------------------
# Don't impute or interpolate missing MPS window values - with
# only 4-5 measured points per subject and a trajectory that
# dips and rises non-linearly, a missing point can't be reliably
# reconstructed. Leave genuine gaps as NA and let each specific
# test use whichever subjects have what it needs. This just
# shows you where the gaps are so you know what's being excluded
# from what, rather than losing track of it silently.
# ============================================================

sim %>%
  group_by(group) %>%
  summarise(across(c(baseline, dip, peak_window, post1, post2),
                    ~ sum(is.na(.)), .names = "missing_{.col}"))

# how many subjects actually have a usable resynthesis_ratio
sim %>%
  group_by(group) %>%
  summarise(n_total = n(),
            n_with_ratio = sum(!is.na(resynthesis_ratio)))

# ============================================================
# 13. POST-HOC PAIRWISE COMPARISONS
# ------------------------------------------------------------
# Using resynthesis_ratio_fold, since that's the version that
# resolved the baseline confound. The continuous IFI regression
# (Section 6) remains your primary test for the overall trend;
# these pairwise contrasts answer the follow-up question of
# which specific groups differ from each other.
# ============================================================

m_group_fold <- lm(resynthesis_ratio_fold ~ group, data = sim)
anova(m_group_fold)

# quick residual check before trusting a parametric post-hoc test
qqnorm(resid(m_group_fold)); qqline(resid(m_group_fold))

# parametric pairwise contrasts (estimated marginal means),
# Tukey-adjusted for the 3 comparisons
emm <- emmeans(m_group_fold, pairwise ~ group, adjust = "tukey")
emm$emmeans     # estimated marginal means + CIs per group
emm$contrasts   # pairwise differences + adjusted p-values

# non-parametric alternative - worth checking these two approaches
# agree, given the small n
dunnTest(resynthesis_ratio_fold ~ group, data = sim, method = "bh")

# ============================================================
# 14. PUBLICATION-QUALITY FIGURES
# ------------------------------------------------------------
# The p-value annotations below are placeholders - update them
# with your own adjusted p-values from emm$contrasts / dunnTest()
# above once you've run them, since exact contrast label order
# can vary depending on factor level coding.
# ============================================================

# ---- Figure: resynthesis ratio by group, with post-hoc brackets ----

p_group <- ggplot(sim, aes(x = group, y = resynthesis_ratio_fold, fill = group)) +
  geom_boxplot(width = 0.5, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2.2, shape = 21, color = "black") +
  geom_signif(comparisons = list(c("Early", "Mid"),
                                  c("Mid", "Late"),
                                  c("Early", "Late")),
              annotations = c("p = ...", "p = ...", "p = ..."),  # fill in
              step_increase = 0.12,
              tip_length = 0.02) +
  scale_fill_viridis_d() +
  labs(x = "Refeed timing group",
       y = "Resynthesis ratio (fold-change corrected)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

p_group
ggsave("resynthesis_ratio_by_group.png", p_group, width = 5, height = 4.5, dpi = 300)

# ---- Figure: dose-response across IFI as a continuous variable ----

p_doseresponse <- ggplot(sim, aes(x = IFI, y = resynthesis_ratio_fold)) +
  geom_smooth(method = "lm", color = "black", fill = "grey70", se = TRUE) +
  geom_point(aes(color = group), size = 2.5) +
  scale_color_viridis_d() +
  labs(x = "Inter-feed interval (min)",
       y = "Resynthesis ratio (fold-change corrected)",
       color = "Group") +
  theme_classic(base_size = 13)

p_doseresponse
ggsave("resynthesis_ratio_doseresponse.png", p_doseresponse,
       width = 5.5, height = 4.5, dpi = 300)

# ---- Figure: mean +/- SEM MPS trajectories (reuses traj_long from Section 5) ----

traj_summary <- traj_long %>%
  group_by(group, time) %>%
  summarise(mean_MPS = mean(MPS, na.rm = TRUE),
            sem_MPS = sd(MPS, na.rm = TRUE) / sqrt(sum(!is.na(MPS))),
            .groups = "drop")

p_traj <- ggplot(traj_summary, aes(x = time, y = mean_MPS, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_MPS - sem_MPS, ymax = mean_MPS + sem_MPS),
                width = 8) +
  scale_color_viridis_d() +
  labs(x = "Time (min)", y = "MPS", color = "Group") +
  theme_classic(base_size = 13)

p_traj
ggsave("mps_trajectories_mean_sem.png", p_traj, width = 6, height = 4.5, dpi = 300)

# ------------------------------------------------------------
# optional: combine panels into one multi-panel figure using
# the patchwork package, e.g.:
#   library(patchwork)
#   (p_traj | p_group) / p_doseresponse
# ------------------------------------------------------------
