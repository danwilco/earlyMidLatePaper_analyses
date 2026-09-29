# ============================================================
# Refeed-timing study: resynthesis index analysis
# One pre-refeed MPS value + one post-refeed MPS value per subject
# Groups defined by inter-feed interval (IFI): 90, 150, 210 min
# ============================================================

library(dplyr)
library(ggplot2)

# ============================================================
# 1. AUC HELPER FUNCTIONS (trapezoidal rule, self-contained)
# same logic as calculate_auc.R - included here so this script
# runs on its own
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
# 2. SIMULATED EXAMPLE DATA
# ------------------------------------------------------------
# Replace with your real data. You'll likely build this from
# two sources:
#   (a) a per-subject summary table: subject, IFI, pre_MPS, post_MPS
#   (b) a long-format AA time series: subject, IFI, time, conc
#       from which you compute pre/post AA AUC using window_auc()
# ============================================================

set.seed(42)
n_per_group <- 6
IFI_levels <- c(90, 150, 210)
fixed_followup <- 150   # common post-refeed observation length (min)

# ---- (a) per-subject MPS summary ----
mps_df <- expand.grid(subj_in_group = 1:n_per_group, IFI = IFI_levels) %>%
  mutate(
    subject = paste0("S", IFI, "_", subj_in_group),
    baseline_MPS = round(rnorm(n(), 0.035, 0.004), 4),
    # pre-refeed MPS decays toward baseline the longer you wait
    pre_MPS = round(baseline_MPS + 0.05 * exp(-IFI / 80) + rnorm(n(), 0, 0.004), 4),
    # post-refeed resynthesis recovers more fully the longer the wait
    post_MPS = round(baseline_MPS + 0.06 * (1 - exp(-IFI / 100)) + rnorm(n(), 0, 0.005), 4)
  )

# resynthesis index: increase attributable to the refeed itself
# (use pre_MPS as the reference point - "how much did MPS jump
# from wherever it had already declined to")
mps_df <- mps_df %>%
  mutate(resynthesis_index = post_MPS - pre_MPS,
         # alternative definition if your question is instead
         # "did MPS get restimulated back to a fully-fed level":
         resynthesis_vs_baseline = post_MPS - baseline_MPS) # use this one!

# ---- (b) simulated AA time series per subject ----
aa_list <- lapply(seq_len(nrow(mps_df)), function(i) {
  ifi <- mps_df$IFI[i]
  end_time <- ifi + fixed_followup
  t <- seq(0, end_time, by = 30)
  
  # simple two-bolus kinetic shape: rise/decay after t=0, then
  # a second rise/decay after the refeed at t = ifi
  base <- 1800
  peak1 <- 3200
  peak2 <- 3000
  conc <- base +
    (peak1 - base) * exp(-((t - 45)^2) / (2 * 40^2)) +
    (peak2 - base) * exp(-((t - (ifi + 45))^2) / (2 * 40^2)) +
    rnorm(length(t), 0, 40)
  
  data.frame(subject = mps_df$subject[i], IFI = ifi, time = t, conc = round(conc, 1))
})
aa_df <- bind_rows(aa_list)

# ---- compute pre- and post-refeed AA AUC per subject ----
auc_summary <- aa_df %>%
  group_by(subject, IFI) %>%
  group_modify(~ {
    baseline_val <- .x$conc[.x$time == 0][1]
    pre  <- window_auc(.x$time, .x$conc, 0, unique(.y$IFI), baseline = baseline_val)
    post <- window_auc(.x$time, .x$conc, unique(.y$IFI), unique(.y$IFI) + fixed_followup,
                       baseline = baseline_val)
    tibble::tibble(
      pre_AA_AUC  = pre["total_auc"],
      pre_AA_iAUC = pre["iauc"],
      post_AA_AUC = post["total_auc"],
      post_AA_iAUC = post["iauc"],
      AA_at_refeed = approx(.x$time, .x$conc, xout = unique(.y$IFI))$y
    )
  }) %>%
  ungroup()

# ---- merge MPS and AA summaries into one analysis-ready table ----
df <- left_join(mps_df, auc_summary, by = c("subject", "IFI"))

# ============================================================
# 3. VISUAL CHECK FIRST - look before you model
# ============================================================

ggplot(df, aes(x = IFI, y = resynthesis_index)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Inter-feed interval (min)", y = "Resynthesis index",
       title = "Resynthesis index vs inter-feed interval") +
  theme_minimal()
# check this plot for a plateau shape before committing to a
# straight-line model - a curved pattern would suggest trying
# a log(IFI) or asymptotic term instead of raw IFI

# ============================================================
# 4. PRIMARY ANALYSIS: does resynthesis scale with IFI?
# ============================================================

# ---- continuous IFI as the primary test (more power than 3-group ANOVA) ----
m_ifi <- lm(resynthesis_index ~ IFI, data = df)
summary(m_ifi)
confint(m_ifi)

# ---- categorical sensitivity check, in case the relationship isn't linear ----
df$IFI_factor <- factor(df$IFI)
summary(aov(resynthesis_index ~ IFI_factor, data = df))
kruskal.test(resynthesis_index ~ IFI_factor, data = df)   # non-parametric alternative

# ============================================================
# 5. LINKING PLASMA AA TO THE RESYNTHESIS INDEX
# ============================================================

# ---- (i) did AA clear back toward baseline before the refeed? ----
cor.test(df$AA_at_refeed, df$resynthesis_index, method = "pearson")
cor.test(df$AA_at_refeed, df$resynthesis_index, method = "spearman")

# ---- (ii) does post-refeed AA exposure track the post-refeed MPS response? ----
cor.test(df$post_AA_iAUC, df$resynthesis_index, method = "pearson")
cor.test(df$post_AA_iAUC, df$resynthesis_index, method = "spearman")

# ---- (optional) does AA explain resynthesis beyond IFI alone? ----
# with n=15-24 total, this 2-predictor model is at the edge of what's
# reasonable - treat as exploratory, not confirmatory
m_combined <- lm(resynthesis_index ~ IFI + AA_at_refeed, data = df)
summary(m_combined)

# ============================================================
# 6. OPTIONAL EXTENSION: if you split pre-refeed MPS into
# multiple sub-windows instead of one summary value, use this
# block to look at the pre-refeed decline trajectory and its
# association with AA decline, within-subject
# ============================================================

# library(rmcorr)
# library(lmerTest)
#
# Expected long-format columns: subject, IFI, window, MPS_window, AA_window_AUC
#
# rmcorr(participant = subject, measure1 = AA_window_AUC, measure2 = MPS_window,
#        dataset = pre_refeed_long_df)
#
# lmer(MPS_window ~ AA_window_AUC + IFI + (1 | subject), data = pre_refeed_long_df)