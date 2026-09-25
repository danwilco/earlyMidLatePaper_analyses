## Linegraphs for plasma concentration in R

# load in libraries
library(tidyverse) # for data loading and manipulation
library(tidyplots) # for publication ready plots
library(lme4) # for linear model analyses
library(garfify) # for post-hoc constrasts


# tidyplots looks the most useful for publication ready visualisations.
# load in some real data - and adjust with additions and changes

early <- read_csv("Desktop/ClaudeProjects/earlyMidLateMPS/plasmaData/plasmaLeu_Early.csv")
timecourseEarly <- tidyplot(early, x = time, y = concentration) |>
  add_mean_dot(size = 3) |>
  add_mean_line(linewidth = 1.25, alpha = 0.5) |>
  add_sem_errorbar(width = 10) |>
  theme_tidyplot() |>
  adjust_size(width = 150, height = 80) |>
  adjust_colors("#ff5c33") |>
  adjust_x_axis_title(title = "Time (mins)") |>
  adjust_y_axis_title("$Arterialised~~Plasma~~Leucine~~Concentration~~(mu*M)$") |>
  adjust_font(fontsize = 12)

mid <- read_csv("Desktop/ClaudeProjects/earlyMidLateMPS/plasmaData/plasmaLeu_Mid.csv")
timecourseMid <- tidyplot(mid, x = time, y = concentration) |>
  add_mean_dot(size = 3) |>
  add_mean_line(linewidth = 1.25, alpha = 0.5) |>
  add_sem_errorbar(width = 10) |>
  theme_tidyplot() |>
  adjust_size(width = 150, height = 80) |>
  adjust_colors("#4da6ff") |>
  adjust_x_axis_title(title = "Time (mins)") |>
  adjust_y_axis_title("$Arterialised~~Plasma~~Leucine~~Concentration~~(mu*M)$") |>
  adjust_font(fontsize = 12)

late <- read_csv("Desktop/ClaudeProjects/earlyMidLateMPS/plasmaData/plasmaLeu_Late.csv")
timecourseLate <- tidyplot(late, x = time, y = concentration) |>
  add_mean_dot(size = 3) |>
  add_mean_line(linewidth = 1.25, alpha = 0.5) |>
  add_sem_errorbar(width = 10) |>
  theme_tidyplot() |>
  adjust_size(width = 150, height = 80) |>
  adjust_colors("#000000") |>
  adjust_x_axis_title(title = "Time (mins)") |>
  adjust_y_axis_title("$Arterialised~~Plasma~~Leucine~~Concentration~~(mu*M)$") |>
  adjust_y_axis(limits = c(0,800)) |>
  adjust_font(fontsize = 12)

## quick statistical testing.

# summarise data

summary <- early |>
  group_by(time) |>
  get_summary_stats(concentration, type = "mean_sd")

sum_early <- data.frame(summary)

  
normality <- early |>
    group_by(time) |>
    shapiro_test(concentration)

norm_early <- data.frame(normality)

res_early <- anova_test(data = early,dv = concentration, wid = subject, within = time) 
get_anova_table(res_early) 


# Random intercept and random slope for time per subject
early$time <- as.factor(early$time)
model <- lmer(concentration ~ time + (1 | subject), data = early)

# View the ANOVA-style summary table
anova(model)

## now run grafify posthoc_vsRef to look at contrasts back to baseline for post-hoc analyses
posthoc_vsRef(Model = model, Fixed_Factor = "time", Ref_Level = 1)




