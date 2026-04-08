# Two distributions: same mean (0), variance (1), skewness (0), different kurtosis
# 1. Standard normal: excess kurtosis = 0 (kurtosis = 3)
# 2. Scaled t-distribution with df=5: excess kurtosis = 6 (kurtosis = 9)
#    Var(t_df) = df/(df-2), so scale by sqrt((df-2)/df) to get variance = 1

library(ggplot2)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork", repos = "https://cloud.r-project.org")
library(patchwork)

set.seed(42)
n <- 50000

# --- Generate samples ---
x_normal <- rnorm(n, mean = 0, sd = 1)

df_t <- 5
x_t <- rt(n, df = df_t) * sqrt((df_t - 2) / df_t)  # scale to variance = 1

# Verify moments
cat("Normal:  mean =", round(mean(x_normal), 3),
    " var =", round(var(x_normal), 3),
    " skew =", round(mean((x_normal - mean(x_normal))^3) / sd(x_normal)^3, 3),
    " kurtosis =", round(mean((x_normal - mean(x_normal))^4) / sd(x_normal)^4, 3), "\n")

cat("Scaled t(5): mean =", round(mean(x_t), 3),
    " var =", round(var(x_t), 3),
    " skew =", round(mean((x_t - mean(x_t))^3) / sd(x_t)^3, 3),
    " kurtosis =", round(mean((x_t - mean(x_t))^4) / sd(x_t)^4, 3), "\n")

dat <- data.frame(
  value = c(x_normal, x_t),
  dist  = rep(c("Normal (kurtosis = 3)", "Scaled t5 (kurtosis = 9)"), each = n)
)

# --- Plot 1: Overlaid densities (full range) ---
p1 <- ggplot(dat, aes(x = value, colour = dist)) +
  geom_density(linewidth = 0.8) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(title = "Same mean, variance, skewness — different kurtosis",
       subtitle = "Full density",
       x = "Value", y = "Density", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# --- Plot 2: Zoom into tails (|x| > 2) ---
p2 <- ggplot(dat, aes(x = value, colour = dist)) +
  geom_density(linewidth = 0.8) +
  coord_cartesian(xlim = c(2.5, 7), ylim = c(0, 0.03)) +
  labs(title = "Right tail zoom (left tail is symmetric)",
       subtitle = "The t5 has visibly more mass beyond 3",
       x = "Value", y = "Density", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# --- Plot 3: Sorted values to show where they differ ---
# Take 200 draws from each, sort them, and plot
set.seed(7)
n_small <- 200
s_normal <- sort(rnorm(n_small))
s_t      <- sort(rt(n_small, df = df_t) * sqrt((df_t - 2) / df_t))

dat_sorted <- data.frame(
  rank = rep(1:n_small, 2),
  value = c(s_normal, s_t),
  dist  = rep(c("Normal (kurtosis = 3)", "Scaled t5 (kurtosis = 9)"), each = n_small)
)

p3 <- ggplot(dat_sorted, aes(x = rank, y = value, colour = dist)) +
  geom_point(size = 0.8, alpha = 0.7) +
  labs(title = "200 sorted draws from each distribution",
       subtitle = "Same in the middle, different at the extremes",
       x = "Rank", y = "Value", colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

combined <- p1 / p2 / p3
ggsave("output/figures/tmp_kurtosis_examples.png", combined,
       width = 8, height = 14, dpi = 150)
cat("\nSaved to output/figures/tmp_kurtosis_examples.png\n")
