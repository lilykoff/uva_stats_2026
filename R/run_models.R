library(svyfosr)
library(tidyverse)

mims_sl = read_rds(here::here("data", "mims_covariates.rds")) |>
  mutate(age_cat = cut(age_in_years_at_screening, breaks = c(0, 17, 29, 39, 49, 59, 69, Inf),
                       labels = c("<18", "18-29", "30-39", "40-49", "50-59", "60-69", "70+")), include.lowest = TRUE)

# get pointwise models from this df

lm_wls_multi <- function(X, Y, w) {
  # X: (n × p), Y: (n × L), w: (n)
  n = nrow(X)
  if (is.null(w)) w = rep(1, n)
  W_half <- sqrt(w)
  Xw <- X * W_half
  Yw <- Y * W_half
  # Solve weighted least squares for all outcomes at once
  coef_mat <- qr.coef(qr(Xw), Yw)
  # returns (p × L) coefficient matrix
  coef_mat
}


X_des = model.matrix(~ gender + age_cat, mims_sl)
Y_mat = mims_sl |> select(starts_with("min")) |> as.matrix()
wt_vec = mims_sl$full_sample_2_year_mec_exam_weight

res = lm_wls_multi(X_des, Y_mat, wt_vec)
dim(res)

male_vec = res[2,]

write_rds(male_vec, here::here("data", "male_effect_pw.rds"))

mims_mat = mims_sl |>
  select(starts_with("min")) |>
  as.matrix()

mims_sl =
  mims_sl |>
  rename(weight = full_sample_2_year_mec_exam_weight)
model_fit = svyfosr::svyfui(mims_mat ~ gender + age_cat,
                            data = mims_sl,
                            weights = weight,
                            family = gaussian(),
                            boot_type = "BRR",
                            num_boots = 100,
                            parallel = TRUE,
                            nknots_min = 20,
                            nknots_min_fpca = 35,
                            n_cores = parallelly::availableCores() - 1,
                            seed = 2213)

model_fit$boots |> dim()

# extract boots for gender

boots = model_fit$boots[2, , ]

write_rds(boots, here::here("data", "bootstrap_coefs.rds"))


plt_df = model_fit$tidy_df |>
  filter(var_name %in% c("(Intercept)", "genderMale"))

write_rds(plt_df, here::here("data", "model_fit_plt.rds"))

plt_df |>
  ggplot(aes(x = l, y = beta_hat)) +
  facet_wrap(.~var_name, scales = "free_y") +
  geom_ribbon(aes(x = l, ymin = lower_pw, ymax = upper_pw, fill = "Pointwise"), alpha = .5) +
  geom_ribbon(aes(x = l, ymin = lower_joint, ymax = upper_joint, fill = "Joint"), alpha = 0.3) +
  geom_line(linewidth = 1.1)  +
  labs(x = "Functional Domain", y = "Coefficient Estimate")
  scale_fill_manual(values = c("Joint" = joint_fill, "Pointwise" =  pw_fill), name = "Confidence Interval")
