library(tidyverse)
I = 10e5 # superpopulation size
num_strata = 30 # number of strata
min_psu = 75
max_psu = 125
seed = 2213
set.seed(seed)
# generate dirichlet probabilities for stratum assignments, using concentration param 4
dirichlet_probs = gtools::rdirichlet(1, rep(4, num_strata))
set.seed(seed)
# generate stratum assignments
stratum_assignments = sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs) #
psu_assignments  = rep(NA, I)

# loop thru strat and assign individuals to PSUs
# between 75 and 125 psus per stratum
for (s in 1:num_strata) {
  set.seed(seed + s)
  num_in_strata = sum(stratum_assignments == s)
  num_psu = round(runif(1, min_psu, max_psu), 0)
  set.seed(seed + s)
  dps = gtools::rdirichlet(1, rep(10, num_psu))
  set.seed(seed + s)
  psu_in_stratum = sample(1:num_psu,
                          num_in_strata,
                          replace = TRUE,
                          prob = dps)
  psu_assignments[stratum_assignments == s]  = paste0(s, "_", psu_in_stratum)
}

# generate X
set.seed(seed)
X_des = cbind(1, rnorm(I, 0, 2))

# generate global intercept and slope functions
L = 50 # length of functional domain
grid  = seq(0, 1, length = L)
beta_fixed  = matrix(NA, 2, L)

beta_fixed[1, ]  = -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
beta_fixed[2, ]  = dnorm(grid, 0.6, 0.15) / 20

rownames(beta_fixed)  = c("Intercept", "x")


btrue =
  beta_fixed %>%
  as_tibble() %>%
  mutate(beta = c("beta_0", "beta_1")) %>%
  pivot_longer(cols = -beta) %>%
  mutate(x = as.numeric(sub(".*V", "", name)),
         type = "betaTrue")

btrue_df =
  btrue %>%
  mutate(name = if_else(beta == "beta_0", "(Intercept)", "X")) %>%
  select(-beta) %>%
  rename(s = x, beta = value)

write_rds(btrue_df, here::here("data", "btrue_df_sim.rds"))

set.seed(seed)
strata_scale = 0.125
stratum_scaling  = rnorm(num_strata, mean = 1, sd = strata_scale)

beta1_by_stratum  = matrix(rep(stratum_scaling, each = L), nrow = num_strata, byrow = TRUE) * matrix(rep(beta_fixed[2, ], times = num_strata),
                                                                                                     nrow = num_strata,
                                                                                                     byrow = TRUE) # num_strata x L matrix
b1_strat_df =
  b1_strat %>%
  as_tibble() %>%
  mutate(strata = row_number(),
         scale = stratum_scaling)
# assign to individuals
beta1_by_indiv  = beta1_by_stratum[stratum_assignments, ]

write_rds(b1_strat_df, here::here("data", "b1_strat.rds"))

strata_sigma = 0.05 # stratum-specific noise
psu_factor = 0.5 # importance of PSU effect relative to strata effect
psu_sigma = sqrt(strata_sigma ^ 2 * psu_factor) # psu-specific noise

# create strata-specific random effects
nbasis  = 5
basis  = fda::create.bspline.basis(c(0, 1), nbasis)
Phi  = fda::eval.basis(grid, basis)



set.seed(seed)
strata_scores  = matrix(rnorm(num_strata * nbasis, 0, strata_sigma),
                        num_strata,
                        nbasis) # num_strata x nbasis matrix

strata_random_effects  = strata_scores %*% t(Phi) # num_strata x nbasis matrix

strata_random_effects %>%
  as_tibble() %>%
  mutate(strata = row_number()) %>%
  pivot_longer(cols = -strata,
               names_to = "l",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  ggplot(aes(x = l, y = value, color = factor(strata))) +
  geom_line(linewidth = 1.1) +
  scale_color_viridis_d(option = "A") +
  labs(x = "Functional Domain", y = "Value", title = "Random effects by stratum") +
  theme(legend.position = "none")

total_psu = length(unique(psu_assignments))

set.seed(seed)
psu_scores  = matrix(rnorm(total_psu * nbasis, 0, psu_sigma), total_psu, nbasis)
psu_random_effects  = psu_scores %*% t(Phi)

psu_random_effects %>%
  as_tibble() %>%
  mutate(psu = row_number()) %>%
  filter(psu %in% 1:50) %>%
  pivot_longer(cols = -psu,
               names_to = "l",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  ggplot(aes(x = l, y = value, color = factor(psu))) +
  geom_line(linewidth = 1.1) +
  scale_color_viridis_d(option = "A") +
  labs(x = "Functional Domain", y = "Value", title = "Random effects by PSU",
       subtitle = "For 50 randomly selected PSU") +
  theme(legend.position = "none")




snr_b = 1

# adjust random effect based on relative importance of random effects
fixef_signal  = matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
  X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

# include stratum-specific slope variation in the random effects
slope_re  = (stratum_scaling[stratum_assignments] - 1) *
  matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

strata_effects_indiv  = strata_random_effects[stratum_assignments, ]
psu_effects_indiv  = psu_random_effects[as.numeric(factor(psu_assignments)), ]
random_effects  = strata_effects_indiv + psu_effects_indiv

ranef  = slope_re + random_effects
ranef  = sd(fixef_signal) / sd(ranef) / snr_b * ranef
rm(random_effects)
lin_pred = fixef_signal + ranef

psu_df = expand_grid(strata = 1:10,
                     psu = 1:5) %>%
  mutate(x = paste0(strata, "_", psu))

re_df =
  ranef %>%
  as_tibble() %>%
  bind_cols(psu = psu_assignments,
            strata = stratum_assignments)

write_rds(re_df, here::here("data", "re_df.rds"))

re_df |>
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  geom_smooth(se = FALSE, linewidth = .7) +
  facet_wrap(~strata, nrow = 2) +
  labs(x = "Functional Domain", y = "Value", title = "Random effects by stratum and PSU") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

fixef_df =
  fixef_signal %>%
  as_tibble() %>%
  bind_cols(psu = psu_assignments,
            strata = stratum_assignments)

write_rds(fixef_df, here::here("data", "fixef_df.rds"))


fixef_df |>
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  geom_smooth(se = FALSE, linewidth = .7) +
  facet_wrap(~strata, nrow = 2) +
  labs(x = "Functional Domain", y = "Value", title = "Fixed effects by stratum and PSU") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

lp_df =
  lin_pred %>%
  as_tibble() %>%
  bind_cols(psu = psu_assignments,
            strata = stratum_assignments)
write_rds(lp_df, here::here("data", "lp_df.rds"))

lp_df |>
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  geom_smooth(se = FALSE, linewidth = .7) +
  facet_wrap(~strata, nrow = 2) +
  labs(x = "Functional Domain", y = "Value", title = "Smoothed linear predictors by strata and PSU") +
  theme_light(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

family = "gaussian"
snr_eps = 1
if (family == "gaussian") {
  sd_lp = sd(lin_pred)
  sigma = sd_lp / snr_eps
  set.seed(seed)
  Y_obs = matrix(
    rnorm(
      n = I * L,
      mean = as.vector(t(lin_pred)),
      sd = sigma
    ),
    # need to use t to put in correct order
    nrow = I,
    ncol = L,
    byrow = TRUE
  )
} else if (family == "binomial") {
  p_true = plogis(as.vector(t(lin_pred)))
  set.seed(seed)
  Y_obs  = matrix(
    rbinom(n = I * L, size = 1, prob = p_true),
    nrow = I,
    ncol = L,
    byrow = TRUE
  )
} else if (family == "poisson") {
  lam_true = exp(as.vector(t(lin_pred)))
  set.seed(seed)
  Y_obs  = matrix(
    rpois(n = I * L, lambda = lin_pred_vec),
    nrow = I,
    ncol = L,
    byrow = TRUE
  )
}

y_obs_samp = Y_obs %>%
  as_tibble() %>%
  bind_cols(psu = psu_assignments,
            strata = stratum_assignments) %>%
  filter(strata %in% c(1:10),
         psu %in% psu_df$x)

write_rds(y_obs_samp, here::here("data", "y_obs.rds"))

y_obs_samp |>
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  facet_wrap(~strata, nrow = 2) +
  geom_smooth(se = FALSE, linewidth = 0.7) +
  labs(x = "Functional Domain", y = "Value", title = "Mean smoothed outcomes by stratum and psu in 10 strata") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

# psu = selected_psus[1]
get_p_i = function(i, probs) probs[i] * (1 + sum((probs[-i]) / (1-probs[-i])))

num_selected_psu = 2 # we want to select 2 PSUs
X1 = X_des[, 2]
compression = 2
I_n = 50 # number of individuals targeted to select from each PSU
inf_level = 10 # strength of "informativeness"

p_strata_design  = dirichlet_probs
final_sample  = c() # to store final sample ids
p1  = rep(NA, I) # stage 1 selection probability
p2  = rep(NA, I) # stage 2 selection probability
p_overall  = rep(NA, I) # overall selection probability
psus  = c() # to store PSUs
num_selected_psu = 2
# within each strata, select PSU with replacement using PPS
for (strata in 1:num_strata) {
  # strata = 1
  # individuals in the strata
  inds_in_stratum  = which(stratum_assignments == strata)

  # Get PSU sizes in this stratum
  psu_sizes = table(psu_assignments[inds_in_stratum])
  psu_ids = names(psu_sizes)

  # Sample PSUs WITH replacement using PPS
  set.seed(strata + seed)
  selected_psus = sample(psu_ids,
                         size = num_selected_psu,
                         replace = FALSE,
                         prob = psu_sizes)

  psu_probs  = psu_sizes / sum(psu_sizes)  # PPS
  # probability of selection is 1 - (p(not selected both times))
  # psu_prob_selected  = 1 - (1 - psu_probs) ^ num_selected_psu
  # names(psu_prob_selected)  = psu_ids

  psu_prob_selected = map_dbl(.x = match(selected_psus, psu_ids),
                              .f = get_p_i,
                              psu_probs)

  names(psu_prob_selected)  = selected_psus
  # within each selected PSU select individuals based on X1
  for (psu in selected_psus) {
    inds_in_psu  = which(psu_assignments == psu &
                           stratum_assignments == strata)
    if (inf_level == 0) {
      # Uniform sampling
      n = length(inds_in_psu)
      inclusion_probs = rep(1 / n, n)
    } else {
      # Compute mean outcome in PSU
      y_mean = rowMeans(Y_obs[inds_in_psu, ])

      # Compute inclusion score depending on family
      incl_score = switch(
        family,
        "gaussian" = y_mean * inf_level,
        "poisson"  = log(y_mean) * inf_level,
        "binomial" = qlogis(pmin(pmax(y_mean, 1e-6), 1 - 1e-6)) * inf_level,
        stop("Unknown family")
      )

      # Apply compression and map to probabilities
      score_compressed = pmax(pmin(incl_score, compression), -compression)
      inclusion_probs = plogis(score_compressed)
    }
    inclusion_probs_adj  = inclusion_probs / sum(inclusion_probs) * I_n
    inclusion_probs_adj[inclusion_probs_adj > 1]  = 1

    set.seed(strata + seed + which(selected_psus == psu)) # ensure reproducibility
    sampled_units  = inds_in_psu[rbinom(length(inds_in_psu), 1, inclusion_probs_adj) == 1]

    final_sample  = c(final_sample, sampled_units)
    psus  = c(psus, rep(psu, length(sampled_units)))
    p_psu  = psu_prob_selected[which(names(psu_prob_selected) == psu)]

    p1[inds_in_psu]  = p_psu
    p2[inds_in_psu]  = inclusion_probs_adj
    p_overall[inds_in_psu]  = p_psu * inclusion_probs_adj
  }
}


survey_weights  = 1 / p_overall


dat.sim  = data.frame(
  ID = final_sample,
  X = X1[final_sample],
  strata = stratum_assignments[final_sample],
  psu = sub(".*\\_", "", psus),
  weight = survey_weights[final_sample],
  p_stage1 = p1[final_sample],
  p_stage2 = p2[final_sample]
)
Y_sample = Y_obs[final_sample, ]

dat.sim %>%
  ggplot(aes(x = weight)) +
  geom_histogram(bins = 50, fill = "#56B4E9FF", color = "black") +
  labs(x = "Survey Weight", y = "Count")  +
  theme_light(base_size = 14)

quants = quantile(dat.sim$weight, c(1/3, 2/3))

Y_mat = Y_sample %>% as.matrix()

mean_df =
  dat.sim %>%
  mutate(Y = I(Y_mat)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

write_rds(mean_df, here::here("data", "mean_df_survwts.rds"))
mean_df %>%
  ggplot() +
  geom_spaghetti(aes(y = Y_mean_strat, color = quant), linewidth = 1.1, alpha = 1) +
  scale_color_manual(name = "Weight tertile",
                     values = c("blue", "orange", "red"),
                     labels = c("1", "2", "3")) +
  labs(x = "Functional domain", y = "Smoothed outcome") +
  theme_light(base_size = 14)

# repeat but w/o informative sampling


inf_level = 0 # strength of "informativeness"

p_strata_design  = dirichlet_probs
final_sample  = c() # to store final sample ids
p1  = rep(NA, I) # stage 1 selection probability
p2  = rep(NA, I) # stage 2 selection probability
p_overall  = rep(NA, I) # overall selection probability
psus  = c() # to store PSUs
num_selected_psu = 2
# within each strata, select PSU with replacement using PPS
for (strata in 1:num_strata) {
  # strata = 1
  # individuals in the strata
  inds_in_stratum  = which(stratum_assignments == strata)

  # Get PSU sizes in this stratum
  psu_sizes = table(psu_assignments[inds_in_stratum])
  psu_ids = names(psu_sizes)

  # Sample PSUs WITH replacement using PPS
  set.seed(strata + seed)
  selected_psus = sample(psu_ids,
                         size = num_selected_psu,
                         replace = FALSE,
                         prob = psu_sizes)

  psu_probs  = psu_sizes / sum(psu_sizes)  # PPS
  # probability of selection is 1 - (p(not selected both times))
  # psu_prob_selected  = 1 - (1 - psu_probs) ^ num_selected_psu
  # names(psu_prob_selected)  = psu_ids

  psu_prob_selected = map_dbl(.x = match(selected_psus, psu_ids),
                              .f = get_p_i,
                              psu_probs)

  names(psu_prob_selected)  = selected_psus
  # within each selected PSU select individuals based on X1
  for (psu in selected_psus) {
    inds_in_psu  = which(psu_assignments == psu &
                           stratum_assignments == strata)
    if (inf_level == 0) {
      # Uniform sampling
      n = length(inds_in_psu)
      inclusion_probs = rep(1 / n, n)
    } else {
      # Compute mean outcome in PSU
      y_mean = rowMeans(Y_obs[inds_in_psu, ])

      # Compute inclusion score depending on family
      incl_score = switch(
        family,
        "gaussian" = y_mean * inf_level,
        "poisson"  = log(y_mean) * inf_level,
        "binomial" = qlogis(pmin(pmax(y_mean, 1e-6), 1 - 1e-6)) * inf_level,
        stop("Unknown family")
      )

      # Apply compression and map to probabilities
      score_compressed = pmax(pmin(incl_score, compression), -compression)
      inclusion_probs = plogis(score_compressed)
    }
    inclusion_probs_adj  = inclusion_probs / sum(inclusion_probs) * I_n
    inclusion_probs_adj[inclusion_probs_adj > 1]  = 1

    set.seed(strata + seed + which(selected_psus == psu)) # ensure reproducibility
    sampled_units  = inds_in_psu[rbinom(length(inds_in_psu), 1, inclusion_probs_adj) == 1]

    final_sample  = c(final_sample, sampled_units)
    psus  = c(psus, rep(psu, length(sampled_units)))
    p_psu  = psu_prob_selected[which(names(psu_prob_selected) == psu)]

    p1[inds_in_psu]  = p_psu
    p2[inds_in_psu]  = inclusion_probs_adj
    p_overall[inds_in_psu]  = p_psu * inclusion_probs_adj
  }
}


survey_weights  = 1 / p_overall


dat.sim  = data.frame(
  ID = final_sample,
  X = X1[final_sample],
  strata = stratum_assignments[final_sample],
  psu = sub(".*\\_", "", psus),
  weight = survey_weights[final_sample],
  p_stage1 = p1[final_sample],
  p_stage2 = p2[final_sample]
)
Y_sample = Y_obs[final_sample, ]

dat.sim %>%
  ggplot(aes(x = weight)) +
  geom_histogram(bins = 50, fill = "#56B4E9FF", color = "black") +
  labs(x = "Survey Weight", y = "Count")  +
  theme_light(base_size = 14)

quants = quantile(dat.sim$weight, c(1/3, 2/3))

Y_mat = Y_sample %>% as.matrix()

mean_df =
  dat.sim %>%
  mutate(Y = I(Y_mat)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

write_rds(mean_df, here::here("data", "mean_df_survwts_random.rds"))
mean_df %>%
  ggplot() +
  geom_spaghetti(aes(y = Y_mean_strat, color = quant), linewidth = 1.1, alpha = 1) +
  scale_color_manual(name = "Weight tertile",
                     values = c("blue", "orange", "red"),
                     labels = c("1", "2", "3")) +
  labs(x = "Functional domain", y = "Smoothed outcome") +
  theme_light(base_size = 14)


