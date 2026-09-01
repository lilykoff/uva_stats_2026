## download data

library(tidyverse)
### process nhanes data to get minute-level steps per subject

if(!dir.exists(here::here("data"))) dir.create(here::here("data"))
# check if data is downloaded
data_files = c(here::here("data", "nhanes_1440_scsslsteps.csv.xz"),
               here::here("data", "nhanes_1440_PAXPREDM.csv.xz"),
               here::here("data", "nhanes_1440_PAXFLGSM.csv.xz"),
               here::here("data", "nhanes_1440_PAXMTSM.csv.xz"))


# fn to check if files exist and download
check_exists_download_file = function(f){
  # f = data_files[1]
  if (!file.exists(f)){
    url = paste0(
      "https://physionet.org/files/minute-level-step-count-nhanes/1.0.2/csv/",
      basename(f)
    )
    curl::curl_download(url, destfile = f)
  }
}

# download covariates file too (from GH)
covar_file = here::here("data", "covariates_accel_mortality_df.rds")

if (!file.exists(covar_file)){
  curl::curl_download(
    "https://raw.githubusercontent.com/lilykoff/nhanes_steps_mortality/main/data/covariates_accel_mortality_df.rds",
    destfile = covar_file
  )
}


# download files
if (!all(file.exists(data_files))) {
  lapply(data_files, check_exists_download_file)
}

# target files
mims_ml = here::here("data", "processed_mims_multilevel.rds")
mims_sl = here::here("data", "mims_covariates.rds")

if (!all(file.exists(c(mims_ml, mims_sl)))) {
  if (!file.exists(here::here("data", "processed_mims_multilevel.rds")) || force) {
    steps = read_csv(here::here("data", "nhanes_1440_scsslsteps.csv.xz"))
    wear = read_csv(here::here("data", "nhanes_1440_PAXPREDM.csv.xz"))
    flags = read_csv(here::here("data", "nhanes_1440_PAXFLGSM.csv.xz"))
    mims = read_csv(here::here("data", "nhanes_1440_PAXMTSM.csv.xz"))

    all.equal(wear %>% select(SEQN, PAXDAYM), flags %>% select(SEQN, PAXDAYM))
    all.equal(mims %>% select(SEQN, PAXDAYM), flags %>% select(SEQN, PAXDAYM))
    all.equal(steps %>% select(SEQN, PAXDAYM), flags %>% select(SEQN, PAXDAYM))


    mims_mat = mims %>% select(starts_with("min")) %>% as.matrix()
    wear_mat = wear %>% select(starts_with("min")) %>% as.matrix()
    flag_mat = flags %>% select(starts_with("min")) %>% as.matrix()
    ## replace entires in wear mat with NA if flag is TRUE
    wear_mat_filt = ifelse(flag_mat, NA, wear_mat)
    ## replace entries in mims mat with NA if flag is TRUE
    mims_mat_filt = ifelse(flag_mat, NA, mims_mat)

    # create logical matrix where wear is 1, 2, or 4
    wear_mat_lgl = matrix(wear_mat_filt %in% c(1, 2, 4),
                          nrow = nrow(wear_mat_filt), ncol = ncol(wear_mat_filt))
    # create logical matrix where true if MIMS > 0
    mims_mat_lgl = matrix(mims_mat_filt > 0,
                          nrow = nrow(mims_mat_filt), ncol = ncol(mims_mat_filt))

    wake_mat_lgl = matrix(wear_mat_filt == 1,
                          nrow = nrow(wear_mat_filt), ncol = ncol(wear_mat_filt))

    # use rowsums to quickly get total number of minute per day with wear
    wear_summary_vec = rowSums(wear_mat_lgl, na.rm = TRUE)
    # use rowsums to get total number of minutes with nonzero MIMS
    mims_summary_vec = rowSums(mims_mat_lgl, na.rm = TRUE)
    # same for wake
    wake_summary_vec = rowSums(wake_mat_lgl, na.rm = TRUE)

    # get wear summary
    wear_summary =
      wear %>%
      select(SEQN, PAXDAYM) %>%
      bind_cols(non_flag_wear = wear_summary_vec) %>%
      bind_cols(non_zero_MIMS = mims_summary_vec) %>%
      bind_cols(wake_wear = wake_summary_vec) %>%
      mutate(include = non_flag_wear >= 1368 & wake_wear >= 420 & non_zero_MIMS >= 420)


    include_days =
      wear_summary %>%
      group_by(SEQN) %>%
      mutate(days_per_sub = sum(include)) %>%
      ungroup() %>%
      mutate(include_final = days_per_sub >= 3 & include) %>%
      filter(include_final) %>%
      select(SEQN, PAXDAYM)

    steps_small =
      steps %>%
      right_join(include_days, by = c("SEQN", "PAXDAYM"))

    mims_small =
      mims %>%
      right_join(include_days, by = c("SEQN", "PAXDAYM"))

    wear_small =
      wear %>%
      right_join(include_days, by = c("SEQN", "PAXDAYM"))

    flags_small =
      flags %>%
      right_join(include_days, by = c("SEQN", "PAXDAYM"))

    # get rid of steps minutes that are nonwear OR flagged
    wear_mat = wear_small %>% select(starts_with("min")) %>% as.matrix()
    flag_mat = flags_small %>% select(starts_with("min")) %>% as.matrix()
    mims_mat = mims_small %>% select(starts_with("min")) %>% as.matrix()
    steps_mat = steps_small %>% select(starts_with("min")) %>% as.matrix()

    wear_mat_lgl = matrix(wear_mat == 3,
                          nrow = nrow(wear_mat), ncol = ncol(wear_mat))
    ## replace entires in wear mat with NA if flag is TRUE
    steps_mat_filt = ifelse(wear_mat_lgl, NA, steps_mat)
    ## replace entries in mims mat with NA if flag is TRUE
    steps_mat_filt2 = ifelse(flag_mat, NA, steps_mat_filt)

    mims_mat_filt = ifelse(wear_mat_lgl, NA, mims_mat)
    mims_mat_filt2 = ifelse(flag_mat, NA, mims_mat_filt)

    steps_final =
      steps_small %>%
      select(-starts_with("min")) %>%
      bind_cols(steps_mat_filt2)

    mims_final =
      mims_small %>%
      select(-starts_with("min")) %>%
      bind_cols(mims_mat_filt2)


    write_rds(steps_final, here::here("data", "processed_steps_multilevel.rds"))
    write_rds(mims_final, mims_ml)

    mims_persub =
      mims_final %>%
      group_by(SEQN) %>%
      summarize(across(starts_with("min"), ~mean(.x, na.rm = TRUE)), .groups = "drop")

    steps_persub =
      steps_final %>%
      group_by(SEQN) %>%
      summarize(across(starts_with("min"), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  } else {
    mims_final = read_rds(here::here("data", "processed_mims_multilevel.rds"))
  }

  covars = read_rds(here::here("data", "covariates_accel_mortality_df.rds"))

  mims_persub =
    mims_final %>%
    group_by(SEQN) %>%
    summarize(across(starts_with("min"), ~mean(.x, na.rm = TRUE)), .groups = "drop")


  pa_df_persub =
    mims_persub %>% mutate(SEQN = as.character(SEQN)) %>%
    left_join(covars, by = "SEQN") %>%
    rename(psu = masked_variance_pseudo_psu,
           strata = masked_variance_pseudo_stratum) %>%
    select(SEQN, data_release_cycle, psu, strata,
           gender, age_in_years_at_screening, full_sample_2_year_interview_weight,
           full_sample_2_year_mec_exam_weight, race_hispanic_origin,
           cat_bmi, starts_with("min_"))

  write_rds(pa_df_persub, mims_sl)
}



