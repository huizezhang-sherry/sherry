library(VGAMdata)
library(VGAM)
library(dplyr)

data("smqP", package = "VGAMdata")

cols <- c('TC30', 'TCquit', 'SCA', 'educ', 'marital',
          'SIA', 'SEQN', 'age', 'gender', 'race2',
          "passiveSmoke_job", "passiveSmoke_rest",
          "passiveSmoke_bar", "passiveSmoke_car",
          "passiveSmoke_home", "passiveSmoke_other",
          "passiveSmoke_ecigarettes")
df <- smqP |> select(all_of(cols))

# Step 1: complete_case_filtering
# Apply complete-case filtering to all columns except TC30, TCquit, and SCA
cols_require_complete <- setdiff(cols, c("TC30", "TCquit", "SCA"))
df_complete <- df |>
  filter(complete.cases(pick(all_of(cols_require_complete))))

# Step 2: age_restriction
# Restrict analysis to respondents aged less than 80 years
df_age_restricted <- df_complete |> filter(age < 80)

# Step 3: data_quality_correction
# Set both SCA and SIA to NA where SCA - SIA < 0
df_age_restricted <- df_age_restricted |>
  mutate(
    SCA_corrected = if_else(!is.na(SCA) & !is.na(SIA) & (SCA - SIA < 0), NA_integer_, SCA),
    SIA_corrected = if_else(!is.na(SCA) & !is.na(SIA) & (SCA - SIA < 0), NA_integer_, SIA)
  )

# Step 4: impute_marital_young
# Recode marital status to "Never married" for everyone age < 20
df_age_restricted <- df_age_restricted |>
  mutate(
    marital_imputed = if_else(age < 20, "Never married", as.character(marital))
  )

# Step 5: impute_education_young
# Recode education based on age for young respondents
df_age_restricted <- df_age_restricted |>
  mutate(
    educ_imputed = case_when(
      age < 14 ~ "< 9th grade",
      age >= 14 & age <= 16 ~ "9-11th grade",
      age >= 17 & age <= 19 ~ "High school",
      TRUE ~ as.character(educ)
    )
  )

# Step 6: define_ever_smokers_sample
# Include if: SIA not missing AND (both SCA and TCquit not missing OR SCA missing and TC30 not missing)
df_ever_smokers <- df_age_restricted |>
  filter(
    !is.na(SIA_corrected) &
    ((!is.na(SCA_corrected) & !is.na(TCquit)) | (is.na(SCA_corrected) & !is.na(TC30)))
  )

df_ever_smokers <- droplevels(df_ever_smokers)

# Step 7: parent_distribution_choice
# Adopt negative binomial distribution as the parent distribution
parent_distribution <- negbinomial

# Step 8: special_value_mid_risk
# Ages 12, 13, 14, 19, 20, 21 as parametrically-inflated group (i.mix)
i_mix_spec <- c(12, 13, 14, 19, 20, 21)

# Step 9: special_value_high_risk
# Ages 15, 16, 17, 18 as nonparametrically-inflated group (i.mlm)
i_mlm_spec <- c(15, 16, 17, 18)

# Step 10: alteration_operator_decision
# Do not use alteration operators
a_mix_spec <- NULL
a_mlm_spec <- NULL

# Step 11: deflation_operator_decision
# Do not use deflation operators
d_mix_spec <- NULL
d_mlm_spec <- NULL

# Step 12: truncation_decision
# Do not use truncation
truncation_spec <- NULL

# Step 13: covariate_selection
# Include race2, gender, marital, and educ as covariates
covariate_formula <- SIA_corrected ~ race2 + gender + marital_imputed + educ_imputed

# Step 14: fit_final_model
# Fit a single vglm with response, parent distribution, special-value structure, and covariates
fitted_model <- vglm(
  covariate_formula,
  family = gaitdnbinomial(
    i.mix = i_mix_spec, i.mlm = i_mlm_spec,
    a.mix = a_mix_spec, a.mlm = a_mlm_spec,
    d.mix = d_mix_spec, d.mlm = d_mlm_spec,
    truncate = truncation_spec
  ),
  data = df_ever_smokers
)
