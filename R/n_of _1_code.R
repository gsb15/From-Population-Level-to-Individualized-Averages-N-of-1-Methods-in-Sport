# Open-Access GARMIN Running Through N-of-1 (HR as outcome; session-level)
# ===================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(broom)
  library(patchwork)
  library(scales)
})


# Helpers
# -------------------------
tidy_gls <- function(mod, conf.int = TRUE, conf.level = 0.95) {
  stopifnot(inherits(mod, "gls"))
  tt <- summary(mod)$tTable
  out <- tibble(
    term      = rownames(tt),
    estimate  = tt[, "Value"],
    std.error = tt[, "Std.Error"],
    statistic = tt[, "t-value"],
    p.value   = tt[, "p-value"]
  )
  if (conf.int) {
    ints <- try(nlme::intervals(mod, level = conf.level)$coef, silent = TRUE)
    if (!inherits(ints, "try-error")) {
      out$conf.low  <- ints[, 1]
      out$conf.high <- ints[, 3]
    } else {
      df_wald <- max(1, mod$dims$N - mod$dims$p)
      q <- stats::qt(1 - (1 - conf.level)/2, df = df_wald)
      out$conf.low  <- out$estimate - q * out$std.error
      out$conf.high <- out$estimate + q * out$std.error
    }
  }
  rownames(out) <- NULL
  out
}

get_glance_vals <- function(mod) {
  g <- broom::glance(mod)
  tibble(AIC = g$AIC, BIC = g$BIC, logLik = g$logLik)
}

metrics_from_fit <- function(obs, pred) {
  res <- obs - pred
  tibble(
    RMSE = sqrt(mean(res^2, na.rm = TRUE)),
    MAE  = mean(abs(res), na.rm = TRUE),
    R2   = 1 - sum((obs - pred)^2, na.rm = TRUE) /
      sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  )
}

add_session_lags <- function(dat) {
  dat %>%
    arrange(t) %>%
    mutate(
      hr_lag1    = lag(hr, 1),
      speed_lag1 = lag(speed, 1),
      dist_lag1  = lag(distance, 1)
    )
}

center_predictors <- function(dat) {
  dat %>%
    mutate(
      hr_c   = as.numeric(scale(hr,   center = TRUE, scale = FALSE)),
      t_c    = as.numeric(scale(t,    center = TRUE, scale = FALSE)),
      speed_c= as.numeric(scale(speed,center = TRUE, scale = FALSE)),
      dist_c = as.numeric(scale(distance, center = TRUE, scale = FALSE)),
      hr_lag1_c    = if ("hr_lag1" %in% names(dat)) as.numeric(scale(hr_lag1, TRUE, FALSE)) else NA_real_,
      dist_lag1_c  = if ("dist_lag1" %in% names(dat)) as.numeric(scale(dist_lag1, TRUE, FALSE)) else NA_real_,
      speed_lag1_c = if ("speed_lag1" %in% names(dat)) as.numeric(scale(speed_lag1, TRUE, FALSE)) else NA_real_
    )
}

prepare_model_frame <- function(dat, rhs_terms, include_outcome = "hr_c") {
  dc <- center_predictors(dat)
  needed <- c(include_outcome, rhs_terms)
  miss <- setdiff(needed, names(dc))
  if (length(miss)) stop("Missing columns for model terms: ", paste(miss, collapse = ", "))
  dc <- dc %>% tidyr::drop_na(all_of(needed))
  if (nrow(dc) < 3) stop("Not enough complete rows after NA drop.")
  dc
}

fit_ols_hr <- function(dat, include_time = TRUE, include_lags = TRUE) {
  rhs <- c("speed_c", "dist_c")
  if (include_time) rhs <- c(rhs, "t_c")
  if (include_lags) rhs <- c(rhs, "hr_lag1_c", "dist_lag1_c", "speed_lag1_c")
  dc <- prepare_model_frame(dat, rhs, include_outcome = "hr_c")
  fm <- as.formula(paste("hr_c ~", paste(rhs, collapse = " + ")))
  lm(fm, data = dc, na.action = na.omit)
}

fit_gls_ar1_hr <- function(dat, include_time = TRUE, include_lags = TRUE) {
  rhs <- c("speed_c", "dist_c")
  if (include_time) rhs <- c(rhs, "t_c")
  if (include_lags) rhs <- c(rhs, "hr_lag1_c", "dist_lag1_c", "speed_lag1_c")
  dc <- prepare_model_frame(dat, rhs, include_outcome = "hr_c")
  fm <- as.formula(paste("hr_c ~", paste(rhs, collapse = " + ")))
  tryCatch(
    nlme::gls(
      fm, data = dc,
      correlation = nlme::corAR1(form = ~ t),
      method = "REML", na.action = na.omit
    ),
    error = function(e) NULL
  )
}

make_panels_by_x <- function(df_pred, xvar = c("speed","distance"),
                             title = NULL,
                             xlab = NULL,
                             ylab = "Heart rate (bpm)",
                             loess_span = 0.9,
                             ncol = 3) {
  xvar <- match.arg(xvar)
  stopifnot(all(c("runner_id","hr","pred","model",xvar) %in% names(df_pred)))
  
  ggplot(df_pred, aes(x = .data[[xvar]], y = hr)) +
    geom_point(alpha = 0.35, size = 1.1) +
    geom_smooth(method = "loess", se = FALSE, span = loess_span, linewidth = 0.9) +
    geom_line(aes(y = pred, color = model), linewidth = 1.05, alpha = 0.9) +
    facet_wrap(~ runner_id, scales = "free_x", ncol = ncol) +
    labs(title = title, x = xlab %||% xvar, y = ylab, color = "Model") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) && nzchar(a[1])) a else b


# 0) CONFIG / IO
# -------------------------
data_dir  <- "D:/Lit Review"
data_file <- "RunningThrough_Reliability_Data.csv"
out_dir   <- file.path(data_dir, "n1_hr3_simplified_outputs")

min_obs_required <- 10
case_ids <- c(173, 673, 686, 1614, 3078, 4118)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(data_dir)


# 1) LOAD & PREP
# -------------------------
df_raw <- read.csv(data_file)

df <- df_raw %>%
  mutate(
    runner_id  = as.factor(redcap_record_id),
    week       = suppressWarnings(as.numeric(week)),
    gender     = as.factor(gender),
    age_strata = as.factor(age_strata),
    speed      = as.numeric(averageSpeedInMetersPerSecond),
    distance   = as.numeric(distance_kilometers),
    hr         = as.numeric(averageHeartRateInBeatsPerMinute)
  )

df_sess <- df %>%
  group_by(runner_id) %>%
  filter(
    n() >= min_obs_required,
    all(!is.na(speed), !is.na(distance), !is.na(hr), !is.na(gender), !is.na(age_strata))
  ) %>%
  arrange(runner_id, week, .by_group = TRUE) %>%
  mutate(t = row_number()) %>%
  ungroup()

cat("Runners kept (>= ", min_obs_required, " complete sessions): ",
    n_distinct(df_sess$runner_id), "\n", sep = "")

df_cases <- df_sess %>% filter(as.integer(as.character(runner_id)) %in% case_ids)


# 2) DESCRIPTIVES
# -------------------------
overall_summary <- df_sess %>%
  summarise(
    n_runners  = n_distinct(runner_id),
    total_obs  = n(),
    mean_speed = mean(speed, na.rm = TRUE),
    sd_speed   = sd(speed, na.rm = TRUE),
    mean_dist  = mean(distance, na.rm = TRUE),
    sd_dist    = sd(distance, na.rm = TRUE),
    mean_hr    = mean(hr, na.rm = TRUE),
    sd_hr      = sd(hr, na.rm = TRUE)
  ) %>%
  transmute(
    n_runners, total_obs,
    `Speed (m/s)`   = sprintf("%.2f ± %.2f", mean_speed, sd_speed),
    `Distance (km)` = sprintf("%.1f ± %.1f", mean_dist, sd_dist),
    `HR (bpm)`      = sprintf("%.0f ± %.0f", mean_hr, sd_hr)
  )

print(overall_summary)


# 3) GROUP MODEL (session-level; HR outcome)
# ============================================================
df_group_sess <- df_sess %>%
  group_by(runner_id) %>%
  mutate(across(c(hr, speed, distance), as.numeric)) %>%
  arrange(t, .by_group = TRUE) %>%
  add_session_lags() %>%
  ungroup() %>%
  mutate(
    t_c = as.numeric(scale(t, center = TRUE, scale = FALSE)),
    speed_c = as.numeric(scale(speed, center = TRUE, scale = FALSE)),
    dist_c  = as.numeric(scale(distance, center = TRUE, scale = FALSE)),
    hr_lag1_c    = as.numeric(scale(hr_lag1, center = TRUE, scale = FALSE)),
    speed_lag1_c = as.numeric(scale(speed_lag1, center = TRUE, scale = FALSE)),
    dist_lag1_c  = as.numeric(scale(dist_lag1, center = TRUE, scale = FALSE))
  )

n_used_nlme <- function(mod, n_total) {
  if (!is.null(mod$na.action)) {
    length(setdiff(seq_len(n_total), as.integer(mod$na.action)))
  } else {
    n_total
  }
}

safe_fit_one_group <- function(dat, fixed_form, random_form = NULL, use_ar1 = TRUE) {
  if (!is.null(random_form)) {
    if (use_ar1) {
      try(nlme::lme(
        fixed = fixed_form,
        random = random_form,
        correlation = nlme::corAR1(form = ~ t | runner_id),
        data = dat, method = "REML", na.action = na.omit
      ), silent = TRUE)
    } else {
      try(nlme::lme(
        fixed = fixed_form,
        random = random_form,
        data = dat, method = "REML", na.action = na.omit
      ), silent = TRUE)
    }
  } else {
    if (use_ar1) {
      try(nlme::gls(
        model = fixed_form,
        correlation = nlme::corAR1(form = ~ t),
        data = dat, method = "REML", na.action = na.omit
      ), silent = TRUE)
    } else {
      try(nlme::gls(
        model = fixed_form,
        data = dat, method = "REML", na.action = na.omit
      ), silent = TRUE)
    }
  }
}

fixed_base <- hr ~ speed_c + dist_c + t_c + gender + age_strata
fixed_lags <- hr ~ speed_c + dist_c + t_c + hr_lag1_c + dist_lag1_c + speed_lag1_c + gender + age_strata

group_candidates <- list(
  list(name = "M0_base_RI_AR1",         fixed = fixed_base, random = ~ 1 | runner_id, use_ar1 = TRUE),
  list(name = "M1_base_RI_noAR",        fixed = fixed_base, random = ~ 1 | runner_id, use_ar1 = FALSE),
  list(name = "M2_lags_RI_AR1",         fixed = fixed_lags, random = ~ 1 | runner_id, use_ar1 = TRUE),
  list(name = "M3_lags_RI_noAR",        fixed = fixed_lags, random = ~ 1 | runner_id, use_ar1 = FALSE),
  list(name = "M4_lags_RSspeed_AR1",    fixed = fixed_lags, random = ~ 1 + speed_c | runner_id, use_ar1 = TRUE),
  list(name = "M5_lags_RSspeed_noAR",   fixed = fixed_lags, random = ~ 1 + speed_c | runner_id, use_ar1 = FALSE),
  list(name = "M6_lags_RSspdDist_AR1",  fixed = fixed_lags, random = ~ 1 + speed_c + dist_c | runner_id, use_ar1 = TRUE),
  list(name = "M7_lags_RSspdDist_noAR", fixed = fixed_lags, random = ~ 1 + speed_c + dist_c | runner_id, use_ar1 = FALSE)
)

group_fits <- lapply(group_candidates, function(spec) {
  fit <- safe_fit_one_group(df_group_sess, spec$fixed, spec$random, spec$use_ar1)
  list(spec = spec, fit = fit)
})

group_model_comp <- purrr::map_dfr(group_fits, function(x) {
  spec <- x$spec; fit <- x$fit
  if (inherits(fit, "try-error")) {
    tibble(model = spec$name, converged = FALSE,
           AIC = NA_real_, BIC = NA_real_, logLik = NA_real_,
           n_used = NA_integer_)
  } else {
    g <- get_glance_vals(fit)
    tibble(model = spec$name, converged = TRUE,
           AIC = g$AIC, BIC = g$BIC, logLik = g$logLik,
           n_used = n_used_nlme(fit, nrow(df_group_sess)))
  }
}) %>%
  arrange(if_else(is.na(AIC), Inf, AIC))

print(group_model_comp)

best_name <- group_model_comp %>% filter(converged) %>% slice(1) %>% pull(model)
group_lme_hr <- purrr::keep(group_fits, ~ !inherits(.x$fit, "try-error") && .x$spec$name == best_name)[[1]]$fit

cat("\nSelected group model by AIC:", best_name, "\n")
print(summary(group_lme_hr))

df_group_sess$pred_group_hr <- NA_real_
idx <- if (!is.null(group_lme_hr$na.action)) {
  setdiff(seq_len(nrow(df_group_sess)), as.integer(group_lme_hr$na.action))
} else seq_len(nrow(df_group_sess))

pred0 <- if (inherits(group_lme_hr, "lme")) predict(group_lme_hr, level = 0) else predict(group_lme_hr)
df_group_sess$pred_group_hr[idx] <- as.numeric(pred0)


# Group figure: observed vs predicted by session index
# -------------------------
group_t_obs <- df_group_sess %>%
  group_by(t) %>%
  summarise(mean_hr = mean(hr, na.rm = TRUE),
            se_hr = sd(hr, na.rm = TRUE)/sqrt(n()),
            .groups = "drop")

group_t_pred <- df_group_sess %>%
  group_by(t) %>%
  summarise(mean_pred = mean(pred_group_hr, na.rm = TRUE),
            se_pred = sd(pred_group_hr, na.rm = TRUE)/sqrt(sum(!is.na(pred_group_hr))),
            .groups = "drop")

fig_group_t <- ggplot() +
  geom_ribbon(data = group_t_obs,
              aes(x = t, ymin = mean_hr - 1.96*se_hr, ymax = mean_hr + 1.96*se_hr),
              alpha = 0.18) +
  geom_line(data = group_t_obs, aes(x = t, y = mean_hr), linewidth = 1.0) +
  geom_line(data = group_t_pred, aes(x = t, y = mean_pred),
            linewidth = 1.0, linetype = "dashed") +
  labs(title = "Group HR: Observed (solid) vs Predicted (dashed) by Session Index",
       x = "Session index (t)", y = "HR (bpm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

print(fig_group_t)


# 4) N-of-1 MODELS (per runner; OLS vs GLS-AR1)
# ============================================================
rhs_full <- c("speed_c","dist_c","t_c","hr_lag1_c","dist_lag1_c","speed_lag1_c")

compare_models_hr <- df_cases %>%
  arrange(runner_id, t) %>%
  group_by(runner_id) %>%
  group_modify(~{
    dat <- .x %>% add_session_lags() %>% mutate(t = t)
    
    # OLS
    m_lm  <- fit_ols_hr(dat, include_time = TRUE, include_lags = TRUE)
    dc_lm <- prepare_model_frame(dat, rhs_full, include_outcome = "hr_c")
    p_lm  <- predict(m_lm, newdata = dc_lm, na.action = na.omit)
    met_lm <- metrics_from_fit(dc_lm$hr_c, p_lm)
    gl_lm  <- get_glance_vals(m_lm)
    
    # GLS-AR1
    m_gls <- fit_gls_ar1_hr(dat, include_time = TRUE, include_lags = TRUE)
    if (!is.null(m_gls)) {
      dc_gls <- prepare_model_frame(dat, rhs_full, include_outcome = "hr_c")
      p_gls  <- as.numeric(predict(m_gls, newdata = dc_gls))
      met_gls <- metrics_from_fit(dc_gls$hr_c, p_gls)
      gl_gls  <- get_glance_vals(m_gls)
    } else {
      met_gls <- tibble(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_)
      gl_gls  <- tibble(AIC = NA_real_, BIC = NA_real_, logLik = NA_real_)
    }
    
    coef_lm  <- broom::tidy(m_lm, conf.int = TRUE) %>% mutate(model = "OLS")
    coef_gls <- if (!is.null(m_gls)) tidy_gls(m_gls, conf.int = TRUE) %>% mutate(model = "GLS_AR1") else tibble()
    
    bind_rows(
      tibble(model = "OLS",    AIC = gl_lm$AIC,  BIC = gl_lm$BIC,  logLik = gl_lm$logLik,
             RMSE = met_lm$RMSE, MAE = met_lm$MAE, R2 = met_lm$R2, coef = list(coef_lm)),
      tibble(model = "GLS_AR1",AIC = gl_gls$AIC, BIC = gl_gls$BIC, logLik = gl_gls$logLik,
             RMSE = met_gls$RMSE, MAE = met_gls$MAE, R2 = met_gls$R2, coef = list(coef_gls))
    )
  }) %>%
  ungroup()

model_comp_table_hr <- compare_models_hr %>%
  select(runner_id, model, AIC, BIC, logLik, RMSE, MAE, R2) %>%
  arrange(runner_id, desc(model))

coef_comp_table_hr <- compare_models_hr %>%
  select(runner_id, model, coef) %>%
  unnest(coef) %>%
  mutate(term = recode(term,
                       "(Intercept)" = "Intercept",
                       "t_c"         = "Observation index (centered)",
                       "speed_c"     = "Speed (centered)",
                       "dist_c"      = "Session distance (centered)",
                       "hr_lag1_c"   = "Lag HR (centered)",
                       "dist_lag1_c" = "Lag Session distance (centered)",
                       "speed_lag1_c"= "Lag Speed (centered)"
  )) %>%
  arrange(runner_id, model, term)

winner_hr <- model_comp_table_hr %>%
  group_by(runner_id) %>%
  summarise(
    best_by_AIC  = model[which.min(AIC)],
    best_by_RMSE = model[which.min(RMSE)],
    .groups = "drop"
  )

cat("\n—— model_comp_table_hr ——\n"); print(model_comp_table_hr)
cat("\n—— winner_hr ——\n"); print(winner_hr)


# N-of-1 Figures: speed & distance panels (observed + model lines)
# -------------------------

pred_df <- df_cases %>%
  arrange(runner_id, t) %>%
  group_by(runner_id) %>%
  group_modify(~{
    dat <- .x %>% add_session_lags()
    dc  <- prepare_model_frame(dat, rhs_full, include_outcome = "hr_c")
    
    m_lm  <- fit_ols_hr(dat, include_time = TRUE, include_lags = TRUE)
    p_lm  <- as.numeric(predict(m_lm, newdata = dc))
    out_lm <- dc %>% transmute(runner_id, hr = hr, speed = speed, distance = distance, model = "OLS", pred = hr_c + p_lm)
    
    m_gls <- fit_gls_ar1_hr(dat, include_time = TRUE, include_lags = TRUE)
    if (!is.null(m_gls)) {
      p_gls <- as.numeric(predict(m_gls, newdata = dc))
      out_gls <- dc %>% transmute(runner_id, hr = hr, speed = speed, distance = distance, model = "GLS_AR1", pred = hr_c + p_gls)
      bind_rows(out_lm, out_gls)
    } else {
      out_lm
    }
  }) %>%
  ungroup()

fig1_speed <- make_panels_by_x(
  pred_df, xvar = "speed",
  title = "Runner-specific HR vs Speed (observed + model predictions)",
  xlab = "Speed (m/s)"
)

fig2_dist <- make_panels_by_x(
  pred_df, xvar = "distance",
  title = "Runner-specific HR vs Distance (observed + model predictions)",
  xlab = "Distance (km)"
)

print(fig1_speed)
print(fig2_dist)


# N-of-1 Figure: coefficient forests (OLS vs GLS-AR1)
# -------------------------
predictor_order_hr <- c(
  "Observation index (centered)",
  "Speed (centered)",
  "Session distance (centered)",
  "Lag HR (centered)",
  "Lag Session distance (centered)",
  "Lag Speed (centered)"
)

coef_plot_df_hr <- coef_comp_table_hr %>%
  filter(term != "Intercept", model %in% c("GLS_AR1","OLS")) %>%
  mutate(
    term  = factor(term, levels = predictor_order_hr),
    model = factor(model, levels = c("GLS_AR1","OLS"))
  )

fig3_forests <- ggplot(coef_plot_df_hr,
                       aes(x = estimate, y = term, color = model, shape = model)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8, stroke = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.6),
                 height = 0.2, linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ runner_id, scales = "free_x") +
  labs(title = "Runner-Specific Coefficient Forests (HR outcome)",
       x = "Coefficient (95% CI)", y = "Predictor",
       color = "Model", shape = "Model") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

print(fig3_forests)


# 5) SAVE (optional)
# -------------------------
# write.csv(overall_summary,       file.path(out_dir, "overall_summary.csv"), row.names = FALSE)
# write.csv(group_model_comp,      file.path(out_dir, "group_model_comparison.csv"), row.names = FALSE)
# write.csv(model_comp_table_hr,   file.path(out_dir, "n1_model_comparison.csv"), row.names = FALSE)
# write.csv(coef_comp_table_hr,    file.path(out_dir, "n1_coefficients.csv"), row.names = FALSE)
# write.csv(winner_hr,             file.path(out_dir, "n1_winners.csv"), row.names = FALSE)
# ggsave(file.path(out_dir, "fig_group_t.png"), fig_group_t, width = 8, height = 5, dpi = 320)
# ggsave(file.path(out_dir, "fig1_speed.png"),  fig1_speed,  width = 8, height = 7, dpi = 320)
# ggsave(file.path(out_dir, "fig2_dist.png"),   fig2_dist,   width = 8, height = 7, dpi = 320)
# ggsave(file.path(out_dir, "fig3_forests.png"),fig3_forests,width = 9, height = 7, dpi = 320)


# End
# ===================================================================
