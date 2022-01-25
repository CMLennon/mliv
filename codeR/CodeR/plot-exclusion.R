

# Notes ----------------------------------------------------------------------------------
#   Goal: Plot 'higher-order' exclusion-restriction results.


# Setup ----------------------------------------------------------------------------------
  # Package
  library(pacman)
  p_load(
    viridis, plotly, tidyverse, ggridges, extrafont, latex2exp, patchwork,
    data.table, magrittr, here
  )


# Load data: Results ---------------------------------------------------------------------
  # Load F1 transformation 3 (exists in two parts)
  f13a_dt = here(
    "Resources-Nonlinear",
    "sim-results-nonlin-FALSE-1-datatrans-3-cf3-dgp-ovb-2sls-lasso-plasso-rfcv-pca-gboost-split-liml-jive.csv"
  ) %>% read_csv()
  f13b_dt = here(
    "Resources-Nonlinear",
    "sim-results-nonlin-FALSE-1-datatrans-3-cf3-dgp-ovb-2sls-split-nnetw.csv"
  ) %>% read_csv() %>% filter(str_detect(model, "Neural Net"))
  # Combine datasets
  f13_dt = rbindlist(list(f13a_dt, f13b_dt), use.names = T, fill = T)
  rm(f13a_dt, f13b_dt)
  # Fix names (normalized RF shows up twice—one is boosted)
  f13_dt[model == "First stage: Random Forest norm, CV" & term == "fst_stage_preds_boost", `:=`(
    model = "First stage: Boosted trees norm, CV"
  )]


# Data work: Cleaning --------------------------------------------------------------------
  # Change names
  f13_dt[, model := str_replace_all(model, "First Stage", "First stage")]
  f13_dt[str_detect(model, "2SLS"), model := "First stage: OLS"]
  f13_dt[str_detect(model, "Naive OLS"), model := "Naive OLS"]
  f13_dt[str_detect(model, "Oracle Model"), model := "\"Oracle\" model"]
  f13_dt[, model := str_replace_all(model, "Random Forest", "RF")]
  f13_dt[, model := str_replace_all(model, "LASSO", "Lasso")]
  f13_dt[, model := str_replace_all(model, "post-", "Post-")]
  f13_dt[, model := str_replace_all(model, "Trees", "trees")]
  f13_dt[, model := str_replace_all(model, "First stage: USSIV", "Split-sample IV")]
  f13_dt[, model := str_replace_all(model, "First stage: JIVE", "Jackknife IV (JIVE)")]
  f13_dt[model == "First stage: LIML (Fuller)", model := "LIML (Fuller)"]
  # Rename neural nets (OOS)
  f13_dt[model == "First stage: Unrestricted Neural Net OOS", model := "First stage: Neural net, unrestricted, OOS"]
  # Rename neural nets (normalized)
  f13_dt[str_detect(model, "Unrestricted Neural Net norm OOS"), model := "First stage: Neural net, unrestricted, OOS, norm"]
  # Abbreviation for unrestricted
  f13_dt[, model := str_replace_all(model, "unrestricted", "unrest.")]
  # Abbreviation for normalized
  f13_dt[, model := str_replace_all(model, "norm", "norm.")]
  # Ensure there's a comma before "norm."
  f13_dt[, model := str_replace_all(model, "(?<=[A-z])\\snorm", ", norm")]
  # Remove 'CV'
  f13_dt[, model := str_remove_all(model, ", CV")]
  # Fix Lasso
  f13_dt[model == "First stage: norm. Lasso selection", model := "First stage: Lasso selection, norm."]
  # Change "First stage" to "2SLS"
  f13_dt[, model := str_replace_all(model, "First stage:", "2SLS:")]
  # Order the models
  all_models = c(
    "2SLS: RF, norm.",
    "2SLS: RF",
    "2SLS: Boosted trees, norm.",
    "2SLS: Boosted trees",
    "2SLS: Neural net, unrest., OOS, norm.",
    "2SLS: Neural net, unrest., OOS",
    "2SLS: Lasso selection, norm.",
    "2SLS: Lasso selection",
    "2SLS: Post-Lasso selection",
    "2SLS: PCA",
    "Jackknife IV (JIVE)",
    "Split-sample IV",
    "LIML (Fuller)",
    "2SLS: OLS",
    "Naive OLS",
    "\"Oracle\" model"
  )
  f13_dt[, model := factor(model, levels = all_models %>% rev(), ordered = T)]
  # # Drop OOS neural net models
  # f13_dt %<>% .[str_detect(model, "OOS", negate = T)]


# Figure: Interaction-genderated bias (line plots) ---------------------------------------
  # The figure: Loop over classes of models
  plot_list = lapply(
    X = list(
      "1-linear" = f13_dt[str_detect(model, "LIML|JIVE|Split|OLS")],
      "2-selection" = f13_dt[str_detect(model, "Lasso|PCA")],
      "3-trees" = f13_dt[str_detect(model, "RF|tree")],
      "4-nets" = f13_dt[str_detect(model, "net")]
    ),
    FUN = function(d) {
      ggplot(
        data = d[, .(
          est_mean = mean(estimate),
          est_median = median(estimate),
          est_hi = quantile(estimate, probs = 0.975),
          est_lo = quantile(estimate, probs = 0.025)
        ), by = .(model, deg = int_degree)],
        aes(x = deg, color = model, fill = model)
      ) +
      geom_hline(yintercept = 1, size = 0.5, linetype = "dashed") +
      geom_ribbon(
        aes(ymin = est_lo, ymax = est_hi),
        alpha = 0.2,
        color = NA
      ) +
      geom_line(
        aes(y = est_mean),
        size = 0.4
      ) +
      geom_point(
        aes(y = est_mean),
        size = 0.7
      ) +
      ylab("Estimate") +
      scale_x_continuous(
        "Level of interaction of instruments in structural error term",
        breaks = 1:7
      ) +
      scale_fill_viridis_d("", option = "magma", end = 0.9, direction = -1) +
      scale_color_viridis_d("", option = "magma", end = 0.9, direction = -1) +
      theme_minimal(
        base_size = 7,
        base_family = "Charter"
      ) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.key.height = unit(1, "mm"),
        legend.key.width = unit(4.5, "mm")
      ) +
      coord_cartesian(ylim = c(0.23, 2.4))
    }
  )
  # Save individual plots
  for (i in seq_along(plot_list)) {
    ggsave(
      plot = plot_list[[i]],
      path = here("Figures"),
      filename = paste0("sim-exclusion-interaction-", names(plot_list)[i], ".pdf"),
      device = cairo_pdf,
      width = 4.5,
      height = 2.5
    )
  }


# Figure: Interaction-genderated bias (density plots) ------------------------------------
  # Create the individual figures
  figure_list = lapply(
    X = list(
      2 = f13_dt[int_degree == 2],
      3 = f13_dt[int_degree == 3],
      4 = f13_dt[int_degree == 4],
      5 = f13_dt[int_degree == 5]
    ),
    FUN = function(d) { 
      ggplot(
        data = d,
        aes(x = estimate, y = model)
      ) +
      geom_density_ridges(
        aes(fill = model, color = model, vline_color = model),
        size = 0.2,
        alpha = 0.85, 
        rel_min_height = 0.0005,
        scale = 5,
        vline_size = 0.5,
        quantile_lines = T,
        quantiles = 2
      ) +
      geom_vline(xintercept = 1, linetype = "longdash", size = 0.2) +
      xlab("Estimate") +
      ylab("") +
      scale_fill_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
      scale_color_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
      scale_discrete_manual(
        "vline_color",
        values = magma(d[,model] %>% uniqueN(), end = 0.95, direction = 1),
        name = "Method"
      ) +
      theme_minimal(
        base_size = 7,
        base_family = "Charter"
      ) +
      theme(legend.position = "none") +
      coord_cartesian(xlim = c(f13_dt[int_degree %in% 2:5,min(estimate)], f13_dt[int_degree %in% 2:5,max(estimate)]))
    }
  )
  # Save individual plots
  for (i in seq_along(figure_list)) {
    ggsave(
      plot = figure_list[[i]],
      path = here("Figures"),
      filename = paste0("sim-exclusion-interaction-density-", (2:5)[i], ".pdf"),
      device = cairo_pdf,
      width = 4.5,
      height = 2.75
    )
  }

