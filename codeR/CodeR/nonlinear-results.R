# Load all results -----------------------------------------------------------------------
# Read in the results, grouping DGPs into separate objects

one_iter = function(iter,shuffle, data_transform ='1', complex = 'T') {
 #This function reads maimaiin a dataset, runs regressions using IV, SSIV, oracle model, neural network predictions
 #and OOS NN predictions, alongside the SSML approach. These get put into a dataframe
  #----------------------------------------
  ## INPUTS ##
  #Complex - always true for this function. See Belloni et al for the meaning here.
  
  #shuffle - 1 or 2. Takes a string that changes how the strength of instruments corresponds to their power
    # 1. - shuffled strength of coefficients (shuffled beta pattern)
    # 2. - strongest instrument at 1 and declining to 100.
  
  #iter - gives iteration of simulation we are on
  
  #data_transform - this can take on one of three values
    #1. - No transformation. There is no nonlinear noise.
    #2. - Slight learnable disturbance. This adds the noise term if(abs(sin(z1)) > .995), sign(sin(z1))*.1, else, 0 twice over to y and once to x
    #3. - All variables changed to binary policy variables. z1-z100 are normalized 0-1, as are the error terms
  
    #numeric - this maps to the dgp where x = f(z) + e => z_1 + z_2 + z_3 +... z_i^(data_transform) + e and y = x + z_i^data_transform + v
   
  #Default is no noise
  
  data_table_is = dir(here("dgp_folder_weird"),paste0("T-",shuffle, '-', iter, '_',data_transform,'.csv'), full.names = T) %>% 
    map_dfr(fread)
  
  data_table_os = dir(here("dgp_folder_weird"),paste0("T-",shuffle, '-', iter, '_', data_transform,'oos.csv'), full.names = T) %>% 
    map_dfr(fread)
  
  iter_dt = data.table()
  
  iter_dt = lm(y ~ fst_stage_preds_is_1, data = data_table_is) %>% tidy(quick = T) %>%
    filter(grepl('fst', term)) %>%
    mutate(model = "First stage: Deep Neural Net norm IS")%>% mutate(true_var = var(data_table_is$x1)) %>% mutate(true_cov_u = cov(data_table_is$x1, data_table_is$y - data_table_is$x1)) %>%
    mutate(var = var(data_table_is$fst_stage_preds_is_1)) %>% mutate(cov_u = cov(data_table_is$fst_stage_preds_is_1, data_table_is$y - data_table_is$x1), cov_e = cov(data_table_is$fst_stage_preds_is_1, data_table_is$x1 - data_table_is$fst_stage_preds_is_1)) %>% bind_rows(.,iter_dt) %>% 
    data.table()
  
  iter_dt = lm(y ~ x1_hat_nnetw, data = data_table_is) %>% tidy(quick = T) %>%
    filter(grepl('x1', term)) %>%
    mutate(model = "First stage: Deep Neural Net IS")%>% mutate(true_var = var(data_table_is$x1)) %>% mutate(true_cov_u = cov(data_table_is$x1, data_table_is$y - data_table_is$x1)) %>%
    mutate(var = var(data_table_is$fst_stage_preds_is_1)) %>% mutate(cov_u = cov(data_table_is$fst_stage_preds_is_1, data_table_is$y - data_table_is$x1), cov_e = cov(data_table_is$fst_stage_preds_is_1, data_table_is$x1 - data_table_is$fst_stage_preds_is_1)) %>% bind_rows(.,iter_dt) %>% 
    data.table()
  
  iter_dt = lm(y ~ fst_stage_preds_oos_2, data = data_table_os) %>% tidy(quick = T) %>%
    filter(grepl('fst', term)) %>%
    mutate(model = "First stage: Deep Neural Net norm OOS")%>% mutate(true_var = var(data_table_os$x1)) %>% mutate(true_cov_u = cov(data_table_os$x1, data_table_os$y - data_table_os$x1)) %>%
    mutate(var = var(data_table_os$fst_stage_preds_oos_2)) %>% mutate(cov_u = cov(data_table_os$fst_stage_preds_oos_2, data_table_os$y - data_table_os$x1), cov_e = cov(data_table_os$fst_stage_preds_oos_2, data_table_os$x1 - data_table_os$fst_stage_preds_oos_2)) %>% bind_rows(.,iter_dt) %>% 
    data.table()
  
  iter_dt = 
  
}
fread_itercount = function(data, iter){
  dt = fread(data) %>% mutate(iter = iter)
  return(dt)
}
#f1_dt = here("dgp_folder_weird") %>% dir("f-1", full.names = T) %>% map_dfr(fread)
t1_1_is_dt = dir(here("dgp_folder_weird"),"T-1", full.names = T)[grepl("(?=1.csv)",dir(here("dgp_folder_weird"),"T-1", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t1_2_is_dt = dir(here("dgp_folder_weird"),"T-1", full.names = T)[grepl("(?=2.csv)",dir(here("dgp_folder_weird"),"T-1", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t1_3_is_dt =  dir(here("dgp_folder_weird"),"T-1", full.names = T)[grepl("(?=3.csv)",dir(here("dgp_folder_weird"),"T-1", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t2_1_is_dt = dir(here("dgp_folder_weird"),"T-2", full.names = T)[grepl("(?=1.csv)",dir(here("dgp_folder_weird"),"T-2", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t2_2_is_dt = dir(here("dgp_folder_weird"),"T-2", full.names = T)[grepl("(?=2.csv)",dir(here("dgp_folder_weird"),"T-2", full.names = T),perl=TRUE)]%>% 
   map_dfr(fread)
t2_3_is_dt = dir(here("dgp_folder_weird"),"T-2", full.names = T)[grepl("(?=3.csv)",dir(here("dgp_folder_weird"),"T-2", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)

t1_1_os_dt = dir(here("dgp_folder_weird"),"T-1", full.names = T)[grepl("(?=1oos.csv)",dir(here("dgp_folder_weird"),"T-1", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t1_2_os_dt = dir(here("dgp_folder_weird"),"T-1", full.names = T)[grepl("(?=2oos.csv)",dir(here("dgp_folder_weird"),"T-1", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t1_3_os_dt = dir(here("dgp_folder_weird"),"T-1", full.names = T)[grepl("(?=3oos.csv)",dir(here("dgp_folder_weird"),"T-1", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t2_1_os_dt =dir(here("dgp_folder_weird"),"T-2", full.names = T)[grepl("(?=1oos.csv)",dir(here("dgp_folder_weird"),"T-2", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t2_2_os_dt =dir(here("dgp_folder_weird"),"T-2", full.names = T)[grepl("(?=2oos.csv)",dir(here("dgp_folder_weird"),"T-2", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)
t2_3_os_dt = t2_1_os_dt =dir(here("dgp_folder_weird"),"T-2", full.names = T)[grepl("(?=3oos.csv)",dir(here("dgp_folder_weird"),"T-2", full.names = T),perl=TRUE)]%>% 
  map_dfr(fread)


#t3_dt = here("dgp_folder_weird") %>% dir("t-3", full.names = T) %>% map_dfr(fread)
# Drop non-CV RF
f1_dt = f1_dt[str_detect(model, "Forest$", negate = T)]
t1_dt = t1_dt[str_detect(model, "Forest$", negate = T)]
t2_dt = t2_dt[str_detect(model, "Forest$", negate = T)]
t3_dt = t3_dt[str_detect(model, "Forest$", negate = T)]
# Change names
f1_dt[, model := str_replace_all(model, "First Stage", "First stage")]
t1_dt[, model := str_replace_all(model, "First Stage", "First stage")]
t2_dt[, model := str_replace_all(model, "First Stage", "First stage")]
t3_dt[, model := str_replace_all(model, "First Stage", "First stage")]
f1_dt[str_detect(model, "2SLS"), model := "First stage: OLS"]
t1_dt[str_detect(model, "2SLS"), model := "First stage: OLS"]
t2_dt[str_detect(model, "2SLS"), model := "First stage: OLS"]
t3_dt[str_detect(model, "2SLS"), model := "First stage: OLS"]
f1_dt[str_detect(model, "Naive OLS"), model := "Naive OLS"]
t1_dt[str_detect(model, "Naive OLS"), model := "Naive OLS"]
t2_dt[str_detect(model, "Naive OLS"), model := "Naive OLS"]
t3_dt[str_detect(model, "Naive OLS"), model := "Naive OLS"]
f1_dt[str_detect(model, "Oracle Model"), model := "\"Oracle\" model"]
t1_dt[str_detect(model, "Oracle Model"), model := "\"Oracle\" model"]
t2_dt[str_detect(model, "Oracle Model"), model := "\"Oracle\" model"]
t3_dt[str_detect(model, "Oracle Model"), model := "\"Oracle\" model"]
f1_dt[, model := str_replace_all(model, "Random Forest", "Random forest")]
t1_dt[, model := str_replace_all(model, "Random Forest", "Random forest")]
t2_dt[, model := str_replace_all(model, "Random Forest", "Random forest")]
t3_dt[, model := str_replace_all(model, "Random Forest", "Random forest")]
f1_dt[, model := str_replace_all(model, "LASSO", "Lasso")]
t1_dt[, model := str_replace_all(model, "LASSO", "Lasso")]
t2_dt[, model := str_replace_all(model, "LASSO", "Lasso")]
t3_dt[, model := str_replace_all(model, "LASSO", "Lasso")]
f1_dt[, model := str_replace_all(model, "post-", "Post-")]
t1_dt[, model := str_replace_all(model, "post-", "Post-")]
t2_dt[, model := str_replace_all(model, "post-", "Post-")]
t3_dt[, model := str_replace_all(model, "post-", "Post-")]
f1_dt[, model := str_replace_all(model, "Trees", "trees")]
t1_dt[, model := str_replace_all(model, "Trees", "trees")]
t2_dt[, model := str_replace_all(model, "Trees", "trees")]
t3_dt[, model := str_replace_all(model, "Trees", "trees")]
f1_dt[, model := str_replace_all(model, "First stage: USSIV", "Split-sample IV")]
t1_dt[, model := str_replace_all(model, "First stage: USSIV", "Split-sample IV")]
t2_dt[, model := str_replace_all(model, "First stage: USSIV", "Split-sample IV")]
t3_dt[, model := str_replace_all(model, "First stage: USSIV", "Split-sample IV")]
f1_dt[, model := str_replace_all(model, "First stage: JIVE", "Jackknife IV (JIVE)")]
t1_dt[, model := str_replace_all(model, "First stage: JIVE", "Jackknife IV (JIVE)")]
t2_dt[, model := str_replace_all(model, "First stage: JIVE", "Jackknife IV (JIVE)")]
t3_dt[, model := str_replace_all(model, "First stage: JIVE", "Jackknife IV (JIVE)")]
f1_dt[model == "First stage: LIML (Fuller)", model := "LIML (Fuller)"]
t1_dt[model == "First stage: LIML (Fuller)", model := "LIML (Fuller)"]
t2_dt[model == "First stage: LIML (Fuller)", model := "LIML (Fuller)"]
t3_dt[model == "First stage: LIML (Fuller)", model := "LIML (Fuller)"]
# Rename neural nets
f1_dt[str_detect(model, "Shallow") & !str_detect(model, "OOS"), model := "First stage: Neural net, shallow"]
t1_dt[str_detect(model, "Shallow") & !str_detect(model, "OOS"), model := "First stage: Neural net, shallow"]
t2_dt[str_detect(model, "Shallow") & !str_detect(model, "OOS"), model := "First stage: Neural net, shallow"]
t3_dt[str_detect(model, "Shallow") & !str_detect(model, "OOS"), model := "First stage: Neural net, shallow"]
f1_dt[str_detect(model, "Narrow") & !str_detect(model, "OOS"), model := "First stage: Neural net, narrow"]
t1_dt[str_detect(model, "Narrow") & !str_detect(model, "OOS"), model := "First stage: Neural net, narrow"]
t2_dt[str_detect(model, "Narrow") & !str_detect(model, "OOS"), model := "First stage: Neural net, narrow"]
t3_dt[str_detect(model, "Narrow") & !str_detect(model, "OOS"), model := "First stage: Neural net, narrow"]
f1_dt[str_detect(model, "Unrestricted") & !str_detect(model, "OOS"), model := "First stage: Neural net, unrestricted"]
t1_dt[str_detect(model, "Unrestricted") & !str_detect(model, "OOS"), model := "First stage: Neural net, unrestricted"]
t2_dt[str_detect(model, "Unrestricted") & !str_detect(model, "OOS"), model := "First stage: Neural net, unrestricted"]
t3_dt[str_detect(model, "Unrestricted") & !str_detect(model, "OOS"), model := "First stage: Neural net, unrestricted"]
# Rename neural nets, OOS
f1_dt[str_detect(model, "Shallow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, shallow"]
t1_dt[str_detect(model, "Shallow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, shallow"]
t2_dt[str_detect(model, "Shallow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, shallow"]
t3_dt[str_detect(model, "Shallow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, shallow"]
f1_dt[str_detect(model, "Narrow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, narrow"]
t1_dt[str_detect(model, "Narrow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, narrow"]
t2_dt[str_detect(model, "Narrow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, narrow"]
t3_dt[str_detect(model, "Narrow") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, narrow"]
f1_dt[str_detect(model, "Unrestricted") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, unrestricted"]
t1_dt[str_detect(model, "Unrestricted") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, unrestricted"]
t2_dt[str_detect(model, "Unrestricted") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, unrestricted"]
t3_dt[str_detect(model, "Unrestricted") & str_detect(model, "OOS"), model := "First stage: Neural net OOS, unrestricted"]
# Abbreviation for unrestricted
f1_dt[, model := str_replace_all(model, "unrestricted", "unrest.")]
t1_dt[, model := str_replace_all(model, "unrestricted", "unrest.")]
t2_dt[, model := str_replace_all(model, "unrestricted", "unrest.")]
t3_dt[, model := str_replace_all(model, "unrestricted", "unrest.")]
# Order the models
all_models = c(
  "First stage: Random forest, CV",
  "First stage: Boosted trees",
  "First stage: Neural net OOS, unrest.",
  "First stage: Neural net OOS, narrow",
  "First stage: Neural net OOS, shallow",
  "First stage: Neural net, unrest.",
  "First stage: Neural net, narrow",
  "First stage: Neural net, shallow",
  "First stage: Lasso selection",
  "First stage: Post-Lasso selection",
  "First stage: PCA",
  "Jackknife IV (JIVE)",
  "Split-sample IV",
  "LIML (Fuller)",
  "First stage: OLS",
  "Naive OLS",
  "\"Oracle\" model"
)
f1_dt[, model := factor(model, levels = all_models %>% rev(), ordered = T)]
t1_dt[, model := factor(model, levels = all_models %>% rev(), ordered = T)]
t2_dt[, model := factor(model, levels = all_models %>% rev(), ordered = T)]
t3_dt[, model := factor(model, levels = all_models %>% rev(), ordered = T)]
# Drop OOS neural net models
f1_dt %<>% .[str_detect(model, "OOS", negate = T)]
t1_dt %<>% .[str_detect(model, "OOS", negate = T)]
t2_dt %<>% .[str_detect(model, "OOS", negate = T)]
t3_dt %<>% .[str_detect(model, "OOS", negate = T)]


# Data work: Drop unwanted columns -------------------------------------------------------
f1_dt[, c("std.error", "statistic", "p.value") := NULL]
t1_dt[, c("std.error", "statistic", "p.value") := NULL]
t2_dt[, c("std.error", "statistic", "p.value") := NULL]
t3_dt[, c("std.error", "statistic", "p.value") := NULL]


# Data work: Covariances -----------------------------------------------------------------
# Rename covariance
setnames(f1_dt, old = "cov", new = "cov_xhat_u")
setnames(t1_dt, old = "cov", new = "cov_xhat_u")
setnames(t2_dt, old = "cov", new = "cov_xhat_u")
setnames(t3_dt, old = "cov", new = "cov_xhat_u")
# Transfer over 'cov_u' to 'cov_xhat_u'
f1_dt[is.na(cov_xhat_u), cov_xhat_u := cov_u]
t1_dt[is.na(cov_xhat_u), cov_xhat_u := cov_u]
t2_dt[is.na(cov_xhat_u), cov_xhat_u := cov_u]
t3_dt[is.na(cov_xhat_u), cov_xhat_u := cov_u]
# Drop 'cov_u'
f1_dt[, cov_u := NULL]
t1_dt[, cov_u := NULL]
t2_dt[, cov_u := NULL]
t3_dt[, cov_u := NULL]
# Calculate the covariance between e and 
beta = 1
f1_dt[, cov_xhat_e := (var * (estimate - beta) - cov_xhat_u) / beta]
t1_dt[, cov_xhat_e := (var * (estimate - beta) - cov_xhat_u) / beta]
t2_dt[, cov_xhat_e := (var * (estimate - beta) - cov_xhat_u) / beta]
t3_dt[, cov_xhat_e := (var * (estimate - beta) - cov_xhat_u) / beta]
# Drop 'cov_e' (missing for most models; calculations match output)
f1_dt[, cov_e := NULL]
t1_dt[, cov_e := NULL]
t2_dt[, cov_e := NULL]
t3_dt[, cov_e := NULL]
# Calculate the bias (error for beta)
f1_dt[, beta_bias := estimate - beta]
t1_dt[, beta_bias := estimate - beta]
t2_dt[, beta_bias := estimate - beta]
t3_dt[, beta_bias := estimate - beta]


# # Figures: Coefficient densities ---------------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d,
#         aes(x = estimate, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_vline(xintercept = 1, linetype = "longdash", size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab("Estimate") +
#       ylab("Density") +
#       scale_fill_viridis_d("Method", option = "magma", end = 0.95, direction = 1) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$estimate %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# Figures: Coefficient densities with 'ggridges' -----------------------------------------
# Create the individual figures
figure_list = lapply(
  X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
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
        vline_size = 0.8,
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
        base_size = 7.5,
        base_family = "Charter"
      ) +
      theme(legend.position = "none")
  }
)
# Save individual plots
for (i in seq_along(figure_list)) {
  ggsave(
    plot = figure_list[[i]],
    path = here("Figures"),
    filename = paste0("sim-density-ridge-", names(figure_list)[i], ".pdf"),
    device = cairo_pdf,
    width = 4.5,
    height = 2.5
  )
}


# # Two-plot figure (stacked): Comparing non-Belloni and Belloni
# two_gg = (
#   (figure_list$f1 + ggtitle("'Low-complexity' case: Seven strong instruments")) /
#   (figure_list$t2 + ggtitle("'High-complexity' case: 100 mixed-strength instruments"))
# ) +
# # plot_annotation(tag_levels = "A", tag_suffix = ".") +
# plot_layout(guides = "collect") & 
# theme(legend.position = "none") &
# guides(fill = guide_legend(ncol = 3, byrow = T))
# # Three-plot figure (stacked): Comparing Belloni
# three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
# plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
# plot_layout(guides = "collect") & 
# theme(legend.position = "none") &
# xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$estimate %>% range()) &
# guides(fill = guide_legend(ncol = 3, byrow = T))
# # Save the two-plot figure
# ggsave(
#   plot = two_gg,
#   path = here("Figures"),
#   filename = "sim-density-ridge-two.pdf",
#   device = cairo_pdf,
#   height = 7,
#   width = 6.5
# )
# # Save the three-plot figure
# ggsave(
#   plot = three_gg,
#   path = here("Figures"),
#   filename = "sim-density-ridge-three.pdf",
#   device = cairo_pdf,
#   height = 7,
#   width = 6.5
# )


# # Figures: Covariance densities with u ---------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d,
#         aes(x = cov_xhat_u, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_vline(xintercept = 0, linetype = "longdash", size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Cov($\\widehat{x}$,u)")) +
#       ylab("Density") +
#       scale_fill_manual(
#         "Method",
#         # values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#         values = magma(d[,model] %>% uniqueN(), end = 0.95, direction = 1)
#       ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$cov %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-u-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-u-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# # Figures: Covariance densities with e ---------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d,
#         aes(x = cov_xhat_e, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_vline(xintercept = 0, linetype = "longdash", size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Cov($\\widehat{x}$,e)")) +
#       ylab("Density") +
#       scale_fill_manual(
#         "Method",
#         values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#       ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme() +
#       coord_cartesian(ylim = c(0,40))
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$cov %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-cov-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# # Figures: Variance densities ------------------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d[str_detect(model, "Oracle", negate = T)],
#         aes(x = var, fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Var(\\widehat{x})")) +
#       ylab("Density") +
#       scale_fill_manual(
#         "Method",
#         values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#       ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$var %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )


# # Figures: Bias-share densities ----------------------------------------------------------
#   # Create the individual figures
#   figure_list = lapply(
#     X = list(f1 = f1_dt, t1 = t1_dt, t2 = t2_dt, t3 = t3_dt),
#     FUN = function(d) { 
#       ggplot(
#         data = d[(
#           str_detect(model, "Lasso|forest|trees") & str_detect(model, "Post", negate = T)
#         )],
#         aes(x = cov_xhat_e / (cov_xhat_e + beta * cov_xhat_u), fill = model)
#       ) +
#       geom_hline(yintercept = 0, size = 0.2) +
#       geom_density(alpha = 0.5, color = NA) +
#       xlab(TeX("Var(\\widehat{x})")) +
#       ylab("Density") +
#       # scale_fill_manual(
#       #   "Method",
#       #   values = viridis::magma(9, end = 0.95, direction = 1) %>% tail(8)
#       # ) +
#       theme_minimal(
#         base_size = 7.5,
#         base_family = "Merriweather"
#       ) +
#       theme()
#     }
#   )
#   # Two-plot figure (stacked): Comparing non-Belloni and Belloni
#   two_gg = (
#     (figure_list$f1 + ggtitle("'Easy' case: Seven strong instruments")) /
#     (figure_list$t2 + ggtitle("'Hard' case: 100 mixed-strength instruments (Belloni et al.)"))
#   ) +
#   # plot_annotation(tag_levels = "A", tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Three-plot figure (stacked): Comparing Belloni
#   three_gg = (figure_list$t1 / figure_list$t2 / figure_list$t3) +
#   plot_annotation(tag_levels = "1", tag_prefix = "B",  tag_suffix = ".") +
#   plot_layout(guides = "collect") & 
#   theme(legend.position = "bottom") &
#   xlim(rbindlist(list(t1_dt, t2_dt, t3_dt))$var %>% range()) &
#   guides(fill = guide_legend(ncol = 3, byrow = T))
#   # Save the two-plot figure
#   ggsave(
#     plot = two_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-two.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
#   # Save the three-plot figure
#   ggsave(
#     plot = three_gg,
#     path = here("Figures"),
#     filename = "sim-density-var-three.pdf",
#     device = cairo_pdf,
#     height = 7,
#     width = 6.5
#   )
