

# Notes ----------------------------------------------------------------------------------
#   Goal: Create tables for results.


# Data notes -----------------------------------------------------------------------------
#   'true_var'        = variance of x1
#   'true_cov'        = covariance between u and x, if u = all of y != x
#   'var'             = variance of x1_hat
#   'cov' or 'cov_u'  = covariance of x1_hat and u (eventually cov_xhat_u)
#   'cov_e'           = covariance of x1_hat and e (eventually cov_xhat_e)


# Setup ----------------------------------------------------------------------------------
  # Load packages
  library(pacman)
  p_load(plotly, tidyverse, scales, huxtable, data.table, magrittr, here)


# Load all results -----------------------------------------------------------------------
  # Read in the results, grouping DGPs into separate objects
  f1_dt = here("Resources") %>% dir("f-1", full.names = T) %>% map_dfr(fread)
  t1_dt = here("Resources") %>% dir("t-1", full.names = T) %>% map_dfr(fread)
  t2_dt = here("Resources") %>% dir("t-2", full.names = T) %>% map_dfr(fread)
  t3_dt = here("Resources") %>% dir("t-3", full.names = T) %>% map_dfr(fread)
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
  f1_dt[str_detect(model, "Oracle Model"), model := "``Oracle\'\' model"]
  t1_dt[str_detect(model, "Oracle Model"), model := "``Oracle\'\' model"]
  t2_dt[str_detect(model, "Oracle Model"), model := "``Oracle\'\' model"]
  t3_dt[str_detect(model, "Oracle Model"), model := "``Oracle\'\' model"]
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
  

# Merge results --------------------------------------------------------------------------
  # Join the results together
  sim_dt = rbindlist(list(
  	"Non-Belloni" = f1_dt,
  	"Belloni, Type 1" = t1_dt,
  	"Belloni, Type 2" = t2_dt,
  	"Belloni, Type 3" = t3_dt
  ), idcol = T)
  # Summarize models
  sum_dt = sim_dt[, .(
  	mean = mean(estimate),
  	sd = sd(estimate)
  ), by = .(.id, model)]
	# Order rows
	setorder(sum_dt, -.id, model)


# Table: Main summary table --------------------------------------------------------------
	# Define the caption
	tmp_caption = "This table provides the means and standard deviations of the distributions illustrated in Figure~\\ref{fig:results-densities}. Each \\textbf{column} contains a separate DGP: (a) contains the \\textit{low-complexity} DGP with 7 (equally) strong instruments; (b)--(d) contain the three \\textit{high-complexity} cases with 100 instruments of mixed strengths. \\textbf{Rows} differ by estimator. For each DGP-estimator combination, we summarize the estimates for the parameter of interest $(\\beta)$ across 1,000 iterations using a mean and standard deviation (the standard deviation is in parentheses)."
	# Find the unique number of methods
	nm = sum_dt[,model] %>% uniqueN()
  # Start the matrix
  sum_mat = matrix(
	  data = c(rbind(
	  	comma(sum_dt[,mean], accuracy = 0.001),
	  	comma(sum_dt[,sd], accuracy = 0.001) %>% paste0("(", ., ")")
	  )),
	  ncol = 4,
	  byrow = F
	)
  # Add method names
  sum_mat %<>% cbind(c(rbind(sum_dt[1:nm,model %>% as.character()], rep("", nm))), .)
  # Add spacing
  sum_mat = cbind(sum_mat[,1], "&", sum_mat[,2], "&", sum_mat[,3], "&", sum_mat[,4], "&", sum_mat[,5])
  sum_mat %<>% cbind(
  	.,
  	c(rbind(rep("\\\\[-0.1em]", nm), rep("\\\\[0.3em]", nm)))
  )
  # Convert to character vector for TeX
  sum_tex = sum_mat %>% knitr::kable(format = "pandoc") %>%
    as.character() %>% tail(-1) %>% head(-1)
  # Add top part
  sum_tex %<>% c(
  	"% NOTE: Do not edit. This file will be overwritten.",
  	"% This file is generated by 'codeR/table-results.R'",
  	"",
  	"\\begin{table}[!htbp] \\centering",
  	"  \\caption{\\textbf{Simulation results}: Mean and standard deviation for methods and DGPs}",
  	"  \\label{tab:main-results}",
  	paste(c("\\begin{tabular}{@{\\extracolsep{4pt}}l", rep("c", 4), "@{}}"), collapse = ""),
  	"	 \\toprule",
    "  \\ra{1.5}",
    paste0(
  	  " & ",
      "\\textit{Low-complexity} case & ",
      "\\multicolumn{3}{c}{\\textit{High-complexity} cases} \\\\ ",
      "\\cmidrule(lr){2-2} \\cmidrule(lr){3-5}"
    ),
    paste0(
      "  & ",
      "7 strong instruments & ",
      "\\multicolumn{3}{c}{100 \\textit{mixed} instruments} \\\\ "
    ),
  	"  & ($A$) & ($B$) & ($C$) & ($D$) \\\\",
  	"  \\midrule",
  	"",
  	.
  )
  # Add bottom part
  sum_tex %<>% c(
  	"",
    "\\bottomrule",
    "\\\\[-0.9em]",
    "\\end{tabular}",
    paste0("\\caption*{\\raggedright\\small ",  tmp_caption," }"),
    "\\end{table}"
  )	
  # Capture to write
  sum_tex = capture.output(expr = cat(paste0(sum_tex, "\n")))
  # Write to file
  write(
  	sum_tex,
  	here("Tables", "main-results-table.tex")
  )


# Table: Bias decomposition --------------------------------------------------------------
  # Define the caption
  tmp_caption = "A cell's value provides the given statistic's mean (\\textbf{column}) in 1,000 iterations of the given combination of DGP (\\textbf{Panel}) and estimator (\\textbf{row}). We omit LIML as it is not a two-stage method and thus does not produce $\\hat{x}$."
  # Join the results together
  bias_dt = rbindlist(list(
    "Non-Belloni" = f1_dt[str_detect(model, "LIML", negate = T)],
    "Belloni, Type 1" = t1_dt[str_detect(model, "LIML", negate = T)],
    "Belloni, Type 2" = t2_dt[str_detect(model, "LIML", negate = T)],
    "Belloni, Type 3" = t3_dt[str_detect(model, "LIML", negate = T)]
  ), idcol = "dgp")
  # Drop unwanted columns
  bias_dt[, c("std.error", "statistic", "p.value") := NULL]
  # Calculate bias components
  bias_dt[, `:=`(b = cov_xhat_u, c = 1 / var)]
  bias_dt[, `:=`(a = (estimate - 1) / c - b)]
  bias_dt[, cov_xxh := a + var]
  # Summarize bias components by model and DGP
  bias_dt %<>% .[, .(
    "N" = .N %>% scales::comma(),
    "Bias" = mean(estimate - 1) %>% round(2),
    "a=Cov(xh,e)" = mean(a) %>% round(2),
    "b=Cov(xh,u)" = mean(b) %>% round(2),
    "c=1/Var(xh)" = mean(c) %>% round(2),
    "Var(xh)" = mean(var) %>% round(2),
    "Var(x)" = mean(true_var) %>% round(2),
    "Cov(x,xh)" = mean(cov_xxh) %>% round(2),
    "Corr(x,xh)" = mean(cov_xxh / sqrt(var * true_var)) %>% round(2),
    "Cov(x,u)" = mean(true_cov) %>% round(2)
  ), by = .(dgp, model)]
  # Order rows
  setorder(bias_dt, dgp, model)
  # The DGPs
  dgp_v = c(
    "Non-Belloni",
    "Belloni, Type 1",
    "Belloni, Type 2",
    "Belloni, Type 3"
  )
  dgp_names = c(
    "\\textit{Low-complexity} case",
    "\\textit{High-complexity} case 1",
    "\\textit{High-complexity} case 2",
    "\\textit{High-complexity} case 3"
  )
  # Create a matrix of results for each DGP
  bias_tex = lapply(
    X = seq_along(dgp_v),
    FUN = function(i) {
      # Grab subset
      m_dt = bias_dt[dgp == dgp_v[i], -c("dgp", "N")]
      # Add "&" column separators
      m_mat = m_dt[[1]] %>% paste0("\\enspace ", .) %>% as.matrix()
      for (j in 2:ncol(m_dt)) {
        m_mat %<>% cbind("&") %>% cbind(m_dt[[j]])
      }
      # Add column endings
      m_mat %<>% cbind("\\\\")
      # Convert to character
      m_chr = m_mat %>% knitr::kable(format = "pandoc") %>% as.character() %>% head(-1) %>% tail(-1)
      # Add midrule and DGP name to the top
      m_chr %<>% c(
        "\\midrule",
        paste0("\\multicolumn{3}{@{}l}{\\textbf{Panel ", LETTERS[i], " DGP:} ", dgp_names[i], "} \\\\"),
        .
      )
    }
  ) %>% unlist()
  # Column names
  col_names = c(
    "\\enspace \\textbf{Estimator}",
    "Bias",
    "$\\text{Cov}(\\hat{x},e)$",
    "$\\text{Cov}(\\hat{x},u)$",
    "$1/\\text{Var}(\\hat{x})$",
    "$\\text{Var}(\\hat{x})$",
    "$\\text{Var}(x)$",
    "$\\text{Cov}(x,\\hat{x})$",
    "$\\text{Corr}(x,\\hat{x})$",
    "$\\text{Cov}(x,u)$"
  )
  # Header, column names
  bias_tex %<>% c(
    " & & \\multicolumn{3}{c}{\\textbf{Bias components}} \\\\ \\cmidrule(lr){3-5}", 
    " & $(a+b)c$ & $a$ & $b$ & $c$ \\\\",
    paste(col_names, collapse = " & ") %>% paste("\\\\"),
    .
  )
  # Add other top table 'stuff'
  bias_tex %<>% c(
    "% NOTE: Do not edit. This file will be overwritten.",
    "% This file is generated by 'codeR/table-results.R'",
    "",
    "\\begin{table}[!htbp] \\centering",
    "  \\caption{\\textbf{Simulation results}: Decomposing bias components (means from simulation)}",
    "  \\label{tab:bias-decomp}",
    "\\resizebox{0.9\\columnwidth}{!}{%",
    paste(
      c("\\begin{tabular}{@{\\extracolsep{-5pt}}l", rep("c", length(col_names) - 1), "@{}}"), 
      collapse = ""
    ),
    "  \\toprule",
    "  \\ra{1.5}",
    "",
    .
  )
  # Add bottom 'stuff'
  bias_tex %<>% c(
    "",
    "\\bottomrule",
    "\\\\[-0.9em]",
    "\\end{tabular}%",
    "}",
    paste0("\\caption*{\\raggedright\\small ",  tmp_caption," }"),
    "\\end{table}"
  ) 
  # Capture to write
  bias_tex = capture.output(expr = cat(paste0(bias_tex, "\n")))
  # Write to file
  write(
    bias_tex,
    here("Tables", "bias-decomp-table.tex")
  )


