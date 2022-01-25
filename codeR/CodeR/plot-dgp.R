
# Notes ----------------------------------------------------------------------------------
#   Goal: Plot DGPs, i.e., build correlation and beta-strength diagrams

# Setup ----------------------------------------------------------------------------------
  # Packages
  library(pacman)
  p_load(
    tidyverse, scales, viridis, latex2exp, extrafont,
    MASS, data.table, here, magrittr
  )

# Generate data: Simple case -------------------------------------------------------------
  # Set seed
  seed = 12345
  set.seed(seed)
  # Number of observations and instruments
  n = 1e3
  k = 7
  # Simple DGP
  source(here("codeR", "simple_data_function.R"))
  i_dt_corr = data_gen(n,k,shuffle)
  betas = rep(1, 7)
  # Names
  z_names = names(i_dt_corr)[grepl("z", names(i_dt_corr))]
  # Dataset of coefficients
  betadata = data.frame(z_names = as.numeric(stringr::str_extract(z_names, "\\d+")), beta = betas)
  # Correlation values
  i_dt_corr[, corr_mat := .(list(cor(.SD))), .SDcols = z_names]
  # Data table of correlation matrix
  corr_mat_dt = i_dt_corr$corr_mat[[1]] %>% data.table()
  corr_mat_dt$corrwith = rownames(corr_mat_dt)
  corr_mat_dt = melt(corr_mat_dt, id.vars = c('corrwith'))
  corr_mat_dt %<>%
    mutate(value = ifelse(is.na(value), 0, value)) %>%
    mutate(value = squish(corr_mat_dt$value, range = c(.00000000001, 1)))
  corr_mat_dt %<>% mutate(
    corrwith = as.numeric(stringr::str_extract(corrwith, "\\d+")),
    variable =  as.numeric(stringr::str_extract(variable, "\\d+")),
    corr = value
# NOTE: Replaced to give full matrix (below gives upper triangle)
    # corr = ifelse(corrwith > variable, value, NA)
  )

# Figure: Coefficient plot (simple) ------------------------------------------------------
  # Coefficient plot (simple case)
  plot_coef_simple = ggplot(
    data = betadata,
    aes(x = z_names, y = beta)
  ) +
  geom_bar(stat = "identity", fill = "grey40") +
  geom_hline(yintercept = 0, size = 1/8) +
  scale_x_continuous(
    TeX("Instrument number ($z_{i}$)"),
    breaks = 1:9
  ) +
  scale_y_continuous(
    TeX("Coefficient ($\\pi_{i}$)"),
    breaks = seq(0, 1, 1)
  ) +
  theme_minimal(base_size = 7.5, base_family = "Merriweather") +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.bottom = element_text(margin = margin(t = -10)),
    axis.text.y = element_text(margin = margin(r = -1)),
    plot.margin = grid::unit(c(0, 0, 0, 0), "in")
  )

# Figure: Correlation matrix (simple) ----------------------------------------------------
  # Correlation plot (simple case)
  plot_corr_matrix_simple = ggplot(
    data = corr_mat_dt,
    aes(x = corrwith, y = variable, fill = corr)
  ) +
  geom_tile() +
  scale_x_continuous(
    TeX("Instrument number ($z_{i}$)"),
    position = "top",
    breaks = 1:9
  ) +
  scale_y_reverse(
    TeX("Instrument number ($z_{i}$)"),
    position = "left",
    breaks = 1:9
  ) +
  theme_minimal(base_size = 7.5, base_family = "Merriweather") +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(margin = margin(b = -10)),
    axis.text.y = element_text(margin = margin(r = -10)),
    plot.margin = grid::unit(c(0, 0, 0, 0), "in")
  ) +
  scale_fill_viridis(
    "Correlation",
    option = "magma",
    na.value = "white",
    begin = 0,
    direction = -1,
    breaks = seq(0, 1, 0.25)
  ) +
  coord_fixed()

# Copmlex case: Generate data and coefficient plot ---------------------------------------
  # Sizes: Sample and instruments
  n = 1000
  k = 100
  # Create empty list for coefficient plot
  coef_plot_list = list()
  # Iterate over the 'shuffle' parameter
  # '1' - shuffle beta pattern
  # '2' - set strongest beta at z_1 (Belloni reproduction)
  # '3' - set strongest beta at z_50, descends from there
  for (shuffle in c("1", "2", "3")) {
    set_C = function(musq=musq, beta_pattern=beta_pattern, Sig_z=Sig_z, n=n) {
      #solving for a C value that changes behavior of lasso selectors.
      C = sqrt(musq/((n+musq)*(crossprod(beta_pattern,Sig_z) %*% beta_pattern)))
      return(C)
    }
    # Generate two errors according to Bellloni (pg 26)
    # First need covariance and variance set, to produce covar matrix
    musq = 1250
    var_er = 1
    cor_er = 0.5
    var_inst = 1
    # Produce Sigma Z. variance is 1, correlation (thus cov) is equal to .5^(|j-h|) for z_j and z_h
    Sig_z = matrix(0,nrow = k,ncol = k)
    for (j in (1:k)) {
      for (h in (1:k)) {
        Sig_z[c(j),c(h)] = .5^abs(j-h)
      }
    }
    z_n1 = mvrnorm(n = n, mu = rep(0,k), Sigma = Sig_z, empirical = TRUE)
    # Use exponential pattern
    beta_pattern = unlist(lapply(c(0:(k-1)), function(x) (.7)^x), use.names = FALSE)
    Sig_z = as.matrix(Sig_z)
    if(shuffle == '1'){
      beta_pattern <- as.vector(beta_pattern)
      beta_pattern = sample(beta_pattern)
    } else if(shuffle == '2'){
      beta_pattern <- as.vector(beta_pattern)
    } else{
      beta_pattern <- as.vector(beta_pattern)
      beta_pattern = c(tail(beta_pattern, -50), beta_pattern[1:50])
    }
    C = set_C(musq=musq, n = n, beta_pattern = beta_pattern,Sig_z = Sig_z)
    C = c(C)
    # Now that C is known, set variance for error term on X
    sigma_v = abs(1 - C*crossprod(beta_pattern,Sig_z) %*% beta_pattern)
    # We can also set our full pi matrix to find our true instrument coefficients
    betas = as.vector(C*beta_pattern)
    cov_er = cor_er*sqrt(sigma_v)
    # Use multivariate distribution given matrix above to get a matrix of values 
    # for sample size of 'n1'
    covarmat = matrix(c(1,cov_er,cov_er, sigma_v), nrow = 2, ncol = 2)
    # We now need to get and separate our error terms
    evmat_n1 = mvrnorm(n = n, mu = c(0,0), Sigma = covarmat, empirical = TRUE)
    errors = as.data.frame(evmat_n1)
    names(errors) <- c("e", "v")
    # Separate errors
    en1 = evmat_n1[,1]
    vn1 = evmat_n1[,2]
    # Construct x's and y's and add them to our z's to create a full data matrix
    X = z_n1 %*% betas + vn1
    # x created
    # x_true
    x_true = X - vn1
    Y = X + en1 + 1
    # y created
    znames = paste(rep("z", k), c(1:k), sep = "_")
    z = as.data.frame(z_n1)
    colnames(z) <- znames
    # X (or d, as in Belloni paper) name
    X = as.data.frame(X)
    colnames(X) <- "X"
    # Target value Y name
    Y = as.data.frame(Y, names = "Y")
    colnames(Y) <- "Y"
    # final matrix
    # true_x = X - vn1
    product<-cbind(z,X[,1],Y,errors, x_true)
    product = data.table(product)
    product = setnames(product, c(znames, "x1", "y", "en1", 'vn1', 'true_x'))
    product[,id := c(1:n)]
    i_dt_corr = product
    z_names = names(i_dt_corr)[grepl("z", names(i_dt_corr))]
    betadata = data.frame(
      z_names = as.numeric(stringr::str_extract(z_names, "\\d+")),
      beta = betas
    )
    i_dt_corr[, corr_mat := .(list(cor(.SD))), .SDcols = z_names]
    corr_mat_dt = i_dt_corr$corr_mat[[1]] %>% data.table()
    corr_mat_dt$corrwith = rownames(corr_mat_dt)
    corr_mat_dt = melt(corr_mat_dt, id.vars = c('corrwith'))
    corr_mat_dt %<>%
      mutate(value = ifelse(is.na(value), 0, value)) %>%
      mutate(value = squish(corr_mat_dt$value, range = c(.00000000001, 1)))
    corr_mat_dt = corr_mat_dt %>% mutate(
      corrwith = as.numeric(stringr::str_extract(corrwith, "\\d+")),
      variable =  as.numeric(stringr::str_extract(variable, "\\d+")),
      corr = value
      # corr = ifelse(corrwith > variable, value, NA)
    )
    # Set color using value of 'shuffle'
    fill_col = magma(3, begin = 0.4, end = 0.8)[as.numeric(shuffle)]
    # Figure: Coefficient plot (complex case)
    plot_coef = ggplot(
      data = betadata,
      aes(x = z_names, y = beta)
    ) +
    geom_bar(stat = "identity", fill = fill_col) +
    geom_hline(yintercept = 0, size = 1/8) +
    scale_x_continuous(
      TeX("Instrument number ($z_{i}$)"),
      breaks = c(1, 25, 50, 75, 100)
    ) +
    scale_y_continuous(
      TeX("Coefficient ($\\pi_{i}$)"),
      breaks = c(0, 0.25, 0.5),
      limits = c(0, 0.53)
    ) +
    theme_minimal(base_size = 7.5, base_family = "Merriweather") +
    theme(
      # legend.position = "bottom",
      legend.position = "none",
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.text.x.bottom = element_text(margin = margin(t = -10)),
      axis.text.y = element_text(margin = margin(r = -1)),
      plot.margin = grid::unit(c(0, 0, 0, 0), "in")
    )
    # Add the plot to the coef_plot_list
    coef_plot_list[[as.numeric(shuffle)]] = plot_coef
  }
    
# Figure: Correlation matrix (complex) ---------------------------------------------------
  # Correlation matrix
  plot_corr_matrix = ggplot(
    data = corr_mat_dt,
    aes(x = corrwith, y = variable, fill = corr)
  ) + 
  geom_tile() +
  scale_x_continuous(
    TeX("Instrument number ($z_{i}$)"),
    position = "top",
    breaks = c(1, 25, 50, 75, 100)
  ) +
  scale_y_reverse(
    TeX("Instrument number ($z_{i}$)"),
    position = "left",
    breaks = c(1, 25, 50, 75, 100)
  ) +
  theme_minimal(base_size = 7.5, base_family = "Merriweather") +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(margin = margin(b = -10)),
    axis.text.y = element_text(margin = margin(r = -10)),
    plot.margin = grid::unit(c(0, 0, 0, 0), "in")
  ) +
  scale_fill_viridis(
    "Correlation",
    option = "magma",
    na.value = "white",
    begin = 0,
    direction = -1,
    breaks = seq(0, 1, 0.25)
  ) +
  coord_fixed()

# Figure: Correlation line (complex) -----------------------------------------------------
  # Correlation plot
  plot_corr = ggplot(
    data = corr_mat_dt %>% filter(variable %in% c(1, 50)),
    aes(
      x = corrwith,
      y = value,
      color = as.character(variable),
      shape = as.character(variable)
    )
  ) +
  geom_hline(yintercept = 0, size = 1/8) +
  geom_line() +
  geom_point() +
  scale_color_manual(
    "Instrument",
    values = magma(2, begin = 0.25, end = 0.85, direction = -1),
    labels = c(expression(z[1]), expression(z[50]))
  ) +
  scale_shape_manual(
    "Instrument",
    values = c(19, 4),
    labels = c(expression(z[1]), expression(z[50]))
  ) +
  scale_x_continuous(
    TeX("Instrument number ($z_i$)"),
    limits = c(0.5, 100.5),
    expand = c(0,0)
  ) +
  scale_y_continuous("Correlation") +
  theme_minimal(base_size = 7.5, base_family = "Merriweather") +
  theme(
    legend.position = "bottom",
    # plot.margin = grid::unit(c(0, 0, 0, 0), "in")
  )

# Save figures ---------------------------------------------------------------------------
  # Simple: Coefficients
  # Simple: Correlation matrix
  # Complex: Coefficients 1
  ggsave(
    plot = coef_plot_list[[1]],
    path = here("Figures"),
    filename = "dgp-complex1-coef.pdf",
    device = cairo_pdf,
    height = 2.5,
    width = 3.25
  )
  # Complex: Coefficients 2
  ggsave(
    plot = coef_plot_list[[2]],
    path = here("Figures"),
    filename = "dgp-complex2-coef.pdf",
    device = cairo_pdf,
    height = 2.5,
    width = 3.25
  )
  # Complex: Coefficients 3
  ggsave(
    plot = coef_plot_list[[3]],
    path = here("Figures"),
    filename = "dgp-complex3-coef.pdf",
    device = cairo_pdf,
    height = 2.5,
    width = 3.25
  )
  # Complex: Correlation matrix
  ggsave(
    plot = plot_corr_matrix,
    path = here("Figures"),
    filename = "dgp-complex-corr-mat.pdf",
    device = cairo_pdf,
    height = 2.5,
    width = 3.25
  )
  # Complex: Correlation lines
  ggsave(
    plot = plot_corr,
    path = here("Figures"),
    filename = "dgp-complex-corr-line.pdf",
    device = cairo_pdf,
    height = 2.5,
    width = 3.25
  )
