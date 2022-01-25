
# Data notes -----------------------------------------------------------------------------
#   - Suffix '.x' the observed data (from the given iteration)
#   - Suffix '.y' comes from the 25 iterations in which we ran all methods (to get ranks)

# Setup ----------------------------------------------------------------------------------
  # Packages
  library(pacman)
  p_load(tidyverse, latex2exp, scales, viridis, data.table, magrittr, here)

# Load data ------------------------------------------------------------------------------  
  cv_dt = here("traininfo", "training_info_nnet_shuf_2_complexT_nnetd.csv") %>% fread()
  # cv_dt = here("traininfo", "training_info_nnet_shuf_3_complexT_nnetd.csv") %>% fread()

# Data work: Find and summarize winners --------------------------------------------------
  # Grab desired variables
  cv_dt %<>% .[, .(
    iter,
    est = beta,
    rank = rank_topx,
    n_param = num_parameters,
    width = flag_width1,
    depth = flag_depth,
    dropout = flag_dropout,
    loss = mean_eval_loss.x,
    loss_avg = mean_eval_loss.y
  )]
  # Order by iteration and loss
  setorder(cv_dt, iter, loss)
  cv_dt[, out := rep(c("winner", "loser"), times = .N/2)]
  # Long to wide
  cv_dt %<>% dcast(
    iter + est ~ out,
    value.var = c("rank", "n_param", "width", "depth", "dropout", "loss", "loss_avg")
  )

# Figure: Bias and model complexity ------------------------------------------------------
  gg_tmp = ggplot(
    data = cv_dt,
    aes(y = est, color = est < 1.2)
  ) +
  geom_segment(
    aes(x = depth_winner, xend = depth_loser, yend = est),
    size = 1/11
  ) +
  geom_point(
    aes(x = depth_winner),
    shape = 19,
    size = 1/2
  ) +
  scale_y_continuous("Estimate") +
  # scale_x_log10("Number of parameters", labels = scales::comma) +
  scale_x_continuous("Depth", breaks = 1:6) +
  scale_color_manual(values = magma(100)[c(15,70)]) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor.x = element_blank()
  )
  # Save
  ggsave(
    plot = gg_tmp,
    path = here("Figures"),
    filename = "nn-depth-bias-t2.pdf",
    device = cairo_pdf,
    width = 3.25,
    height = 7
  )


# Figure: Bias and loss ------------------------------------------------------------------
  gg_tmp = ggplot(
    data = cv_dt,
    aes(x = loss_winner, y = est, color = est < 1.2)
  ) +
  geom_hline(
    yintercept = 1,
    size = 1/3,
    color = "grey40",
    linetype = "solid"
    ) +
  geom_hline(
    yintercept = cv_dt[est < 1.2, mean(est)],
    size = 1/4,
    color = magma(100)[70],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = cv_dt[est >= 1.2, mean(est)],
    size = 1/4,
    color = magma(100)[15],
    linetype = "dashed"
  ) +
  geom_point(
    shape = 19,
    size = 1/4,
  ) +
  scale_y_continuous("Estimate") +
  scale_x_continuous("Loss (MSE)", labels = comma) +
  scale_color_manual(values = magma(100)[c(15,70)]) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  ) +
  theme(legend.position = "none")
  # Save
  ggsave(
    plot = gg_tmp,
    path = here("Figures"),
    filename = "nn-depth-loss-t2.pdf",
    device = cairo_pdf,
    width = 3.25,
    height = 3.25
  )


# Figure: Probability of choosing shallower option ---------------------------------------
  # Find winner depth, min depth, and max depth for each iteration
  depth_dt = cv_dt[, .(
    depth_winner,
    depth_max = max(depth_winner, depth_loser),
    depth_min = min(depth_winner, depth_loser)
  ), by = iter]
  # Find the probability that CV chose the shallower option for each pair
  depth_dt %<>% .[, .(
    prob_shallower = mean(depth_winner == depth_min),
    n_shallower = sum(depth_winner == depth_min),
    .N
  ), by = .(depth_min, depth_max)]
  # Plot
  gg_tmp = ggplot(
    data = depth_dt,
    aes(x = depth_min, y = depth_max, fill = prob_shallower)
  ) +
  geom_tile() +
  geom_text(
    aes(
      # label = paste0("frac(", n_shallower, ",", N, ")"),
      label = paste0(n_shallower, "/", N),
      color = (prob_shallower > 0.75)
    ),
    # parse = T,
    family = "Charter",
    size = 1.5,
    show.legend = F
  ) +
  scale_x_continuous("Minimum depth option", breaks = 1:6) +
  scale_y_continuous("Maximum depth option", breaks = 1:6) +
  scale_fill_viridis(
    "Probability CV chose the more shallow depth",
    option = "magma",
    end = 0.95
  ) +
  scale_color_manual(values = c("white", "black")) +
  theme_minimal(
    base_size = 7.5,
    base_family = "Charter"
  ) +
  theme(
    legend.position = "bottom",
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.bottom = element_text(margin = margin(t = -5)),
    axis.text.y = element_text(margin = margin(r = -5)),
    plot.margin = grid::unit(c(0, 0, 0, 0), "in")
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 12,
      barheight = 0.75
    )
  ) +
  coord_equal()
  # Save
  ggsave(
    plot = gg_tmp,
    path = here("Figures"),
    filename = "nn-depth-choice-t2.pdf",
    device = cairo_pdf,
    width = 3.25,
    height = 3.25
  )
