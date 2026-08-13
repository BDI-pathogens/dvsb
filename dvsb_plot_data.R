# Plot y expected by plate
ggplot(df_cal %>% 
         filter(cal == 1)) +
  geom_line(aes(xlog, y_mean, group = label, col = label)) +
  labs(x = "x = Ab concentration (log)",
       y = "y = ELISA OD",
       col = "") +
  coord_cartesian(expand = F) +
  ylim(0, NA)

# Plot observed y
ggplot(df_cal) +
  geom_point(aes(jitter(x), y, col = label)) +
  scale_x_log10(breaks = xs) +
  labs(x = "Concentration",
       y = "OD",
       col = "")
ggplot(df_cal %>% filter(plate <= 25)) +
  geom_point(aes(x, y)) +
  scale_x_log10(breaks = xs) +
  labs(x = "Concentration",
       y = "OD",
       col = "plate") +
  facet_wrap(~label)

# Calculate indepent ML 4PL curves per plate and add to the plot
if (FALSE) {
  library(dr4pl)
  df_ml <- df_cal %>%
    split(.$plate) %>%
    map(function(df_) {
      fit <- dr4pl(data = df_, dose = x, response = y)$parameters
      tibble(param = names(fit), value = as.numeric(fit))
    }) %>%
    bind_rows(.id = "plate") %>%
    mutate(plate = as.integer(plate)) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    inner_join(df_plate, by = "plate") %>%
    expand_grid(xlog = seq(min(xlogs), max(xlogs), length.out = 100)) %>%
    mutate(x = exp(xlog)) %>%
    mutate(`independent\nmax-likelihood` = theta_1 + (theta_4 - theta_1) / (1 + (x / theta_2)^theta_3),
           truth = PL4(xlog, f_1, f_2, f_3, f_4)) %>%
    pivot_longer(c("truth", "independent\nmax-likelihood"), names_to = "y", values_to = "value") %>%
    mutate(y = factor(y, levels = c("truth", "independent\nmax-likelihood")))
  ggplot() +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "plate",
         linetype = "") +
    geom_line(data = df_ml %>% filter(y == "truth"),
              aes(x, value, col = as.factor(plate))) 
  ggsave("~/foo_1.pdf", height = 5.5, width = 6)
  ggplot(df_cal) +
    geom_point(aes(x, y, col = as.factor(plate))) +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "plate",
         linetype = "") +
    facet_wrap(~plate) +
    geom_line(data = df_ml %>% filter(y == "truth"),
              aes(x, value, col = as.factor(plate))) 
  ggsave("~/foo_2.pdf", height = 8, width = 3.2)
  ggplot(df_cal %>% filter(plate == 18)) +
    geom_point(aes(x, y)) +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD") +
    geom_line(data = df_ml %>% filter(y == "truth") %>% filter(plate == 18),
              aes(x, value)) 
  ggsave("~/foo_2b.pdf", height = 4, width = 4)
  ggplot(df_cal) +
    geom_point(aes(x, y, col = as.factor(plate))) +
    scale_x_log10(breaks = xs) +
    labs(x = "x = Ab concentration",
         y = "y = OD",
         col = "plate",
         linetype = "") +
    facet_wrap(~plate) +
    geom_line(data = df_ml, aes(x, value, col = as.factor(plate), linetype = y)) 
  ggsave("~/foo_3.pdf", height = 8, width = 4)
}

ggplot(df_sam %>% mutate(x_mix_group = paste("x mix group =", x_mix_group))) +
  geom_histogram(aes(x), fill = "grey", bins = 30) +
  scale_x_log10(limits = c(NA, NA)) +
  coord_cartesian(expand = F) +
  labs(x = "x = Ab concentration",
       y = "Number of samples") +
  facet_wrap(~x_mix_group, ncol = 8) +
  geom_vline(xintercept = exp(mu_pos)) +
  geom_vline(xintercept = exp(mu_neg))
if (FALSE) {
  ggsave("~/foo_4.pdf", height = 4, width = 5)
}


ggplot(df_sam %>% mutate(x_mix_group = paste("x mix group =", x_mix_group))) +
  geom_histogram(aes(y), fill = "grey", bins = 30) +
  scale_x_log10(limits = c(NA, NA)) +
  coord_cartesian(expand = F) +
  labs(x = "x = Ab concentration",
       y = "Number of samples") +
  facet_wrap(~x_mix_group, ncol = 8) 

# Plot calibrators and samples by plate
if (FALSE) {
  bind_rows(df_cal %>% mutate(type = "calibrator") ,
            df_sam %>% mutate(type = if_else(pos, "+ sample", "- sample"))) %>%
    ggplot() +
    geom_point(aes(jitter(x), y, col = type)) +
    scale_x_log10(breaks = xs) +
    facet_wrap(~label) +
    labs(x = "Concentration",
         y = "Observed optical density",
         col = "") +
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)) 
  ggsave("~/lassa_serology_cross-sectional_data.pdf", height = 9, width = 12)  
}