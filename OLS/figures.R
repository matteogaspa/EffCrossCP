#load data -----
rm(list = ls())
load("result_ols_n.RData")
library(ggplot2)
library(patchwork)
library(showtext)
showtext_auto()

ps    <- seq(5, 200, by = 5)
dim(res_cov); dim(res_len)
cov_p <- apply(res_cov, c(1,3), mean)
len_p <- apply(res_len, c(1,3), mean)
dim(cov_p); dim(len_p)

len_plot <- data.frame(
  p = ps,
  mod_cc = len_p[,1],
  e_mod_cc = len_p[,2],
  u_mod_cc = len_p[,3],
  eu_mod_cc = len_p[,4],
  cross = len_p[,5],
  split = len_p[,6],
  full = len_p[,7],
  jack = len_p[,8],
  splt2 = len_p[,9]
)

# plot 1 - set width + cov
p1 <- ggplot(len_plot, aes(x = p)) +
  geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = split, color = "split"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "gold2"),
    name = "") +
  labs(title = "Averaged prediction set size", x = "p", y = "Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(len_plot, aes(x = p)) +
  geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = split, color = "split"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "gold2"),
    name = "") +
  labs(title = "(zoomed in)", x = "p", y = "Size") + 
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

#combined_plot_size <- (p1 + p2) + plot_layout(guides = "collect") 

cov_plot <- data.frame(
  p = ps,
  mod_cc = cov_p[,1],
  e_mod_cc = cov_p[,2],
  u_mod_cc = cov_p[,3],
  eu_mod_cc = cov_p[,4],
  cross = cov_p[,5],
  split = cov_p[,6],
  full = cov_p[,7],
  jack = cov_p[,8],
  splt2 = cov_p[,9]
)

p3 <- ggplot(cov_plot, aes(x = p)) +
  geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = split, color = "split"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "gray") +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "gray") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "gold2"),
    name = "") +
  labs(title = "Empirical coverage", x = "p", y = "Coverage") + 
  coord_cartesian(ylim = c(0.70, 1)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))


combined_plot <- (p1 + p2 + p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 12))
combined_plot
ggsave("ols_plot1.pdf", combined_plot, width = 8, height = 5)

# plot 2 
p1n <- ggplot(len_plot, aes(x = p)) +
  geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = split, color = "split"), linewidth = 1.3, linetype = "solid") +
  #geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = splt2, color = "split (2α)"), linewidth = 1, linetype = "solid") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "gold2", "split (2α)" = "dodgerblue2"),
    name = "") +
  labs(title = "Averaged prediction set size", x = "p", y = "Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2n <- ggplot(len_plot, aes(x = p)) +
  geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = split, color = "split"), linewidth = 1.3, linetype = "solid") +
  #geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = splt2, color = "split (2α)"), linewidth = 1, linetype = "solid") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "gold2", "split (2\u03B1)" = "dodgerblue2"),
    name = "") +
  labs(title = "(zoomed in)", x = "p", y = "Size") + 
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

p3n <- ggplot(cov_plot, aes(x = p)) +
  geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = split, color = "split"), linewidth = 1.3, linetype = "solid") +
  #geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = splt2, color = "split (2α)"), linewidth = 1, linetype = "solid") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "gray") +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "gray") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "gold2", "split (2α)" = "dodgerblue2"),
    name = "") +
  labs(title = "Empirical coverage", x = "p", y = "Coverage") + 
  coord_cartesian(ylim = c(0.70, 1)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))


combined_plot2 <- (p1n + p2n + p3n) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 12))
combined_plot2
ggsave("ols_plot2.pdf", combined_plot2, width = 8, height = 5)



# plot 3 
p1b <- ggplot(len_plot, aes(x = p)) +
  #geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = split, color = "split"), linewidth = 1.3, linetype = "solid") +
  geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = splt2, color = "split (2α)"), linewidth = 1, linetype = "solid") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "forestgreen", "split (2α)" = "dodgerblue2"),
    name = "") +
  labs(title = "Averaged prediction set size", x = "p", y = "Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2b <- ggplot(len_plot, aes(x = p)) +
  #geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = split, color = "split"), linewidth = 1.3, linetype = "solid") +
  geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = splt2, color = "split (2α)"), linewidth = 1, linetype = "solid") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "forestgreen", "split (2\u03B1)" = "dodgerblue2"),
    name = "") +
  labs(title = "(zoomed in)", x = "p", y = "Size") + 
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

p3b <- ggplot(cov_plot, aes(x = p)) +
  #geom_line(aes(y = mod_cc, color = "mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = e_mod_cc, color = "e-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = u_mod_cc, color = "u-mod-cross"), linewidth = 1.25, linetype = "solid") +
  geom_line(aes(y = eu_mod_cc, color = "eu-mod-cross"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = cross, color = "cross"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = split, color = "split"), linewidth = 1.3, linetype = "solid") +
  geom_line(aes(y = full, color = "full"), linewidth = 1, linetype = "solid") +
  geom_line(aes(y = jack, color = "jackknife+"), linewidth = 1, linetype = "solid") +
  #geom_line(aes(y = splt2, color = "split (2α)"), linewidth = 1, linetype = "solid") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "gray") +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "gray") +
  scale_color_manual(
    values = c("mod-cross" = "black", "e-mod-cross" = "blue", "u-mod-cross" = "green3", "eu-mod-cross" = "red",
               "cross" = "orange", "split" = "pink2", "full" = "purple", "jackknife+" = "forestgreen", "split (2α)" = "dodgerblue2"),
    name = "") +
  labs(title = "Empirical coverage", x = "p", y = "Coverage") + 
  coord_cartesian(ylim = c(0.70, 1)) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))


combined_plot3 <- (p1b + p2b + p3b) + plot_layout(guides = "collect") 
combined_plot3
ggsave("ols_plot3.pdf", combined_plot3, width = 12, height = 5)

