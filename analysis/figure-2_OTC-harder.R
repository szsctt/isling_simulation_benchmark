#### import data ####


library(tidyverse)
library(cowplot)


source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/sim_functions.R")
source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/figure-2_functions.R")


results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment1_OTC_chr1"
folders <- c("condition-breakdown_OTC-harder", "condition-breakdown_OTC-harder-2")

score_window_plot <- 5
coords_score_type_plot <- "coords_mean"
dist_plot_offset <- 0.5
facet_scales <- "free_y"

plot_tool_order <- c("isling", "Seeksv", "Polyidus", "ViFi")

# https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
offset <- dist_plot_offset
num_y_breaks <- 2

df_list <- list()
results <- importAllIntScoreExperiments(results_dir, folders, score_window_plot, coords_score_type)


df_list[['scores']] <- results$int_scores
incomplete <- results$incomplete

df_list[['scores']] <- df_list[['scores']] %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window_plot)

df_list[['found']] <- importAllFoundScores(results_dir, folders, incomplete, 
                                     score_window_plot, coords_score_type_plot)

df_list[['sim']] <- importAllSimScores(results_dir, folders, incomplete, 
                                   score_window_plot, coords_score_type_plot)

df_list[['scores']] %>% 
  select(batch, experiment) %>% 
  distinct()

# take note of the experiments in each batch
exp_names <- list()
for (b in unique(df_list[['scores']]$batch)) {
  exp_names[[b]] <- df_list[['scores']] %>% 
    filter(batch == b) %>% 
    pull(experiment) %>%
    unique()
}
exp_names

for (e in unique(df_list$scores$experiment)) {
  vars <- getExpVars(df_list$score, e)
  print(glue("experiment: {e}, variables: {vars}"))
}



plot_exps <- list(
  "OTC-fcov" = "fold coverage",
  "chromosomes" = "chromosome"
)
batch_plot <- "condition-breakdown_OTC-harder"
makeFigure2(plot_exps, df_list, "plots/figure2_v3_otc-harder.pdf", 12, 12)


# filter coverage == 6 from OTC-fcov experiment
for (name in names(df_list)) {
  df_list[[name]] <- df_list[[name]] %>% 
    filter(case_when(
      experiment == "OTC-fcov" ~ fcov != 6,
      TRUE ~ TRUE
    )) 
}

p <- makeFigure2(plot_exps, df_list, "plots/figure2_v4_otc-harder.pdf", 8, 11)

legend_plot <- df_list$scores %>% 
  filter(batch == batch_plot) %>% 
  ggplot(aes(x = PPV, y = TPR, color = tool)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom")
#print(legend_plot)


p2 <- cowplot::plot_grid(p, get_legend(legend_plot), ncol = 1, rel_heights = c(1, 0.05))
cowplot::save_plot("plots/figure2_v4_otc-harder.pdf", p2, base_width = 12, base_height = 6)


plotlist <- list()

for (e in names(plot_exps)) {
  for (y_var in c("TPR", "PPV")) {
    exp_var <- getExpVars(df_list[['scores']], e)[1]
  plotlist[[paste0(e, "_", y_var)]] <- df_list[['scores']] %>% 
    filter(batch == batch_plot) %>% 
    filter(experiment == e) %>% 
    ggplot(aes(x = !!sym(exp_var), y = !!sym(y_var), color = analysis_tool)) + 
    geom_line() 
  }
}
names(plotlist)
print(plotlist[["OTC-fcov_TPR"]])
print(plotlist[["OTC-fcov_PPV"]])
print(plotlist[["chromosomes_TPR"]])
print(plotlist[["chromosomes_PPV"]])



#### attempt 3 ####
results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment1_OTC_chr1/condition-breakdown_OTC-harder-2"
folders <- c("chr-fcov")
results <- importAllIntScoreExperiments(results_dir, folders, score_window_plot, coords_score_type)

unique(results$int_scores$experiment)

results$int_scores <- results$int_scores %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window_plot)


batch_plot <- "condition-breakdown_OTC-harder-2"
exp_plot <- "chr-fcov"
pal <- scales::hue_pal()(4)[1:3]

new_df_plot <- results$int_scores %>% 
  mutate(tool = str_replace(analysis_condition, "\\d+", "")) %>% 
  mutate(tool = case_when(
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "Seeksv",
    tool == "vifi" ~"ViFi",
    TRUE ~ tool
  )) %>% 
  ungroup() %>% 
  pivot_longer(TPR:PPV, names_to = "score_type", values_to = "score") %>% 
  group_by(host_name, fcov, analysis_condition, score_type, tool) %>% 
  summarise(mean_score = mean(score),
            sd_score = sd(score),
            sem_score = sd(score)/sqrt(n()),
            n_score = n(),
            mean_minus_sem = mean_score - sem_score,
            mean_plus_sem = mean_score + sem_score) 


PPV_plot <- new_df_plot %>% 
  filter(score_type == "PPV") %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, plot_tool_order)) %>% 
  ggplot(aes(x = fcov, y = mean_score, color = tool)) + 
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = mean_minus_sem, ymax = mean_plus_sem)) +
  facet_wrap(~host_name, nrow=1) +
  scale_x_log10() + 
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  ylim(0, 1) +
  xlab('Fold coverage') +
  ylab('PPV')

print(PPV_plot)

TPR_plot <- new_df_plot %>% 
  filter(score_type == "TPR") %>% 
  mutate(tooln = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, plot_tool_order)) %>% 
  ggplot(aes(x = fcov, y = mean_score, color = tool)) + 
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = mean_minus_sem, ymax = mean_plus_sem)) +
  facet_wrap(~host_name, nrow=1) +
  scale_x_log10() + 
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylim(0, 1) +
  xlab('Fold coverage') +
  ylab('TPR')


print(TPR_plot)

p <- cowplot::plot_grid(PPV_plot, TPR_plot, ncol = 1, labels = "AUTO")
print(p)

legend <- cowplot::get_legend(TPR_plot + theme(legend.position = "bottom") + guides(color = guide_legend("tool")))

p_legend <- cowplot::plot_grid(p, legend, rel_heights = c(1, 0.1), ncol = 1)
print(p_legend)
cowplot::save_plot("plots/figure2_v5_otc-harder.pdf", p_legend, base_height = 4.2, base_width = 5.5)


#### try with new data - condition-breakdown 3 ####

results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment1_OTC_chr1/condition-breakdown_OTC-harder-3"
folders <- c("chr-fcov")

df_list <- list()
results <- importAllIntScoreExperiments(results_dir, folders, score_window_plot, coords_score_type)

# import all data so that we can better look at the difference between analysis conditions
df_list[['scores']] <- results$int_scores
incomplete <- results$incomplete

df_list[['scores']] <- df_list[['scores']] %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window_plot)

df_list[['found']] <- importAllFoundScores(results_dir, folders, incomplete, 
                                           score_window_plot, coords_score_type_plot)

df_list[['sim']] <- importAllSimScores(results_dir, folders, incomplete, 
                                       score_window_plot, coords_score_type_plot)


unique(results$int_scores$experiment)

results$int_scores <- results$int_scores %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window_plot)

plot_exps <- list(
  "chr-fcov" = "chr/fcov"
)

#makeFigure2(plot_exps, df_list, "plots/figure2_v6_otc-harder.pdf", 12, 12)





batch_plot <- "condition-breakdown_OTC-harder-3"
exp_plot <- "chr-fcov"
pal <- scales::hue_pal()(4)[1:3]
min_mapq_keep <- 25

new_df_plot <- results$int_scores %>% 
  filter(analysis_condition != "vifi") %>% 
  mutate(tool = str_replace(analysis_condition, "\\d+", "")) %>% 
  mutate(tool = case_when(
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "Seeksv",
    tool == "vifi" ~ "ViFi",
    TRUE ~ tool
  )) %>% 
  ungroup() %>% 
  pivot_longer(TPR:PPV, names_to = "score_type", values_to = "score") %>% 
  group_by(host_name, fcov, tool, score_type, min_mapq) %>% 
  summarise(mean_score = mean(score),
            sd_score = sd(score),
            sem_score = sd(score)/sqrt(n()),
            n_score = n(),
            mean_minus_sem = mean_score - sem_score,
            mean_plus_sem = mean_score + sem_score) %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, plot_tool_order)) %>% 
  filter(ifelse(tool == "isling", min_mapq==min_mapq_keep, TRUE))  

unique(new_df_plot$min_mapq)

PPV_plot <- new_df_plot %>% 
  filter(score_type == "PPV") %>% 
  ggplot(aes(x = fcov, y = mean_score, color =tool)) + 
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = mean_minus_sem, ymax = mean_plus_sem)) +
  facet_wrap(~host_name, nrow=1) +
  scale_x_log10() + 
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  ylim(0, 1) +
  xlab('Fold coverage') +
  ylab('PPV') 

#+
#  scale_color_manual(values = pal)

print(PPV_plot + theme(legend.position = "bottom"))

TPR_plot <- new_df_plot %>% 
  filter(score_type == "TPR") %>% 
  ggplot(aes(x = fcov, y = mean_score, color = tool)) + 
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = mean_minus_sem, ymax = mean_plus_sem)) +
  facet_wrap(~host_name, nrow=1) +
  scale_x_log10() + 
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
#        strip.text.x = element_blank()
) +
  ylim(0, 1) +
  xlab('Fold coverage') +
  ylab('TPR')

#+
#  scale_color_manual(values = pal)


print(TPR_plot + theme(legend.position = "bottom"))

p <- cowplot::plot_grid(PPV_plot, TPR_plot, ncol = 1, labels = "AUTO")
print(p)

legend <- cowplot::get_legend(TPR_plot + theme(legend.position = "bottom") + guides(color = guide_legend("tool")))

p_legend <- cowplot::plot_grid(p, legend, rel_heights = c(1, 0.1), ncol = 1)
print(p_legend)

cowplot::save_plot("plots/figure2_v7_otc-harder.pdf", p_legend, base_height = 4.2, base_width = 5.5)

p <- cowplot::plot_grid(PPV_plot, TPR_plot, ncol = 1)
print(p)
p_legend <- cowplot::plot_grid(p, legend, rel_heights = c(1, 0.1), ncol = 1)
print(p_legend)
cowplot::save_plot("plots/otc-harder-fcov.pdf", p_legend, base_height = 4.2, base_width = 5.5)
