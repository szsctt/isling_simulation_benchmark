library(tidyverse)
library(cowplot)

source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/sim_functions.R")

results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment1_OTC_chr1"

score_window <- 5
coords_score_type_plot <- "coords_mean"
dist_plot_offset <- 0.5
facet_scales <- "free_y"

# https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
offset <- dist_plot_offset
num_y_breaks <- 2


#### figure 1 ####

# import scored integrations
exp_dir <- file.path(results_dir, "easier-harder")
int_scores <- importIntScoreExperiment(exp_dir) 

# only use postprocessed data from our pipeline
int_scores <- int_scores %>% 
  filter(case_when((analysis_tool != "pipeline") ~ TRUE,
                   post ~ TRUE,
                   TRUE ~ FALSE
  )) %>% 
  mutate(analysis_tool = str_replace(analysis_tool, 'analysis|pipeline', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis|pipeline', 'isling')) 

#import distances from each found integration to nearest simulated integration
found_scores <- importNearestSimToFound(exp_dir, score_window, coords_score_type_plot) 

# only use postprocessed data from our pipeline
found_scores <- found_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling'))

# import distances from each simulated integration to the nearest found integration
sim_scores <- importNearestFoundToSim(exp_dir, score_window, coords_score_type_plot) 

sim_scores <- sim_scores%>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>%  
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling'))


# we had some different analysis and simulation conditions - figure out which ones to keep
isling_conditions <- int_scores 



int_scores %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window) %>% 
  ggplot(aes(x = PPV, y = TPR, color = analysis_condition)) +
  geom_point() +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_wrap(experiment~condition) +
  theme_classic() +
  theme(
    strip.background = element_rect(
      color = 'white', fill = "#f2f2f2", linetype = NULL,
    ),
    legend.position = "bottom"
  )



found_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>%  
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition)) +
  geom_freqpoly(bins = 20) +
  scale_x_log10() +
  xlab("Distance") +
  ylab("Count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
  facet_wrap(experiment~condition)  +
  theme(legend.position = "bottom")


found_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>% 
  filter(str_detect(analysis_condition, "isling")) %>% 
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition)) +
  geom_freqpoly(bins = 20) +
  scale_x_log10() +
  xlab("Distance") +
  ylab("Count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
  facet_wrap(experiment~condition)  +
  theme(legend.position = "bottom")

sim_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>%  
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition)) +
  geom_freqpoly(bins = 20) +
  scale_x_log10() +
  xlab("Distance") +
  ylab("Count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
  facet_wrap(experiment~condition)  +
  theme(legend.position = "bottom")

sim_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>% 
  filter(str_detect(analysis_condition, "isling")) %>% 
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition)) +
  geom_freqpoly(bins = 20) +
  scale_x_log10() +
  xlab("Distance") +
  ylab("Count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
  facet_wrap(experiment~condition)  +
  theme(legend.position = "bottom")


found_scores %>% 
  mutate(correct_chr = (chr == "chr1")) %>% 
  ggplot(aes(x = analysis_condition, fill = correct_chr)) +
  geom_bar() +
  ylab("Count") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90)
  )  +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
  scale_fill_grey()  +
  guides(fill=guide_legend("correct chr"))  +
  facet_wrap(experiment~condition) 




# filter data for plotting
isling_to_keep <- "isling2"
condition_to_keep <- "cond5"

plot_found_scores <- found_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE)) %>%  # keep isling2 only
  filter(condition == condition_to_keep) %>% 
  mutate(tool = str_match(analysis_condition, "(.+)\\d+")[,2]) %>% 
  mutate(tool = str_replace(tool, "polyidus", "Polyidus")) %>% 
  mutate(tool = str_replace(tool, "seeksv", "Seeksv")) %>% 
  mutate(tool = str_replace(tool, "vifi", "ViFi"))


# double check filtering
plot_found_scores %>% 
  group_by(experiment, tool, score_dist, score_type, condition, replicate) %>% 
  summarise(tp = sum(score == "tp"),
            fp = sum(score == "fp"))

plot_sim_scores <- sim_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>% 
  filter(condition == condition_to_keep) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE)) %>%   # keep isling2 only
  mutate(tool = str_match(analysis_condition, "(.+)\\d+")[,2])  %>% 
  mutate(tool = str_replace(tool, "polyidus", "Polyidus")) %>% 
  mutate(tool = str_replace(tool, "seeksv", "Seeksv")) %>% 
  mutate(tool = str_replace(tool , "vifi", "ViFi"))

# double check filtering
plot_sim_scores %>% 
  group_by(experiment, tool, score_dist, score_type, condition, replicate) %>% 
  summarise(tp = sum(score == "tp"),
            fn = sum(score == "fn"),
            total_with_reads = n())

# scored integrations
plot_int_scores <- int_scores %>% 
  filter(condition == condition_to_keep) %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE))  

plot_int_scores <- plot_int_scores %>% 
  select(experiment, analysis_tool, window, coords_score_type, condition, replicate, tp, fp, fn, PPV, TPR)  %>% 
  mutate(analysis_tool = str_replace(analysis_tool, "polyidus", "Polyidus")) %>% 
  mutate(analysis_tool = str_replace(analysis_tool, "seeksv", "Seeksv")) %>% 
  mutate(analysis_tool = str_replace(analysis_tool, "vifi", "ViFi")) 


OTC_found <- plot_found_scores %>% 
  filter(experiment == "OTC-harder") %>% 
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = tool)) +
  geom_freqpoly(bins = 20, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(tool), scales = "free_y") +
  xlab("Distance") +
  ylab("Count") +
  theme_classic()  + 
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks))

print(OTC_found)


OTC_sim <- plot_sim_scores %>% 
  filter(experiment == "OTC-harder") %>% 
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = tool)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(tool), scales = "free_y") +
  xlab("Distance") +
  ylab("Count") +
  theme_classic()  + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks))
print(OTC_sim)

OTC_scores <- plot_int_scores %>% 
  filter(experiment == "OTC-harder") %>% 
  ggplot(aes(x = PPV, y = TPR, colour = analysis_tool)) +
  geom_point(show.legend = FALSE, alpha = 0.5) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Positive predictive value") +
  ylab("True positive rate") +
  theme_classic()
print(OTC_scores)


AAV_found <- plot_found_scores %>% 
  filter(experiment == "AAV-harder") %>% 
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = tool)) +
  geom_freqpoly(bins = 30, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(tool), scales = "free_y") +
  xlab("Distance") +
  ylab("Count") +
  theme_classic()  + 
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank()
  )  +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) 
print(AAV_found)


AAV_sim <- plot_sim_scores %>% 
  filter(experiment == "AAV-harder") %>% 
  mutate(dist = dist+offset) %>% 
  ggplot(aes(x = dist, colour = tool)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(tool), scales = "free_y") +
  xlab("Distance") +
  ylab("Count") +
  theme_classic() + 
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom"
  )  +
  scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks))
print(AAV_sim)

AAV_scores <- plot_int_scores %>% 
  filter(experiment == "AAV-harder") %>% 
  ggplot(aes(x = PPV, y = TPR, colour = analysis_tool)) +
  geom_point(show.legend = FALSE, alpha = 0.5) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Positive predictive value") +
  ylab("True positive rate") +
  theme_classic()
print(AAV_scores)

dir.create("plots")
figure_1 <- cowplot::plot_grid(AAV_found, AAV_sim, AAV_scores, OTC_found, OTC_sim, OTC_scores, labels = "AUTO")
print(figure_1)
cowplot::save_plot("plots/figure_1.pdf", figure_1, base_height = 4.2)

#a, b, c : AAV_harder
#d, e, f: OTC_harder
#a, d: distance from each output integration to nearest simulated integration
#b, e: distance from each simulated integration to nearest output integration
#c, f: PPV and TPR with cutoff for distance = 5 (n = 3)
#red - isling, green - polyidus, blue - seeksv, purple - vifi


#### figure 1, version 2 ####

foundplots <- list()
for (e in unique(int_scores$experiment)) {
  foundplots[[e]]  <- plot_found_scores %>% 
    filter(experiment == e) %>% 
    mutate(tool = forcats::fct_reorder(tool, desc(tool))) %>% 
    mutate(dist = dist+offset) %>% 
    ggplot(aes(x = dist, colour = tool)) +
    geom_freqpoly(bins = 50) +
    scale_x_log10() +
    facet_grid(rows = vars(tool), scales = "free_y") +
    xlab("Distance") +
    ylab("Count") +
    theme_classic()  + 
    theme(
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks))
}
print(foundplots[['OTC-harder']] + theme(legend.position = "bottom"))
print(foundplots[['AAV-harder']] + theme(legend.position = "bottom"))

simplots <- list()
for (e in unique(int_scores$experiment)) {
  simplots[[e]] <- plot_sim_scores  %>% 
    filter(experiment == e) %>% 
    mutate(tool = forcats::fct_reorder(tool, desc(tool))) %>% 
    mutate(dist = dist+offset) %>% 
    ggplot(aes(x = dist, colour = tool)) +
    geom_freqpoly(bins = 100) +
    scale_x_log10() +
    facet_grid(rows = vars(tool), scales = "free_y") +
    xlab("Distance") +
    ylab("Count") +
    theme_classic() + 
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "none"
    )  +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks))
}

print(simplots[['OTC-harder']] + theme(legend.position = "bottom"))

print(simplots[['AAV-harder']] + theme(legend.position = "bottom"))

chrplots <- list()
for (e in unique(int_scores$experiment)) {
  chrplots[[e]] <- plot_found_scores %>% 
    filter(experiment == e) %>% 
    mutate(correct_chr = (chr == "chr1")) %>% 
    mutate(tool = forcats::fct_reorder(tool, desc(tool))) %>% 
    ggplot(aes(x = tool, fill = correct_chr)) +
    geom_bar() +
    ylab("Count") +
    theme_classic() + 
    theme(
      axis.title.x = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "none"
    )  +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
    scale_fill_grey()  +
    guides(fill=guide_legend("correct chr")) 
} 
chrplots[['OTC-easier']] <- chrplots[['OTC-easier']] + theme(legend.position = c(0.8, 0.8))

print(chrplots[['OTC-harder']])
print(chrplots[['AAV-harder']])


scoreplots <- list()
for (e in unique(int_scores$experiment)) {
  scoreplots[[e]] <- plot_int_scores %>% 
    filter(experiment == e) %>% 
    mutate(analysis_tool = forcats::fct_reorder(analysis_tool, desc(analysis_tool))) %>% 
    ggplot(aes(x = PPV, y = TPR, colour = analysis_tool)) +
    geom_point(alpha = 0.7) +
    xlim(0, 1) +
    ylim(0, 1) +
    xlab("PPV") +
    ylab("TPR") +
    theme_classic() +
    theme(legend.position = "none") +
    guides(color=guide_legend("tool")) 
} 

print(scoreplots[['OTC-harder']] + theme(legend.position = "bottom"))
print(scoreplots[['AAV-harder']] + theme(legend.position = "bottom"))



figure_1 <- cowplot::plot_grid(simplots[['OTC-harder']], foundplots[['OTC-harder']], 
                               chrplots[['OTC-harder']], scoreplots[['OTC-harder']],
                               simplots[['AAV-harder']], foundplots[['AAV-harder']], 
                               chrplots[['AAV-harder']], scoreplots[['AAV-harder']],
                               labels = "AUTO",  ncol=4)
print(figure_1)
cowplot::save_plot("plots/figure_1_v2.pdf", figure_1, base_height = 8)

#### figure 1 version 3 ####

simchrplots <- list()
for (e in unique(int_scores$experiment)) {
  simchrplots[[e]] <- plot_sim_scores %>% 
    filter(experiment == e) %>% 
    mutate(correct_chr = (chr == "chr1")) %>% 
    ggplot(aes(x = tool, fill = correct_chr)) +
    geom_bar() +
    ylab("Count") +
    theme_classic() + 
    theme(
      axis.title.x = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "none"
    )  +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
    scale_fill_grey()  +
    guides(fill=guide_legend("correct chr")) 
} 
simchrplots[['OTC-easier']] <- simchrplots[['OTC-easier']] + theme(legend.position = c(0.8, 0.8))

print(simchrplots[['OTC-harder']])

print(simchrplots[['AAV-harder']])

cowplotlist <- list()
for (e in unique(int_scores$experiment)) {
  cowplotlist[[paste0(e, "_", "sim")]] <- (simplots[[e]])
  cowplotlist[[paste0(e, "_", "simchr")]] <- (simchrplots[[e]]) 
  cowplotlist[[paste0(e, "_", "found")]] <- (foundplots[[e]])
  cowplotlist[[paste0(e, "_", "foundchr")]] <- (chrplots[[e]])   
}

figure_1 <- cowplot::plot_grid(plotlist = cowplotlist, labels = "AUTO", ncol = 4)
print(figure_1)
cowplot::save_plot("plots/figure_1_v3.pdf", figure_1, base_height = 8)

#### figure 1 version 4 ####

plotexps <- c("AAV-harder", "OTC-harder")
cowplotlist <- list()
for (e in plotexps) {
  cowplotlist[[paste0(e, "_", "sim")]] <- (simplots[[e]]) + 
    facet_grid(rows = vars(tool)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  cowplotlist[[paste0(e, "_", "found")]] <- (foundplots[[e]])  + 
    facet_grid(rows = vars(tool)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  cowplotlist[[paste0(e, "_", "foundchr")]] <- (chrplots[[e]])  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

figure_1 <- cowplot::plot_grid(plotlist = cowplotlist, labels = "AUTO", ncol=3)
print(figure_1)
cowplot::save_plot("plots/figure_1_v4.pdf", figure_1, base_height = 4.2)

#### figure 1 version 5 ####

plotexps <- c("AAV-harder", "OTC-harder")
cowplotlist <- list()
for (e in plotexps) {
  cowplotlist[[paste0(e, "_", "sim")]] <- (simplots[[e]]) + 
    facet_grid(rows = vars(tool)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  cowplotlist[[paste0(e, "_", "found")]] <- (foundplots[[e]])  + 
    facet_grid(rows = vars(tool)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  cowplotlist[[paste0(e, "_", "foundchr")]] <- (chrplots[[e]])  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  cowplotlist[[paste0(e, "_", "score")]] <- (scoreplots[[e]])  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
legend <- get_legend(
  cowplotlist[[paste0(e, "_", "sim")]] + theme(legend.position = "bottom")
)

figure_1 <- cowplot::plot_grid(plotlist = cowplotlist, labels = "AUTO", ncol=4)
figure_1 <- cowplot::plot_grid(figure_1, legend, ncol=1, rel_heights = c(1, 0.1))
print(figure_1)
cowplot::save_plot("plots/figure_1_v5.pdf", figure_1, base_height = 4.2)


#### figure 1 version 6 ####

plotexps <- c("AAV-harder", "OTC-harder")
rows <- list()
for (e in plotexps) {
  cowplotlist <- list()
  cowplotlist[[paste0(e, "_", "sim")]] <- (simplots[[e]]) + 
    facet_grid(rows = vars(tool)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  cowplotlist[[paste0(e, "_", "found")]] <- (foundplots[[e]])  + 
    facet_grid(rows = vars(tool)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  cowplotlist[[paste0(e, "_", "foundchr")]] <- (chrplots[[e]])  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  cowplotlist[[paste0(e, "_", "scores")]] <- scoreplots[[e]] +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  rows[[paste0(e, "null")]] <- NULL
  rows[[e]] <- cowplot::plot_grid(plotlist = cowplotlist, labels = c("i", "ii", "iii", "iv"), ncol=4)
}
legend <- get_legend(
  cowplotlist[[paste0(e, "_", "scores")]] + theme(legend.position = "bottom")
)

figure_1 <- cowplot::plot_grid(NULL, rows[['AAV-harder']], NULL, rows[['OTC-harder']], labels = c("A", "", "B", ""), ncol=2, nrow=length(plotexps), rel_widths=c(0.05, 1))
figure_1 <- cowplot::plot_grid(figure_1, legend, ncol=1, rel_heights = c(1, 0.1))
print(figure_1)
cowplot::save_plot("plots/figure_1_v6.pdf", figure_1, base_height = 4.2)
