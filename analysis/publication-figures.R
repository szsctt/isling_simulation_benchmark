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
  mutate(analysis_tool = str_replace(analysis_tool, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling'))

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


# filter data for plotting
isling_to_keep <- "isling2"

plot_found_scores <- found_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE)) %>%  # keep isling2 only
  mutate(tool = str_match(analysis_condition, "(.+)\\d+")[,2])


# double check filtering
plot_found_scores %>% 
  group_by(experiment, tool, score_dist, score_type, condition, replicate) %>% 
  summarise(tp = sum(score == "tp"),
            fp = sum(score == "fp"))

plot_sim_scores <- sim_scores %>% 
  filter(score_type == coords_score_type_plot) %>% 
  filter(score_dist == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE)) %>%   # keep isling2 only
  mutate(tool = str_match(analysis_condition, "(.+)\\d+")[,2])

# double check filtering
plot_sim_scores %>% 
  group_by(experiment, tool, score_dist, score_type, condition, replicate) %>% 
  summarise(tp = sum(score == "tp"),
            fn = sum(score == "fn"),
            total_with_reads = n())

# scored integrations
plot_int_scores <- int_scores %>% 
  filter(coords_score_type == coords_score_type_plot) %>% 
  filter(window == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE))

plot_int_scores %>% 
  select(experiment, analysis_tool, window, coords_score_type, condition, replicate, tp, fp, fn, PPV, TPR)




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
print(foundplots[['OTC-easier']])
print(foundplots[['OTC-harder']])
print(foundplots[['AAV-easier']])
print(foundplots[['AAV-harder']])

simplots <- list()
for (e in unique(int_scores$experiment)) {
  simplots[[e]] <- plot_sim_scores  %>% 
    filter(experiment == e) %>% 
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
print(simplots[['OTC-easier']])
print(simplots[['OTC-harder']])
print(simplots[['AAV-easier']])
print(simplots[['AAV-harder']])

chrplots <- list()
for (e in unique(int_scores$experiment)) {
chrplots[[e]] <- plot_found_scores %>% 
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
chrplots[['OTC-easier']] <- chrplots[['OTC-easier']] + theme(legend.position = c(0.8, 0.8))
print(chrplots[['OTC-easier']])
print(chrplots[['OTC-harder']])
print(chrplots[['AAV-easier']])
print(chrplots[['AAV-harder']])


scoreplots <- list()
for (e in unique(int_scores$experiment)) {
  scoreplots[[e]] <- plot_int_scores %>% 
  filter(experiment == e) %>% 
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
scoreplots[['OTC-easier']] <- scoreplots[['OTC-easier']] + theme(legend.position = c(0.2, 0.8))
print(scoreplots[['OTC-easier']])
print(scoreplots[['OTC-harder']])
print(scoreplots[['AAV-easier']])
print(scoreplots[['AAV-harder']])



figure_1 <- cowplot::plot_grid(simplots[['OTC-easier']], foundplots[['OTC-easier']], 
                               chrplots[['OTC-easier']], scoreplots[['OTC-easier']],
                               simplots[['OTC-harder']], foundplots[['OTC-harder']], 
                               chrplots[['OTC-harder']], scoreplots[['OTC-harder']],
                               simplots[['AAV-easier']], foundplots[['AAV-easier']], 
                               chrplots[['AAV-easier']], scoreplots[['AAV-easier']],
                               simplots[['AAV-harder']], foundplots[['AAV-harder']], 
                               chrplots[['AAV-harder']], scoreplots[['AAV-harder']],
                               labels = "AUTO")
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
print(simchrplots[['OTC-easier']])
print(simchrplots[['OTC-harder']])
print(simchrplots[['AAV-easier']])
print(simchrplots[['AAV-harder']])

cowplotlist <- list()
for (e in unique(int_scores$experiment)) {
  cowplotlist[[paste0(e, "_", "sim")]] <- (simplots[[e]])
  cowplotlist[[paste0(e, "_", "simchr")]] <- (simchrplots[[e]]) 
  cowplotlist[[paste0(e, "_", "found")]] <- (foundplots[[e]])
  cowplotlist[[paste0(e, "_", "foundchr")]] <- (chrplots[[e]])   
}

figure_1 <- cowplot::plot_grid(plotlist = cowplotlist, labels = "AUTO")
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
  cowplotlist[[paste0(e, "_", "score")]] + theme(legend.position = "bottom")
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

#### figure 2 ####

# import scored integration sites for condition-breakdown'
exp_dir <- file.path(results_dir, "condition-breakdown-1")
int_scores <- importIntScoreExperiment(exp_dir) 

# only use postprocessed data from our pipeline
int_scores <- int_scores %>% 
  filter(case_when((analysis_tool != "pipeline") ~ TRUE,
                   post ~ TRUE,
                   TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '\\d', ''))
  
# for 'AAV' experiment, filter only for results where we used the correct viral reference
int_scores <- int_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
  ungroup()

# check the filtering is correct
int_scores %>% 
  filter(experiment == "AAV") %>% 
  select(virus_name, analysis_virus) %>% 
  distinct()

# intial plot of scored ints for each condition
unique(int_scores$experiment)
exp_vars <- list(
  "AAV" = "virus_name",
  "chromosomes" = "host_name",
  "OTC-deletion" = "p_delete",
    "OTC-epi" = "epi_num",
    "OTC-fcov" = "fcov",
    "OTC-gap" = "p_gap",
    "OTC-min_len" = "min_len",
    "OTC-min_sep" = "min_sep",
    "OTC-overlap" = "p_overlap",
    "OTC-rearrange" = "p_rearrange",
    "OTC-whole" = "p_whole"
)
scoreplots <- list()
for (e in unique(int_scores$experiment)) {
scoreplots[[e]] <- int_scores %>% 
  filter(experiment == e) %>% 
    mutate(!!exp_vars[[e]] := as.factor(!!sym(exp_vars[[e]]))) %>% 
    ggplot(aes(x = PPV, y = TPR, color = analysis_tool, shape = !!sym(exp_vars[[e]]))) +
    geom_point() +
    xlim(0, 1) +
    ylim(0, 1) +
    facet_grid(vars(window), vars(coords_score_type)) +
  ggtitle(e)
}

print(scoreplots['AAV'])
print(scoreplots['chromosomes'])
print(scoreplots['OTC-deletion'])
print(scoreplots['OTC-rearrange'])
print(scoreplots['OTC-fcov'])
print(scoreplots['OTC-min_len'])
print(scoreplots['OTC-min_sep'])
print(scoreplots['OTC-gap'])
print(scoreplots['OTC-whole'])
print(scoreplots['OTC-epi'])


int_scores <- int_scores %>% 
  filter(window == score_window) %>% 
  filter(coords_score_type == coords_score_type_plot)

scoreplots <- list()
for (e in unique(int_scores$experiment)) {
  scoreplots[[e]] <- int_scores %>% 
    filter(experiment == e) %>% 
    mutate(!!exp_vars[[e]] := as.factor(!!sym(exp_vars[[e]]))) %>% 
    ggplot(aes(x = PPV, y = TPR, color = analysis_tool, shape = !!sym(exp_vars[[e]]))) +
    geom_point() +
    xlim(0, 1) +
    ylim(0, 1) +
    ggtitle(e) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.8))
}

figure_2 <- cowplot::plot_grid(plotlist = scoreplots)
print(figure_2)
cowplot::save_plot("plots/condition-breakdown-1_scores_all.pdf", figure_2, base_height = 8)


#import distances from each found integration to nearest simulated integration
found_scores <- importNearestSimToFound(exp_dir, score_window, coords_score_type_plot) 

# only use postprocessed data from our pipeline
found_scores <- found_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  filter(score_dist == score_window) %>% 
  filter(score_type == coords_score_type_plot) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '\\d', ''))


# add analysis conditions
conds <- importSimulationConditions(exp_dir)

found_scores <- found_scores %>% 
  left_join(conds, by=c("experiment", "condition", "replicate"))

found_scores <- found_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
  ungroup()


found_scores %>% 
  mutate(dist = dist+dist_plot_offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition, linetype = analysis_host)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(analysis_condition), cols = vars(experiment), scales = "free_y")
ggsave("plots/found_scores_2.pdf")

foundplots <- list()
for (e in names(exp_vars)) {
  foundplots[[e]] <- found_scores %>% 
    filter(experiment == e) %>% 
    mutate(!!exp_vars[[e]] := as.factor(!!sym(exp_vars[[e]]))) %>% 
    mutate(dist = dist+offset) %>% 
    ggplot(aes(x = dist, color = analysis_condition, linetype = !!sym(exp_vars[[e]]))) +
    geom_freqpoly(bins = 100) +
    scale_x_log10() +
    facet_grid(rows = vars(analysis_condition)) +
    xlab("Distance") +
    ylab("Count") +
    theme_classic() + 
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      legend.position = "none"
    )  +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
    labs(title = e, x="Distance", y="Count", color="tool", linetype=exp_vars[[e]])

}
foundplots[[names(exp_vars)[1]]]

figure_2 <- cowplot::plot_grid(plotlist = foundplots)
print(figure_2)

cowplot::save_plot("plots/condition-breakdown-1_found-dist-all_coords-mean.pdf", figure_2, base_height = 8)

chrplots <- list()
for (e in names(exp_vars)) {
  chrplots[[e]] <- found_scores %>% 
    filter(experiment == e) %>% 
    mutate(correct_chr = (chr == host_name)) %>% 
    mutate(!!exp_vars[[e]] := as.factor(!!sym(exp_vars[[e]]))) %>% 
    ggplot(aes(x = !!sym(exp_vars[[e]]), fill = correct_chr)) +
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
    guides(fill=guide_legend("correct chr")) +
    facet_grid(cols = vars(analysis_condition)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(e)
} 

chrplots[[names(exp_vars)[1]]]

figure_2 <- cowplot::plot_grid(plotlist = chrplots)
print(figure_2)

cowplot::save_plot("plots/condition-breakdown-1_found-chr-all_coords-mean.pdf", figure_2, base_height = 8)


# import distances from each simulated integration to the nearest found integration
sim_scores <- importNearestFoundToSim(exp_dir, score_window, coords_score_type_plot) 

sim_scores <- sim_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '\\d', ''))  %>% 
  filter(score_dist == score_window) %>% 
  filter(score_type == "coords_mean")

sim_scores <- sim_scores %>% 
  left_join(conds, by=c("experiment", "condition", "replicate"))


sim_scores <- sim_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
  ungroup()

sim_scores %>% 
  mutate(dist = dist+dist_plot_offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition, linetype = analysis_host)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(analysis_condition), cols = vars(experiment), scales = "free_y")
ggsave("plots/sim_scores_2.pdf")


simplots <- list()
for (e in names(exp_vars)) {
  simplots[[e]] <- sim_scores %>% 
    filter(experiment == e) %>% 
    mutate(!!exp_vars[[e]] := as.factor(!!sym(exp_vars[[e]]))) %>% 
    mutate(dist = dist+offset) %>% 
    ggplot(aes(x = dist, color = analysis_condition, linetype = !!sym(exp_vars[[e]]))) +
    geom_freqpoly(bins = 100) +
    scale_x_log10() +
    facet_grid(rows = vars(analysis_condition)) +
    theme_classic() + 
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
#      legend.position = "none"
    )  +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
    labs(title = e, x="Distance", y="Count", color="tool", linetype=exp_vars[[e]]) +
    guides(color = FALSE) 
  
}
simplots[[names(exp_vars)[1]]] + theme(legend.position = "right")

figure_2 <- cowplot::plot_grid(plotlist = simplots)
print(figure_2)

cowplot::save_plot("plots/condition-breakdown-1_sim-dist-all_coords-mean.pdf", figure_2, base_height = 8)



#### figure 2, version 0 ####
# create new list with interleaved plots
cowplotlist < list()
for (e in names(exp_vars)) {
  cowplotlist[[paste0(e, "sim")]] <- simplots[[e]] + theme(legend.position = "none")
  cowplotlist[[paste0(e, "found")]] <- foundplots[[e]] + theme(legend.position = "none")
  cowplotlist[[paste0(e, "found_chr")]] <- chrplots[[e]] + theme(legend.position = "none")
  cowplotlist[[paste0(e, "scores")]] <- scoreplots[[e]] + theme(legend.position = "none")
}

figure_2 <- cowplot::plot_grid(plotlist = cowplotlist, ncol=4)
print(figure_2)

cowplot::save_plot("plots/figure2_v0.pdf", figure_2, base_height = 12)

#### figure 2, version 1 ####
# the plot only a few conditions
exp_plot <- c("OTC-fcov", "chromosomes", "AAV", "OTC-rearrange", "OTC-deletion", "OTC-whole")
cowplotlist <- list()
for (e in exp_plot) {
    cat("adding experiment", e, "\n")
    cowplotlist[[paste0(e, "sim")]] <- simplots[[e]] + 
      theme(
        legend.position = "right",
        plot.title = element_blank()
            ) +
      guides(color="none") 
    cowplotlist[[paste0(e, "found")]] <- foundplots[[e]] + 
      theme(legend.position = "none",
            plot.title = element_blank()) 
    cowplotlist[[paste0(e, "found_chr")]] <- chrplots[[e]] + 
      theme(legend.position = "none",
            plot.title = element_blank())
    cowplotlist[[paste0(e, "scores")]] <- scoreplots[[e]] + 
      theme(legend.position = c(0.25, 0.6),
            plot.title = element_blank()) +
      guides(color="none") 
}
cowplotlist[[paste0(e, "sim")]]
cowplotlist[[paste0(e, "found")]] 
cowplotlist[[paste0(e, "found_chr")]]
cowplotlist[[paste0(e, "scores")]]

p <- int_scores %>% 
  ggplot(aes(x = PPV, y = TPR, colour = tool)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom")

legend <- get_legend(p)

figure_2 <- cowplot::plot_grid(plotlist = cowplotlist, ncol=4, nrow=length(exp_plot), labels="AUTO")
figure_2 <- cowplot::plot_grid(figure_2, legend, ncol=1, rel_heights = c(1, 0.05))
#print(figure_2)

cowplot::save_plot("plots/figure2_v1.pdf", figure_2, base_height = 12, base_width = 12)
