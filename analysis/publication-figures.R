library(tidyverse)
library(lubridate)
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

exp_dir <- file.path(results_dir, "condition-breakdown-2")
int_scores2 <- importIntScoreExperiment(exp_dir) 

int_scores <- bind_rows(int_scores, int_scores2) %>% 
  filter(window == 5) %>% 
  filter(coords_score_type == "coords_mean")

# only use postprocessed data from our pipeline
int_scores <- int_scores %>% 
  filter(case_when((analysis_tool != "pipeline") ~ TRUE,
                   post ~ TRUE,
                   TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_tool, 'pipeline', 'isling')) %>% 
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

# check that we've got three replicates for each condition in each experiment, for each tool
n_per_cond <- int_scores %>% 
  filter(window == 5) %>% 
  filter(coords_score_type == "coords_mean") %>% 
  mutate(batch = str_split(filename, "/", simplify = TRUE)[,10]) %>% 
  select(window, coords_score_type, experiment, condition, replicate, tool, batch) %>% 
  distinct() %>% 
  group_by(batch, window, coords_score_type, experiment, condition, tool) %>% 
  summarise(count = n())

incomplete_exps <- n_per_cond %>% filter(count != 3) %>% pull(experiment)
incomplete_conds <- n_per_cond %>% filter(count != 3) %>% pull(condition)
incomplete_batches <- n_per_cond %>% filter(count != 3) %>% pull(batch)
incomplete_tools <- n_per_cond %>% filter(count != 3) %>% pull(tool)

# if we don't have three replicates for each tool, drop that condition
int_scores <- int_scores %>% 
  mutate(batch = str_split(filename, "/", simplify = TRUE)[,10]) %>% 
  ungroup() %>% 
  rowwise() %>% 
  filter(!((experiment %in% incomplete_exps) & (batch %in% incomplete_batches) & (condition %in% incomplete_conds))) %>% 
  ungroup()  

# what's left after this
int_scores %>% 
  select(batch, experiment, condition, replicate) %>% 
  distinct() %>% 
  group_by(batch, experiment) %>% 
  summarise(conditions = n_distinct(condition),
            observations = n())


# we had an OTC-fcov experiment in both -1 and -2, so only keep one
int_scores %>% 
  filter(experiment == "OTC-fcov") %>% 
  select(batch, condition, fcov) %>% 
  distinct()

int_scores <- int_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "OTC-fcov", batch == "condition-breakdown-2", TRUE)) %>% 
  ungroup()


# initial plot with just conditions - no info about manipulated variables
scoreplots <- list()
for (e in unique(int_scores$experiment)) {
  scoreplots[[e]] <- int_scores %>% 
    filter(experiment == e) %>% 
    ggplot(aes(x = PPV, y = TPR, shape = analysis_tool, color = condition)) +
    geom_point() +
    xlim(0, 1) +
    ylim(0, 1) +
    theme(legend.position = "none") +
    ggtitle(e)
}

unique(int_scores$experiment)
print(scoreplots[['AAV']] + theme(legend.position = "right"))
print(scoreplots[['OTC-deletion']] + theme(legend.position = "right"))
print(scoreplots[['OTC-epi']] + theme(legend.position = "right"))
print(scoreplots[['OTC-fcov']] + theme(legend.position = "right"))
print(scoreplots[['OTC-gap']] + theme(legend.position = "right"))
print(scoreplots[['OTC-min_len']] + theme(legend.position = "right"))
print(scoreplots[['OTC-min_sep']] + theme(legend.position = "right"))
print(scoreplots[['OTC-overlap']] + theme(legend.position = "right"))
print(scoreplots[['OTC-rearrange']] + theme(legend.position = "right"))
print(scoreplots[['OTC-whole']] + theme(legend.position = "right"))
print(scoreplots[['chromosomes']] + theme(legend.position = "right"))
print(scoreplots[['OTC-frag_len']] + theme(legend.position = "right"))
print(scoreplots[['OTC-host_deletion']] + theme(legend.position = "right"))
print(scoreplots[['OTC-juncs']]) + theme(legend.position = "right")
print(scoreplots[['OTC-rearrange-epi']] + theme(legend.position = "right"))
print(scoreplots[['OTC-rearrange-frags']] + theme(legend.position = "right"))

figure_2 <- cowplot::plot_grid(plotlist = scoreplots)
print(figure_2)
cowplot::save_plot("plots/condition-breakdown-1_scores_all.pdf", figure_2, base_height = 8)

#### figure 2 plots ####

# list of experiments and which variables were manipulated in each
unique(int_scores$experiment)
exp_vars <- list(
  "AAV" = c("virus_name"),
  "chromosomes" = c("host_name"),
  "OTC-deletion" = c("p_delete"),
  "OTC-epi" = c("epi_num"),
  "OTC-fcov" = c("fcov"),
  "OTC-gap" = c("p_gap"),
  "OTC-min_len" = c("min_len"),
  "OTC-min_sep" = c("min_sep"),
  "OTC-overlap" = c("p_overlap"),
  "OTC-rearrange" = c("p_rearrange"),
  "OTC-whole" = c("p_whole"),
  "OTC-rearrange-frags" = c("p_rearrange", "p_delete", "lambda_split"),
  "OTC-rearrange-epi" = c("p_rearrange", "p_delete", "epi_num"),
  "OTC-juncs" = c("p_overlap", "p_gap", "lambda_junction"),
  'OTC-host_deletion' = c("p_host_deletion", "lambda_host_deletion"),
  "OTC-frag_len" = c("frag_len")
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
    ggtitle(e) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.8))
}

figure_2 <- cowplot::plot_grid(plotlist = scoreplots)
print(figure_2)
cowplot::save_plot("plots/condition-breakdown-1_scores_all.pdf", figure_2, base_height = 8)


#import distances from each found integration to nearest simulated integration
found_scores <- tibble(
  batch = c("condition-breakdown-1","condition-breakdown-2"),
  exp_dir = map(batch, ~file.path(results_dir, .))
) %>% 
  mutate(data = map(exp_dir, ~importNearestSimToFound(., score_window, coords_score_type_plot))) %>% 
  unnest(data)

unique(found_scores$experiment)
found_scores %>% 
  group_by(batch, experiment) %>% 
  summarise(n_found = n())


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

# remove any data from conditions in which all three replicates weren't completed
found_scores <- found_scores %>% 
  ungroup() %>% 
  rowwise() %>% 
  filter(!((experiment %in% incomplete_exps) & (batch %in% incomplete_batches) & (condition %in% incomplete_conds))) %>% 
  ungroup()  

# check filtering
found_scores %>% 
  mutate(batch = str_split(filename, "/", simplify = TRUE)[,10]) %>% 
  select(experiment, condition, replicate, analysis_condition, batch) %>% 
  distinct() %>% 
  group_by(batch, experiment, condition, analysis_condition) %>% 
  summarise(count = n()) %>% 
  filter(count != 3)
# there's only results from one replicate of seeksv in the OTC-fcov experiment, condition0 - 
# but I checked the files and they're present, and it's because there weren't any found integrations for two replicates in this condition

# add analysis conditions
conds  <- tibble(
  batch = c("condition-breakdown-1","condition-breakdown-2"),
  exp_dir = map(batch, ~file.path(results_dir, .))
)  %>% 
  mutate(conds = map(exp_dir, ~importSimulationConditions(.))) %>% 
  unnest(conds)
  
found_scores <- found_scores %>% 
  left_join(conds, by=c("experiment", "condition", "replicate", "batch"))

# filter out any conditions were different virus was integrated and analysed
found_scores <- found_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
  ungroup()

found_scores %>% 
  filter(experiment == "AAV") %>% 
  select(virus_name, analysis_virus) %>% 
  distinct()

# we had a fcov experiment in batch 1 and 2, so just keep batch 1
found_scores <- found_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "OTC-fcov", batch == "condition-breakdown-2", TRUE)) %>% 
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
    mutate(dist = dist+offset) %>% 
    ggplot(aes(x = dist, color = condition)) +
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
    ggplot(aes(x = condition, fill = correct_chr)) +
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
#import distances from each found integration to nearest simulated integration
sim_scores <- tibble(
  batch = c("condition-breakdown-1","condition-breakdown-2"),
  exp_dir = map(batch, ~file.path(results_dir, .))
) %>% 
  mutate(data = map(exp_dir, ~importNearestFoundToSim(., score_window, coords_score_type_plot))) %>% 
  unnest(data)

unique(sim_scores$experiment)
sim_scores %>% 
  group_by(batch, experiment) %>% 
  summarise(n_sim = n())


# only use postprocessed data from our pipeline
sim_scores <- sim_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  filter(score_dist == score_window) %>% 
  filter(score_type == coords_score_type_plot) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '\\d', ''))

# remove any data from conditions in which all three replicates weren't completed
sim_scores <- sim_scores %>% 
  ungroup() %>% 
  rowwise() %>% 
  filter(!((experiment %in% incomplete_exps) & (batch %in% incomplete_batches) & (condition %in% incomplete_conds))) %>% 
  ungroup()  

# check filtering
sim_scores %>% 
  mutate(batch = str_split(filename, "/", simplify = TRUE)[,10]) %>% 
  select(experiment, condition, replicate, analysis_condition, batch) %>% 
  distinct() %>% 
  group_by(batch, experiment, condition, analysis_condition) %>% 
  summarise(count = n()) %>% 
  filter(count != 3)

# add analysis conditions
sim_scores <- sim_scores %>% 
  left_join(conds, by=c("experiment", "condition", "replicate", "batch"))

# filter out any conditions were different virus was integrated and analysed
sim_scores <- sim_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
  ungroup()

sim_scores %>% 
  filter(experiment == "AAV") %>% 
  select(virus_name, analysis_virus) %>% 
  distinct()

# we had a fcov experiment in batch 1 and 2, so just keep batch 1
sim_scores <- sim_scores %>% 
  rowwise() %>% 
  filter(ifelse(experiment == "OTC-fcov", batch == "condition-breakdown-2", TRUE)) %>% 
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
    mutate(dist = dist+offset) %>% 
    ggplot(aes(x = dist, color = condition)) +
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





#### figure 2 function - make a row of four panels ####
figure2Row <- function(exp_name, exp_vars, combined_var_name) {
  
  
  # first part
  found <- found_scores %>% 
    filter(experiment == exp_name) %>% 
    rowwise() %>% 
    mutate(!!combined_var_name:= )
}

exp_name <-"OTC-juncs"
exp_vars <- c("p_overlap", "p_gap", "lambda_junction")
combined_var_name <- "junction properties"

names_df <- int_scores %>% 
  filter(experiment == exp_name) %>% 
  select(all_of(exp_vars))

for (var in colnames(names_df)) {
  names_df <- names_df %>% 
    mutate(!!var := paste0(var, ": ", !!sym(var))) 
}

mutate(!!exp_vars[[e]] := as.factor(!!sym(exp_vars[[e]])))



combineVarNames <- function(df, exp_name, exp_vars) {
  names_df <- df %>% 
    filter(experiment == exp_name) %>% 
    select(all_of(exp_vars))
  
  # include variable names in each row of names_df
  for (var in colnames(names_df)) {
    names_df <- names_df %>% 
      mutate(!!var := paste0(var, ": ", !!sym(var))) 
  }
  
  names_df <- names_df %>% 
    rowwise() %>% 
    mutate(label = paste0(across(everything()), collapse=", "))
  
  return(names_df)
}

combineVarNames(int_scores, "OTC-juncs", c("p_overlap", "p_gap", "lambda_junction"))

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


#### figure 3 ####
# runtime stuff

results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment2_time/time"
files <- list.files(results_dir, recursive = TRUE, pattern="_runtime.tsv")

coltypes <- cols(
  tool = col_character(),
  dataset = col_character(),
  sample = col_character(),
  replicate = col_double(),
  exit_value = col_double(),
  user_time = col_double(),
  system_time = col_double(),
  elapsed_time = col_character(),
  CPU = col_character(),
  shared_text = col_double(),
  unshared_data = col_double(),
  max_rss = col_double(),
  fs_inputs = col_double(),
  fs_outputs = col_double(),
  major_page_faults = col_double(),
  minor_page_faults = col_double(),
  swaps = col_double(),
  command = col_character()
)

times <- tibble(
  filename = file.path(results_dir, files),
  experiment = dirname(files),
  data = map(filename, ~read_tsv(., col_types=coltypes))
) %>% 
  unnest(data) %>% 
  mutate(elapsed_time = case_when(
    str_detect(elapsed_time, "\\d{0,2}:\\d{2}\\.\\d{2}") ~ lubridate::ms(elapsed_time),
    str_detect(elapsed_time, "\\d{0,2}:\\d{2}\\:\\d{2}") ~ lubridate::hms(elapsed_time),
    TRUE ~ lubridate::as.period(NA)
  )) %>% 
  mutate(elapsed_time = as.duration(elapsed_time)) %>% 
  rename(time_rep = replicate)

conds <- importSimulationConditions(results_dir)

times <- left_join(times, conds, by=c("experiment", "sample"))

times %>% 
  ggplot(aes(x = condition, y = elapsed_time, color=tool)) +
  geom_point() +
  facet_wrap(vars(experiment)) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

plots <- list()
for (e in unique(times$experiment)) {
  plots[[e]] <- times %>% 
    filter(experiment == e) %>% 
    ggplot(aes(x = condition, y = elapsed_time, color=tool)) +
    geom_point() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
          )  
}
plots['coverage']

p <- times %>% 
  ggplot(aes(x = condition, y = elapsed_time, color=tool)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom")

legend <- get_legend(p)

figure_3 <- cowplot::plot_grid(plotlist = plots, labels="AUTO")
figure_3 <- cowplot::plot_grid(figure_3, legend, ncol=1, rel_heights = c(1, 0.05))
print(figure_3)

cowplot::save_plot("plots/figure3_v1.pdf", figure_3)

