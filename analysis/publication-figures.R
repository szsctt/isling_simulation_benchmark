library(tidyverse)
library(cowplot)

source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/sim_functions.R")

results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment1_OTC_chr1"

score_window <- 5
coords_score_type <- "coords_mean"
dist_plot_offset <- 0.5


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
found_scores <- importNearestSimToFound(exp_dir) 

# only use postprocessed data from our pipeline
found_scores <- found_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling'))

# import distances from each simulated integration to the nearest found integration
sim_scores <- importNearestFoundToSim(exp_dir) 

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
  filter(score_type == coords_score_type) %>% 
  filter(experiment == "OTC-harder" | experiment == "AAV-harder") %>% 
  filter(score_dist == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE)) %>%  # keep isling2 only
  mutate(tool = str_match(analysis_condition, "(.+)\\d+")[,2])



# double check filtering
plot_found_scores %>% 
  group_by(experiment, tool, score_dist, score_type, condition, replicate) %>% 
  summarise(tp = sum(score == "tp"),
            fp = sum(score == "fp"))

plot_sim_scores <- sim_scores %>% 
  filter(score_type == coords_score_type) %>% 
  filter(experiment == "OTC-harder" | experiment == "AAV-harder") %>% 
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
  filter(coords_score_type == coords_score_type) %>% 
  filter(experiment == "OTC-harder" | experiment == "AAV-harder") %>% 
  filter(window == score_window) %>% 
  filter(ifelse(str_detect(analysis_condition, "isling"), str_detect(analysis_condition, isling_to_keep), TRUE))

plot_int_scores %>% 
  select(experiment, analysis_tool, window, coords_score_type, condition, replicate, tp, fp, fn, PPV, TPR)


# https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
offset <- dist_plot_offset
num_y_breaks <- 2

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


OTC_scores <- plot_int_scores %>% 
  filter(experiment == "OTC-harder") %>% 
  ggplot(aes(x = PPV, y = TPR, colour = analysis_tool)) +
  geom_point(show.legend = FALSE, alpha = 0.5) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Positive predictive value") +
  ylab("True positive rate") +
  theme_classic()


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

AAV_scores <- plot_int_scores %>% 
  filter(experiment == "AAV-harder") %>% 
  ggplot(aes(x = PPV, y = TPR, colour = analysis_tool)) +
  geom_point(show.legend = FALSE, alpha = 0.5) +
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Positive predictive value") +
  ylab("True positive rate") +
  theme_classic()

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

#### figure 2 ####

# import scored integration sites for condition-breakdown'
exp_dir <- file.path(results_dir, "condition-breakdown")
int_scores <- importIntScoreExperiment(exp_dir) 

# only use postprocessed data from our pipeline
int_scores <- int_scores %>% 
  filter(case_when((analysis_tool != "pipeline") ~ TRUE,
                   post ~ TRUE,
                   TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '0', ''))  

# intial plot of scored ints for each condition
int_scores %>% 
  filter(str_detect(experiment, "harder")) %>% 
    ggplot(aes(x = PPV, y = TPR, color = analysis_tool, shape = host_name)) +
    geom_point() +
    xlim(0, 1) +
    ylim(0, 1) +
    facet_grid(vars(experiment)) 


int_scores <- int_scores %>% 
  filter(window == score_window) %>% 
  filter(coords_score_type == "coords_mean")

#import distances from each found integration to nearest simulated integration
found_scores <- importNearestSimToFound(exp_dir) 

# only use postprocessed data from our pipeline
found_scores <- found_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '0', ''))   %>% 
  filter(score_dist == score_window) %>% 
  filter(score_type == "coords_mean")


# add analysis conditions
conds <- importSimulationConditions(exp_dir)

found_scores <- found_scores %>% 
  left_join(conds, by=c("experiment", "condition", "replicate"))

found_scores %>% 
  filter(str_detect(experiment, "harder")) %>% 
  mutate(dist = dist+dist_plot_offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition, linetype = analysis_host)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(analysis_condition), cols = vars(experiment), scales = "free_y")

found_scores %>% 
  filter(str_detect(experiment, "harder")) %>% 
  mutate(dist = dist+dist_plot_offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition, linetype = analysis_host)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(cols = vars(experiment), rows = vars(host_name), scales = "free_y")


# import distances from each simulated integration to the nearest found integration
sim_scores <- importNearestFoundToSim(exp_dir) 

sim_scores <- sim_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, '0', ''))  %>% 
  filter(score_dist == score_window) %>% 
  filter(score_type == "coords_mean")

sim_scores <- sim_scores %>% 
  left_join(conds, by=c("experiment", "condition", "replicate"))

sim_scores %>% 
  filter(str_detect(experiment, "harder")) %>% 
  mutate(dist = dist+dist_plot_offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition, linetype = analysis_host)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(rows = vars(analysis_condition), cols = vars(experiment), scales = "free_y")

sim_scores %>% 
  filter(str_detect(experiment, "harder")) %>% 
  mutate(dist = dist+dist_plot_offset) %>% 
  ggplot(aes(x = dist, colour = analysis_condition, linetype = analysis_host)) +
  geom_freqpoly(bins = 100, show.legend = FALSE) +
  scale_x_log10() +
  facet_grid(cols = vars(experiment), rows = vars(host_name), scales = "free_y")
