#!/usr/bin/env R

#### scoring functions ####
score_dists <- function(dists, group_cols, above_threshold = "fp", score_col = "d_coords_mean", threshold = 5) {
  print(glue("checking column {score_col} with threshold {threshold}"))
  return(dists %>% 
           mutate(score = case_when(
             !!ensym(score_col) == -1 ~ above_threshold,
             !!ensym(score_col) > threshold ~ above_threshold,
             !!ensym(score_col) <= threshold ~ "tp",
             TRUE ~ as.character(NA)
           )) %>% 
           group_by(across(one_of(c(group_cols, "score")))) %>% 
           tally() %>% 
           mutate(window = threshold) %>% 
           mutate(score_type = score_col) %>% 
           ungroup())
  
}

all_scores <- function(found_dists, sim_dists, group_cols, threshold, score_col) {
  
  found <- score_dists(found_dists, group_cols, above_threshold ="fp", score_col, threshold) %>% 
    filter(score == "fp")
  sim <- score_dists(sim_dists, group_cols, above_threshold ="fn", score_col, threshold)
  
  
  return(
    bind_rows(found, sim) %>% 
      pivot_wider(names_from = score, values_from = n, values_fill = 0) %>% 
      ungroup() %>% 
      rowwise() %>% 
      mutate(PPV = tp / (tp + fp)) %>% 
      mutate(TPR = tp / (tp + fn)) %>% 
      ungroup()
  )
}



#### import data for AAV-OTC experiment ####

import_AAV_OTC <- function() {
  
  scoring_dist <- 5
  jaccard_pad <- 0
  dist_add <- 0.5
  plot_score_type <- "d_shortest"
  analysis_condition_plot <- "isling1"
  tool_order <- c("isling", "Polyidus", "Seeksv", "ViFi", "VSeq-Toolkit")
  condition_plot <- "cond0"

AAV_OTC_dir <- "../out/experiment1_OTC_chr1/AAV-OTC"
score_types <- c("shortest")

sim_conds <- importSimulationConditions(AAV_OTC_dir)

analysis_conds <- importAnalysisConditions(AAV_OTC_dir) %>% 
  select(-one_of("adapter_1", "adapter_2"))


jac <- importJaccardExperiment(AAV_OTC_dir) %>% 
  mutate(pad = as.integer(str_extract(sim_file, "(?<=pad)\\d+"))) %>% 
  mutate(analysis_tool = str_replace(analysis_tool, "analysis", "isling")) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, "analysis", "isling")) %>% 
  mutate(analysis_condition_short = str_replace(analysis_condition_short, "analysis", "isling"))

jac <- jac %>% 
  filter(!is.na(host_name))

#import distances from each found integration to nearest simulated integration
found_scores <- importNearestSimToFound(AAV_OTC_dir) 


# only use postprocessed data from our pipeline
found_scores <- found_scores %>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(tool = str_replace(analysis_condition, "\\d+", "")) %>% 
  filter(experiment == "AAV" | experiment == "OTC")

# import distances from each simulated integration to the nearest found integration
sim_scores <- importNearestFoundToSim(AAV_OTC_dir) 

sim_scores <- sim_scores%>% 
  filter(case_when(
    !str_detect(analysis_condition, "pipeline") ~ TRUE,
    post ~ TRUE,
    TRUE ~ FALSE
  )) %>%  
  mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
  mutate(tool = str_replace(analysis_condition, "\\d+", "")) %>% 
  filter(experiment == "AAV" | experiment == "OTC")


plot_cols <- paste0("d_", score_types)

group_cols<- c("unique", "analysis_condition", "experiment", "condition", "replicate",  "analysis_host", "analysis_virus", "post", "tool")

dists = c(0, 2, 5, 10, 30, 60, 100)


int_scores <- tibble()

for (s in plot_cols) {
  
  for (t in dists) {
    
    int_scores <- bind_rows(
      int_scores, all_scores(found_scores, sim_scores, 
                             group_cols = group_cols, threshold = t, score_col = s) 
    )
  }
  
}

jac_plot <- jac %>% 
  filter(case_when(
    str_detect(analysis_condition_short, "analysis|pipeline|isling") ~ analysis_condition_short == analysis_condition_plot,
    TRUE ~ TRUE
  )) %>% 
  filter(pad == jaccard_pad) %>% 
  filter(condition == condition_plot) %>% 
  mutate(tool = case_when(
    analysis_tool == "polyidus" ~ "Polyidus",
    analysis_tool == "seeksv" ~ "Seeksv",
    analysis_tool == "vifi" ~ "ViFi",
    analysis_tool == "vseq-toolkit" ~ "VSeq-Toolkit",
    TRUE ~ analysis_tool
  )) %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, tool_order)) %>% 
  mutate(experiment = str_replace(experiment, "-harder", ""))

sim_plot <- sim_scores %>% 
  filter(condition == condition_plot) %>% 
  filter(case_when(
    str_detect(analysis_condition, "analysis|pipeline|isling") ~ analysis_condition == analysis_condition_plot,
    TRUE ~ TRUE
  )) %>% 
  mutate(tool = case_when(
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "Seeksv",
    tool == "vifi" ~ "ViFi",
    tool == "vseq-toolkit" ~ "VSeq-Toolkit",
    TRUE ~ tool
  ))  %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, tool_order)) %>% 
  mutate(experiment = str_replace(experiment, "-harder", ""))

found_plot <- found_scores %>% 
  filter(condition == condition_plot) %>% 
  filter(case_when(
    str_detect(analysis_condition, "analysis|pipeline|isling") ~ analysis_condition == analysis_condition_plot,
    TRUE ~ TRUE
  )) %>% 
  mutate(tool = case_when(
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "Seeksv",
    tool == "vifi" ~ "ViFi",
    tool == "vseq-toolkit" ~ "VSeq-Toolkit",
    TRUE ~ tool
  ))  %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, tool_order)) %>% 
  mutate(experiment = str_replace(experiment, "-harder", ""))



scores_plot <- int_scores %>% 
  filter(condition == condition_plot) %>% 
  filter(window == scoring_dist) %>% 
  filter(score_type == plot_score_type) %>% 
  filter(case_when(
    str_detect(analysis_condition, "analysis|pipeline|isling") ~ analysis_condition == analysis_condition_plot,
    TRUE ~ TRUE
  )) %>% 
  mutate(tool = case_when(
    tool == "polyidus" ~ "Polyidus",
    tool == "seeksv" ~ "Seeksv",
    tool == "vifi" ~ "ViFi",
    tool == "vseq-toolkit" ~ "VSeq-Toolkit",
    TRUE ~ tool
  )) %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, tool_order)) %>% 
  mutate(experiment = str_replace(experiment, "-harder", "")) 



return(list(jac, found_scores, sim_scores, jac_plot, sim_plot, found_plot, scores_plot))
}


import_OTC_conds <- function() {
  
  
  OTC_conds_dir <-"../out/experiment1_OTC_chr1/OTC-condition-breakdown/"
  
  jac <- importJaccardExperiment(OTC_conds_dir) %>% 
    mutate(pad = as.integer(str_extract(sim_file, "(?<=pad)\\d+"))) %>% 
    mutate(analysis_tool = str_replace(analysis_tool, "analysis", "isling"))
  
  jac <- jac %>% 
    filter(pad != 2)
  
  
  #import distances from each found integration to nearest simulated integration
  found_scores <- importNearestSimToFound(OTC_conds_dir) 
  
  
  # only use postprocessed data from our pipeline
  found_scores <- found_scores %>% 
    filter(case_when(
      !str_detect(analysis_condition, "pipeline") ~ TRUE,
      post ~ TRUE,
      TRUE ~ FALSE
    )) %>% 
    mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
    mutate(tool = str_replace(analysis_condition, "\\d+", "")) 
  
  # import distances from each simulated integration to the nearest found integration
  sim_scores <- importNearestFoundToSim(OTC_conds_dir) 
  
  sim_scores <- sim_scores%>% 
    filter(case_when(
      !str_detect(analysis_condition, "pipeline") ~ TRUE,
      post ~ TRUE,
      TRUE ~ FALSE
    )) %>%  
    mutate(analysis_condition = str_replace(analysis_condition, 'analysis', 'isling')) %>% 
    mutate(tool = str_replace(analysis_condition, "\\d+", "")) 
  
  group_cols<- c("unique", "analysis_condition", "experiment", "condition", "replicate",  "analysis_host", "analysis_virus", "post", "tool")
  
  dists = c(5)
  
  int_scores <- tibble()
  int_scores <- bind_rows(
    int_scores, all_scores(found_scores, sim_scores, 
                           group_cols = group_cols, threshold = 5, score_col = "d_shortest") 
  )
  
  
  int_scores <- int_scores %>% 
    mutate(tool = str_replace(tool, "pipeline", "isling")) %>% 
    mutate(tool = str_replace(tool, "seeksv", "Seeksv")) %>% 
    mutate(tool = str_replace(tool, "polyidus", "Polyidus")) %>% 
    mutate(tool = str_replace(tool, "vifi", "ViFi")) %>% 
    mutate(tool = str_replace(tool, "vseq-toolkit", "VSeq-Toolkit"))
  
  conds <- jac %>% 
    distinct(condition, fcov, host_name)
  
  int_scores <- int_scores %>% 
    left_join(conds, by="condition")
  
  return(list(jac, int_scores, found_scores, sim_scores))
  
}

makeFigue2 <- function(isling_condition) {
  
  jaccard_pad <- 0
  
  
  jac_summary <- jac %>% 
    filter(pad == jaccard_pad) %>% 
    filter(case_when(
      str_detect(analysis_condition, "analysis") ~ str_detect(analysis_condition, isling_condition),
      TRUE ~ TRUE
    )) %>% 
    mutate(analysis_tool = case_when(
      analysis_tool == "polyidus" ~ "Polyidus",
      analysis_tool == "seeksv" ~ "Seeksv",
      analysis_tool == "vifi" ~ "ViFi",
      analysis_tool == "vseq-toolkit" ~ "VSeq-Toolkit",
      TRUE ~ analysis_tool
    )) %>% 
    mutate(analysis_tool = as.factor(analysis_tool)) %>% 
    mutate(analysis_tool = forcats::fct_relevel(analysis_tool, tool_order))  %>% 
    group_by(condition, fcov, analysis_tool, analysis_condition, pad, host_name) %>% 
    summarise(n = n(),
              jac_mean = mean(jaccard),
              jac_sem = sd(jaccard)/sqrt(n),
              jac_min = jac_mean - jac_sem,
              jac_max = jac_mean + jac_sem) 
  
  jac_plot <- jac_summary %>% 
    ggplot(aes(x = fcov, y = jac_mean, color = analysis_tool)) +
    geom_point() +
    geom_line() +
    geom_linerange(aes(ymin = jac_min, ymax = jac_max)) +
    facet_grid(cols = vars(host_name)) +
    scale_x_log10() +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    ylim(0, 1) +
    ylab("Jaccard")
  
  
  PPV_summary <- int_scores %>% 
    filter(case_when(
      str_detect(analysis_condition, "isling") ~ str_detect(analysis_condition, isling_condition),
      TRUE ~ TRUE
    )) %>% 
    mutate(tool = as.factor(tool)) %>% 
    mutate(tool = forcats::fct_relevel(tool, tool_order)) %>% 
    group_by(condition, fcov, tool, analysis_condition, host_name) %>% 
    summarise(n = n(),
              PPV_mean = mean(PPV),
              PPV_sem = sd(PPV)/sqrt(n),
              PPV_min = PPV_mean - PPV_sem,
              PPV_max = PPV_mean + PPV_sem) 
  
  PPV_plot <- PPV_summary %>% 
    ggplot(aes(x = fcov, y = PPV_mean, color = tool)) +
    geom_point() +
    geom_line() +
    geom_linerange(aes(ymin = PPV_min, ymax = PPV_max)) +
    facet_grid(cols = vars(host_name)) +
    scale_x_log10() +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_blank()) +
    ylim(0, 1) +
    ylab("PPV")
  
  
  TPR_summary <- int_scores %>% 
    filter(case_when(
      str_detect(analysis_condition, "isling") ~ str_detect(analysis_condition, isling_condition),
      TRUE ~ TRUE
    )) %>% 
    mutate(tool = as.factor(tool)) %>% 
    mutate(tool = forcats::fct_relevel(tool, tool_order)) %>% 
    group_by(condition, fcov, tool, analysis_condition, host_name) %>% 
    summarise(n = n(),
              TPR_mean = mean(TPR),
              TPR_sem = sd(TPR)/sqrt(n),
              TPR_min = TPR_mean - TPR_sem,
              TPR_max = TPR_mean + TPR_sem) 
  
  TPR_plot <- TPR_summary %>% 
    ggplot(aes(x = fcov, y = TPR_mean, color = tool)) +
    geom_point() +
    geom_line() +
    geom_linerange(aes(ymin = TPR_min, ymax = TPR_max)) +
    facet_grid(cols = vars(host_name)) +
    scale_x_log10() +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_blank()) +
    ylim(0, 1) +
    ylab("TPR") +
    xlab("Fold coverage")
  
  
  p <- cowplot::plot_grid(jac_plot, PPV_plot, TPR_plot, labels="AUTO", ncol = 1)
  legend <- cowplot::get_legend(jac_plot + theme(legend.position = "bottom"))
  p2 <- cowplot::plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.05))
  print(p2)
  cowplot::save_plot(glue::glue("../figures/figure3.{scoring_dist}dist.pdf"), p2)
  cowplot::save_plot(glue::glue("../figures/figure3.{scoring_dist}dist.png"), p2)
  
  return(list(p2, jac_summary, PPV_summary, TPR_summary))
  
  
}
