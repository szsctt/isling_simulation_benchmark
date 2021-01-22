# functions for making figure 2

score_window <- 5
coords_score_type_plot <- "coords_mean"
dist_plot_offset <- 0.5
facet_scales <- "free_y"
# https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
offset <- dist_plot_offset
num_y_breaks <- 2

#### import int scores ####

importAllIntScoreExperiments <- function(results_dir, folders, score_window, coords_score_type) {
  
  # import data from all folders
  int_scores <- tibble(
    batch = folders,
    f = file.path(results_dir, folders)
  ) %>% 
    mutate(data = map(f, importIntScoreExperiment)) %>% 
    unnest(data)
  
  # filter for desired window for scoring and coordinate type
  int_scores <- int_scores %>% 
    filter(window == score_window) %>% 
    filter(coords_score_type == coords_score_type)
  
  # only use postprocessed data from our pipeline, and rename
  int_scores <- int_scores %>% 
    filter(case_when((analysis_tool != "pipeline") ~ TRUE,
                     post ~ TRUE,
                     TRUE ~ FALSE
    )) %>% 
    mutate(analysis_condition = str_replace(analysis_condition, 'analysis|pipeline', 'isling')) %>% 
    mutate(analysis_condition = str_replace(analysis_tool, 'analysis|pipeline', 'isling')) 
  
  
  # for 'AAV' experiment, filter only for results where we used the correct viral reference
  int_scores <- int_scores %>% 
    rowwise() %>% 
    filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
    ungroup()
  
  # check that we've got three replicates for each condition in each experiment, for each tool
  n_per_cond <- int_scores %>% 
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
  
  incomplete <- list('exps' = incomplete_exps, 'conds' = incomplete_conds, 'batches' = incomplete_batches)
  
  return(list('int_scores' = int_scores, 'incomplete' = incomplete))
}

importAllFoundScores <- function(results_dir, folders, incomplete, score_window, coords_score_type) {
  
  # import data
  found_scores <- tibble(
    batch = folders,
    f = map_chr(batch, ~file.path(results_dir, .)),
    data = map(f, ~importNearestSimToFound(., score_window, coords_score_type))
  ) %>% 
    unnest(data)
  
  # only use postprocessed data from our pipeline
  found_scores <- found_scores %>% 
    filter(case_when(
      !str_detect(analysis_condition, "pipeline") ~ TRUE,
      post ~ TRUE,
      TRUE ~ FALSE
    )) %>% 
    filter(score_dist == score_window) %>% 
    filter(score_type == coords_score_type_plot) %>% 
    mutate(analysis_condition = str_replace(analysis_condition, 'analysis|pipeline', 'isling'))
  
  
  # remove any data from conditions in which all three replicates weren't completed
  found_scores <- found_scores %>% 
    ungroup() %>% 
    rowwise() %>% 
    filter(!((experiment %in% incomplete[['exps']]) & (batch %in% incomplete[['batches']]) & (condition %in% incomplete[['conds']]))) %>% 
    ungroup()  
  
  # import analysis conditions and add to tibble
  conds  <- tibble(
    batch = folders,
    exp_dir = map(batch, ~file.path(results_dir, .))
  )  %>% 
    mutate(conds = map(exp_dir, ~importSimulationConditions(.))) %>% 
    unnest(conds)
  
  found_scores <- found_scores %>% 
    left_join(conds, by=c("experiment", "condition", "replicate", "batch"))
  
  # filter out any conditions where different virus was integrated and analysed
  found_scores <- found_scores %>% 
    rowwise() %>% 
    filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
    ungroup()
  
  return(found_scores)
  
}

importAllSimScores <- function(results_dir, folders, incomplete, score_window, coords_score_type) {
  
  # import data
  sim_scores <- tibble(
    batch = folders,
    exp_dir = map_chr(batch, ~file.path(results_dir, .))
  ) %>% 
    mutate(data = map(exp_dir, ~importNearestFoundToSim(., score_window, coords_score_type_plot))) %>% 
    unnest(data)
  
  # for our data, only use postprocessed data, and change names
  sim_scores <- sim_scores %>% 
    filter(case_when(
      !str_detect(analysis_condition, "pipeline") ~ TRUE,
      post ~ TRUE,
      TRUE ~ FALSE
    )) %>% 
    filter(score_dist == score_window) %>% 
    filter(score_type == coords_score_type_plot) %>% 
    mutate(analysis_condition = str_replace(analysis_condition, 'analysis|pipeline', 'isling')) 
  
  # remove any data from conditions in which all three replicates weren't completed
  sim_scores <- sim_scores %>% 
    ungroup() %>% 
    rowwise() %>% 
    filter(!((experiment %in% incomplete[['exps']]) & (batch %in% incomplete[['batches']]) & (condition %in% incomplete[['conds']]))) %>% 
    ungroup()  

  # import analysis conditions and add to tibble
  conds  <- tibble(
    batch = folders,
    exp_dir = map(batch, ~file.path(results_dir, .))
  )  %>% 
    mutate(conds = map(exp_dir, ~importSimulationConditions(.))) %>% 
    unnest(conds)
    
  # add analysis conditions
  sim_scores <- sim_scores %>% 
    left_join(conds, by=c("experiment", "condition", "replicate", "batch"))
  
  
  #filter out any conditions were different virus was integrated and analysed
  sim_scores <- sim_scores %>% 
    rowwise() %>% 
    filter(ifelse(experiment == "AAV", virus_name == analysis_virus, TRUE)) %>% 
    ungroup()
  
  return(sim_scores)
  
}



#### plotting functions ####


# function to make plot of found scores
distPlot <- function(df,  exp_name, combined_var_name){
  filt_df <- df  %>% 
    filter(experiment == exp_name) %>% 
    mutate(dist = dist+offset) 
  
  if (nrow(filt_df) == 0) {
    stop(glue("no rows found for experiment {exp_name}"))
  }
  
  short_var_name <- paste0(combined_var_name, "_short")
  filt_df <- filt_df %>% 
    mutate(!!combined_var_name := combineVarNames(filt_df, exp_name)) %>% 
    mutate(!!short_var_name := shortCombinedVarNames(filt_df, exp_name))
  
  p <- filt_df %>% 
    mutate(!!combined_var_name := as.factor(!!sym(combined_var_name))) %>% 
    mutate(!!short_var_name := as.factor(!!sym(short_var_name))) %>% 
    ggplot(aes(x = dist, color = analysis_condition)) +
    geom_freqpoly(bins = 500) +
    scale_x_log10() +
    facet_grid(rows = vars(!!sym(short_var_name))) +
    theme_classic() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      legend.position = "none"
    )  +
    scale_y_continuous(breaks = scales::pretty_breaks(n = num_y_breaks)) +
    labs(x="Distance", y="Count", color="tool", linetype=combined_var_name)
  
  return(p)
}


scorePlot <- function(int_score_df,  exp_name, combined_var_name) {
  filt_df <- int_score_df  %>% 
    filter(experiment == exp_name)
  
  if (nrow(filt_df) == 0) {
    stop(glue("no rows found for experiment {exp_name}"))
  }
  
  # add names for plotting
  short_var_name <- paste0(combined_var_name, "_short")
  filt_df <- filt_df %>% 
    mutate(!!combined_var_name := combineVarNames(filt_df, exp_name))  %>% 
    mutate(!!short_var_name := shortCombinedVarNames(filt_df, exp_name))
  
  p <- filt_df %>% 
    mutate(!!combined_var_name := as.factor(!!sym(combined_var_name))) %>% 
    mutate(!!short_var_name := as.factor(!!sym(short_var_name))) %>% 
    ggplot(aes(x = PPV, y = TPR, shape = !!sym(short_var_name), color = analysis_tool)) +
    geom_point(alpha = 0.5) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.3, 0.6))  +
    guides(color = FALSE) +
    labs(shape=combined_var_name)
}

chrPlot <- function(found_score_df, exp_name, combined_var_name){
  filt_df <- found_score_df  %>% 
    filter(experiment == exp_name)
  if (nrow(filt_df) == 0) {
    stop(glue("no rows found for experiment {exp_name}"))
  }
  
  # add names for plotting
  short_var_name <- paste0(combined_var_name, "_short")
  filt_df <- filt_df %>% 
    mutate(!!combined_var_name := combineVarNames(filt_df, exp_name)) %>% 
    mutate(!!short_var_name := shortCombinedVarNames(filt_df, exp_name))
  
  p <- filt_df %>% 
    mutate(correct_chr = (chr == host_name)) %>% 
    mutate(!!combined_var_name := as.factor(!!sym(combined_var_name))) %>% 
    mutate(!!short_var_name := as.factor(!!sym(short_var_name))) %>% 
    ggplot(aes(x = !!sym(short_var_name), fill = correct_chr)) +
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
} 

# only make plots we want to look at
makeFigure2 <- function(to_plot_list, df_list, output_filename, bheight, bwidth) {
  selected_plots <- list()
  for (e in names(to_plot_list)) {
    print(glue("generating plots for {e}"))
    selected_plots[[paste0(e, "_sim")]] <- distPlot(df_list$sim, e,  to_plot_list[[e]])
    selected_plots[[paste0(e, "_found")]] <- distPlot(df_list$found,  e,  to_plot_list[[e]]) 
    selected_plots[[paste0(e, "_chr")]] <- chrPlot(df_list$found,  e,  to_plot_list[[e]])
    selected_plots[[paste0(e, "_scores")]] <- scorePlot(df_list$scores,  e,  to_plot_list[[e]])
  }
  
  p <- cowplot::plot_grid(plotlist=selected_plots, ncol=4, labels="AUTO")
  
  cowplot::save_plot(output_filename, p, base_width=bwidth, base_height = bheight)
  
  return(p)
}
