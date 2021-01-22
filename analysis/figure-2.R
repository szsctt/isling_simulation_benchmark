#### imports and constants ####

library(tidyverse)
library(cowplot)


source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/sim_functions.R")

results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment1_OTC_chr1"
folders <- c("condition-breakdown-1", "condition-breakdown-2")

score_window <- 5
coords_score_type_plot <- "coords_mean"
dist_plot_offset <- 0.5
facet_scales <- "free_y"

# https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
offset <- dist_plot_offset
num_y_breaks <- 2

#### import scores ####

# import scored integration sites for condition-breakdown'
int_scores <- tibble(
    batch = folders,
    f = file.path(results_dir, folders)
  ) %>% 
  mutate(data = map(f, importIntScoreExperiment)) %>% 
  unnest(data)


int_scores <- int_scores %>% 
  filter(window == score_window) %>% 
  filter(coords_score_type == coords_score_type_plot)

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

int_scores %>% 
  filter(experiment == "OTC-fcov") %>% 
  select(batch, condition, fcov) %>% 
  distinct()

# for consistency with other data frames, add virus_name and host_name columns
#int_scores <- int_scores %>% 
#  mutate(virus_name = analysis_virus) %>% 
#  mutate(host_name = analysis_host)


#### import found scores ####
#import distances from each found integration to nearest simulated integration
found_scores <- tibble(
  batch = folders,
  f = map_chr(batch, ~file.path(results_dir, .)),
  data = map(f, ~importNearestSimToFound(., score_window, coords_score_type_plot))
) %>% 
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

# filter out any conditions where different virus was integrated and analysed
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


#### import sim scores ####
sim_scores <- tibble(
  batch = folders,
  exp_dir = map_chr(batch, ~file.path(results_dir, .))
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

#### filtering data for display ####

# our initial plots (one per experiment) show that some experiments have 
# too many conditions for an informative plot

for (e in unique(int_scores$experiment)) {
  vars <- getExpVars(int_scores, e)
  print(glue("experiment: {e}, variables: {vars}"))
}


# for OTC-host_deletion
int_scores %>% 
  filter(experiment == "OTC-host_deletion") %>% 
  select(condition, p_host_deletion, lambda_host_deletion) %>% 
  distinct()

# we have p_host_deletion = (0, 0.25, 0.5, 0.7, 1) and lambda_host_deleltion (10, 100, 1000, 10000).
# keep only p_host_deletion = (0,  0.5,  1) and lambda_host_deleltion (100, 10000)
# we can also remove lambda_host_deleltion =  100 if p_host_deletion = 0, since if there is no host deletion
# it doesn't matter what the lambda was
filter_OTC_host_deletion <- function(df) {
  return(df %>% 
           ungroup() %>% 
           rowwise() %>% 
           filter(case_when(
             experiment != 'OTC-host_deletion' ~ TRUE,
             p_host_deletion %in% c(0.5, 1) ~ TRUE,
             lambda_host_deletion %in% c(100, 10000) ~ TRUE,
             p_host_deletion == 0 & lambda_host_deletion == 100 ~ TRUE,
             TRUE ~ FALSE
           )) %>% 
           ungroup()
  )
}
nrow(int_scores)
nrow(filter_OTC_host_deletion(int_scores))

#nrow(sim_scores)
#nrow(filter_OTC_host_deletion(sim_scores))

plot_int_scores <- filter_OTC_host_deletion(int_scores)
plot_sim_scores <- filter_OTC_host_deletion(sim_scores)
plot_found_scores <- filter_OTC_host_deletion(found_scores)

plot_int_scores %>% 
  filter(experiment == "OTC-host_deletion") %>% 
  select(condition, p_host_deletion, lambda_host_deletion) %>% 
  distinct()



# OTC-juncs: p_overlap,  p_gap, lambda_junction
int_scores %>% 
  filter(experiment == "OTC-juncs") %>% 
  select(condition, p_overlap, p_gap, lambda_junction) %>% 
  distinct()

# we have p_overlap = (0, 0.25, 0.5), p_gap = (0, 0.25, 0.5), and lambda_junction (1, 4, 10).
# split each level of lambda_junction into a different level
int_scores %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(experiment = case_when(
    experiment != "OTC-juncs" ~ experiment,
    experiment == "OTC-juncs" ~ paste0(experiment, "_", lambda_junction)
  )) %>% 
  head()

label_OTC_juncs <- function(df) {
  return(df %>% 
           ungroup() %>% 
           rowwise() %>% 
           mutate(experiment = case_when(
             experiment != "OTC-juncs" ~ experiment,
             experiment == "OTC-juncs" ~ paste0(experiment, "_", lambda_junction)
           )) %>% 
           ungroup()
  )
}
plot_int_scores <- label_OTC_juncs(plot_int_scores)
plot_sim_scores <- label_OTC_juncs(plot_sim_scores)
plot_found_scores <- label_OTC_juncs(plot_found_scores)


# OTC-rearrange-epi:p_rearrange,  p_delete, epi_num
int_scores %>% 
  filter(experiment == "OTC-rearrange-epi") %>% 
  select(condition, p_rearrange, p_delete, epi_num) %>% 
  distinct()

int_scores %>% 
  filter(experiment == "OTC-rearrange-epi") %>% 
  pull(p_rearrange) %>% 
  unique()

int_scores %>% 
  filter(experiment == "OTC-rearrange-epi") %>% 
  pull(p_delete) %>% 
  unique()

int_scores %>% 
  filter(experiment == "OTC-rearrange-epi") %>% 
  pull(epi_num) %>% 
  unique()

# similarly, split into separate eperiments based on episome number
label_OTC_rearrange_epi <- function(df) {
  return(df %>% 
           ungroup() %>% 
           rowwise() %>% 
           mutate(experiment = case_when(
             experiment != "OTC-rearrange-epi" ~ experiment,
             experiment == "OTC-rearrange-epi" ~ paste0(experiment, "_", epi_num)
           )) %>% 
           ungroup()
  )
}

plot_int_scores <- label_OTC_rearrange_epi(plot_int_scores)
plot_sim_scores <- label_OTC_rearrange_epi(plot_sim_scores)
plot_found_scores <- label_OTC_rearrange_epi(plot_found_scores)

unique(plot_int_scores$experiment)
unique(plot_sim_scores$experiment)
unique(plot_found_scores$experiment)


### OTC-rearrange-frags: p_rearrange, p_delete, lambda_split
# again, split into different eperiments based on lambda_split value
label_OTC_rearrange_frags <- function(df) {
  return(df %>% 
           ungroup() %>% 
           rowwise() %>% 
           mutate(experiment = case_when(
             experiment != "OTC-rearrange-frags" ~ experiment,
             experiment == "OTC-rearrange-frags" ~ paste0(experiment, "_", lambda_split)
           )) %>% 
           ungroup()
  )
}

plot_int_scores <- label_OTC_rearrange_frags(plot_int_scores)
plot_sim_scores <- label_OTC_rearrange_frags(plot_sim_scores)
plot_found_scores <- label_OTC_rearrange_frags(plot_found_scores)


nrow(plot_sim_scores)
plot_sim_scores %>% 
  pull(experiment) %>% 
  unique()

plot_sim_scores %>% 
  filter(experiment == "OTC-rearrange-epi_100") %>% 
  select(experiment, condition, epi_num, p_rearrange, p_delete) %>% 
  distinct()

getExpVars(plot_sim_scores, "OTC-rearrange-epi_100")


#### plotting functions ####



#getExpVars(int_scores, "OTC-juncs")
#getExpVars(sim_scores, "OTC-juncs")
#getExpVars(found_scores, "OTC-juncs")


#combineVarNames(int_scores, "OTC-juncs")
#combineVarNames(int_scores, "OTC-epi")
#combineVarNames(sim_scores, "OTC-epi")
#combineVarNames(found_scores, "OTC-epi")

#shortCombinedVarNames(int_scores, "OTC-juncs")


combineVarNames(int_scores, "chromosomes")
shortCombinedVarNames(int_scores, "chromosomes")

conds %>% 
  filter(experiment == "chromosomes") %>% 
  distinct()


int_scores %>% 
  pull(host) %>% 
  unique()

int_scores %>% 
  pull(host_name) %>% 
  unique()



#print(distPlot(found_scores,"OTC-juncs", "junction_properties"))
#print(distPlot(sim_scores,"OTC-juncs", "junction_properties"))
#print(distPlot(found_scores, "AAV", "virus"))
#print(distPlot(sim_scores, "AAV", "virus"))




#print(scorePlot(int_scores, "OTC-juncs", "junction properties"))
#print(scorePlot(int_scores,"AAV", "viruses"))
#print(scorePlot(int_scores,"chromosomes", "chromosome"))



#print(chrPlot(found_scores, "OTC-juncs", "junction properties"))

figureRow <- function(found_df, sim_df, int_df, exp_name, combined_var_name) {
  
  # make individual plots
  found <- distPlot(found_df,  exp_name,  combined_var_name)
  sim <- distPlot(found_df, exp_name,  combined_var_name)
  
  # make legend for found and sim
  
  chr <- chrPlot(found_df,  exp_name,  combined_var_name)
  scores <- scorePlot(int_df,  exp_name,  combined_var_name)
  
  p <- cowplot::plot_grid(sim, found, chr, scores, nrow=2)
  
  return(p)
  
  # add title
  title <- ggdraw() +
    draw_label(combined_var_name, fontface = "bold", x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  return(p)
  #return(plot_grid(title, p, ncol = 1, rel_heights = 0.01, 1))
  
}


row <- figureRow(found_scores, sim_scores, int_scores, "AAV", "viruses")
#print(row)
cowplot::save_plot(glue("plots/2_viruses.pdf"), row, base_width = 7, base_height=7)


# re-plot to see how data looks now
allplots <- list()
for (e in unique(plot_int_scores$experiment)) {
  print(glue("plotting {e}"))
  row <- figureRow(plot_found_scores, plot_sim_scores, plot_int_scores,  e, e)
  cowplot::save_plot(glue("plots/2_{e}.pdf"), row, base_width = 7, base_height=7)
  allplots[[e]] <- row
}


df_list <- list(
  'sim' = sim_scores,
  'found' = found_scores,
  'scores' = int_scores
)
to_plot <- list("chromosomes" = "chromosome",
                "OTC-fcov" = "fold coverage",
                "OTC-whole" = "p(whole)",
                "OTC-rearrange-epi_100"= 'p(rearrange)/p(delete)',
                "OTC-juncs_4" = "p(overlap)/p(gap)"
)
makeFigure2(to_plot, df_list, "plots/figure2_v1.pdf", 12, 12)

# we still have too many conditions in OTC-rearrange-epi_100 and OTC-juncs_4
to_plot <- list("chromosomes" = "chromosome",
                "OTC-fcov" = "fold coverage",
                "OTC-whole" = "p(whole)",
                "OTC-rearrange"= 'p(rearrange)'
)
p <- makeFigure2(to_plot, df_list, "plots/figure2_v2.pdf", 12, 12)

