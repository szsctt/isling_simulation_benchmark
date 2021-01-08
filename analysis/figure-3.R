#### figure 3 ####
# runtime stuff
library(tidyverse)
library(cowplot)
library(lubridate)

source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/sim_functions.R")


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
    str_detect(elapsed_time, "\\d{0,2}:\\d{2}\\.\\d{2}") ~ lubridate::as.duration(lubridate::ms(elapsed_time)),
    str_detect(elapsed_time, "\\d{0,2}:\\d{2}\\:\\d{2}") ~ lubridate::as.duration(lubridate::hms(elapsed_time)),
    TRUE ~ lubridate::as.duration(NA)
  )) %>% 
  rename(time_rep = replicate)

conds <- importSimulationConditions(results_dir)

times <- left_join(times, conds, by=c("experiment", "sample"))

times %>% 
  ggplot(aes(x = condition, y = elapsed_time, color=tool)) +
  geom_point() +
  facet_grid(rows = vars(experiment), scales = 'free_y') +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


plots <- list()
for (e in unique(times$experiment)) {
  filt_df <- times %>% 
    filter(experiment == e)
  
  
  plots[[e]] <- filt_df %>% 
    mutate(label = combineVarNames(filt_df, conds, e)) %>% 
    ggplot(aes(x = label, y = elapsed_time, color=tool)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )  
}
plots['coverage']
plots['viral_load']


p <- times %>% 
  ggplot(aes(x = condition, y = elapsed_time, color=tool)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "bottom")

legend <- get_legend(p)

figure_3 <- cowplot::plot_grid(plotlist = plots, labels="AUTO")
figure_3 <- cowplot::plot_grid(figure_3, legend, ncol=1, rel_heights = c(1, 0.05))
print(figure_3)

cowplot::save_plot("plots/figure3_v1.pdf", figure_3)



coverage <- times %>% 
  filter(experiment == "coverage") %>% 
  ggplot(aes(x = fcov, y = elapsed_time, color=tool)) +
  geom_point(alpha = 0.5)  +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  xlab("Fold coverage") +
  ylab("Runtime (s)") +
  scale_y_log10() +
  scale_x_log10()

print(coverage)

viral_load <- times %>% 
  filter(experiment == "viral_load") %>% 
  mutate(int_num = int_num + 1) %>% 
  ggplot(aes(x = int_num, y = elapsed_time, color=tool)) +
  geom_point(alpha = 0.5)  +
  facet_wrap(vars(epi_num)) +
  theme_classic()  + 
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
#    axis.title.y = element_blank(),
  )  +
  xlab("Number of integrations") +
  ylab("Runtime (s)") +
  scale_y_log10() +
  scale_x_log10(breaks=c(1, 100, 10000))
  

print(viral_load)


figure_3 <- cowplot::plot_grid(coverage, viral_load, labels="AUTO")
figure_3 <- cowplot::plot_grid(figure_3, legend, ncol=1, rel_heights = c(1, 0.05))
print(figure_3)

cowplot::save_plot("plots/figure3_v2.pdf", figure_3)
