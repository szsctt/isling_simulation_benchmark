#### figure 3 ####
# runtime stuff
library(tidyverse)
library(cowplot)
library(lubridate)

source("/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/analysis/sim_functions.R")

plot_tool_order <- c("isling", "Seeksv", "Polyidus", "ViFi", 'VSeq-Toolkit')

experiments <- c("coverage", 'viral_load')

results_dir <- "/datasets/work/hb-viralint/work/simulations/intvi_simulation-experiments/out/experiment2_time/"
files <- list.files(results_dir, recursive = TRUE, pattern="_runtime.tsv")

files <- files[ basename(dirname(files)) %in% experiments]

files <- files[dirname(dirname(files)) == "."]


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
  mutate(duration_elapsed_time = case_when(
    str_detect(elapsed_time, "\\d{0,2}:\\d{2}\\.\\d{2}") ~ lubridate::as.duration(lubridate::ms(elapsed_time)),
    str_detect(elapsed_time, "\\d{0,2}:\\d{2}\\:\\d{2}") ~ lubridate::as.duration(lubridate::hms(elapsed_time)),
    TRUE ~ lubridate::as.duration(NA)
  )) %>% 
  rename(time_rep = replicate)

# how many times does each exit value occur?
table(times$exit_value)
times %>% 
  group_by(tool, exit_value) %>% 
  summarise(n = n())

# check number of replicates for each pair of fastq files
times %>% 
  group_by(experiment, sample, tool) %>% 
  summarise(n = n())

times %>% 
  group_by(experiment, sample, tool) %>% 
  summarise(n = n()) %>% 
  filter(n != 3)

times %>% 
  group_by(experiment, sample) %>% 
  summarise(n = n()) %>% 
  filter(n != 15)

# add simulation conditions
conds <- importSimulationConditions(results_dir)

times <- left_join(times, conds, by=c("experiment", "sample")) %>% 
  mutate(tool = str_replace(tool, "polyidus", "Polyidus")) %>% 
  mutate(tool = str_replace(tool, "seeksv", "Seeksv")) %>% 
  mutate(tool = str_replace(tool, "vifi", "ViFi")) %>% 
  mutate(tool = str_replace(tool, "vseq-toolkit", "VSeq-Toolkit"))  


times %>% 
  mutate(tool = as.factor(tool)) %>% 
  ggplot(aes(x = condition, y = duration_elapsed_time, color=tool)) +
  geom_boxplot() +  
  facet_grid(rows = vars(experiment), scales = 'free') +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name='tool')
ggsave("plots/runtime_elapsed-time.pdf")

# memory usage - max_rss
times %>% 
  mutate(tool = as.factor(tool)) %>% 
  ggplot(aes(x = condition, y = swaps, color=tool)) +
  geom_boxplot() +  
  facet_grid(rows = vars(experiment), scales = 'free') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name='tool')
ggsave("plots/runtime_max-rss.pdf")

# swaps
times %>% 
  mutate(tool = as.factor(tool)) %>% 
  ggplot(aes(x = condition, y = max_rss, color=tool)) +
  geom_boxplot() +  
  facet_grid(rows = vars(experiment), scales = 'free') +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name='tool')
ggsave("plots/runtime_swaps.pdf")

# filesystem output
times %>% 
  mutate(tool = as.factor(tool)) %>% 
  ggplot(aes(x = condition, y = fs_outputs, color=tool)) +
  geom_boxplot() +  
  facet_grid(rows = vars(experiment), scales = 'free') +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name='tool')
ggsave("plots/runtime_fs-output.pdf")

# CPU
times %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(CPU = str_replace(CPU, "%", "")) %>% 
  mutate(CPU = as.integer(CPU)) %>% 
  ggplot(aes(x = condition, y = CPU, color=tool)) +
  geom_boxplot() +  
  facet_grid(rows = vars(experiment), scales = 'free') +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name='tool')
ggsave("plots/runtime_CPU.pdf")

plots <- list()
for (e in unique(times$experiment)) {
  filt_df <- times %>% 
    filter(experiment == e)
  
  
  plots[[e]] <- filt_df %>% 
    mutate(label = combineVarNames(filt_df, e)) %>% 
    ggplot(aes(x = label, y = duration_elapsed_time, color=tool)) +
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
  ggplot(aes(x = condition, y = duration_elapsed_time, color=tool)) +
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
  ggplot(aes(x = fcov, y = duration_elapsed_time, color=tool)) +
  geom_point()  +
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
  ggplot(aes(x = int_num, y = duration_elapsed_time, color=tool)) +
  geom_point()  +
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


viral_load <- times %>% 
  filter(experiment == "viral_load") %>% 
  mutate(int_num = as_factor(int_num)) %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, plot_tool_order)) %>% 
   ggplot(aes(x = int_num, y = duration_elapsed_time, color=tool)) +
  geom_boxplot()  +
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
  scale_y_log10()


print(viral_load + theme(legend.position = "bottom"))

coverage <- times %>% 
  filter(experiment == "coverage") %>% 
  mutate(fcov = as.factor(fcov)) %>% 
  mutate(tool = as.factor(tool)) %>% 
  mutate(tool = forcats::fct_relevel(tool, plot_tool_order)) %>% 
  ggplot(aes(x = fcov, y = duration_elapsed_time, color=tool)) +
  geom_boxplot()  +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  xlab("Fold coverage") +
  ylab("Runtime (s)") +
  scale_y_log10() 

print(coverage + theme(legend.position = "bottom"))

legend <- get_legend(coverage + theme(legend.position = "bottom"))

figure_3 <- cowplot::plot_grid(coverage, viral_load, labels="AUTO")
figure_3 <- cowplot::plot_grid(figure_3, legend, ncol=1, rel_heights = c(1, 0.05))
print(figure_3)

cowplot::save_plot("plots/figure3_v3.pdf", figure_3)
