importData <- function(sim_sum_path, analy_sum_path, results_sum_path) {
  # function to import summaries from an experiment
  
  sim_cols <- c('experiment', 'condition', 'replicate', 'sim_host', 'sim_virus',
                'int_num', 'epi_num', 'p_whole', 'p_rearrange', 
                'p_delete', 'lambda_split', 'p_overlap', 'p_gap',
                'lambda_junction', 'p_host_deletion', 'lambda_host_deletion',
                'read_len', 'fcov', 'frag_len', 'frag_std',
                'seq_sys', 'out_directory', 'host_fasta', 'virus_fasta',
                'random_seed', 'sim_fa_filename', 'sim_int_info_filename',
                'sim_epi_info_filename', 'read_sam_filename',
                'sorted_bam_filename', 'annotated_info_filename',
                'sample', 'unique_sim')
  
  sim <- read_tsv(sim_sum_path, col_names = sim_cols, skip=1)
  
  analy_cols <- c('dataset', 'config_dataset', 'sample', 'analysis_host', 
                  'host_fasta_analysis', 'analysis_virus', 'virus_fasta_analysis', 
                  'merge', 'dedup', 'unique_analysis', 'outdir', 'bwa_mem_params',
                  'R1_file', 'R2_file', 'bam_file', 'adapter_1',
                  'adapter_2', 'postargs')
  
  # import analysis table and add experiment and analysis_condition columns
  analysis <- read_tsv(analy_sum_path, col_names = analy_cols, skip=1) %>% 
    mutate(experiment = str_split(dataset, "_", simplify=TRUE)[,1]) %>% 
    mutate(analysis_condition = str_split(dataset, "_", simplify=TRUE)[,2])
  
    
  # join simulation and analysis tables
  joined <- left_join(analysis, sim, by=c('sample', 'experiment'))
  
  # import scored reads summary
  scored_cols <- c('sim_info_file', 'sim_bam_file', 'analysis_info',
                   'results_file', 'true_positive', 'true_negative',
                   'false_positive', 'false_negative')
  
  # import and extract info for joining
  scored <- read_tsv(results_sum_path, col_names = scored_cols, skip = 1)  %>% 
    mutate(results_info = basename(results_file)) %>% 
    mutate(experiment = str_split(results_info, "_", simplify=TRUE)[,1]) %>% 
    mutate(results_info_2 = str_split(results_info, "_", simplify=TRUE)[,2]) %>% 
    mutate(analysis_condition = str_split(results_info_2, "\\.", simplify = TRUE)[,1]) %>% 
    mutate(sample = paste0(
                            str_split(results_info_2, "\\.", simplify = TRUE)[,2], ".",  
                            str_split(results_info_2, "\\.", simplify = TRUE)[,3]
                            )
    ) %>% 
    mutate(analysis_host = str_split(results_info_2, "\\.", simplify = TRUE)[,4]) %>% 
    mutate(analysis_virus = str_split(results_info_2, "\\.", simplify = TRUE)[,5]) %>% 
    mutate(post = str_detect(results_info_2, "post")) %>% 
    mutate(config_dataset = paste0(experiment, "_", analysis_condition)) %>% 
    arrange(sample, analysis_host, analysis_virus, post) 
  
  scored <- scored %>% 
    mutate(TPR = true_positive / (true_positive + false_negative)) %>% 
    mutate(TNR = true_negative / (true_negative + false_positive)) %>% 
    mutate(PPV = true_positive / (true_positive + false_positive)) %>% 
    mutate(accuracy = (true_positive + true_negative) / 
             (true_positive + true_negative + false_positive + false_negative)) %>% 
    mutate(balanced_accuracy = (TPR + TNR)/2)
    
    #https://en.wikipedia.org/wiki/Precision_and_recall
  # join scored with simulation and analysis tables
  return(left_join(scored, joined, 
                   by = c('sample', 'analysis_host', 'analysis_virus', 
                          'config_dataset', 'analysis_condition', 'experiment')))
  
}