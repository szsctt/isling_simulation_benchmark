
importReadScoreExperiment <- function(experiment_path) {
  
  # expected filenames of anaysis, simulated and read score files
  analysis_filename <- "^analysis_conditions\\.tsv$"
  sim_filename <- "^simulation_summary\\.tsv$"
  results_filename <- "scored_reads_summary\\.tsv$"
  
  # get directories in the specified directories to look for
  search_dirs <- list.dirs(experiment_path, recursive = FALSE)
  
  # check that we found at least one directory
  if (length(search_dirs) == 0) {
    stop("no files found in specified directory")
  }
  
  # to keep track of which directories we imported data from
  imported_dirs <- c()
  
  read_scores <- tibble()
  
  for (dir in search_dirs) {
    
    # look for a file with simulation info
    sim_file <- list.files(dir, pattern = sim_filename)
    
    # look for a file with analysis pipeline paramters
    analysis_file <- list.files(dir, pattern = analysis_filename)
    
    # look for a summary of read scores
    results_file <- list.files(dir, pattern = results_filename)
    
    # check we found exactly one file for each
    if (length(sim_file) != 1) {
      next
    }
    else if (length(analysis_file) != 1) {
      next
    }
    else if (length(results_file) != 1) {
      next
    }
    
    # import data
    tmp <- importReadData(
      paste0(dir, "/", sim_file), 
      paste0(dir, "/", analysis_file), 
      paste0(dir, "/", results_file)
      )
    
    read_scores <- bind_rows(read_scores, tmp)
    
  }
  
  return(read_scores)
  
}


#importReadScoreExperiment("../out/experiment0_short-refs")

importReadData <- function(sim_sum_path, analy_sum_path, results_sum_path) {
  # function to import read score summaries from an experiment
  
  sim_cols <- c('experiment', 'condition', 'replicate', 'sim_host', 'sim_virus',
                'int_num', 'epi_num', 'min_sep', 'p_whole', 'p_rearrange', 
                'p_delete', 'lambda_split', 'p_overlap', 'p_gap',
                'lambda_junction', 'p_host_deletion', 'lambda_host_deletion',
                'read_len', 'fcov', 'frag_len', 'frag_std',
                'seq_sys', 'out_directory', 'host_fasta', 'virus_fasta',
                'random_seed', 'sim_fa_filename', 'sim_int_info_filename',
                'sim_epi_info_filename', 'read_sam_filename',
                'sorted_bam_filename', 'annotated_info_filename',
                'sample', 'unique_sim')
  
  sim <- read_tsv(sim_sum_path, col_names = sim_cols, skip=1)
  
  analy_cols <- c('experiment', 'experimet_dup', 'analysis_condition', 'tool', 'analysis_host',
                  'host_fasta', 'analysis_virus', 'virus_fasta', 'bam_suffix', 'read_folder',
                  'R1_suffix', 'R2_suffix', 'outdir', 'merge', 'dedup', 'postargs',
                  'seq_sys', 'adapter_1', 'adapter_2', 'bwa_mem_params', 'score_ints',
                   'score_ints_window', 'score_ints_tool')
  
  # import analysis table and add experiment and analysis_condition columns
  analysis <- read_tsv(analy_sum_path, col_names = analy_cols, skip=1)
    
  # join simulation and analysis tables
  joined <- left_join(analysis, sim, by=c('experiment'))
  
  # import scored reads summary
  scored_cols <- c('sim_info_file', 'sim_bam_file', 'analysis_info',
                   'results_file', 'junc_type', 'score_type',
                   'true_positive', 'true_negative',
                   'false_positive', 'false_negative')
  
  # import and extract info for joining
  scored <- read_tsv(results_sum_path, col_names = scored_cols, skip = 1)  
    
  scored <- scored %>% 
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
                          'config_dataset'='analysis_condition', 'experiment')))
  
}

# import reads from scored file with column 'col' (numbered starting from 1) with entry matching regex 'str'
get_reads_matching_condition <- function(col, str, scored_reads_file) {
  
  # filter scored_reads_file into temporary file
  tmp <- tempfile(tmpdir = ".")
  command <- paste0("awk 'match($", col, ", /", str, "/)' ", scored_reads_file, " | uniq > ", tmp)
  system(command, intern=TRUE)
  
  # read temporary file using read_tsv
  colnames <- c("readID", "intID", "found_score", "host_score", "virus_score", "found", "n_found", "side", "correct_side",
                "type",  "corect_type", "correct_host_chr", "host_start_dist", "host_stop_dist", "host_coords_overlap",
                "correct_virus", "virus_coords_overlap", "virus_start_dist", "virus_stop_dist", "ambig_diff")
  
  scored <- read_tsv(tmp, col_names = colnames)
  file.remove(tmp)
  
  return(scored)
}


# import integrations from scored file with column 'col_num' (numbered starting from 1) with entry matching regex 'str'
# ie for integrations with reads that are have 'found_score' (column 3) as false negatives, 
# use get_ints_matching_contidion("3", "fn", ...), and number of reads for each integration matching this conditiona
annotate_ints_matching_condition <- function(col_num, regex_str, scored_reads_file, annotated_int_file) {
  
  # get integration IDs and types in temp file
  tmp <- tempfile(tmpdir = ".")
  command <- paste0("awk 'match($", col_num, ", /", regex_str, "/)' ", scored_reads_file, " | sort | uniq | cut -f2,10 | sort | uniq -c > ", tmp)
  system(command, intern=TRUE)

  # read temp file
  int_ids <- read_table2(tmp, col_names = c("n", "int_id", "type")) %>% 
    pivot_wider(names_from = 'type', values_from = 'n', names_prefix=paste0("n_", regex_str, "_"))
  file.remove(tmp)
  
  # import integrations
  ints <- read_tsv(annotated_int_file)
  
  # join info from temp file with ints
  ints <- left_join(ints, int_ids, by = c("id" = "int_id")) %>% 
    ungroup()
  
  # replace NA counts with 0
  # https://github.com/sparklyr/sparklyr/issues/1062

  for (read_type in c("discord", "chimeric")) {
    colname <- paste0("n_", regex_str, "_", read_type)
    ints <- ints %>% 
      mutate(!!colname := replace_na(.[[colname]], 0))
  }
  
  # also add counts of how many left and right discordant and chimeric reads were annotated
  for (read_type in c("discord", "chimeric")) {
    for (side in c("left", "right")) {
      new_colname <- paste0("annotated_", read_type, "_", side)
      old_colname <- paste0(side, "_", read_type)
      ints <- ints %>% 
        mutate(!!new_colname := lengths(str_split(.[[old_colname]], ";"))) %>% 
        mutate(!!new_colname := ifelse(is.na(.[[old_colname]]), 0, .[[new_colname]]))
    }
  }

  # add total numbers for discordant and chimeric reads for each integration
  ints <- ints %>% 
    rowwise() %>% 
    mutate(annotated_chimeric = annotated_chimeric_left + annotated_chimeric_right) %>% 
    mutate(annotated_discord = annotated_discord_left + annotated_discord_right) %>% 
    ungroup()
  
  return(ints)
}


#scored_file <- "../out/experiment0_prelim/rep68-easier/scored_reads/rep68-easier_analysis0.cond0.rep0.chr1.rep68.tsv"
#int_file <- "../out/experiment0_prelim/rep68-easier/sim_ints/cond0.rep0.int-info.annotated.tsv"
#test <- annotate_ints_matching_condition("3", "fn", scored_file, int_file)

# bin integrations annotated with number of missing reads, etc, using `annotate_ints_matching_condition` function above.
# note that this assumes that there's only one chromosome!
add_bins_to_annotated_ints <- function(ints, width, chr_len) {
  chrs <- unique(ints$chr)
  if (length(chrs) > 1) {
    stop("found more than one chromosome in dataset")
  }
  n_bins <- ceiling(chr_len/width)
  ints <- ints %>% 
    mutate(bin = cut(hPos, breaks=seq(0, width*n_bins, width), labels= seq(0, n_bins-1))) %>% 
    mutate(bin = as.numeric(as.character(bin)))
  
  return(ints)
}

#scored_file <- "../out/experiment0_prelim/rep68-easier/scored_reads/rep68-easier_analysis0.cond0.rep0.chr1.rep68.tsv"
#int_file <- "../out/experiment0_prelim/rep68-easier/sim_ints/cond0.rep0.int-info.annotated.tsv"
#test <- annotate_ints_matching_condition("3", "fn", scored_file, int_file)
#add_bins_to_annotated_ints(test, 10000, 248956422)

missing_read_info <- function(scored_reads_glob, int_file_string) {
  scored_files <- Sys.glob(scored_reads_glob)
  scored_files <- scored_files[str_detect(scored_files, "summary", negate=TRUE)]
  
  return(scored_files)
  
  cat('processing', length(scored_files), "files\n")
  
  df <- tibble(
    scored_reads = scored_files,
    scored_reads_basename = basename(scored_files),
    experiment = basename(dirname(dirname(scored_reads))),
    analysis_condition = str_split(basename(scored_reads_basename), "\\.", simplify=TRUE)[,1],
    condition = str_split(basename(scored_reads_basename), "\\.", simplify=TRUE)[,2],
    replicate = str_split(basename(scored_reads_basename), "\\.", simplify=TRUE)[,3],
    analysis_host =  str_split(basename(scored_reads_basename), "\\.", simplify=TRUE)[,4],
    analysis_virus =  str_split(basename(scored_reads_basename), "\\.", simplify=TRUE)[,5],
    post = str_detect(scored_reads_basename, "post"),
    sample=paste0(condition, replicate),
    analysis_condition_number = str_extract(analysis_condition, "\\d+")
  )
  df <- df %>% 
    mutate(replicate = as.integer(str_extract(replicate, "\\d+"))) %>% 
    mutate(int_file = str_replace(int_file_string, "\\{experiment\\}", experiment)) %>% 
    mutate(int_file = str_replace(int_file, "\\{condition\\}", condition)) %>% 
    mutate(int_file = str_replace(int_file, "\\{replicate\\}", as.character(replicate)))
  
  return(df)
}



#test <-missing_read_info("../out/experiment0_prelim/*/scored_reads/*.tsv",
#                         "../out/experiment0_prelim/{experiment}/sim_ints/{condition}.rep{replicate}.int-info.annotated.tsv")
#test

fn_chimeric_reads_explained_by_short_integrations <- function(read_scores) {
  # import the read scores and look for false negatives (found_score) that can be explained 
  # by short (undetectable) integrations
  
  # these integrations are indicated by reads that cross both the left and right junctions of the integration
  # or reads that cross more than one integration
  
  if (!file.exists(read_scores)) {
    stop(paste0("file ", read_scores, " does not exist"))
  }
  
  read_scores <- readr::read_tsv(read_scores) 
  
  fn_reads <- read_scores %>% 
    dplyr::filter(found_score == "fn") %>% 
    dplyr::filter(type == "chimeric") 
  
  if (nrow(fn_reads) == 0) {
    return(NA)
  }
  
  fn_reads <- fn_reads %>% 
    distinct() %>% 
    group_by(readID) %>% 
    summarise(number_sides = n_distinct(side),
              number_ints = n_distinct(intID)) %>% 
    rowwise() %>% 
    mutate(multiple_sides_or_ints = number_sides > 1 | number_ints > 1) %>% 
    ungroup() %>% 
    summarise(n_explained = sum(multiple_sides_or_ints),
              total_fn = n())
  
  return(fn_reads$n_explained / fn_reads$total_fn)
  
  
}

test <- fn_chimeric_reads_explained_by_short_integrations("../out/experiment0_short-refs/test-harder/scored_reads/test-harder_analysis0.cond1.rep3.human.AAV.tsv")
