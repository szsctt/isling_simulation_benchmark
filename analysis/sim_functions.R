library(glue)
library(tidyverse)



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
  
  sim_cols <- cols(experiment = col_character(), 
                condition = col_character(), 
                replicate = col_character(), 
                host_name = col_character(), 
                virus_name = col_character(),
                int_num = col_integer(), 
                epi_num = col_integer(), 
                min_sep = col_integer(), 
                min_len = col_integer(), 
                p_whole = col_double(), 
                p_rearrange = col_double(), 
                p_delete = col_double(), 
                lambda_split = col_double(), 
                p_overlap = col_double(), 
                p_gap = col_double(),
                lambda_junction = col_double(), 
                p_host_deletion = col_double(), 
                lambda_host_deletion = col_double(),
                read_len = col_integer(), 
                fcov = col_double(), 
                frag_len = col_double(), 
                frag_std = col_double(),
                seq_sys = col_character(), 
                out_directory = col_character(), 
                host_fasta = col_character(), 
                virus_fasta = col_character(),
                random_seed = col_integer(), 
                sim_fa_filename = col_character(), 
                sim_int_info_filename = col_character(),
                sim_epi_info_filename = col_character(), 
                read_sam_filename = col_character(),
                sorted_bam_filename = col_character(), 
                annotated_info_filename = col_character(),
                sample = col_character(), 
                unique = col_character())
  
  sim <- read_tsv(sim_sum_path, col_types = sim_cols) 
    
  sim <- sim %>% 
    rename(sim_host = host_name) %>% 
    rename(sim_virus = virus_name) %>% 
    rename(sim_unique = unique)
  
  # import analysis table and add experiment and analysis_condition columns
  analy_cols <- cols(
    experiment= col_character(),   
    host = col_character(),    
    host_fasta = col_character(),      
    virus = col_character(),   
    virus_fasta = col_character(),     
    analysis_condition = col_character(),      
    merge = col_double(),  
    trim = col_logical(), 
    dedup = col_double(),   
    postargs = col_character(),
    adapter_1 = col_character(),        
    adapter_2 = col_character(),       
    bwa_mem_params = col_character(),  
    clip_cutoff = col_integer(),  
    cigar_tol = col_integer(),       
    min_mapq = col_integer(),        
    tool = col_character(),    
    merge_dist = col_integer(),    
    merge_n_min = col_integer(),      
    score_ints_window = col_integer(),       
    score_ints = col_logical(),     
    score_merged_ints = col_logical(),        
    score_reads = col_logical(),      
    exp = col_character(),     
    read_folder = col_character(),    
    R1_suffix = col_character(),       
    R2_suffix = col_character(),       
    outdir = col_character(),  
    bam_suffix = col_character(),      
    host_mappability = col_character(),        
    host_mappability_exclude = col_character(),        
    host_genes = col_character(),     
    host_exons = col_character(),      
    host_oncogenes = col_character(),  
    host_centromeres = col_character(),       
    host_conserved_regions = col_character(), 
    host_segdup = col_character(),     
    detection_mode = col_character(), 
    flank_region_size = col_integer(),     
    sensitivity_level = col_integer(),     
    min_contig_length = col_integer(),      
    blastn_evalue_thrd = col_double(),      
    similarity_thrd = col_double(),
    chop_read_length = col_double(),        
    minIdentity = col_double()
  )
  analysis <- read_tsv(analy_sum_path, col_types = analy_cols) %>% 
    rename(analysis_host = host,
           analysis_virus = virus)
  
  # join simulation and analysis tables
  joined <- left_join(analysis, sim, by=c('experiment'))
  
  # import scored reads summary
  scored_cols <-cols(
    sim_info_file = col_character(),
    sim_sam_file = col_character(),
    analysis_info_file = col_character(),
    results_file = col_character(),
    junc_type = col_character(),
    score_type = col_character(),
    true_positives = col_double(),
    true_negatives = col_double(),
    false_positives = col_double(),
    false_negatives = col_double()
  )
  
  # import and extract info for joining
  scored <- read_tsv(results_sum_path, col_types = scored_cols)  
    
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
    mutate(TPR = true_positives / (true_positives + false_negatives)) %>% 
    mutate(TNR = true_negatives / (true_negatives + false_positives)) %>% 
    mutate(PPV = true_positives / (true_positives + false_positives)) %>% 
    mutate(accuracy = (true_positives + true_negatives) / 
             (true_positives + true_negatives + false_positives + false_negatives)) %>% 
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

#test <- fn_chimeric_reads_explained_by_short_integrations("../out/experiment0_short-refs/test-harder/scored_reads/test-harder_analysis0.cond1.rep3.human.AAV.tsv")

import_fp_reads <- function(read_scores_path, score_type = 'host_score', junction_type = 'chimeric') {
  # given a path to a scored reads file, import the false positives 
  # of a specified score type and junction type
  
  if (!file.exists(read_scores_path)) {
    stop(glue("file {read_scores} does not exist"))
  }
  
  read_scores <- readr::read_tsv(read_scores_path) 
  
  if (!(score_type %in% colnames(read_scores))) {
    stop(glue("coudln't find column {score_type}"))
  }
  if (!(junction_type %in% c('chimeric', 'discord'))) {
    stop(glue("junction type {junction_type} is invalid (must be either 'chimeric' or 'discord')"))
  }
  
  fp_reads <- read_scores %>% 
    filter(across(matches(score_type), ~ . == "fp")) %>% 
    dplyr::filter(type == junction_type) %>% 
    distinct() 
  
  return(fp_reads)
  
}

fp_reads_explained_by_wrong_location <- function(read_scores_path, score_type = 'host_score', junction_type = 'chimeric') {
  # import the read scores and look for false negatives (found_score) that can be explained 
  # by short (undetectable) integrations
  
  # these integrations are indicated by reads that cross both the left and right junctions of the integration
  # or reads that cross more than one integration

  fp_reads <- import_fp_reads(read_scores_path, score_type, junction_type)
  
  if (nrow(fp_reads) == 0) {
    return(NA)
  }
  
  return(sum(fp_reads$found) / nrow(fp_reads))
  
}
#test_file <- "../out/experiment0_short-refs/test-easier/scored_reads/test-easier_analysis0.cond1.rep1.human.AAV.uniq.tsv"
#fp_reads_explained_by_wrong_location(test_file, 'virus_score', 'discord')

#import_fp_reads(test_file, 'virus_score', 'discord')


importIntScoreExperiment <- function(exp_path) {
  
  # import scored integrations
  int_scores <- importIntScoresFromSummaries(exp_path)
  
  # import analysis conditions
  analysis_conditions <- importAnalysisConditions(exp_path)
  
  # add analysis condition info to int scores
  int_scores <- int_scores %>% left_join(analysis_conditions, by="unique") %>% 
    mutate(merge = case_when(
      merge == 1 ~ "merged",
      merge == 0 ~ "unmerged",
      TRUE ~ "unmerged"
    ))
  
  # import simulation conditions
  sim_conditions <- importSimulationConditions(exp_path)

  return(left_join(int_scores, sim_conditions, by=c("experiment", "condition", "replicate")))
  
}

importSimulationConditions <- function(exp_path) {
  
  sim_conditions_file <- "simulation_summary.tsv"
  
  # get files with simulation condtions
  cond_files <- list.files(exp_path, pattern = sim_conditions_file, recursive = TRUE)
  
  if (length(cond_files) == 0) {
    stop(glue("didn't find any files in {exp_path}"))
  }
  
  # import simulation conditions
  sim_conditions <- tibble(
    file = file.path(exp_path, cond_files),
    conds = map(file, ~read_tsv(.))
  ) %>% 
    unnest(conds) %>% 
    select(-file) %>% 
    distinct()   
  
  return(sim_conditions)
}

importAnalysisConditions <- function(exp_path) {
  analysis_condtions_file <- "analysis_conditions.tsv" 
  
  # get files with analysis conditions
  cond_files <- list.files(exp_path, pattern = analysis_condtions_file, recursive = TRUE)
  cond_files <- cond_files[!str_detect(cond_files, "pipeline_analysis_conditions.tsv")]
  
   # import analysis conditions
  analysis_conditions <- tibble(
    file = file.path(exp_path, cond_files),
    conds = map(file, ~read_tsv(.))
  ) %>% 
    unnest(conds) %>% 
    select(-file) %>% 
    distinct() 
  
  return(analysis_conditions)
  }

importIntScoresFromSummaries <- function(exp_path) {
  
  # list files
  scored_ints_folder <- "scored_ints"
  summary_suffix <- "_summary.tsv"
  
  folders <- list.dirs(exp_path)
  folders <- folders[str_detect(folders, scored_ints_folder)]
  
  if (length(folders) == 0) {
    stop(glue("no results found for experiment path {exp_path}"))
  }
  
  # get files for each folder
  filenames <- c()
  for (dir in folders) {
    files <- list.files(dir, pattern = summary_suffix)
    files <- file.path(dir, files)
    filenames <- c(filenames, files)
  }
  
  # import each file
  int_scores <- tibble(
    filename = filenames,
    scores = map(filenames, ~read_tsv(.))
  ) %>% 
    unnest(scores)
  
  # add extra info to scored integrations
  int_scores <- int_scores %>% 
    mutate(results_file = basename(filename)) %>% 
    mutate(unique = str_split(results_file, "\\.", simplify=TRUE)[,1]) %>% 
    mutate(condition = str_split(results_file, "\\.", simplify=TRUE)[,2]) %>%   
    mutate(replicate = str_split(results_file, "\\.", simplify=TRUE)[,3]) %>% 
    mutate(replicate = as.double(str_extract(replicate, "\\d+"))) %>% 
    mutate(analysis_host = str_split(results_file, "\\.", simplify=TRUE)[,4]) %>% 
    mutate(analysis_virus = str_split(results_file, "\\.", simplify=TRUE)[,5]) %>% 
    mutate(post = str_detect(results_file, "post")) %>% 
    rowwise() %>% 
    mutate(TPR = tp/(tp+fn)) %>% 
    mutate(PPV = tp/(tp+fp)) %>% 
    ungroup()
  
  return(int_scores)
}


importJaccardExperiment <- function(exp_path) {
  
  # list files
  jaccard_folder <- "jaccard"
  
  folders <- list.dirs(exp_path)
  folders <- folders[str_detect(folders, jaccard_folder)]

  # get files for each folder
  filenames <- c()
  for (dir in folders) {
    files <- list.files(dir)
    files <- file.path(dir, files)
    filenames <- c(filenames, files)
  }
  
  column_names = c("result_file", "sim_file", "intersection", "union", "jaccard", "n_intersections")
  
  column_types <- cols(
    result_file = col_character(),
    sim_file = col_character(),
    intersection = col_integer(),
    union = col_integer(),
    jaccard = col_double(),
    n_intersections = col_integer()
  )
  
  # import each file
  jaccard <- tibble(
    filename = filenames,
    scores = map(filenames, ~read_tsv(., col_types = column_types, col_names = column_names))
  ) %>% 
    unnest(scores)
  
  jaccard <- jaccard %>% 
    mutate(file_name = basename(filename)) %>% 
    mutate(analysis_condition = str_split(file_name, "\\.", simplify=TRUE)[,1]) %>% 
    mutate(analysis_condition_short = str_match(analysis_condition, "(analysis|seeksv|polyidus|vifi|vseq-toolkit)\\d+")[,1] ) %>% 
    mutate(r = paste0("(.+)_", analysis_condition_short)) %>% 
    mutate(experiment = str_match(analysis_condition, r)[,2]) %>% 
    select(-r) %>% 
    mutate(analysis_tool = str_match(analysis_condition_short, "(.+)\\d+")[,2]) %>% 
    mutate(condition = str_split(file_name, "\\.", simplify=TRUE)[,2]) %>% 
    mutate(replicate = str_split(file_name, "\\.", simplify=TRUE)[,3]) %>%  
    mutate(replicate = as.double(str_extract(replicate, "\\d+"))) %>%
    mutate(analysis_host = str_split(file_name, "\\.", simplify=TRUE)[,4]) %>%     
    mutate(analysis_virus = str_split(file_name, "\\.", simplify=TRUE)[,5]) %>%   
    mutate(post = str_split(file_name, "\\.", simplify=TRUE)[,6] == "post") %>% 
    select(-file_name)
    
  # import simulation conditions
  sim_conditions <- importSimulationConditions(exp_path)
  
  return(left_join(jaccard, sim_conditions, by=c("experiment", "condition", "replicate")))  
    
}

importDistScoreExperiment <- function(exp_path, type) {
  if (type == "found") {
    summary_suffix <- ".found-results.tsv"
    score_type_sub <- '.found-results'
  }
  else {
    summary_suffix <- ".sim-results.tsv"
    score_type_sub <- '.sim-results'
  }
  scored_ints_folder <- "bedtools_closest"
  
  # get folders to read files from
  folders <- list.dirs(exp_path)
  folders <- folders[str_detect(folders, scored_ints_folder)]
  
  # get files for each folder
  filenames <- c()
  for (dir in folders) {
    files <- list.files(dir, pattern = summary_suffix)
    files <- file.path(dir, files)
    filenames <- c(filenames, files)
  }
  
  # make tibble
  int_scores <- tibble(filename = filenames) 
  
  # add extra info to scored integrations
  int_scores <- int_scores %>% 
    rowwise() %>% 
    mutate(results_file = basename(filename)) %>% 
    mutate(unique = str_split(results_file, "\\.", simplify=TRUE)[,1]) %>% 
    mutate(analysis_condition = str_match(unique(unique), "(analysis|seeksv|polyidus|vifi|vseq-toolkit)\\d+")[,1] ) %>% 
    mutate(r = paste0("(.+)_", analysis_condition)) %>% 
    mutate(experiment = str_match(unique, r)[,2]) %>% 
    select(-r) %>% 
    mutate(condition = str_split(results_file, "\\.", simplify=TRUE)[,2]) %>%   
    mutate(replicate = str_split(results_file, "\\.", simplify=TRUE)[,3]) %>% 
    mutate(replicate = as.double(str_extract(replicate, "\\d+"))) %>% 
    mutate(analysis_host = str_split(results_file, "\\.", simplify=TRUE)[,4]) %>% 
    mutate(analysis_virus = str_split(results_file, "\\.", simplify=TRUE)[,5]) %>% 
    mutate(post = str_detect(results_file, "post")) %>% 
    ungroup()
  
  # import data
  int_scores <- int_scores %>% 
    mutate(data =  map(filename, ~importDistFile(.))) %>% 
    unnest(data)
  
  return(int_scores)
  
}


importDistFile <- function(filename) {
   colspec <- cols(
     chr_1 = col_character(),
     start_1 = col_integer(),
     stop_1 = col_integer(),
     chr_2 = col_character(),
     start_2 = col_integer(),
     stop_2 = col_integer(),
     d_shortest = col_integer()
   )
  #   coords_mean = col_double(),
  #   coords_min = col_double(),
  #   midpoint = col_double()
  # )
  cat("importing file ", filename, "\n")
  
  return(read_tsv(filename, col_names = names(colspec$cols), col_types = colspec))
  #return(read_tsv(filename))
  
}


importNearestSimToFound <- function(exp_path) {
  return(importDistScoreExperiment(exp_path, "found"))
}

importNearestFoundToSim <- function(exp_path) {
  
  return(importDistScoreExperiment(exp_path, "sim"))
}

#### plotting functions ####

# get the variables that are different in an experiment
getExpVars <- function(conds, exp_name) {
  select_cols <- c("host_name", "virus_name", "int_num", "epi_num", "min_sep", "min_len",
              "p_whole", "p_rearrange", "p_delete", "lambda_split", "p_overlap", "p_gap",
              "lambda_junction", "p_host_deletion", "lambda_host_deletion", "read_len", "fcov", 
              "frag_len", "frag_std", "seq_sys")

  cols <- conds %>% 
    filter(experiment == exp_name) 
  
  if (nrow(conds) == 0) {
    stop(glue("no rows found for experiment {exp_name}"))
  }
  
  cols <- cols %>% 
    select(one_of(select_cols)) %>% 
    summarise(across(everything(), n_distinct),) %>% 
    select(where(~sum(.) > 1)) %>% 
    select(one_of(select_cols)) %>% 
    colnames()
  
  cols <- cols[!str_detect(cols, "fasta")]
  
  return(cols)
}

#function to combine variables manipulated in an experiment to give a label for each condition
combineVarNames <- function(df, exp_name) {
  exp_vars <- getExpVars(df, exp_name)
  
  if (length(exp_vars) == 1) {
    return(df[[exp_vars[1]]])
  }
  
  names_df <- df %>% 
    filter(experiment == exp_name) %>% 
    select(all_of(exp_vars))
  
  if (nrow(names_df) == 0) {
    stop(glue("no rows found for experiment {exp_name}"))
  } 
  
  # include variable names in each row of names_df
  for (var in colnames(names_df)) {
    names_df <- names_df %>% 
      mutate(!!var := paste0(var, ": ", !!sym(var))) 
  }
  
  names_df <- names_df %>% 
    rowwise() %>% 
    mutate(label = paste0(across(everything()), collapse=", "))
  
  return(names_df$label)
}

shortCombinedVarNames <- function(df, exp_name) {
  exp_vars <- getExpVars(df, exp_name)
  
  if (length(exp_vars) == 1) {
    return(df[[exp_vars[1]]])
  }
  
  names_df <- df %>% 
    filter(experiment == exp_name) %>% 
    select(all_of(exp_vars))
  
  names_df <- names_df %>% 
    rowwise() %>% 
    mutate(short_label = paste0(across(everything()), collapse = "/"))
  
  
  return(names_df$short_label)
}

get_readids_column <- function(column) {
  column <- column[!is.na(column)]
  return(str_split(column, ";") %>% flatten_chr())
}

# compare readIDs of output integrations to simulated ones
# include all integraitons, not just uniquely localised ones
compareIslingOutputSim <- function(output_file, sim_file) {
  
  out <- read_tsv(output_file)
  sim <- read_tsv(sim_file)
  
  sim_chimeric <- c(
    get_readids_column(sim$left_chimeric),
    get_readids_column(sim$right_chimeric)
  )
  
  sim_chimeric
  
  sim_discordant <- c(
    get_readids_column(sim$left_discord),
    get_readids_column(sim$right_discord)   
  )
  

  out <- out %>% 
    rowwise() %>% 
    mutate(out_in_sim = case_when(
      Type == "chimeric" ~ ReadID %in% sim_chimeric,
      Type == "short" ~ ReadID %in% sim_chimeric,
      Type == "discordant" ~ ReadID %in% sim_discordant
    )) %>% 
    mutate(sim_in_out = case_when (
      Type == "chimeric" ~ sim_chimeric %in% ReadID,
      Type == "short" ~ sim_chimeric %in% ReadID,
      Type == "discordant" ~ sim_discordant %in% ReadID      
    ))
  
  return(out)
  
}

compareAllIslingOutputSim <- function(exp_dir) {
  
  results <- tibble(
    out = list.files(exp_dir, pattern="integrations.post.txt", full.names = T, recursive=T),
    experiment = basename(dirname(dirname(dirname(out)))),
    analysis_condition = basename(dirname(dirname(out))),
    sim_experiment = str_split(analysis_condition, "_", simplify=T)[,1],
    condition = str_extract(basename(out), "cond\\d+"),
    replicate = as.double(str_extract(str_extract(basename(out), "rep\\d+"), "\\d")),
    sim = file.path(exp_dir, sim_experiment, "sim_ints", glue::glue("{condition}.rep{replicate}.int-info.annotated.tsv")),
    data = map2(out, sim, ~compareIslingOutputSim(.x, .y)),
    tp = map_dbl(data, ~sum(.$out_in_sim)),
    fp = map_dbl(data, ~sum(!.$out_in_sim)),
    PPV = tp / (tp + fp)
  )
  
  
  sim_conds <- importSimulationConditions(exp_dir)
  
  return(left_join(results, sim_conds,by=c("condition", "replicate", "sim_experiment"="experiment")))
  
}
