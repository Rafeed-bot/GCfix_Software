#!/bin/Rscript
#
# NDR Range
#
#
# Written by POH Zhong Wee - 2024 Oct 9

# Packages -------
if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if (!require("pacman", quietly = TRUE)){
    install.packages("pacman")
}

required_packages = c("tidyverse", "data.table", "svglite", "argparse")
pacman::p_load(char = required_packages)

# Arguments -------
parser <- ArgumentParser(description = "Example script using argparse in R")

parser$add_argument("-m", "--maindir", type = "character",
                    default = paste0(getwd(), "/temp"),
                    help = "Input file path")

parser$add_argument("-c", "--covlist", type = "character",
                    help = "List of coverage files (one per line)",
                    required = TRUE)

parser$add_argument("-g", "--geneset", type = "character",
                    help = "Geneset file",
                    required = TRUE)

parser$add_argument("-l", "--label", type = "character",
                    default = "run_1",
                    help = "Output label")

args = parser$parse_args()

# Functions -------
check_dir = function(indir){
  
  print(paste0("Directory: ", indir))
  
  if(!dir.exists(indir)){
    
    print(paste0("Directory does not exist. Creating directory..."))
    dir.create(indir, recursive = T, mode = "775")
    
  } else {
    
    print("Directory Exists.")
    
  }
  
}
  
  
# Directories ------
output_dir <- paste0(args$maindir, "/", "output")
check_dir(output_dir)
setwd(output_dir)


# Data ------
coverage_list <- read_lines(args$covlist)
geneset_df <- read_csv(args$geneset)
geneset_order <- geneset_df$group %>% unique


# Run -----
label_dir <- paste0(output_dir, "/", args$label)
check_dir(label_dir)
setwd(label_dir)

total_cov_df <- lapply(coverage_list, function(x){
  
  t_sid <- basename(x) %>%
    str_remove_all("(-ndr|_coverage|.csv|.gz)")
  
  temp_df <- read_csv(x, show_col_types = FALSE) %>%
    mutate(sid = t_sid) %>%
    dplyr::select(sid, everything())
  
  return(temp_df)
  
}) %>% do.call(rbind, .)


select_cols <- colnames(total_cov_df)[str_detect(colnames(total_cov_df), "_mean$")] %>%
  str_remove("_mean$")

base_cols <- colnames(total_cov_df)
for (i in select_cols) {
  base_cols <- base_cols[!str_detect(base_cols, i)]
}


# combine_mean <- function(input_mean, input_size) {weighted.mean(input_mean, input_size)}
combine_sd <- function(input_sd, input_mean, input_size){
  sqrt(weighted.mean(input_sd^2 + input_mean^2, input_size) - weighted.mean(input_mean, input_size)^2)
}

total_cov_sum_df <- lapply(select_cols, function(t_cols){
  
  print(t_cols)

  filter_cols <- paste0(t_cols, c("_mean", "_sd", "_count"))
  
  select_cov <- total_cov_df %>%
    dplyr::select(all_of(c(base_cols, filter_cols)))
  
  colnames(select_cov)[colnames(select_cov) %in% filter_cols] <- c("mean", "sd", "count")
  
  
  select_cov_sum <- select_cov %>%
    group_by(group, rel_pos) %>%
    summarize(sid_mean = mean(mean, na.rm = T),
              sid_sd = sd(mean, na.rm = T),
              total_mean = weighted.mean(mean, count),
              total_sd = combine_sd(sd, mean, count),
              .groups = "drop")
  select_cov_sum$type <- t_cols
  
  select_cov_sum <- select_cov_sum %>%
    dplyr::select(type, everything())
  
  return(select_cov_sum)
  
}) %>% do.call(rbind, .)

fwrite(total_cov_sum_df, file = paste0(args$label, "_ndrcov_summary.csv.gz"), sep = ",", compress = "gzip")

# Plotting -----
plot_dir <- paste0(label_dir, "/plot")
check_dir(plot_dir)

for (t_cols in select_cols) {
  
  print(t_cols)
  
  cov_filter <- total_cov_sum_df %>%
    filter(type == t_cols)
  
  cov_filter <- cov_filter %>%
    mutate(sid_upper = sid_mean + (2*sid_sd),
           sid_lower = sid_mean - (2*sid_sd),
           total_upper = total_mean + (2*total_sd),
           total_lower = total_mean - (2*total_sd))
  
  total_upper_limit <- cov_filter$total_upper %>% quantile(probs = c(0.95))
  total_lower_limit <- cov_filter$total_lower %>% quantile(probs = c(0.05))
  
  cov_filter <- cov_filter %>%
    mutate(total_upper = ifelse(total_upper >= total_upper_limit, total_upper_limit, total_upper),
           total_lower = ifelse(total_lower <= total_lower_limit, total_lower_limit, total_lower))
  
  cov_filter$group <- cov_filter$group %>% factor(levels = geneset_order)
  
  total_plot <- ggplot(cov_filter) +
    
    geom_ribbon(mapping = aes(x = rel_pos, ymin = total_lower, ymax = total_upper,
                              fill = group),
                alpha = 0.1) +
    geom_line(mapping = aes(x = rel_pos, y = total_mean, color = group)) +
    
    labs(x = "Relative Position", y = "Relative Coverage",
         color = "Group", fill = "Group",
         title = t_cols) +
    
    scale_y_continuous(limits = c(total_lower_limit, total_upper_limit)) + 
    
    theme_bw() +
    theme()
  
  outfile <- paste0(args$label, "_", t_cols, "-ndr_coverage-total.svg")
  ggsave(plot = total_plot, filename = outfile, path = plot_dir,
         width = 8, height = 6)
  ggsave(plot = total_plot, filename = str_replace(outfile, ".svg$", ".png"), path = plot_dir,
         width = 8, height = 6)
  
  
  sid_plot <- ggplot(cov_filter) +
    
    geom_ribbon(mapping = aes(x = rel_pos, ymin = sid_lower, ymax = sid_upper,
                              fill = group),
                alpha = 0.1) +
    geom_line(mapping = aes(x = rel_pos, y = sid_mean, color = group)) +
    
    labs(x = "Relative Position", y = "Relative Coverage",
         color = "Group", fill = "Group",
         title = t_cols) +
    
    theme_bw() +
    theme()
  
  sid_outfile <- paste0(args$label, "_", t_cols, "-ndr_coverage-sid.svg")
  ggsave(plot = sid_plot, filename = sid_outfile, path = plot_dir,
         width = 8, height = 6)
  ggsave(plot = sid_plot, filename = str_replace(sid_outfile, ".svg$", ".png"), path = plot_dir,
         width = 8, height = 6)
  
}






















