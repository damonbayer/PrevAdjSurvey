library(tidyverse)
library(fs)

tibble(file_path = dir_ls("~/PrevAdjSurvey/code/simulations/confidence_interval_perfect_assay"))


dir_ls("~/PrevAdjSurvey/code/simulations/confidence_interval_perfect_assay") %>%
  enframe(name = NULL) %>%
  filter(str_starts(path_file(value), "slurm")) %>%
  deframe() %>%
  file_delete()


perfect_assay_raw_vec <-
  combined_imperfect_assay_raw <-
  tibble(file_path = dir_ls("//data/bayerdm/confidence_interval_perfect_assay")) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  filter(str_detect(file_name, "raw")) %>%
  arrange(file_name) %>%
  pull(file_path)

combined_perfect_assay_raw <-
  map_dfr(perfect_assay_raw_vec[350:400],
          ~read_rds(.) %>%
            mutate(width = conf_int_2 - conf_int_1) %>%
            group_by(design, method) %>%
            summarize(conf_int_1 = mean(conf_int_1),
                      conf_int_2 = mean(conf_int_2),
                      width = mean(width),
                      lower_error_freq = mean(lower_error),
                      upper_error_freq = mean(upper_error),
                      coverage = mean(covered), .groups = "drop"))

write_csv(combined_perfect_assay_raw, "~/PrevAdjSurvey/code/simulations/confidence_interval_perfect_assay/combined_perfect_assay_raw.csv")

cat(setdiff(1:400, fnisihed_files), sep = ", ")

fnisihed_files <- dir_ls("/data/bayerdm/confidence_interval_perfect_assay/") %>%
  enframe(name = NULL) %>%
  mutate(file_name = path_file(value)) %>%
  filter(str_starts(file_name, "results_raw")) %>%
  mutate(file_num = as.integer(str_sub(file_name, 13, 16))) %>% pull(file_num)
range(fnisihed_files)
tmp <-
tibble(file_path = dir_ls("//data/bayerdm/fixed_weights_many")) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  filter(str_detect(file_name, "summary")) %>%
  arrange(file_name) %>%
  mutate(data = map(file_path, read_rds)) %>%
  pull(data) %>%
  bind_rows()

write_csv(tmp, "//data/bayerdm/fixed_weights_many/combined_fixed_weights_many_summary.csv")
