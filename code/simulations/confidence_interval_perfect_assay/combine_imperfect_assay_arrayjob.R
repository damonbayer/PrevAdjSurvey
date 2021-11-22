library(tidyverse)
library(fs)

imperfect_assay_raw_vec <-
  combined_imperfect_assay_raw <-
  tibble(file_path = dir_ls("//data/bayerdm/confidence_interval_imperfect_assay")) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  filter(str_detect(file_name, "raw")) %>%
  arrange(file_name) %>%
  pull(file_path)

combined_imperfect_assay_raw <-
  map_dfr(imperfect_assay_raw_vec,
          ~read_rds(.) %>%
            mutate(width = conf_int_2 - conf_int_1) %>%
            group_by(design, method) %>%
            summarize(conf_int_1 = mean(conf_int_1),
                      conf_int_2 = mean(conf_int_2),
                      width = mean(width),
                      lower_error_freq = mean(lower_error),
                      upper_error_freq = mean(upper_error),
                      coverage = mean(covered), .groups = "drop"))
combined_imperfect_assay_raw

write_csv(combined_imperfect_assay_raw, "~/PrevAdjSurvey/code/simulations/confidence_interval_imperfect_assay/combined_imperfect_assay_raw.csv")


combined_imperfect_assay_summary <-
  tibble(file_path = dir_ls("//data/bayerdm/confidence_interval_imperfect_assay")) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  filter(str_detect(file_name, "summary")) %>%
  arrange(file_name) %>%
  mutate(data = map(file_path, read_rds)) %>%
  pull(data) %>%
  bind_rows()

combined_imperfect_assay_summary

write_csv(combined_imperfect_assay_summary, "~/PrevAdjSurvey/code/simulations/confidence_interval_imperfect_assay/combined_imperfect_assay_summary.csv")

experimental_design <-read_rds("//data/bayerdm/confidence_interval_imperfect_assay/experimental_design.rds")
experimental_design %>% select_if(~!is.list(.)) %>% write_csv("~/PrevAdjSurvey/code/simulations/confidence_interval_imperfect_assay/experimental_design.csv")


dir_ls("~/PrevAdjSurvey/code/simulations/confidence_interval_imperfect_assay") %>%
  enframe(name = NULL) %>%
  filter(str_starts(path_file(value), "slurm")) %>%
  deframe() %>%
  file_delete()
combined_imperfect_assay_raw <-
  tibble(file_path = dir_ls("//data/bayerdm/confidence_interval_imperfect_assay")) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  filter(str_detect(file_name, "raw")) %>%
  arrange(file_name) %>%
  mutate(data = map(file_path, read_rds)) %>%
  pull(data) %>%
  bind_rows()

combined_imperfect_assay_raw %>% group_by(method, design) %>% summarize(max(conf_int_2))

