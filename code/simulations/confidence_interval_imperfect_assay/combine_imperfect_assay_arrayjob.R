library(tidyverse)
library(fs)


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

