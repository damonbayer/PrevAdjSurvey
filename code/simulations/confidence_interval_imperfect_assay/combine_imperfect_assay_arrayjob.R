library(tidyverse)
library(fs)


tmp <-
  tibble(file_path = dir_ls("//data/bayerdm/confidence_interval_imperfect_assay")) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  filter(str_detect(file_name, "summary")) %>%
  arrange(file_name) %>%
  mutate(data = map(file_path, read_rds)) %>%
  pull(data) %>%
  bind_rows()

write_csv(tmp, "//data/bayerdm/fixed_weights_many/combined_fixed_weights_many_summary.csv")
