
temp <- readRDS( file.path(project_dir, "data/output/kmeans_clusters/20231024_0947__11clusters_100nstart_scaled.rds"))
temp_old <- readRDS(file.path(project_dir, "data/output/kmeans_clusters/res_v0/20221110_1148__11clusters_100nstart_scaled.rds"))
temp$old_kmeans <- temp_old[rownames(temp),]$kmeans
temp$old_dnam <- temp_old[rownames(temp),]$dnam_quot
temp <- mutate(temp, changed_kmeans = case_when(kmeans == old_kmeans ~ "No",
                                                kmeans ==7 & old_kmeans == 8 ~ "No",
                                                kmeans ==8 & old_kmeans == 7 ~ "No",
                                                TRUE ~ "Yes"))

filter(temp, changed_kmeans == "Yes") %>% pull(changed_kmeans) %>% table()
temp %>% pull(changed_kmeans) %>% table()
filter(temp, old_dnam == 0 & changed_kmeans == "Yes") %>% pull(changed_kmeans) %>% table()
