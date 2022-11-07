library(dplyr)
library(ggplot2)
library(gprofiler2)



mir_res <- readRDS("01-Exploratory_analysis_and_differential_expression/outputs/RDS/mir_res.rds")
targ_mir_mRNA <- readRDS("02-Regulation_targets/outputs/RDS/targ_mir_mRNA.rds")
cor_mir <- readRDS("02-Regulation_targets/outputs/RDS/cor_mir.rds")


mir_mrna_hist_data1 <- targ_mir_mRNA$mature_mirna_id %>% 
  table %>% 
  as.data.frame() %>%
  select(mir=".", Freq)

ggplot2::ggplot(mir_mrna_hist_data1, aes(x=Freq)) +
  geom_histogram(binwidth = 3, color="gray", fill="black") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab("MicroRNA")+
  ylab("Frequency")+
  ggtitle("Histogram of number of target genes per microRNA")+
  ggsave(filename = "Analyses/targ_mir_mRNA_hist1.jpeg")


  
mir_mrna_hist_data2 <- targ_mir_mRNA$target_ensembl %>% 
  table %>% 
  as.data.frame() %>%
  select(mir=".", Freq)


ggplot2::ggplot(mir_mrna_hist_data2, aes(x=Freq)) +
  geom_histogram(binwidth = 0.99, color="black", fill="darkred") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,19)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  xlab("Genes")+
  ylab("Frequency")+
  ggtitle("Histogram of number of microRNA per target gene")+
  ggsave(filename = "Analyses/targ_mir_mRNA_hist2.png")




mir93_gost <- 
  gprofiler2::gost(query = filter(targ_mir_mRNA, mature_mirna_id == "hsa-miR-93-5p") %>%
                     pull(target_ensembl))

gostplot(mir93_gost)

mir26b_gost <- 
  gprofiler2::gost(query = filter(targ_mir_mRNA, mature_mirna_id == "hsa-miR-26b-5p") %>%
                     pull(target_ensembl))

gostplot(mir26b_gost)
