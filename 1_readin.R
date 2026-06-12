library(here)
library(scuttle)


if (!file.exists(here("output", "sample.list"))){
  file.list <- paste0("Y:/Sharmilla/human_sc_data/calvanese_aorta/", c("GSM4968831_Aorta-4wk-658.csv.gz", "GSM4968832_Aorta-5wk-555.csv.gz", # 555 is 5 weeks
                                                                            "GSM4968833_Aorta-5wk-575.csv.gz"))  # 575 is 5.5 week
                                                                       #"GSM4968834_Aorta-6wk-563.csv.gz")) 
  names(file.list) <- c("CS14_wk4", "CS15_wk5", "CS15_wk5.5")#, "week6")                  
  
  sample.list <- lapply(file.list, readSparseCounts, sep = ",", row.names = TRUE)
  
  saveRDS(sample.list, file = here("output", "sample.list.RDS"))
} else {
  sample.list <- readRDS(here("output", "sample.list"))
}

