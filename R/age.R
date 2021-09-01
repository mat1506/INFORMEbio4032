library(here)
library(tidyverse)
library(rbacon)

#Bacon(coredir = here::here("analysis/cores_bacon"))


rbacon::Bacon(core = "TL", coredir = here::here("analysis/cores_bacon"),
              ask = FALSE, plot.pdf = T,depths.file = T,
              thick = 2,ssize=1000,
              hiatus.depths=c(175),
              acc.mean=c(50,80),
              acc.shape=c(1.4,1),
              slump=c(175,183))

X0 <-read_csv(here::here("analysis/data/raw_data/2021-08-31_informeMF2021_data_V01.csv"))

X1 <- read_delim(here::here("analysis/cores_bacon/TL/TL_111_ages.txt"))



if(!file.exists(here::here("analysis/data/derived_data/datagemodel.csv"))){
  datagemodel <- bind_cols(X1[2:5],X0)
  path_out <- here::here("analysis/data/derived_data/","datagemodel.csv")
  write_csv(datagemodel,path_out)
}
