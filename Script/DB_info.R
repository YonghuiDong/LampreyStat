library(MSbox)
library(tidyverse)
## pos =========================================================================
DF1 <- read_csv("Data/Fish_pos.csv") %>%
  filter(CV_QC <= 30) %>%
  filter(MS2 == "DDA for preferred ion") %>%
  filter(`mzCloud Best Match` > 80)


## neg  ========================================================================
DF2 <- read_csv("Data/Fish_neg.csv") %>%
  filter(CV_QC <= 30) %>%
  filter(MS2 == "DDA for preferred ion") %>%
  filter(`mzCloud Best Match` > 80)
 

DF1 <- read_csv("Data/Fish_pos.csv") %>%
  filter(CV_QC <= 30) %>%
  filter(MS2 != "No MS2")

DF2 <- read_csv("Data/Fish_neg.csv") %>%
  filter(CV_QC <= 30) %>%
  filter(MS2 != "No MS2")
