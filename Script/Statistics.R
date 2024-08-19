library(MSbox)
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(RColorBrewer)

#(1) load the data. Data has been clean by CD by QC CV% <= 30
PL <- read_csv("Data/Fish_pos.csv") %>%
  filter(CV_QC <= 30)

DF <- PL %>%
  select(starts_with("Area:")) %>%
  t()

## get Group information
myGroup1 <- gsub("Area: (.+).raw.*", "\\1", rownames(DF))
myGroup2 <- as.factor(gsub('[0-9]+', '', myGroup1))
rownames(DF) <- myGroup1


#(2) Initial Data Quality
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 7))

##(1) heatmap
annoRow <- as.data.frame(myGroup2)
rownames(annoRow) <- rownames(DF)
colnames(annoRow) <- "Tissue"

pheatmap(scale(DF, center = T, scale = T), 
         cluster_cols = T,
         cluster_rows = T, 
         annotation_row = annoRow,
         show_rownames = T,
         show_colnames = F)

##(2) PCA
viewPCA(log10(DF), Group = myGroup2, centering = T, scaling = "auto", frame = T)


##(3) Ion intensity distribution
DF %>%
  log10() %>%
  viewTIC(Group = myGroup2)

##(4) Internal standard check

###(4.1) Type-I plot
ISTD <- PL %>%
  filter(Formula == "C9 H10 Cl N O2") %>%
  select(starts_with("Area:")) %>%
  log10()

ISTD2 <- cbind.data.frame(Int = t(ISTD), Tissue = myGroup2)
ggbarplot(ISTD2, x = "Tissue", y = "Int", 
          color = "Black", 
          fill = "grey80", 
          add = c("mean_sd", "jitter"), 
          add.params = list(alpha = 0.5),
          palette = mycolors) + 
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Log10 Peak Area") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 24)) 

###(4.2) Type-II plot
ggline(ISTD2, x = "Tissue", y = "Int", 
       color = "Blue", 
       add = c("mean_sd", "jitter"), 
       add.params = list(alpha = 0.3, size = 3),
       ylim = c(8, 11)) + 
  ylab("Log10 Peak Area") +
  geom_hline(yintercept = mean(ISTD2$Int), size = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = mean(ISTD2$Int) * 1.05, size = 1, color = "grey", linetype = "dotted") +
  geom_hline(yintercept = mean(ISTD2$Int) * 0.95, size = 1, color = "grey", linetype = "dotted") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 24)) 


#(3) Statistics
##(3.1) compounds highly abundant in Bu

myStat <- PL %>%
  select(-matches('AdjP.*QC')) %>%
  select(-matches('FC.*QC')) %>%
  filter(if_all(matches('AdjP.*Bu') , ~ . < 0.05)) %>%
  ## just by chance the FC here are all tissues/Bu
  filter(if_all(matches('FC.*Bu') , ~ . <= 0.1)) %>%
  select(Name, MW, RT, Formula, MS2, "mzCloud Best Match", matches('FC.*Bu'), starts_with('Area:')) %>%
  filter(MS2 == "DDA for preferred ion")

## Neg
PL <- read_csv("Data/Fish_Neg.csv") %>%
  filter(CV_QC <= 30)

a = doStat(DF, Group = myGroup2)

myStat <- a %>%
  mutate(Name = PL$Name, 
         Formula = PL$Formula,
         MW = PL$MW,
         RT = PL$RT) %>%
  bind_cols(select(PL, starts_with("Area: "))) %>%
  select(-matches('pValue_.*QC')) %>%
  select(-matches('Fold_.*QC')) %>%
  filter(if_all(matches('pValue_.*Bu') , ~ . < 0.05)) %>%
  ## just by chance the FC here are all tissues/Bu
  filter(if_all(matches('Fold_Bu.*') , ~ . >= 100) | if_all(matches('Fold_.*Bu') , ~ . <= 0.01)) %>%
  select(Name, MW, RT, Formula, matches('Fold_Bu.*'), matches('Fold_.*Bu'), starts_with('Area:'))


# Negative ----------------------------------------------------------------
PL2 <- read_csv("Data/Fish_Neg.csv") %>%
  filter(CV_QC <= 30)

DF <- PL2 %>%
  select(starts_with("Area:")) %>%
  t()

## get Group information
myGroup1 <- gsub("Area: (.+).raw.*", "\\1", rownames(DF))
myGroup2 <- as.factor(gsub('[0-9]+', '', myGroup1))
rownames(DF) <- myGroup1

ISTD <- PL2 %>%
  filter(Formula == "C9 H10 Cl N O2") %>%
  select(starts_with("Area:")) %>%
  log10()

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 7))

ISTD2 <- cbind.data.frame(Int = t(ISTD), Tissue = myGroup2)
ggline(ISTD2, x = "Tissue", y = "Int", 
       color = "Blue", 
       add = c("mean_sd", "jitter"), 
       add.params = list(alpha = 0.3, size = 3)) + 
  ylab("Log10 Peak Area") +
  geom_hline(yintercept = mean(ISTD2$Int), size = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = mean(ISTD2$Int) * 1.1, size = 1, color = "grey", linetype = "dotted") +
  geom_hline(yintercept = mean(ISTD2$Int) * 0.90, size = 1, color = "grey", linetype = "dotted")
