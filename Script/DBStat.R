library(tidyverse)
library(MSbox)
library(ggplot2)
library(RColorBrewer)
library(extrafont)
font_import() 
fonts()
loadfonts()

DB <- read_csv("Data/Progenesis.csv")
DB2<- DB %>%
  filter(!is.na(Class)) %>%
  select(Formula, Class) %>%
  add_count(Class, name = "classCount") %>%
  mutate(Class = ifelse(classCount >= 22, as.character(Class), "Other")) %>%
  mutate(Formula = str_remove_all(Formula, "[+]")) %>%
  mutate(Mass = MSbox::mass(Formula, caseSensitive = T)) %>%
  mutate(MRange = cut(Mass, breaks = seq(from = 60, to = 1000, by = 60))) %>%
  add_count(Class, name = "class2Count")

## barchart
p1 <- ggplot(DB2, aes(x = MRange, fill = Class)) + 
  geom_bar(colour = "black", size=0.1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +
  xlab("Mass Range") +
  ylab("Count") +
  theme(text = element_text(size = 12, family ="Times New Roman"),
        legend.text = element_text(size = 12)) 

ggsave("Writing/Figures/PDF/font_ggplot.pdf", plot = p1,  width = 8, height = 4)

## Pie chart
DB3 <- DB %>%
  filter(!is.na(Class)) %>%
  select(Formula, Class) %>%
  add_count(Class, name = "classCount") %>%
  mutate(Class = ifelse(classCount >= 22, as.character(Class), "Other")) %>%
  count(Class) %>%
  mutate(Percent = round(n/sum(n) * 100, 2)) %>%
  mutate(lab.ypos = cumsum(Percent) - 0.5 * Percent)

ggplot(DB3, aes(x = 2, y = Percent, fill = Class))+
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(palette = "Set3") +
  theme_void() + 
  xlim(0.5, 2.5)

## total mass features
DBinfo <- read_csv("Data/DBInfo.csv") %>%
  mutate(Mode = factor(Mode, levels = c("Positive", "Negative"))) 
ggplot(DBinfo, aes(x = Feature, y = Count, fill = Mode)) +
  geom_bar(stat = "identity", position="dodge", color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw()
  
  
