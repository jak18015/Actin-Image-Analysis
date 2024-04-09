# Actin area analysis

# Packages ----
library(tidyverse)
library(ggpubr)

# FUNCTIONS ----
## create columns
treatmentColumn <- function(x) {
  ifelse(grepl(".*cytod.*", x, ignore.case = TRUE),
         "CytoD",
         ifelse(
           grepl(".*jas.*", x, ignore.case = TRUE),
           "Jas",
           ifelse(
             grepl(".*minus.*", x, ignore.case = TRUE),
             "Control",
             ifelse(
               grepl(".*plus.*", x, ignore.case = TRUE),
               "MyoF-KD",
               "ERROR NO TREATMENT"))))
}
trialColumn <- function(x) {
  ifelse(grepl(".*20211028.*", x), "trial1",
         ifelse(grepl(".*20221024.*", x), "trial2",
                ifelse(grepl(".*20221113.*", x), "trial3",
                       ifelse(grepl(".*20220627.*", x), "trial4",
                              ifelse(grepl(".*20230925.*", x), "trial5",
                                     ifelse(grepl(".*20230927.*", x), "trial6",
                                            ifelse(grepl(".*20231003.*", x), "trial7",
                                            "ERROR NO TRIAL")))))))
}
shapeColumn <- function(x) {
  ifelse(grepl("trial1", x), 23,
         ifelse(grepl("trial2", x), NA,
                ifelse(grepl("trial3", x), 21,
                       ifelse(grepl("trial4", x), NA,
                              ifelse(grepl("trial5", x), NA,
                                     ifelse(grepl("trial6", x), 22,
                                            ifelse(grepl("trial7", x), 24,
                                                   "ERROR NO TRIAL")))))))
}

frameCut <- function(x) {
  ifelse(x > 61, NA, x)
}
## significance values
symnum.args <- 
  list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
    symbols = c("****", "***", "**", "*", "ns"))

# PROCESSING ----
a <- setNames(
  read_csv(
    list.files(path = "results/08_analysis/segmentation/csv/",
               recursive = FALSE,
               pattern = "*.csv",
               full.names = TRUE)),
  c("frame", "img", "actArea"))

a <- a %>%
  group_by(img) %>%
  mutate(frame = frameCut(frame),
         img = gsub("_", "-", img))
a <- na.omit(a)
a <- a %>%
  summarize(actArea = mean(actArea))

v <- read_csv(
  list.files(path = "results/08_analysis/mask/csv/",
             recursive = FALSE,
             pattern = "*.csv",
             full.names = TRUE))
v <- setNames(subset(v,select=c("img", "area")),c("img", "vacArea"))
v <- v %>%
  group_by(img) %>%
  mutate(img = gsub("_", "-", img))
           

avArea <- merge(a, v, by = "img")
avArea <- na.omit(avArea)
rm(a,v)

avArea <-
  avArea %>%
  mutate(treatment = treatmentColumn(img),
         trial = trialColumn(img),
         shapeID = shapeColumn(trial),
         normArea = actArea/vacArea)
avArea <- na.omit(avArea)

parasiteArea <- avArea %>%
  group_by(trial, treatment, shapeID, img) %>%
  summarize(actArea = mean(actArea),
            vacArea = mean(vacArea),
            normArea = mean(normArea))


trialArea <- avArea %>%
  group_by(trial, treatment, shapeID) %>%
  summarize(actArea = mean(actArea),
            vacArea = mean(vacArea),
            normArea = mean(normArea),
            n = n())

parasiteCounts <- parasiteArea %>%
  group_by(treatment) %>%
  summarize(n = n())
# STATISTICS ----
t.test(normArea ~ treatment,
       data = trialArea,
       paired = TRUE,
       var.equal = TRUE)
# PLOTS ----
normAreaPlot <-
    ggplot(data = trialArea,
           aes(
             x = treatment,
             y = normArea,
             fill = treatment)
           ) +
    geom_boxplot(data = parasiteArea,
                 width = 0.4,
                 alpha = 0.5,
                 linewidth = 0.5,
                 outlier.alpha = 0
                 ) +
    geom_jitter(data = subset(parasiteArea, treatment == "Control"),
                aes(shape = as.numeric(shapeID)),
                color = "#000000",
                fill = "#555555",
                show.legend = FALSE,
                size = 1,
                stroke = 0.5,
                width = 0.05) +
  geom_jitter(data = subset(parasiteArea, treatment == "MyoF-KD"),
              aes(shape = as.numeric(shapeID)),
              color = "#000000",
              fill = "#E69F00",
              show.legend = FALSE,
              size = 1,
              stroke = 0.5,
              width = 0.05) +
    geom_point(data = subset(trialArea, treatment == "Control"),
               aes(shape = as.numeric(shapeID)),
               color = "#000000",
               fill = "#555555",
               size = 4,
               alpha = 0.75) +
  geom_point(data = subset(trialArea, treatment == "MyoF-KD"),
             aes(shape = as.numeric(shapeID)),
             color = "#000000",
             fill = "#E69F0055",
             size = 4,
             alpha = 0.75) +
  xlab("") +
  ylab("Normalized Area (μm^2)") +
    scale_shape_identity("Replicate") +
    theme_classic2() +
    scale_fill_colorblind() +
  theme(legend.position = "none",
        axis.text = element_text(color = "#000000", size = 10),
        axis.ticks = element_line(color = "#000000")) +
    stat_compare_means(data = trialArea,
                       aes(x = treatment,
                           y = normArea),
                       comparisons = list(c("Control", "MyoF-KD")),
                       method = "t.test",
                       label = "p.signif",
                       symnum.args = symnum.args,
                       family = "Arial"
                       )
normAreaPlot

stop()
