# Actin line scan analysis

# Libraries ====
library(tidyverse)
library(ggpubr)
# Functions ====
# parses the treatment condition from img column
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
# parses the trial number from the img column
trialColumn <- function(x) {
  ifelse(grepl(".*20211028.*", x), "trial1",
    ifelse(grepl(".*20221024.*", x), "trial2",
      ifelse(grepl(".*20221113.*", x), "trial3",
        ifelse(grepl(".*20220627.*", x), "trial4",
          "ERROR NO TRIAL"))))
}
# adds a shape corresponding to trial
shapeColumn <- function(x) {
  ifelse(grepl("trial1", x), 21,
    ifelse(grepl("trial2", x), 22,
      ifelse(grepl("trial3", x), 24,
        ifelse(grepl("trial4", x), 23,
          "ERROR NO TRIAL"))))
}
# PROCESSING ====
d0_rawCSV <- read_csv(list.files(
  path = "results/csv/",
  recursive = TRUE,
  pattern = "*.csv",
  full.names = TRUE
))

intDF <- d0_rawCSV %>%
  mutate(treatment = treatmentColumn(img),
         trial = trialColumn(img)) %>%
  pivot_longer(
    cols = !c(
      "trial",
      "treatment",
      "img",
      "actinCount",
      "px",
      "fStart",
      "fEnd"
    ),
    names_to = "frame",
    values_to = "intensity") 
jasCytoD <- read_csv("results/cyto-jas-intensities.csv")
jasCytoD <- jasCytoD %>%
  group_by(img, actinCount) %>%
  mutate(fStart = min(frame),
         fEnd = max(frame))
intDFPlusDrug <- rbind(intDF, jasCytoD)
intDF <- intDFPlusDrug %>%
  mutate(trial = trialColumn(img)) %>%
  mutate(frame = as.numeric(frame),
         trial = as.factor(trial),
         treatment = as.factor(treatment)) %>%
  group_by(img, actinCount) %>%
  filter(between(frame, fStart, fEnd))

rm(intDFPlusDrug, jasCytoD)

perPixelAlpha <-
  intDF %>%
  group_by(trial, treatment,img, actinCount, px) %>%
  summarize(alphaValue = (max(intensity) / mean(intensity))-1)

perParasiteAlpha <-
  perPixelAlpha %>%
  group_by(trial, treatment, img, actinCount) %>%
  summarize(meanAlphaValue = mean(alphaValue))
  
perTrialAlpha <-
  perParasiteAlpha %>%
  group_by(trial, treatment) %>%
  summarize(meanTrialAlphaValue = mean(meanAlphaValue))

## reorder
treatmentOrder <- c("Control", "MyoF-KD", "CytoD", "Jas")
trialOrder <- c("trial1", "trial2", "trial3", "trial4")
intDF$treatment <- factor(intDF$treatment,
                                 levels = treatmentOrder)
perPixelAlpha$treatment <- factor(perPixelAlpha$treatment,
                                  levels = treatmentOrder)
perParasiteAlpha$treatment <- factor(perParasiteAlpha$treatment,
                                  levels = treatmentOrder)
perTrialAlpha$treatment <- factor(perTrialAlpha$treatment,
                                  levels = treatmentOrder)
# shape identifiers
perTrialAlpha <- 
  perTrialAlpha %>%
  mutate(shapeID = shapeColumn(trial))

# counts
parasiteCounts <- perParasiteAlpha %>%
  mutate(ID = paste0(img, actinCount)) %>%
  group_by(treatment) %>%
  count()


# px vs intensity plot====
m34 <- subset(intDF, grepl("20211028-minusiaa_34", img))
m34 <- subset(m34, actinCount == 2)
for (i in 1:length(m34$px)) {
  m34$px[i] = m34$px[i] - 1
}

p07 <- subset(intDF, grepl("20211028-plusiaa_07", img))
p07 <- subset(p07, px > 4)
p07 <- subset(p07, actinCount == 2)
for (i in 1:length(p07$px)) {
  p07$px[i] = p07$px[i] - 5
}
c06 <- subset(intDF, grepl("CytoD_06", img))
c06 <- subset(c06, actinCount == 1)
for (i in 1:length(c06$px)) {
  c06$px[i] = c06$px[i] - 1
}
j11 <- subset(intDF, grepl("Jas_11", img))
j11 <- subset(j11, px > 3)
for (i in 1:length(j11$px)) {
  j11$px[i] = j11$px[i] - 4
}

pvi <- rbind(m34, p07, c06, j11)

pviPlot <- 
  ggplot(pvi)+
  geom_boxplot(aes(group = interaction(px, treatment),
                   x = px,
                   y = intensity/1000,
                   fill = treatment),
               lwd = 0.1,
               width = 0.9,
               outlier.size = 0.5,
               outlier.stroke = 0.1,
               alpha = 0.75,
               show.legend = FALSE)+
  scale_fill_colorblind()+
  scale_color_colorblind()+
  theme_classic2()+
  labs(fill = "Treatment")+
  xlab("Distance (px)")+
  ylab("Intensity (AU * 10^3)")+
  scale_x_continuous(breaks = c(1,10,20),
                     minor_breaks = c(1:20))+
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     labels = c("0", "5", "10", "15"),
                     limits = c(0,15),
                     expand = c(0,0))+
  facet_wrap(~treatment, 
             ncol = 4,
             scales = "free_y")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.ticks = element_line(color = "#000000"),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title =  element_text(size = 8),
        legend.position = "none",
        aspect.ratio = 0.75,
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.y = element_blank())

pviPlot

#Px Vs Alpha Plot ====
m34 <- subset(perPixelAlpha, grepl("20211028-minusiaa_34", img))
m34 <- subset(m34, actinCount == 2)
for (i in 1:length(m34$px)) {
  m34$px[i] = m34$px[i] - 1
}

p07 <- subset(perPixelAlpha, grepl("20211028-plusiaa_07", img))
p07 <- subset(p07, px > 4)
p07 <- subset(p07, actinCount == 2)
for (i in 1:length(p07$px)) {
  p07$px[i] = p07$px[i] - 5
}
c06 <- subset(perPixelAlpha, grepl("CytoD_06", img))
c06 <- subset(c06, actinCount == 1)
for (i in 1:length(c06$px)) {
  c06$px[i] = c06$px[i] - 1
}
j11 <- subset(perPixelAlpha, grepl("Jas_11", img))
j11 <- subset(j11, px > 3)
for (i in 1:length(j11$px)) {
  j11$px[i] = j11$px[i] - 4
}


pva <- rbind(m34, p07, c06, j11)

pva <- pva %>%
  group_by(img, treatment) %>%
  mutate(mean = mean(alphaValue))

pvaPlot <-
  ggplot(pva)+
  geom_line(aes(x = px,
                y = alphaValue,
                color = treatment))+
  geom_line(aes(x = px,
                y = mean,
                color = treatment),
            linetype = "dashed")+
  scale_fill_colorblind()+
  scale_color_colorblind()+
  theme_classic2()+
  labs(color = "Treatment")+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 10),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title =  element_text(size = 10),
        legend.position = "none")+
  xlab("Distance (px)")+
  ylab("Alpha")+
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3),
                     labels = c(0, 0.1, 0.2, 0.3),
                     limits = c(0,0.35),
                     expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))

pvaPlot
  
  
#Mean Image Treatment Boxplot ====
symnum.args <- 
  list(
    cutpoints = c(
      0, 0.0001, 0.001, 0.01, 0.05, Inf), 
    symbols = c(
      "****", "***", "**", "*", "ns"))

maPlot <-
  ggplot(data = perTrialAlpha,
         mapping = aes(x = treatment,
                       y = meanTrialAlphaValue))+
  geom_boxplot(
    data = perParasiteAlpha,
    aes(x = as.factor(factor(treatment, levels = treatmentOrder)),
        y = meanAlphaValue,
        fill = treatment),
    lwd = 0.25,
    alpha = 0.2,
    outlier.alpha = 0,
    width = 0.9,
    position = "identity"
  )+
  geom_jitter(
    data = perParasiteAlpha,
    aes(x = treatment,
        y = meanAlphaValue,
        fill = treatment),
    width = 0.15,
    shape = 21,
    size = 0.4,
    stroke = 0.3
  )+
  geom_point(
    data = perTrialAlpha %>% filter(treatment == "Control"),
    aes(fill = treatment,
        shape = shapeID),
    size = 2,
    alpha = 0.9,
    fill = "#3a3a3a",
    stroke = 0.5,
    show.legend = TRUE
  ) +
  geom_point(
    data = perTrialAlpha %>% filter(treatment == "MyoF-KD"),
    aes(fill = treatment,
        shape = shapeID),
    size = 2,
    alpha = 0.9,
    fill = "#e69f00ff",
    stroke = 0.5,
    show.legend = TRUE
  ) +
  scale_shape_identity("Replicate") +
  scale_fill_colorblind() +
  scale_color_colorblind() +
  theme_classic2() +
  xlab("Treatment") +
  ylab("Alpha") +
  scale_y_continuous(limits = c(0,0.35),
                     expand = c(0,0),
                     breaks = c(0, 0.1, 0.2, 0.3)
                     )+
  stat_compare_means(data = perTrialAlpha,
                     aes(x = treatment,
                         y = meanTrialAlphaValue),
                     method = "t.test",
                     paired = TRUE,
                     comparisons = list(c("Control", "MyoF-KD")),
                     symnum.args = symnum.args)+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(size = 10, angle = -30, color = "black", hjust = 0),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        legend.text = element_text(size = 10),
        legend.title =  element_text(size = 10),
        legend.position = "none",
        aspect.ratio = 1
  )

maPlot



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
