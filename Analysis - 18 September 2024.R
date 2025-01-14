# Installing Packages ----
install.packages("corrplot") # Used to visualize correlation matrices and create correlograms to identify data patterns.
install.packages("ggrepel") #  Enhances ggplot2 by adding labels that repel each other, avoiding overlap in plots.
install.packages("Hmisc") # Provides tools for data manipulation, descriptive statistics, and creating correlation matrices with significance tests.
install.packages("VennDiagram") # Creates Venn diagrams to visualize intersections, unions, and unique elements among sets.
install.packages("devtools") # Facilitates package development and installation, including GitHub-based packages.
install.github("vqv/ggbiplot")# Extends ggplot2 to visualize PCA results through biplots for exploring variable and observation relationships.
install_github("kassambara/factoextra") # For Grouping in PCA
install.packages("ggeffects") # Simplifies the creation of effect plots for regression models to interpret model predictions.
install.packages("igraph") # Focused on network analysis and visualization, allowing the creation and analysis of graph structures.
install.packages("qgraph") # Builds visualizations for networks and integrates with igraph, often used for correlation networks and psychological studies.
install.packages("vegan")# Provides tools for ecological data analysis, including ordination, diversity analysis, and multivariate statistics.
install.packages("factoextra") # For Grouping in PCA
library(readr) # Offers fast and user-friendly functions for importing datasets (e.g., CSV, TSV) into R.
library(lme4) # Fits linear and generalized linear mixed-effects models, useful for hierarchical data or random effects.
library(lmerTest) # Extends lme4 by adding p-values and other summaries for mixed models.
library(ggplot2) # A fundamental package for creating customizable and advanced visualizations based on the grammar of graphics.
library(vegan)
library(dplyr) # A data manipulation package with easy-to-use functions and piping (%>%) for filtering and summarizing data.
library(corrplot) # For Correlation Matrix and Corellogram
library(igraph) # For Correlation Matrix
library(qgraph) # For Correlation Matrix
library(qgraph)
library(igraph)
library(Hmisc)
library(ggrepel)
library(tidyr) # Simplifies reshaping data (wide to long format) and managing missing values, often used with dplyr.
library(ggbiplot)
library(devtools)
library(ggbiplot)
library(plyr) # Provides functions to split, apply, and combine operations on datasets, though less commonly used due to dplyr.
library(ggeffects)
library(car) # Provides tools for regression analysis, ANOVA, and diagnostic plots, including variance partitioning.
library(reshape2) # Offers functions for reshaping data frames, such as melting and casting, to prepare tidy datasets.
library(vegan)
library(factoextra)
View(AusC)

# Loading Files into R ----
setwd("/Users/murugash/Desktop/PhD/Austria Experiment/R Analysis/Crepis and Lotus Samples - R Analysis")
AusL <- read_csv("Austria Experiment - Lotus Data - Core Response Variables - 15 August 2024.csv")
AusC <- read_csv("Austria Experiment - Crepis Data - Core Response Variables - 15 August 2024.csv")
AusC <- AusC[-c(63:163),]
AusPC <- read_csv("Austria Experiment - Crepis and Leaf Data - RDA - 15 August 2024.csv")
AusPC <- AusPC[-c(85:117),]
AusRD <- read_csv("Austria Experiment - Crepis and Leaf Data - RDA - 15 August 2024.csv")
AusRD <- AusRD[-c(85:117),]
AusL1 <- AusL[AusL$`Plot No.`!= "P16",]
View(AusRD)
dir()
AusC

# Function to Log and Square Root transform Response Variables ----
log_sqrt_transform <- function(df, cols) {
  # Check if cols are numeric (indices)
  if (is.numeric(cols)) {
    cols <- names(df)[cols]
  }
  
  for (col in cols) {
    log_col_name <- paste0("log_", col)
    sqrt_col_name <- paste0("sqrt_", col)
    df[[log_col_name]] <- log(df[[col]])
    df[[sqrt_col_name]] <- sqrt(df[[col]])
  }
  
  return(df)
}

AusL <- log_sqrt_transform(AusL, c(9:23))
AusC <- log_sqrt_transform(AusC, c(9:17))
View(AusC)
# Function to create Histogram of Response Variables ----

hist_df <- function(df, cols) {
  # Check if cols are numeric (indices)
  if (is.numeric(cols)) {
    cols <- names(df)[cols]
  }
  for(col in cols) {
    hist(df[[col]], main = paste0("Histogram of ", col), xlab = col)
  }
}

# Histogram of Response Variables ----

hist_df(AusL, c(9:53))
hist_df(AusC, c(9:35))
hist_df(AusPC, c(9:25))

# Shapiro-Wilk test for Response Variables ----
Lotus_shapiro_results <- apply(AusL[,c(9:53)], 2, shapiro.test)
Lotus_shapiro_results
Crepis_shapiro_results <- apply(AusC[,c(9:35)], 2, shapiro.test)
Crepis_shapiro_results
colnames(AusC)

# Converting Grouping variables, fixed effect predictor variables and response variables to Factors ----
## In Lotus Dataset ----
AusL$CO2 <- as.factor(AusL$CO2)
AusL$Temperature <- as.factor(AusL$Temperature)
AusL$Drought <- as.factor(AusL$Drought)
AusL$`Plot No.` <- as.factor(AusL$`Plot No.`)

## In Crepis Dataset ----
AusC$CO2 <- as.factor(AusC$CO2)
AusC$Temperature <- as.factor(AusC$Temperature)
AusC$Drought <- as.factor(AusC$Drought)
AusC$`Plot No.` <- as.factor(AusC$`Plot No.`)

## In PCA Dataset ----
AusPC$CO2 <- as.factor(AusPC$CO2)
AusPC$Temperature <- as.factor(AusPC$Temperature)
AusPC$Drought <- as.factor(AusPC$Drought)
AusPC$`Plot No.` <- as.factor(AusPC$`Plot No.`)
AusPC$Species <- as.factor(AusPC$Species)

## In RDA Dataset----
AusRD$`Plot No.` <- as.factor(AusRD$`Plot No.`)
AusRD$CO2 <- as.factor(AusRD$CO2)
AusRD$Temperature <- as.factor(AusRD$Temperature)
AusRD$Drought <- as.factor(AusRD$Drought)
AusRD$Species <- as.factor(AusRD$Species)

names(AusC)

# Preparation of Mixed Effect Model ----
## Mixed-effect Models for Response Variables of Lotus ----
l_DAmodel <- lmer(`Display Area (DA)(cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
log_SSPAmodel <- lmer(`log_Specific Standard Petal Area (SSA)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
log_SW1PAmodel <- lmer(`log_Specific Wing 1 Petal Area (SW1A)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
SW2PAmodel <- lmer(`Specific Wing 2 Petal Area (SW2A)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
sq_SKPAmodel <- lmer(`sqrt_Specific Keel Petal Area (SKA)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
log_SPMAmodel <- lmer(`log_Standard Petal Mass per Area (SPMA)(g/cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
log_W1PMAmodel <- lmer(`log_Wing 1 Petal Mass per Area (W1PMA)(g/cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
log_W2PMAmodel <- lmer(`log_Wing 2 Petal Mass per Area (W2PMA)(g/cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
log_KPMAmodel <- lmer(`log_Keel Petal Mass per Area (KPMA)(g/cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
l_sq_SPAmodel <- lmer(`sqrt_Specific Petal Area (SPA)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
l_log_PMAmodel <- lmer(`log_Petal Mass per Area (PMA)(g/cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
l_PDMCmodel <- lmer(`Petal Dry Matter Content (PDMC)(g/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
l_LAmodel <- lmer(`Leaf Area (LA)(cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
l_SLAmodel <- lmer(`Specific Leaf Area (SLA)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)
l_log_LDMCmodel <- lmer(`log_Leaf Dry Matter Content (LDMC)(g/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusL)

## Mixed Effect Models for Response Variables of Crepis ----
c_sq_DAmodel <- lmer(`Display Area (DA)(cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
c_sq_SPAmodel <- lmer(`sqrt_Specific Petal Area (SPA)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
c_log_PMAmodel <- lmer(`log_Petal Mass Per Area (PMA)(g/cm2)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
c_sq_PDMCmodel <- lmer(`sqrt_Petal Dry Matter Content (PDMC)(g/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
c_log_LAmodel <- lmer(`log_Leaf Area (LA)(cm2)`~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
c_log_SLAmodel <- lmer(`log_Specific Leaf Area (SLA)(cm2/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
c_sq_LDMCmodel <- lmer(`sqrt_Leaf Dry Matter Content (LDMC)(g/g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
sq_SNmodel <- lmer(`sqrt_Number of Seeds (SN)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)
log_SMmodel <- lmer(`log_Seed Mass (SM)(g)` ~ CO2 * Temperature * Drought + (1|`Plot No.`), data = AusC)


# Preparation of Linear Models ----
## Linear Models for Response Variables of Lotus ----
l_DAmodel <- lm(`Display Area (DA)(cm2)` ~ CO2 * Temperature * Drought, data = AusL)
log_SSPAmodel <- lm(`log_Specific Standard Petal Area (SSA)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusL)
log_SW1PAmodel <- lm(`log_Specific Wing 1 Petal Area (SW1A)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusL)
SW2PAmodel <- lm(`Specific Wing 2 Petal Area (SW2A)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusL)
sq_SKPAmodel <- lm(`sqrt_Specific Keel Petal Area (SKA)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusL)
log_SPMAmodel <- lm(`log_Standard Petal Mass per Area (SPMA)(g/cm2)` ~ CO2 * Temperature * Drought, data = AusL)
log_W1PMAmodel <- lm(`log_Wing 1 Petal Mass per Area (W1PMA)(g/cm2)` ~ CO2 * Temperature * Drought, data = AusL)
log_W2PMAmodel <- lm(`log_Wing 2 Petal Mass per Area (W2PMA)(g/cm2)` ~ CO2 * Temperature * Drought, data = AusL)
log_KPMAmodel <- lm(`log_Keel Petal Mass per Area (KPMA)(g/cm2)` ~ CO2 * Temperature * Drought, data = AusL)
l_sq_SPAmodel <- lm(`sqrt_Specific Petal Area (SPA)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusL)
l_log_PMAmodel <- lm(`log_Petal Mass per Area (PMA)(g/cm2)` ~ CO2 * Temperature * Drought, data = AusL)
l_PDMCmodel <- lm(`Petal Dry Matter Content (PDMC)(g/g)` ~ CO2 * Temperature * Drought, data = AusL)
l_LAmodel <- lm(`Leaf Area (LA)(cm2)` ~ CO2 * Temperature * Drought, data = AusL)
l_SLAmodel <- lm(`Specific Leaf Area (SLA)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusL)
l_log_LDMCmodel <- lm(`log_Leaf Dry Matter Content (LDMC)(g/g)` ~ CO2 * Temperature * Drought, data = AusL)

## Linear Models for Response Variables of Crepis ----
c_sq_DAmodel <- lm(`Display Area (DA)(cm2)`~ CO2 * Temperature * Drought, data = AusC)
c_sq_SPAmodel <- lm(`sqrt_Specific Petal Area (SPA)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusC)
c_log_PMAmodel <- lm(`log_Petal Mass Per Area (PMA)(g/cm2)` ~ CO2 * Temperature * Drought, data = AusC)
c_sq_PDMCmodel <- lm(`sqrt_Petal Dry Matter Content (PDMC)(g/g)` ~ CO2 * Temperature * Drought, data = AusC)
c_log_LAmodel <- lm(`log_Leaf Area (LA)(cm2)`~ CO2 * Temperature * Drought, data = AusC)
c_log_SLAmodel <- lm(`log_Specific Leaf Area (SLA)(cm2/g)` ~ CO2 * Temperature * Drought, data = AusC)
c_sq_LDMCmodel <- lm(`sqrt_Leaf Dry Matter Content (LDMC)(g/g)` ~ CO2 * Temperature * Drought, data = AusC)
sq_SNmodel <- lm(`sqrt_Number of Seeds (SN)` ~ CO2 * Temperature * Drought, data = AusC)
log_SMmodel <- lm(`log_Seed Mass (SM)(g)` ~ CO2 * Temperature * Drought, data = AusC)



# ANOVA Analysis ----
## Anova test for Reponse Variables in Lotus ----
anova(l_DAmodel) # Significance: Temperature and CO2:Drought in ME; Temperature, CO2:Temperature and CO2:Drought in LM
anova(log_SSPAmodel) # Significance: CO2: Drought marginally significant in LM
anova(log_SW1PAmodel)
anova(SW2PAmodel) # Significance: Drought in LM
anova(sq_SKPAmodel)
anova(log_SPMAmodel) # Significance: CO2: Drought marginally significant in LM
anova(log_W1PMAmodel)
anova(log_W2PMAmodel) # Significance: Drought and CO2: Temperature in LM (Latter marginally significant)
anova(log_KPMAmodel)
anova(l_sq_SPAmodel) # Significance: Drought marginally significant in LM
anova(l_log_PMAmodel) # Significance: Drought marginally significant in LM
anova(l_PDMCmodel)
anova(l_LAmodel) # Significance: CO2 and CO2:Drought in ME; Same in LM
anova(l_SLAmodel) # Significance: Drought and CO2: Temperature in ME (Latter marginally significant); CO2, Drought and CO2: Temperature in LM
anova(l_log_LDMCmodel) # Significance: CO2, Drought and CO2:Drought in ME (Latter marginally significant); CO2, Drought and CO2: Drought in LM

## Anova test for Response Variables in Crepis ----
anova(c_sq_DAmodel) # Significance: Drought in ME; Drought in LM
anova(c_sq_SPAmodel) # Significance: Temperature and CO2: Temperature both marginally significant in ME; Drought and CO2: Temperature in LM (latter marginally significant) 
anova(c_log_PMAmodel) # Significance: Temperature marginally significant in ME; Drought in LM
anova(c_sq_PDMCmodel) # Significance: Temperature marginally significant in ME; Drought and CO2: Temperature in LM(Latter marginally significant)
anova(c_log_LAmodel) # Significance: Temperature marginally significant in LM
anova(c_log_SLAmodel)
anova(c_sq_LDMCmodel) # Significance: Drought and CO2: Drought in LM (Latter marginally significant)
anova(sq_SNmodel) # Significance: Temperature, Drought and CO2: Temperature in ME (Latter marginally significant); Same in LM
anova(log_SMmodel)


# Plotting Significant Response Variables against Predictor Variables ----
## Plots in Lotus ----
plot(AusL$`Display Area (DA)(cm2)`~AusL$Temperature, main = "Temperature on DA", xlab = "Temperature Level", ylab = "Display Area (DA)(cm2)")
plot(AusL$`Specific Petal Area (SPA)(cm2/g)`~AusL$Drought, main= "Drought on SPA", xlab = "Drought Level", ylab = "Specific Petal Area (SPA) (cm2/g)")
plot(AusL$`Leaf Area (LA)(cm2)`~AusL$CO2, main = "CO2 on LA", xlab = "CO2 Level", ylab = "Leaf Area (LA)(cm2)")
plot(AusL$`Specific Leaf Area (SLA)(cm2/g)`~AusL$CO2, main = "CO2 on SLA", xlab = "CO2 Level", ylab = "Specific Leaf Area (SLA)(cm2/g)")
plot(AusL$`Specific Leaf Area (SLA)(cm2/g)`~AusL$Drought, main = "Drought on SLA", xlab = "Drought Level", ylab = "Specific Leaf Area (SLA)(cm2/g)")
plot(AusL$`Leaf Dry Matter Content (LDMC)(g/g)`~AusL$CO2, main = "CO2 on LDMC", xlab = "CO2 Level", ylab = "Leaf Dry Matter Content (LDMC)(g/g)")
plot(AusL$`Leaf Dry Matter Content (LDMC)(g/g)`~AusL$Drought, main = "Drought on LDMC", xlab = "Drought Level", ylab = "Leaf Dry Matter Content (LDMC)(g/g)")

## Plots in Lotus (Using ggplot) ----
### CO2 ----
ggplot(AusL, aes(x = AusL$CO2, y = AusL$`Leaf Area (LA)(cm2)`))+
  geom_boxplot(fill = "#f7766d") +
  labs(
    title = "Lotus - Effect of CO2 on LA",
    x = "CO2 Level",
    y = "Leaf Area (LA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

ggplot(AusL, aes(x = CO2,y= `Specific Leaf Area (SLA)(cm2/g)`)) +
  geom_boxplot(fill = "#f7766d") +
  labs(
    title = "Lotus - Effect of CO2 on SLA",
    x = "CO2 Level",
    y = "Specific Leaf Area (SLA) (cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

ggplot(AusL, aes(x = CO2,y= `Leaf Dry Matter Content (LDMC)(g/g)`)) +
  geom_boxplot(fill = "#f7766d") +
  labs(
    title = "Lotus - Effect of CO2 on LDMC",
    x = "CO2 Level",
    y = "Leaf Dry Matter Content (LDMC) (g/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

### Temperature ----
ggplot(AusL, aes(x = AusL$Temperature, y = AusL$`Display Area (DA)(cm2)`))+
  geom_boxplot(fill = "#f7766d") +
  labs(
    title = "Lotus - Effect of Temperature on DA",
    x = "Temperature Level",
    y = "Display Area (DA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

### Drought ----
ggplot(AusL, aes(x = AusL$Drought, y= AusL$`Specific Leaf Area (SLA)(cm2/g)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Lotus - Effect of Drought on SLA",
    x = "Drought Level",
    y = "Specific Leaf Area (SLA) (cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

ggplot(AusL, aes(x = Drought, y= AusL$`Leaf Dry Matter Content (LDMC)(g/g)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Lotus - Effect of Drought on LDMC",
    x = "Drought Level",
    y = "Leaf Dry Matter Content (LDMC) (g/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

### Interactions ----

ggplot(AusL, aes(x = CO2, y = `Display Area (DA)(cm2)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Lotus - Effect of CO2 and Temperature on DA",
    x = "CO2 Level",
    y = "Display Area (DA)(cm2)",
    fill = "Temperature"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )


ggplot(AusL, aes(x = CO2, y = `Specific Leaf Area (SLA)(cm2/g)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Lotus - Effect of CO2 and Temperature on SLA",
    x = "CO2 Level",
    y = "Specific Leaf Area (SLA) (cm2/g)",
    fill = "Temperature"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )
ggplot(AusL, aes(x = CO2, y = AusL$`Display Area (DA)(cm2)`, fill = Drought)) +
  geom_boxplot() +
  facet_wrap(~ Drought) +
  labs(
    title = "Lotus - Effect of Combined Interaction on DA",
    x = "CO2 + Temperature Level",
    y = "Display Area (DA)(cm2)",
    fill = "Drought"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )

ggplot(AusL, aes(x = CO2, y = AusL$`Leaf Area (LA)(cm2)`, fill = Drought)) +
  geom_boxplot() +
  facet_wrap(~ Drought) +
  labs(
    title = "Lotus - Effect of Combined Interaction on LA",
    x = "CO2 + Temperature Level",
    y = "Leaf Area (LA)(cm2)",
    fill = "Drought"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )

ggplot(AusL, aes(x = CO2, y = `Leaf Dry Matter Content (LDMC)(g/g)`, fill = Drought)) +
  geom_boxplot() +
  facet_wrap(~ Drought) +
  labs(
    title = "Lotus - Effect of Combined Interaction on LDMC",
    x = "CO2 + Temperature Level",
    y = "Leaf Dry Matter Content (LDMC)(g/g)",
    fill = "Drought"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )

anova(lv_LAmodel)
## Plots in Crepis ----
plot(AusC$`Display Area (DA)(cm2)`~AusC$Drought, main = "Drought on DA", xlab = "Drought Level", ylab = "Display Area (DA)(cm2)")
plot(AusC$`Specific Petal Area (SPA)(cm2/g)`~AusC$Drought, main = "Drought on SPA", xlab = "Drought Level", ylab = "Specific Petal Area (SPA) (cm2/g)")
plot(AusC$`Petal Dry Matter Content (PDMC)(g/g)`~AusC$Drought, main = "Drought on PDMC", xlab = "Drought Level", ylab = "Petal Dry Matter Content (PDMC) (g/g)")
plot(AusC$`Leaf Area (LA)(cm2)`~AusC$Temperature, main = "Temperature on Leaf Area", xlab = "Temperature Level", ylab = "Leaf Area (LA)(cm2)")
plot(AusC$`Leaf Dry Matter Content (LDMC)(g/g)`~AusC$Drought, main = "Drought on LDMC", xlab = "Drought Level", ylab = "Leaf Dry Matter Content (PDMC) (g/g)")
plot(AusC$`Number of Seeds (SN)`~AusC$Temperature, main = "Temperature on No. of Seeds", xlab = "Temperature Level", ylab = "No. of Seeds")
plot(AusC$`Number of Seeds (SN)`~AusC$Drought, main = "Drought on No. of Seeds", xlab = "Drought", ylab = "No. of Seeds")

## Plots in Crepis (Using ggplot) ----
View(AusC)

### Drought ----
ggplot(AusC, aes(x = AusC$Drought, y= `Display Area (DA)(cm2)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Crepis - Effect of Drought on DA",
    x = "Drought Level",
    y = "Display Area (DA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

ggplot(AusC, aes(x = AusC$Drought, y= `Specific Petal Area (SPA)(cm2/g)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Crepis - Effect of Drought on SPA",
    x = "Drought Level",
    y = "Specific Petal Area (SPA)(cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

ggplot(AusC, aes(x = AusC$Drought, y= `Petal Dry Matter Content (PDMC)(g/g)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Crepis - Effect of Drought on PDMC",
    x = "Drought Level",
    y = "Petal Dry Matter Content (PDMC)(g/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

ggplot(AusC, aes(x = AusC$Drought, y= `Leaf Dry Matter Content (LDMC)(g/g)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Crepis - Effect of Drought on LDMC",
    x = "Drought Level",
    y = "Leaf Dry Matter Content (g/g)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  ) 

ggplot(AusC, aes(x = AusC$Drought, y= `Number of Seeds (SN)`)) +
  geom_boxplot(fill = "#f7766d")+
  labs(
    title = "Crepis - Effect of Drought on SN",
    x = "Drought Level",
    y = "Number of Seeds (SN)"
  ) +
  theme(
    plot.title = element_text(size = 36, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 32),  # Increase x-axis title size
    axis.title.y = element_text(size = 32),  # Increase y-axis title size
    axis.text.x = element_text(size = 28),   # Increase x-axis text size
    axis.text.y = element_text(size = 28)    # Increase y-axis text size
  )

### Interactions ----

ggplot(AusC, aes(x = CO2, y = `Specific Petal Area (SPA)(cm2/g)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Crepis - Effect of CO2 and Temperature on SPA",
    x = "CO2 Level",
    y = "Specific Petal Area (SPA)(cm2/g)",
    fill = "Temperature"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )

ggplot(AusC, aes(x = CO2, y = `Number of Seeds (SN)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Crepis - Effect of CO2 and Temperature on SN",
    x = "CO2 Level",
    y = "Number of Seeds (SN)",
    fill = "Temperature"
  ) +
  theme(
    plot.title = element_text(size = 30, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 28),  # Increase x-axis title size
    axis.title.y = element_text(size = 28),  # Increase y-axis title size
    axis.text.x = element_text(size = 26),   # Increase x-axis text size
    axis.text.y = element_text(size = 26),   # Increase y-axis text size
    legend.title = element_text(size = 24),  # Increase legend title size
    legend.text = element_text(size = 22)    # Increase legend text size
  )

# Principal Component Analysis ----
## Creation of PCA models ----
View(AusPC)
colnames(AusPC)
lmv <- prcomp(AusPC[AusPC$Species == "Lotus", c(9,18,20:23)], scale = TRUE)
cmv <- prcomp(AusPC[AusPC$Species == "Crepis", c(9,18,20:23)], scale = TRUE)
cav <- prcomp(AusPC[AusPC$Species == "Crepis", c(9,18,20:25)], scale =  TRUE)

lcav <- prcomp(AusPC[,c(16, 18,19,20)], center = TRUE, scale = TRUE)
names(AusPC)

colnames(AusPC)
View(AusPC)
## Plotting of PCA models ---- 
names(AusPC)[names(AusPC) == "Specific Petal Area (SPA)(cm2/g)"] <- "SPA"
names(AusPC)[names(AusPC) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "PDMC"
names(AusPC)[names(AusPC) == "Specific Leaf Area (SLA)(cm2/g)"] <- "SLA"
names(AusPC)[names(AusPC) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "LDMC"
names(AusPC)[names(AusPC) == "Number of Seeds (SN)"] <- "SN"
names(AusPC)[names(AusPC) == "Seed Mass (SM)(g)"] <- "SM"
names(AusPC)[names(AusPC) == "Display Area (DA)(cm2)"] <- "DA"
names(AusPC)[names(AusPC) == "Leaf Area (LA)(cm2)"] <- "LA"

biplot(lmv)
biplot(cmv, scale = 1, alpha = 0)
biplot(lcav, scale = 0)
scores(lmv, display = "species")
?biplot
lmv_scores <- lmv$x
lmv_scores[,3]
summary(lmv)
summary(cmv)
summary(cav)
lmv

### PCA plots using ggbiplot ----

ggbiplot(lmv, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(x = "PC1 (32.7%)", y = "PC2 (28.3%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16)   # Increase size of y-axis title
  )

ggbiplot(cmv, varname.size = 6 ,alpha = 0) + 
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(x = "PC1 (33.6%)", y = "PC2 (21.6%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16)   # Increase size of y-axis title
  )

ggbiplot(cav, varname.size = 6 ,alpha = 0) + 
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(x = "PC1 (28%)", y = "PC2 (18.6%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16)   # Increase size of y-axis title
  )


## Grouping based on Climatic Variables ----
## Using Ordispider Package ----
### For Lotus ----
l_pca_scores <- as.data.frame(lmv$x)   # Extract PCA individual scores
l_groups <- AusPC$Treatment[AusPC$Species == "Lotus"]  # Grouping variable

# Define a consistent color palette for the polygons
unique_groups <- unique(l_groups)
group_colors <- setNames(rainbow(length(unique_groups)), unique_groups)

# Base plot
plot(l_pca_scores$PC1, l_pca_scores$PC2, 
     xlab = "PC1 (32.7%)", ylab = "PC2 (28.3%)", main = "Lotus - PCA",
     col = group_colors[l_groups], pch = 19, cex = 1.5)

# Add spider diagram (ordihull polygons)
ordispider(l_pca_scores, l_groups, col = rainbow(length(unique(l_groups))), label = TRUE)
ordihull(l_pca_scores, l_groups, col = group_colors, draw = "polygon", alpha = 0.4)

# Add legend matching polygon colors
legend("topright",
       legend = unique_groups,
       fill = group_colors,
       border = "black",
       title = "Treatment Groups")

# Adding arrows
l_pca_loadings <- as.data.frame(lmv$rotation[, 1:2])  
arrows(0, 0, 
       l_pca_loadings$PC1 * 2,  
       l_pca_loadings$PC2 * 2,  
       col = "blue", length = 0.1)
text(l_pca_loadings$PC1 * 2.2, l_pca_loadings$PC2 * 2.2, labels = rownames(l_pca_loadings), col = "blue", cex = 0.8)

### For Crepis (Main Variables) ----
cm_pca_scores <- as.data.frame(cmv$x)
c_groups <- AusPC$Treatment[AusPC$Species == "Crepis"]

#Base plot
plot(cm_pca_scores$PC1, cm_pca_scores$PC2,
     xlab = "PC1 (33.6%)", ylab = "PC2 (21.6%)", main = "Crepis - Main Variables - PCA",
     col=as.factor(c_groups), pch = 19, cex = 1.5)

# Add spider diagram
ordispider(cm_pca_scores, c_groups, col = rainbow(length(unique(c_groups))), label = TRUE)
ordihull(cm_pca_scores, c_groups, col = rainbow(length(unique(c_groups))),draw = "polygon")

# Optional: Add legend
legend("topright", legend = unique(c_groups), 
       col = rainbow(length(unique(c_groups))), 
       pch = 19, bty = "o")

# Adding arrows
cm_pca_loadings <- as.data.frame(cmv$rotation[, 1:2])  
arrows(0, 0, 
       cm_pca_loadings$PC1 * 2,  
       cm_pca_loadings$PC2 * 2,  
       col = "blue", length = 0.1)
text(cm_pca_loadings$PC1 * 2.2, cm_pca_loadings$PC2 * 2.2, labels = rownames(cm_pca_loadings), col = "blue", cex = 0.8)

### For Crepis (All Variables) ----
ca_pca_scores <- as.data.frame(cav$x)
c_groups <- AusPC$Treatment[AusPC$Species == "Crepis"]

#Base plot
plot(ca_pca_scores$PC1, ca_pca_scores$PC2,
     xlab = "PC1 (28%)", ylab = "PC2 (18.6%)", main = "Crepis - All Variables - PCA",
     col=as.factor(c_groups), pch = 19, cex = 1.5)

# Add spider diagram
ordispider(ca_pca_scores, c_groups, col = rainbow(length(unique(c_groups))), label = TRUE)
ordihull(ca_pca_scores, c_groups, col = rainbow(length(unique(c_groups))),draw = "polygon")

# Optional: Add legend
legend("topright", legend = unique(c_groups), 
       col = rainbow(length(unique(c_groups))),
       border = "black",
       pch = 19, bty = "o", title = "Treatment Groups")

# Adding arrows
ca_pca_loadings <- as.data.frame(cav$rotation[, 1:2])  
arrows(0, 0, 
       ca_pca_loadings$PC1 * 2,  
       ca_pca_loadings$PC2 * 2,  
       col = "blue", length = 0.1)
text(ca_pca_loadings$PC1 * 2.2, ca_pca_loadings$PC2 * 2.2, labels = rownames(ca_pca_loadings), col = "blue", cex = 0.8)

## Rank Abundance Curves ----
### Creating PCA for each treatment ----
### Lotus ----
lmv <- prcomp(AusPC[AusPC$Species == "Lotus", c(9,18,20:23)], scale = TRUE)
lmv_control <- prcomp(AusPC[AusPC$Species == "Lotus" & AusPC$Treatment == "C0T0D0", c(9,18,20:23)], scale = TRUE)
lmv_drought <- prcomp(AusPC[AusPC$Species == "Lotus" & AusPC$Treatment == "C0T0D1", c(9,18,20:23)], scale = TRUE)
lmv_temperature <- prcomp(AusPC[AusPC$Species == "Lotus" & AusPC$Treatment == "C0T2D0", c(9,18,20:23)], scale = TRUE)
lmv_co2 <- prcomp(AusPC[AusPC$Species == "Lotus" & AusPC$Treatment == "C2T0D0", c(9,18,20:23)], scale = TRUE)
lmv_ct <- prcomp(AusPC[AusPC$Species == "Lotus" & AusPC$Treatment == "C2T2D0", c(9,18,20:23)], scale = TRUE)
lmv_ctd <- prcomp(AusPC[AusPC$Species == "Lotus" & AusPC$Treatment == "C2T2D1", c(9,18,20:23)], scale = TRUE)
lmv_ctd
lmv_control
### Extracting Variance Explained for each treatment
variance_lmv_control <- lmv_control$sdev^2 / sum(lmv_control$sdev^2)
variance_lmv_drought <- lmv_drought$sdev^2 / sum(lmv_drought$sdev^2)
variance_lmv_temperature <- lmv_temperature$sdev^2 / sum(lmv_temperature$sdev^2)
variance_lmv_co2 <- lmv_co2$sdev^2 / sum(lmv_co2$sdev^2)
variance_lmv_ct <- lmv_ct$sdev^2 / sum(lmv_ct$sdev^2)
variance_lmv_ctd <- lmv_ctd$sdev^2 / sum(lmv_ctd$sdev^2)

### Combining Variance Explained Results to a dataframe
AusL_PC <- data.frame(
  PCA_axis = rep(1:6,6),
  Variance_Explained = c(variance_lmv_control, variance_lmv_co2, variance_lmv_temperature, variance_lmv_drought, variance_lmv_ct, variance_lmv_ctd),
  Treatment = rep(c("C0T0D0", "C2T0D0", "C0T2D0", "C0T0D2", "C2T2D0", "C2T2D1"), each = 6)
)

### Adding zero to treatments with insufficient PCA Axis 
if (length(variance_lmv_drought) < 6) {
  variance_lmv_drought <- c(variance_lmv_drought, rep(0, 6 - length(variance_lmv_drought)))
}
if (length(variance_lmv_ct) < 6) {
  variance_lmv_ct <- c(variance_lmv_ct, rep(0, 6 - length(variance_lmv_ct)))
}

ggplot(AusL_PC, aes(x = AusL_PC$PCA_axis, y = AusL_PC$Variance_Explained, color = AusL_PC$Treatment)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Lotus - Effect of Climatic Factors on Trait Covariation",
    x = "PCA Axis",
    y = "Proportion of Variance Explained",
    color = "Treatment"
  ) +
  theme_minimal()

View(AusL_PC)

### Crepis (Main Variables) ----
cmv <- prcomp(AusPC[AusPC$Species == "Crepis", c(9,18,20:23)], scale = TRUE)
cmv_control <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C0T0D0", c(9,18,20:23)], scale = TRUE)
cmv_drought <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C0T0D1", c(9,18,20:23)], scale = TRUE)
cmv_temperature <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C0T2D0", c(9,18,20:23)], scale = TRUE)
cmv_co2 <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C2T0D0", c(9,18,20:23)], scale = TRUE)
cmv_ct <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C2T2D0", c(9,18,20:23)], scale = TRUE)
cmv_ctd <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C2T2D1", c(9,18,20:23)], scale = TRUE)

### Extracting variance explained for each treatment 
variance_cmv_control <- cmv_control$sdev^2 / sum(cmv_control$sdev^2)
variance_cmv_drought <- cmv_drought$sdev^2 / sum(cmv_drought$sdev^2)
variance_cmv_temperature <- cmv_temperature$sdev^2 / sum(cmv_temperature$sdev^2)
variance_cmv_co2 <- cmv_co2$sdev^2 / sum(cmv_co2$sdev^2)
variance_cmv_ct <- cmv_ct$sdev^2 / sum(cmv_ct$sdev^2)
variance_cmv_ctd <- cmv_ctd$sdev^2 / sum(cmv_ctd$sdev^2)

lengths <- sapply(list(variance_cmv_control, variance_cmv_co2, variance_cmv_temperature, variance_cmv_drought, variance_cmv_ct, variance_cmv_ctd), length)
print(lengths)

### Combining variance explained results to a dataframe

AusCM_PC <- data.frame(
  PCA_axis = rep(1:6, 6),
  Variance_Explained = c(variance_cmv_control, variance_cmv_co2, variance_cmv_temperature, variance_cmv_drought, variance_cmv_ct, variance_cmv_ctd),
  Treatment = rep(c("C0T0D0", "C2T0D0", "C0T2D0", "C0T0D2", "C2T2D0", "C2T2D1"), each = 6)
)

ggplot(AusCM_PC, aes(x = AusCM_PC$PCA_axis, y = AusCM_PC$Variance_Explained, color = AusCM_PC$Treatment)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Crepis (Main Variables) - Effect of Climatic Factors on Trait Covariation",
    x = "PCA Axis",
    y = "Proportion of Variance Explained",
    color = "Treatment"
  ) +
  theme_minimal()

### Crepis (All Variables) ----
cav <- prcomp(AusPC[AusPC$Species == "Crepis", c(9,18,20:25)], scale = TRUE)

cav_control <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C0T0D0", c(9,18,20:25)], scale = TRUE)
cav_drought <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C0T0D1", c(9,18,20:25)], scale = TRUE)
cav_temperature <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C0T2D0", c(9,18,20:25)], scale = TRUE)
cav_co2 <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C2T0D0", c(9,18,20:25)], scale = TRUE)
cav_ct <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C2T2D0", c(9,18,20:25)], scale = TRUE)
cav_ctd <- prcomp(AusPC[AusPC$Species == "Crepis" & AusPC$Treatment == "C2T2D1", c(9,18,20:25)], scale = TRUE)

### Extracting variance explained for each treatment 
variance_cav_control <- cav_control$sdev^2 / sum(cav_control$sdev^2)
variance_cav_drought <- cav_drought$sdev^2 / sum(cav_drought$sdev^2)
variance_cav_temperature <- cav_temperature$sdev^2 / sum(cav_temperature$sdev^2)
variance_cav_co2 <- cav_co2$sdev^2 / sum(cav_co2$sdev^2)
variance_cav_ct <- cav_ct$sdev^2 / sum(cav_ct$sdev^2)
variance_cav_ctd <- cav_ctd$sdev^2 / sum(cav_ctd$sdev^2)

### Checking the lengths of variance vectors
lengths <- sapply(list(variance_cav_control, variance_cav_co2, variance_cav_temperature, variance_cav_drought, variance_cav_ct, variance_cav_ctd), length)
print(lengths)

### Adding zero to treatments with insufficient PCA Axis 
if (length(variance_cav_co2) < 8) {
  variance_cav_co2 <- c(variance_cav_co2, rep(0, 8 - length(variance_cav_co2)))
}

if (length(variance_cav_ct) < 8) {
  variance_cav_ct <- c(variance_cav_ct, rep(0, 8 - length(variance_cav_ct)))
}

if (length(variance_cav_ctd) < 8) {
  variance_cav_ctd <- c(variance_cav_ctd, rep(0, 8 - length(variance_cav_ctd)))
}

### Combining variance explained results to a dataframe

AusCA_PC <- data.frame(
  PCA_axis = rep(1:8, 6),
  Variance_Explained = c(variance_cav_control, variance_cav_co2, variance_cav_temperature, variance_cav_drought, variance_cav_ct, variance_cav_ctd),
  Treatment = rep(c("C0T0D0", "C2T0D0", "C0T2D0", "C0T0D2", "C2T2D0", "C2T2D1"), each = 8)
)

ggplot(AusCA_PC, aes(x = AusCA_PC$PCA_axis, y = AusCA_PC$Variance_Explained, color = AusCA_PC$Treatment)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Crepis (All Variables) - Effect of Climatic Factors on Trait Covariation",
    x = "PCA Axis",
    y = "Proportion of Variance Explained",
    color = "Treatment"
  ) +
  theme_minimal()

## Extracting PCA Results----
### For Lotus ----

lotus_pca_trait_loadings <- lmv$rotation*lmv$sdev
lotus_pca_eigenvalues <- lmv$sdev^2
lotus_pca_variance <- lotus_pca_eigenvalues/sum(lotus_pca_eigenvalues)
lotus_pca_cum_variance <- cumsum(lotus_pca_variance)

trait_loadings_df <- as.data.frame(lotus_pca_trait_loadings)
variance_percent <- lotus_pca_variance * 100
cumulative_variance_percent <- lotus_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_pca_combined.csv", row.names = TRUE)

### For Treatments
#### Control
lotus_control_pca_trait_loadings <- lmv_control$rotation
lotus_control_pca_eigenvalues <- lmv_control$sdev^2
lotus_control_pca_variance <- lotus_control_pca_eigenvalues / sum(lotus_control_pca_eigenvalues)
lotus_control_pca_cum_variance <- cumsum(lotus_control_pca_variance)

trait_loadings_df <- as.data.frame(lotus_control_pca_trait_loadings)
variance_percent <- lotus_control_pca_variance * 100
cumulative_variance_percent <- lotus_control_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_control_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_control_pca_combined.csv", row.names = TRUE)

#### CO2
lotus_co2_pca_trait_loadings <- lmv_co2$rotation
lotus_co2_pca_eigenvalues <- lmv_co2$sdev^2
lotus_co2_pca_variance <- lotus_co2_pca_eigenvalues / sum(lotus_co2_pca_eigenvalues)
lotus_co2_pca_cum_variance <- cumsum(lotus_co2_pca_variance)

trait_loadings_df <- as.data.frame(lotus_co2_pca_trait_loadings)
variance_percent <- lotus_co2_pca_variance * 100
cumulative_variance_percent <- lotus_co2_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_co2_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_co2_pca_combined.csv", row.names = TRUE)

#### Temperature
lotus_temperature_pca_trait_loadings <- lmv_temperature$rotation
lotus_temperature_pca_eigenvalues <- lmv_temperature$sdev^2
lotus_temperature_pca_variance <- lotus_temperature_pca_eigenvalues / sum(lotus_temperature_pca_eigenvalues)
lotus_temperature_pca_cum_variance <- cumsum(lotus_temperature_pca_variance)

trait_loadings_df <- as.data.frame(lotus_temperature_pca_trait_loadings)
variance_percent <- lotus_temperature_pca_variance * 100
cumulative_variance_percent <- lotus_temperature_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_temperature_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_temperature_pca_combined.csv", row.names = TRUE)

#### Drought
lotus_drought_pca_trait_loadings <- lmv_drought$rotation
lotus_drought_pca_eigenvalues <- lmv_drought$sdev^2
lotus_drought_pca_variance <- lotus_drought_pca_eigenvalues / sum(lotus_drought_pca_eigenvalues)
lotus_drought_pca_cum_variance <- cumsum(lotus_drought_pca_variance)

trait_loadings_df <- as.data.frame(lotus_drought_pca_trait_loadings)
variance_percent <- lotus_drought_pca_variance * 100
cumulative_variance_percent <- lotus_drought_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_drought_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_drought_pca_combined.csv", row.names = TRUE)

#### CO2 and Temperature
lotus_ct_pca_trait_loadings <- lmv_ct$rotation
lotus_ct_pca_eigenvalues <- lmv_ct$sdev^2
lotus_ct_pca_variance <- lotus_ct_pca_eigenvalues / sum(lotus_ct_pca_eigenvalues)
lotus_ct_pca_cum_variance <- cumsum(lotus_ct_pca_variance)

trait_loadings_df <- as.data.frame(lotus_ct_pca_trait_loadings)
variance_percent <- lotus_ct_pca_variance * 100
cumulative_variance_percent <- lotus_ct_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_ct_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_ct_pca_combined.csv", row.names = TRUE)

#### Combined Interaction
lotus_ctd_pca_trait_loadings <- lmv_ctd$rotation
lotus_ctd_pca_eigenvalues <- lmv_ctd$sdev^2
lotus_ctd_pca_variance <- lotus_ctd_pca_eigenvalues / sum(lotus_ctd_pca_eigenvalues)
lotus_ctd_pca_cum_variance <- cumsum(lotus_ctd_pca_variance)

trait_loadings_df <- as.data.frame(lotus_ctd_pca_trait_loadings)
variance_percent <- lotus_ctd_pca_variance * 100
cumulative_variance_percent <- lotus_ctd_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = lotus_ctd_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "lotus_ctd_pca_combined.csv", row.names = TRUE)


### For Crepis (Main Variables) ----

crepis_m_pca_trait_loadings <- cmv$rotation * cmv$sdev
crepis_m_pca_eigenvalues <- cmv$sdev^2
crepis_m_pca_variance <- crepis_m_pca_eigenvalues / sum(crepis_m_pca_eigenvalues)
crepis_m_pca_cum_variance <- cumsum(crepis_m_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_pca_trait_loadings)
variance_percent <- crepis_m_pca_variance * 100
cumulative_variance_percent <- crepis_m_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_pca_combined.csv", row.names = TRUE)

### For Treatments
#### Control
crepis_m_control_pca_trait_loadings <- cmv_control$rotation
crepis_m_control_pca_eigenvalues <- cmv_control$sdev^2
crepis_m_control_pca_variance <- crepis_m_control_pca_eigenvalues / sum(crepis_m_control_pca_eigenvalues)
crepis_m_control_pca_cum_variance <- cumsum(crepis_m_control_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_control_pca_trait_loadings)
variance_percent <- crepis_m_control_pca_variance * 100
cumulative_variance_percent <- crepis_m_control_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_control_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_control_pca_combined.csv", row.names = TRUE)

#### CO2
crepis_m_co2_pca_trait_loadings <- cmv_co2$rotation
crepis_m_co2_pca_eigenvalues <- cmv_co2$sdev^2
crepis_m_co2_pca_variance <- crepis_m_co2_pca_eigenvalues / sum(crepis_m_co2_pca_eigenvalues)
crepis_m_co2_pca_cum_variance <- cumsum(crepis_m_co2_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_co2_pca_trait_loadings)
variance_percent <- crepis_m_co2_pca_variance * 100
cumulative_variance_percent <- crepis_m_co2_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_co2_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_co2_pca_combined.csv", row.names = TRUE)

#### Temperature
crepis_m_temperature_pca_trait_loadings <- cmv_temperature$rotation
crepis_m_temperature_pca_eigenvalues <- cmv_temperature$sdev^2
crepis_m_temperature_pca_variance <- crepis_m_temperature_pca_eigenvalues / sum(crepis_m_temperature_pca_eigenvalues)
crepis_m_temperature_pca_cum_variance <- cumsum(crepis_m_temperature_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_temperature_pca_trait_loadings)
variance_percent <- crepis_m_temperature_pca_variance * 100
cumulative_variance_percent <- crepis_m_temperature_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_temperature_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_temperature_pca_combined.csv", row.names = TRUE)

#### Drought
crepis_m_drought_pca_trait_loadings <- cmv_drought$rotation
crepis_m_drought_pca_eigenvalues <- cmv_drought$sdev^2
crepis_m_drought_pca_variance <- crepis_m_drought_pca_eigenvalues / sum(crepis_m_drought_pca_eigenvalues)
crepis_m_drought_pca_cum_variance <- cumsum(crepis_m_drought_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_drought_pca_trait_loadings)
variance_percent <- crepis_m_drought_pca_variance * 100
cumulative_variance_percent <- crepis_m_drought_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_drought_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_drought_pca_combined.csv", row.names = TRUE)

#### CO2 and Temperature
crepis_m_ct_pca_trait_loadings <- cmv_ct$rotation
crepis_m_ct_pca_eigenvalues <- cmv_ct$sdev^2
crepis_m_ct_pca_variance <- crepis_m_ct_pca_eigenvalues / sum(crepis_m_ct_pca_eigenvalues)
crepis_m_ct_pca_cum_variance <- cumsum(crepis_m_ct_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_ct_pca_trait_loadings)
variance_percent <- crepis_m_ct_pca_variance * 100
cumulative_variance_percent <- crepis_m_ct_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_ct_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_ct_pca_combined.csv", row.names = TRUE)

#### Combined Interaction
crepis_m_ctd_pca_trait_loadings <- cmv_ctd$rotation
crepis_m_ctd_pca_eigenvalues <- cmv_ctd$sdev^2
crepis_m_ctd_pca_variance <- crepis_m_ctd_pca_eigenvalues / sum(crepis_m_ctd_pca_eigenvalues)
crepis_m_ctd_pca_cum_variance <- cumsum(crepis_m_ctd_pca_variance)

trait_loadings_df <- as.data.frame(crepis_m_ctd_pca_trait_loadings)
variance_percent <- crepis_m_ctd_pca_variance * 100
cumulative_variance_percent <- crepis_m_ctd_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_m_ctd_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_m_ctd_pca_combined.csv", row.names = TRUE)

### For Crepis (All Variables) ----

crepis_a_pca_trait_loadings <- cav$rotation * cav$sdev
crepis_a_pca_eigenvalues <- cav$sdev^2
crepis_a_pca_variance <- crepis_a_pca_eigenvalues / sum(crepis_a_pca_eigenvalues)
crepis_a_pca_cum_variance <- cumsum(crepis_a_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_pca_trait_loadings)
variance_percent <- crepis_a_pca_variance * 100
cumulative_variance_percent <- crepis_a_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_pca_combined.csv", row.names = TRUE)

### For Treatments
#### Control
crepis_a_control_pca_trait_loadings <- cav_control$rotation
crepis_a_control_pca_eigenvalues <- cav_control$sdev^2
crepis_a_control_pca_variance <- crepis_a_control_pca_eigenvalues / sum(crepis_a_control_pca_eigenvalues)
crepis_a_control_pca_cum_variance <- cumsum(crepis_a_control_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_control_pca_trait_loadings)
variance_percent <- crepis_a_control_pca_variance * 100
cumulative_variance_percent <- crepis_a_control_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_control_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_control_pca_combined.csv", row.names = TRUE)

#### CO2
crepis_a_co2_pca_trait_loadings <- cav_co2$rotation
crepis_a_co2_pca_eigenvalues <- cav_co2$sdev^2
crepis_a_co2_pca_variance <- crepis_a_co2_pca_eigenvalues / sum(crepis_a_co2_pca_eigenvalues)
crepis_a_co2_pca_cum_variance <- cumsum(crepis_a_co2_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_co2_pca_trait_loadings)
variance_percent <- crepis_a_co2_pca_variance * 100
cumulative_variance_percent <- crepis_a_co2_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_co2_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_co2_pca_combined.csv", row.names = TRUE)

#### Temperature
crepis_a_temperature_pca_trait_loadings <- cav_temperature$rotation
crepis_a_temperature_pca_eigenvalues <- cav_temperature$sdev^2
crepis_a_temperature_pca_variance <- crepis_a_temperature_pca_eigenvalues / sum(crepis_a_temperature_pca_eigenvalues)
crepis_a_temperature_pca_cum_variance <- cumsum(crepis_a_temperature_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_temperature_pca_trait_loadings)
variance_percent <- crepis_a_temperature_pca_variance * 100
cumulative_variance_percent <- crepis_a_temperature_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_temperature_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_temperature_pca_combined.csv", row.names = TRUE)

#### Drought
crepis_a_drought_pca_trait_loadings <- cav_drought$rotation
crepis_a_drought_pca_eigenvalues <- cav_drought$sdev^2
crepis_a_drought_pca_variance <- crepis_a_drought_pca_eigenvalues / sum(crepis_a_drought_pca_eigenvalues)
crepis_a_drought_pca_cum_variance <- cumsum(crepis_a_drought_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_drought_pca_trait_loadings)
variance_percent <- crepis_a_drought_pca_variance * 100
cumulative_variance_percent <- crepis_a_drought_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_drought_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_drought_pca_combined.csv", row.names = TRUE)

#### CO2 and Temperature
crepis_a_ct_pca_trait_loadings <- cav_ct$rotation
crepis_a_ct_pca_eigenvalues <- cav_ct$sdev^2
crepis_a_ct_pca_variance <- crepis_a_ct_pca_eigenvalues / sum(crepis_a_ct_pca_eigenvalues)
crepis_a_ct_pca_cum_variance <- cumsum(crepis_a_ct_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_ct_pca_trait_loadings)
variance_percent <- crepis_a_ct_pca_variance * 100
cumulative_variance_percent <- crepis_a_ct_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_ct_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_ct_pca_combined.csv", row.names = TRUE)

#### Combined Interaction
crepis_a_ctd_pca_trait_loadings <- cav_ctd$rotation
crepis_a_ctd_pca_eigenvalues <- cav_ctd$sdev^2
crepis_a_ctd_pca_variance <- crepis_a_ctd_pca_eigenvalues / sum(crepis_a_ctd_pca_eigenvalues)
crepis_a_ctd_pca_cum_variance <- cumsum(crepis_a_ctd_pca_variance)

trait_loadings_df <- as.data.frame(crepis_a_ctd_pca_trait_loadings)
variance_percent <- crepis_a_ctd_pca_variance * 100
cumulative_variance_percent <- crepis_a_ctd_pca_cum_variance * 100
combined_df <- rbind(
  trait_loadings_df,
  Eigenvalues = crepis_a_ctd_pca_eigenvalues,
  Variance_Percent = variance_percent,
  Cumulative_Variance_Percent = cumulative_variance_percent
)
write.csv(combined_df, "crepis_a_ctd_pca_combined.csv", row.names = TRUE)



## PCA Plots for Each Treatment ----
### For Lotus ----

ggbiplot(lmv_control, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T0D0", x = "PC1 (42.4%)", y = "PC2 (22.8%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Center and enlarge the title
  )

ggbiplot(lmv_co2, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T0D0", x = "PC1 (45.4%)", y = "PC2 (31.8%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Center and enlarge title
  )

ggbiplot(lmv_temperature, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T2D0", x = "PC1 (44.2%)", y = "PC2 (29.7%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(lmv_drought, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T0D1", x = "PC1 (39.5%)", y = "PC2 (34.9%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(lmv_ct, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T2D0", x = "PC1 (44.2%)", y = "PC2 (31.2%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(lmv_ctd, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T2D1", x = "PC1 (66.8%)", y = "PC2 (24.7%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

### For Crepis (Main Variables) ----
ggbiplot(cmv_control, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T0D0", x = "PC1 (32.1%)", y = "PC2 (30.6%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Center and enlarge the title
  )

ggbiplot(cmv_co2, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T0D0", x = "PC1 (51.2%)", y = "PC2 (41.5%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Center and enlarge title
  )

ggbiplot(cmv_temperature, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T2D0", x = "PC1 (35.5%)", y = "PC2 (27.8%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(cmv_drought, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T0D1", x = "PC1 (42.7%)", y = "PC2 (25.3%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(cmv_ct, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T2D0", x = "PC1 (69.7%)", y = "PC2 (15.7%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(cmv_ctd, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T2D1", x = "PC1 (49.2%)", y = "PC2 (28.7%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

### For Crepis (All Variables) ----
ggbiplot(cav_control, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T0D0", x = "PC1 (33.4%)", y = "PC2 (23.8%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Center and enlarge the title
  )

ggbiplot(cav_co2, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T0D0", x = "PC1 (52.0%)", y = "PC2 (36.4%)") +
  theme(
    axis.title.x = element_text(size = 16),  # Increase size of x-axis title
    axis.title.y = element_text(size = 16),  # Increase size of y-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Center and enlarge title
  )

ggbiplot(cav_temperature, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T2D0", x = "PC1 (31.5%)", y = "PC2 (30.3%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(cav_drought, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C0T0D1", x = "PC1 (41.6%)", y = "PC2 (27.6%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(cav_ct, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T2D0", x = "PC1 (53.9%)", y = "PC2 (20.6%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

ggbiplot(cav_ctd, varname.size = 6, alpha = 0) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(title = "C2T2D1", x = "PC1 (42.9%)", y = "PC2 (24.3%)") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

# Redundancy Analysis ----
### Sub-setting Data frames for Lotus and Crepis ----
lotus_rows <- AusRD[AusRD$Species == "Lotus", ]
names(lotus_rows)
lotus_rows <- lotus_rows[-c(40:72), ]
AusRDL <- lotus_rows[, -c(24,25)]
View(crepis_rows)
names(AusRDL)
crepis_rows <- AusRD[AusRD$Species == "Crepis",]
crepis_rows <- crepis_rows[-c(46:78),]
AusRDC <- crepis_rows[,-c(10:17)]
## Redundancy Analysis for Lotus ----
names(AusRDL)[names(AusRDL) == "Specific Petal Area (SPA)(cm2/g)"] <- "F-SPA"
names(AusRDL)[names(AusRDL) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "F-PDMC"
names(AusRDL)[names(AusRDL) == "Specific Leaf Area (SLA)(cm2/g)"] <- "L-SLA"
names(AusRDL)[names(AusRDL) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "L-LDMC"
names(AusRDL)[names(AusRDL) == "Display Area (DA)(cm2)"] <- "F-DA"
names(AusRDL)[names(AusRDL) == "Leaf Area (LA)(cm2)"] <- "L-LA"
names(AusRDL)[names(AusRDL) == "CO2"] <- "C"
names(AusRDL)[names(AusRDL) == "Temperature"] <- "T"
names(AusRDL)[names(AusRDL) == "Drought"] <- "D"


rda.l <- rda(AusRDL[,c(9,18,20:23)] ~ C*T*D, data = AusRDL, scale = TRUE)
Lotus_RDA <- anova.cca(rda.l, permutations = 9999, by = "terms")
Lotus_RDA <- as.data.frame(Lotus_RDA)
View(Lotus_RDA) # CO2, Temperature, Drought and CO2:Drought are significant
anova(rda.l)
summary(rda.l)
names(AusRDL)
### Extracting Values for Lotus RDA ----
summary_text <- capture.output(summary(rda.l))
writeLines(summary_text, "lotus_rda_summary.txt")
l_rda_eigenvalues <- eigenvals(rda.l, constrained = TRUE)
l_rda_proportion_explained <- l_rda_eigenvalues/sum(l_rda_eigenvalues) 
l_rda_proportion_explained

summary(rda.ln)

# Changing Name of Interaction ----

bp_scores <- scores(rda.ln, display = "bp")
rownames(bp_scores)
# Rename the fourth row
rownames(bp_scores)[4] <- "(C+T)xD"

# View the updated biplot scores
rownames(bp_scores)

rownames(scores(rda.ln, display = "bp"))[4] <- "(CT):D"

### Plot for Lotus RDA ----

rda.ln <- rda(AusRDL[,c(9,18,20:23)] ~ C + T + D + C:D, data = AusRDL, scale = TRUE)
plot(rda.ln, type = "n",xlim = c(-1,1), ylim = c(-1,1), xlab = "RDA1 (60.15%)", ylab = "RDA2 (20.84%)")
site_scores <- scores(rda.ln, display = "sites")
plot(rda.ln, display = "both", cex = 0.7)
points(site_scores, pch = 16, col = "red", cex = 0.5) # Ignore this to not get points
### Adding Arrows and Text for Response Variables
arrows(0, 0, scores(rda.ln, display = "species")[,1], scores(rda.ln, display = "species")[,2], col = 'blue', length = 0.1)
text(scores(rda.ln, display = "species")[,1], scores(rda.ln, display = "species")[,2], labels = rownames(scores(rda.ln, display = "species")), col = 'blue', pos = 3, cex = 1)
### Adding Arrows and Text for Predictor Variables
arrows(0, 0, scores(rda.ln, display = "bp")[,1], scores(rda.ln, display = "bp")[,2], col = 'red', length = 0.1)
text(scores(rda.ln, display = "bp")[,1], scores(rda.ln, display = "bp")[,2], labels = rownames(bp_scores), col = 'red', pos = 3, cex = 1)

## Redundancy Analysis for Crepis ----
names(AusRDC)[names(AusRDC) == "Specific Petal Area (SPA)(cm2/g)"] <- "F-SPA"
names(AusRDC)[names(AusRDC) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "F-PDMC"
names(AusRDC)[names(AusRDC) == "Specific Leaf Area (SLA)(cm2/g)"] <- "L-SLA"
names(AusRDC)[names(AusRDC) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "L-LDMC"
names(AusRDC)[names(AusRDC) == "Display Area (DA)(cm2)"] <- "F-DA"
names(AusRDC)[names(AusRDC) == "Leaf Area (LA)(cm2)"] <- "L-LA"
names(AusRDC)[names(AusRDC) == "Number of Seeds (SN)"] <- "S-SN"
names(AusRDC)[names(AusRDC) == "Seed Mass (SM)(g)"] <- "S-SM"
names(AusRDC)[names(AusRDC) == "CO2"] <- "C"
names(AusRDC)[names(AusRDC) == "Temperature"] <- "T"
names(AusRDC)[names(AusRDC) == "Drought"] <- "D"

rda.c <- rda(AusRDC[,c(9:10,12:17)] ~ C*T*D, data = AusRDC, scale = TRUE)
Crepis_RDA <- anova.cca(rda.c, permutations = 9999, by = "terms")
Crepis_RDA # Drought and CO2: Temperature are significant
Crepis_RDA <- as.data.frame(Crepis_RDA)
View(Crepis_RDA)
summary(rda.c)
plot(AusC$`Leaf Area (LA)(cm2)`~AusC$Temperature)
### Extracting Values for Crepis RDA ----
summary_text <- capture.output(summary(rda.c))
writeLines(summary_text, "crepis_rda_summary.txt")
c_rda_eigenvalues <- eigenvals(rda.c, constrained = TRUE)
c_rda_proportion_explained <- c_rda_eigenvalues/sum(c_rda_eigenvalues) 
c_rda_proportion_explained

### Plot for Crepis RDA ----
rda.cn <- rda(AusRDC[,c(9:10,12:17)] ~ D + C:T + Condition(C+T), data = AusRDC, scale = TRUE)

plot(rda.cn, type = "n", xlim = c(-1,1), ylim = c(-0.5,0.5), xlab = "RDA1 (65.15%)", ylab = "RDA2 (20.54%)")
site_scores <- scores(rda.cn, display = "sites")
plot(rda.cn, display = "both", cex = 0.7)
points(site_scores, pch = 16, col = "red", cex = 0.5) # Ignore this to not get row names 
### Adding Arrows and Text for Response Variables
arrows(0, 0, scores(rda.cn, display = "species")[,1], scores(rda.cn, display = "species")[,2], col = 'blue', length = 0.1)
text(scores(rda.cn, display = "species")[,1], scores(rda.cn, display = "species")[,2], labels = rownames(scores(rda.cn, display = "species")), col = 'blue', pos = 3, cex = 1)
### Adding Arrows and Text for Predictor Variables
arrows(0, 0, scores(rda.cn, display = "bp")[,1], scores(rda.cn, display = "bp")[,2], col = 'red', length = 0.1)
text(scores(rda.cn, display = "bp")[,1], scores(rda.cn, display = "bp")[,2], labels = rownames(scores(rda.cn, display = "bp")), col = 'red', pos = 3, cex = 1)

rda.cn <- rda(AusRDC[,c(18:23)] ~ Drought + Condition(AusRDC$`Temperature Level`*AusRDC$`CO2 Level`), data = AusRDC, scale = TRUE)
View(AusRD)
## Redundancy Analysis for both Crepis and Lotus ----

rda.1 <- rda(AusRD[,c(18,20:22)] ~ CO2*Temperature*Drought, data = AusRD, scale = TRUE)
View(AusRD)
ordiplot(rda.1, type = "text", scale = TRUE)
anova.cca(rda.1, permutations = 9999)
anova.cca(rda.1, permutations = 9999, by = "axis")
anova.cca(rda.1, permutations = 9999, by = "terms")
Combined_RDA <- anova.cca(rda.1, permutations = 9999, by = "terms")
Combined_RDA # Only Drought is Significant
writeLines(capture.output(Combined_RDA), "Combined_rda_anova.txt")
Combined_RDA$F
Combined_RDA$`Pr(>F)`

summary_rda <- summary(rda.1)
summary_rda


# Preparing function which exports anova results to CSV ----
anova_to_csv <- function(model, csv_file_path) {
  # Perform the ANOVA test
  anova_results <- anova(model)
  
  # Convert the results to a data frame
  anova_df <- as.data.frame(anova_results)
  
  # Add a column with the model formula to distinguish different models
  model_formula <- paste(deparse(formula(model)), collapse = " ")
  anova_df$model <- rep(model_formula, nrow(anova_df))
  
  # Check if the CSV file already exists
  if (file.exists(csv_file_path)) {
    # Read the existing data
    existing_df <- read_csv(csv_file_path)
    
    # Append the new results
    combined_df <- bind_rows(existing_df, anova_df)
  } else {
    # If the file does not exist, the combined data frame is the new ANOVA results
    combined_df <- anova_df
  }
  
  # Write the combined DataFrame to the CSV file
  write_csv(combined_df, file = csv_file_path)
  
  cat("ANOVA results have been written to", csv_file_path, "\n")
}

## Output ANOVA results to CSV file for Mixed Effect Linear models in Lotus ----
anova_output_csv <- 'anova_LMLo - 18 August 2024.csv'
anova_to_csv(l_DAmodel, anova_output_csv)
anova_to_csv(log_SSPAmodel, anova_output_csv)
anova_to_csv(log_SW1PAmodel, anova_output_csv)
anova_to_csv(SW2PAmodel, anova_output_csv)
anova_to_csv(sq_SKPAmodel, anova_output_csv)
anova_to_csv(l_LAmodel, anova_output_csv)
anova_to_csv(l_sq_SPAmodel, anova_output_csv)
anova_to_csv(l_PDMCmodel, anova_output_csv)
anova_to_csv(l_SLAmodel, anova_output_csv)
anova_to_csv(l_log_LDMCmodel, anova_output_csv)

## Output ANOVA results to CSV file for Mixed Effect Linear models in Crepis ----
anova_output_csv <- 'anova_LMCr - 18 August 2024.csv'
anova_to_csv(c_sq_DAmodel, anova_output_csv)
anova_to_csv(c_sq_SPAmodel, anova_output_csv)
anova_to_csv(c_sq_PDMCmodel, anova_output_csv)
anova_to_csv(c_log_LAmodel, anova_output_csv)
anova_to_csv(c_log_SLAmodel, anova_output_csv)
anova_to_csv(c_sq_LDMCmodel, anova_output_csv)
anova_to_csv(sq_SNmodel, anova_output_csv)
anova_to_csv(log_SMmodel, anova_output_csv)

## Output ANOVA results to CSV file for Linear models in Lotus ----
anova_output_csv <- 'anova_LLo - 18 August 2024.csv'
anova_to_csv(l_DAmodel, anova_output_csv)
anova_to_csv(log_SSPAmodel, anova_output_csv)
anova_to_csv(log_SW1PAmodel, anova_output_csv)
anova_to_csv(SW2PAmodel, anova_output_csv)
anova_to_csv(sq_SKPAmodel, anova_output_csv)
anova_to_csv(l_sq_SPAmodel, anova_output_csv)
anova_to_csv(l_LAmodel, anova_output_csv)
anova_to_csv(l_PDMCmodel, anova_output_csv)
anova_to_csv(l_SLAmodel, anova_output_csv)
anova_to_csv(l_log_LDMCmodel, anova_output_csv)

## Output ANOVA results to CSV file for Linear models in Crepis ----
anova_output_csv <- 'anova_LCr - 18 August 2024.csv'
anova_to_csv(c_sq_DAmodel, anova_output_csv)
anova_to_csv(c_sq_SPAmodel, anova_output_csv)
anova_to_csv(c_sq_PDMCmodel, anova_output_csv)
anova_to_csv(c_log_LAmodel, anova_output_csv)
anova_to_csv(c_log_SLAmodel, anova_output_csv)
anova_to_csv(c_sq_LDMCmodel, anova_output_csv)
anova_to_csv(sq_SNmodel, anova_output_csv)
anova_to_csv(log_SMmodel, anova_output_csv)


# Correlogram ----
## For Lotus ----
View(AusL)
Lotcor_Data <- AusL[,c(9,18,20:23)]
names(Lotcor_Data)
names(Lotcor_Data)[names(Lotcor_Data) == "Display Area (DA)(cm2)"] <- "DA"
names(Lotcor_Data)[names(Lotcor_Data) == "Specific Petal Area (SPA)(cm2/g)"] <- "SPA"
names(Lotcor_Data)[names(Lotcor_Data) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "PDMC"
names(Lotcor_Data)[names(Lotcor_Data) == "Specific Leaf Area (SLA)(cm2/g)"] <- "SLA"
names(Lotcor_Data)[names(Lotcor_Data) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "LDMC"
names(Lotcor_Data)[names(Lotcor_Data) == "Leaf Area (LA)(cm2)"] <- "LA"
Lotcor_Matrix <- cor(Lotcor_Data, use = "complete.obs", method = "pearson")
corrplot(Lotcor_Matrix, method = "color", type = "full", 
         addCoef.col = "black", # Add black text with correlation values
         tl.col = "black", tl.srt = 45, # Rotate labels and set their color
         number.cex = 0.8) 
corrplot.mixed(Lotcor_Matrix, order = 'AOE')

### Correlogram which Includes p-values ----
l_corr_matrix <- rcorr(as.matrix(Lotcor_Data))
l_corr_matrix
l_r_values <- l_corr_matrix$r
l_p_values <- l_corr_matrix$P
l_r_values[upper.tri(l_r_values)] <- NA
l_p_values[lower.tri(l_p_values)] <- NA
l_r_melt <- melt(l_r_values, na.rm = TRUE)
l_p_melt <- melt(l_p_values, na.rm = TRUE)
l_merged_data <- merge(l_r_melt, l_p_melt, by = c("Var1", "Var2"), all = TRUE)
colnames(l_merged_data) <- c("Var1", "Var2", "correlation", "p_value")
l_merged_data$label <- with(l_merged_data, 
                            ifelse(is.na(correlation), 
                                   ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), 
                                   round(correlation, 2)))
l_merged_data$fontface <- with(l_merged_data, 
                               ifelse(is.na(correlation) & p_value <= 0.05, "bold", "plain"))
ggplot(data = l_merged_data, aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue4", mid = "#ffffff", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  geom_text(aes(label = label, fontface = fontface), color = "black", size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust size as needed
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed()

### Correlogram which includes p-values and the central diagonal is white ----

# Calculate correlation matrix
l_corr_matrix <- rcorr(as.matrix(Lotcor_Data))
l_r_values <- l_corr_matrix$r
l_p_values <- l_corr_matrix$P

# Set upper triangle and diagonal elements to NA for r values
l_r_values[upper.tri(l_r_values)] <- NA
diag(l_r_values) <- NA

# Set lower triangle elements to NA for p values
l_p_values[lower.tri(l_p_values)] <- NA

# Melt data into long format
l_r_melt <- melt(l_r_values, na.rm = TRUE)
l_p_melt <- melt(l_p_values, na.rm = TRUE)

# Merge correlation and p-value data
l_merged_data <- merge(l_r_melt, l_p_melt, by = c("Var1", "Var2"), all = TRUE)
colnames(l_merged_data) <- c("Var1", "Var2", "correlation", "p_value")

# Create labels and fontface based on correlation and p-values
l_merged_data$label <- with(l_merged_data, 
                            ifelse(is.na(correlation), 
                                   ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), 
                                   round(correlation, 2)))
l_merged_data$fontface <- with(l_merged_data, 
                               ifelse(is.na(correlation) & p_value <= 0.05, "bold", "plain"))

# Plot correlogram
ggplot(data = l_merged_data, aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue4", mid = "#ffffff", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  geom_text(aes(label = label, fontface = fontface), color = "black", size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed()

## For Crepis ----
names(AusC)
Crecor_Data <- AusC[,c(9,10,12:17)]
names(Crecor_Data)
names(Crecor_Data)[names(Crecor_Data) == "Display Area (DA)(cm2)"] <- "DA"
names(Crecor_Data)[names(Crecor_Data) == "Specific Petal Area (SPA)(cm2/g)"] <- "SPA"
names(Crecor_Data)[names(Crecor_Data) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "PDMC"
names(Crecor_Data)[names(Crecor_Data) == "Specific Leaf Area (SLA)(cm2/g)"] <-"SLA"
names(Crecor_Data)[names(Crecor_Data) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "LDMC"
names(Crecor_Data)[names(Crecor_Data) == "Number of Seeds (SN)"] <- "SN"
names(Crecor_Data)[names(Crecor_Data) == "Seed Mass (SM)(g)"] <- "SM"
names(Crecor_Data)[names(Crecor_Data) == "Leaf Area (LA)(cm2)"] <- "LA"
Crecor_Matrix <- cor(Crecor_Data, use = "complete.obs", method = "pearson")
corrplot(Crecor_Matrix, method = "color", type = "full", 
         addCoef.col = "black", # Add black text with correlation values
         tl.col = "black", tl.srt = 45, # Rotate labels and set their color
         number.cex = 0.8) 
corrplot.mixed(Crecor_Matrix, order = 'AOE')

### Correlogram which Includes p-values ----
c_corr_matrix <- rcorr(as.matrix(Crecor_Data))
c_corr_matrix
c_r_values <- c_corr_matrix$r
c_p_values <- c_corr_matrix$P
c_r_values[upper.tri(c_r_values)] <- NA
c_p_values[lower.tri(c_p_values)] <- NA
c_r_melt <- melt(c_r_values, na.rm = TRUE)
c_p_melt <- melt(c_p_values, na.rm = TRUE)
c_merged_data <- merge(c_r_melt, c_p_melt, by = c("Var1", "Var2"), all = TRUE)
colnames(c_merged_data) <- c("Var1", "Var2", "correlation", "p_value")
c_merged_data$label <- with(c_merged_data, 
                            ifelse(is.na(correlation), 
                                   ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), 
                                   round(correlation, 2)))
c_merged_data$fontface <- with(c_merged_data, 
                               ifelse(is.na(correlation) & p_value <= 0.05, "bold", "plain"))


ggplot(data = c_merged_data, aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue4", mid = "#ffffff", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation",) +
  geom_text(aes(label = label, fontface = fontface), color = "black", size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed()

### Correlogram which includes p-values and the central diagonal is white ----

c_corr_matrix <- rcorr(as.matrix(Crecor_Data))
c_r_values <- c_corr_matrix$r
c_p_values <- c_corr_matrix$P

# Set upper triangle and diagonal elements to NA for r values
c_r_values[upper.tri(c_r_values)] <- NA
diag(c_r_values) <- NA

# Set lower triangle elements to NA for p values
c_p_values[lower.tri(c_p_values)] <- NA

# Melt data into long format
c_r_melt <- melt(c_r_values, na.rm = TRUE)
c_p_melt <- melt(c_p_values, na.rm = TRUE)

# Merge correlation and p-value data
c_merged_data <- merge(c_r_melt, c_p_melt, by = c("Var1", "Var2"), all = TRUE)
colnames(c_merged_data) <- c("Var1", "Var2", "correlation", "p_value")

# Create labels and fontface based on correlation and p-values
c_merged_data$label <- with(c_merged_data, 
                            ifelse(is.na(correlation), 
                                   ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)), 
                                   round(correlation, 2)))
c_merged_data$fontface <- with(c_merged_data, 
                               ifelse(is.na(correlation) & p_value <= 0.05, "bold", "plain"))

# Plot correlogram
ggplot(data = c_merged_data, aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue4", mid = "#ffffff", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  geom_text(aes(label = label, fontface = fontface), color = "black", size = 5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed()

# Variance Partitioning using Barplot ----

## Creation of Separate Models ----
### For Lotus ----
lv_DAmodel <- lm(AusL$`Display Area (DA)(cm2)`~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusL)
lv_sq_SPAmodel <- lm(AusL$`sqrt_Specific Petal Area (SPA)(cm2/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusL)
lv_PDMCmodel <- lm(AusL$`Petal Dry Matter Content (PDMC)(g/g)` ~CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusL)
lv_LAmodel <- lm(AusL$`Leaf Area (LA)(cm2)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusL)
lv_SLAmodel <- lm(AusL$`Specific Leaf Area (SLA)(cm2/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusL)
lv_log_LDMCmodel <- lm(AusL$`log_Leaf Dry Matter Content (LDMC)(g/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusL)

### For Crepis ----
cv_sq_DAmodel <- lm(AusC$`Display Area (DA)(cm2)`~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_sq_SPAmodel <- lm(AusC$`Specific Petal Area (SPA)(cm2/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_sq_PDMCmodel <- lm(AusC$`Petal Dry Matter Content (PDMC)(g/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_log_LAmodel <- lm(AusC$`log_Leaf Area (LA)(cm2)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_log_SLAmodel <- lm(AusC$`log_Specific Leaf Area (SLA)(cm2/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_sq_LDMCmodel <- lm(AusC$`sqrt_Leaf Dry Matter Content (LDMC)(g/g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_sq_SNmodel <- lm(AusC$`sqrt_Number of Seeds (SN)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
cv_log_SMmodel <- lm(AusC$`log_Seed Mass (SM)(g)` ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusC)
anova(cv_sq_DAmodel)

### Creation of Plot for Lotus ----
l_model <- list(
  DA = lv_DAmodel,
  SPA = lv_sq_SPAmodel,
  PDMC = lv_PDMCmodel,
  LA = lv_LAmodel,
  SLA = lv_SLAmodel,
  LDMC = lv_log_LDMCmodel
)

## For both Significant and Marginally Significant Variables ----
all_results_df <- data.frame()

# Loop through each model in the l_model list
for (model_name in names(l_model)) {
  model <- l_model[[model_name]]
  
  # Perform ANOVA using the Anova function from the car package
  anova_results <- anova(model)
  
  # Convert ANOVA results to a data frame
  anova_results_df <- as.data.frame(anova_results)
  
  # Calculate the proportion of variance explained for each factor
  anova_results_df <- anova_results_df %>%
    mutate(PropVar = `Sum Sq` / sum(`Sum Sq`),
           Factor = rownames(anova_results_df),
           Trait = model_name) %>% # Add trait name
    filter(`Pr(>F)` < 0.1)%>%
    filter(Factor != "(Intercept)" & Factor != "Residuals") # Keep only significant variables # Exclude intercept and residuals
  
  # Combine results into the all_results_df
  all_results_df <- rbind(all_results_df, anova_results_df)
}


desired_order <- c("LDMC", "SLA", "LA", "PDMC", "SPA", "DA")

# Convert 'Trait' to a factor with the desired order
all_results_df$Trait <- factor(all_results_df$Trait, levels = desired_order)

# Create a stacked bar plot for all models 
ggplot(all_results_df, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained", title = "L. corniculatus - Variance Explained by Climatic Variables") +
  scale_fill_manual(values = c("#7e57c2", "#00bfc4","#fbc02d","#f7756d", "#ff9800"), 
                    labels = c("CO2", "(CO2 + Temperature) : Drought", "CO2 : Temperature", "Drought", "Temperature")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())


## For Significant Variables Only ----

all_results_df <- data.frame()

# Loop through each model in the l_model list
for (model_name in names(l_model)) {
  model <- l_model[[model_name]]
  
  # Perform ANOVA using the Anova function from the car package
  anova_results <- anova(model)
  
  # Convert ANOVA results to a data frame
  anova_results_df <- as.data.frame(anova_results)
  
  # Calculate the proportion of variance explained for each factor
  anova_results_df <- anova_results_df %>%
    mutate(PropVar = `Sum Sq` / sum(`Sum Sq`),
           Factor = rownames(anova_results_df),
           Trait = model_name) %>% # Add trait name
    filter(`Pr(>F)` < 0.05)%>%
    filter(Factor != "(Intercept)" & Factor != "Residuals") # Keep only significant variables # Exclude intercept and residuals
  
  # Combine results into the all_results_df
  all_results_df <- rbind(all_results_df, anova_results_df)
}

desired_order <- c("LDMC", "SLA", "LA", "PDMC", "SPA", "DA")

# Convert 'Trait' to a factor with the desired order
all_results_df$Trait <- factor(all_results_df$Trait, levels = desired_order)

# Create a stacked bar plot for all models 
ggplot(all_results_df, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained", title = "L. corniculatus - Variance Explained by Climatic Variables") +
  scale_fill_manual(values = c("#7e57c2", "#00bfc4","#fbc02d","#f7756d", "#ff9800"), 
                    labels = c("CO2", "(CO2 + Temperature) : Drought", "CO2 : Temperature", "Drought", "Temperature")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())


#### Including RDA ----

rda.lv <- rda(AusRDL[,c(9,18,20:23)] ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusRDL, scale = TRUE)

lv_RDA <- anova.cca(rda.lv, permutations = 9999, by = "terms")

lv_rda_results_df <- lv_RDA %>% 
  as.data.frame() %>%  # Ensure it is a data frame
  mutate(
    PropVar = Variance / sum(Variance),  # Calculate proportion of variance explained
    Factor = rownames(lv_RDA),           # Add Factor names as a new column
    Trait = "RDA"                        # Label this as RDA to differentiate in the plot
  ) %>%
  filter(Factor != "Residual")  # Correct filter condition to exclude residuals

all_results_selected <- all_results_df[, c("PropVar", "Factor", "Trait")]

# Select necessary columns from lv_rda_results_df
lv_rda_results_selected <- lv_rda_results_df[, c("PropVar", "Factor", "Trait")]

# Combine the selected columns
combined_results_df <- rbind(all_results_selected, lv_rda_results_selected)

combined_results_df
desired_order <- c("RDA", "LDMC", "SLA", "LA", "PDMC", "SPA", "DA")

# Convert 'Trait' to a factor with the desired order
combined_results_df$Trait <- factor(combined_results_df$Trait, levels = desired_order)
combined_results_df
# Create a stacked bar plot for all models
ggplot(combined_results_df, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = " Traits", y = "Proportion of Variance Explained", title = "L. corniculatus - Variance Explained by Climatic Variables") +
  scale_fill_manual(values = c("seagreen3", "mediumpurple", "royalblue3", "goldenrod2", "darkorange2"), 
                    labels = c("CO2", "Temperature", "Drought", "CO2:Temperature", "(CO2 + Temperature) :Drought")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())

### Creation of plot for Crepis ----

c_model <- list(
  DA = cv_sq_DAmodel,
  SPA = cv_sq_SPAmodel,
  PDMC = cv_sq_PDMCmodel,
  LA = cv_log_LAmodel,
  SLA = cv_log_SLAmodel,
  LDMC = cv_sq_LDMCmodel,
  SN = cv_sq_SNmodel,
  SM = cv_log_SMmodel
)


## For both Significant and Marginally Significant Variables ----

all_results_df <- data.frame()

# Loop through each model in the l_model list
for (model_name in names(c_model)) {
  model <- c_model[[model_name]]
  
  # Perform ANOVA using the Anova function from the car package
  anova_results <- anova(model)
  
  # Convert ANOVA results to a data frame
  anova_results_df <- as.data.frame(anova_results)
  
  # Calculate the proportion of variance explained for each factor
  anova_results_df <- anova_results_df %>%
    mutate(PropVar = `Sum Sq` / sum(`Sum Sq`),
           Factor = rownames(anova_results_df),
           Trait = model_name) %>% # Add trait name
    filter(`Pr(>F)` < 0.1)%>%
    filter(Factor != "(Intercept)" & Factor != "Residuals") # Keep only significant variables # Exclude intercept and residuals
  
  # Combine results into the all_results_df
  all_results_df <- rbind(all_results_df, anova_results_df)
}

desired_order <- c("SM", "SN", "LDMC", "SLA", "LA", "PDMC", "SPA", "DA")

# Convert 'Trait' to a factor with the desired order
all_results_df$Trait <- factor(all_results_df$Trait, levels = desired_order)

# Create a stacked bar plot for all models 
ggplot(all_results_df, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = " Traits", y = "Proportion of Variance Explained", title = "C. capillaris - Variance Explained by Climatic Variables") +
  scale_fill_manual(values = c("#7e57c2", "#00bfc4", "#f7756d","#fbc02d"), 
                    labels = c("(CO2 + Temperature) : Drought", "CO2 : Temperature", "Drought", "Temperature")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())

## For Significant Variables Alone ----

all_results_df <- data.frame()

# Loop through each model in the l_model list
for (model_name in names(c_model)) {
  model <- c_model[[model_name]]
  
  # Perform ANOVA using the Anova function from the car package
  anova_results <- anova(model)
  
  # Convert ANOVA results to a data frame
  anova_results_df <- as.data.frame(anova_results)
  
  # Calculate the proportion of variance explained for each factor
  anova_results_df <- anova_results_df %>%
    mutate(PropVar = `Sum Sq` / sum(`Sum Sq`),
           Factor = rownames(anova_results_df),
           Trait = model_name) %>% # Add trait name
    filter(`Pr(>F)` < 0.05)%>%
    filter(Factor != "(Intercept)" & Factor != "Residuals") # Keep only significant variables # Exclude intercept and residuals
  
  # Combine results into the all_results_df
  all_results_df <- rbind(all_results_df, anova_results_df)
}

desired_order <- c("SM", "SN", "LDMC", "SLA", "LA", "PDMC", "SPA", "DA")

# Convert 'Trait' to a factor with the desired order
all_results_df$Trait <- factor(all_results_df$Trait, levels = desired_order)

# Create a stacked bar plot for all models 
ggplot(all_results_df, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained", title = "C. capillaris - Variance Explained by Climatic Variables") +
  scale_fill_manual(values = c("#f7756d", "#00bfc4")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())


#### Including RDA ----
rda.cv <- rda(AusRDC[,c(9:10,12:17)] ~ CO2 + Temperature + Drought + CO2:Temperature + CO2:Drought, data = AusRDC, scale = TRUE)

cv_RDA <- anova.cca(rda.cv, permutations = 9999, by = "terms")

cv_rda_results_df <- lv_RDA %>% 
  as.data.frame() %>%  # Ensure it is a data frame
  mutate(
    PropVar = Variance / sum(Variance),  # Calculate proportion of variance explained
    Factor = rownames(lv_RDA),           # Add Factor names as a new column
    Trait = "RDA"                        # Label this as RDA to differentiate in the plot
  ) %>%
  filter(Factor != "Residual")  # Correct filter condition to exclude residuals

all_results_selected <- all_results_df[, c("PropVar", "Factor", "Trait")]

# Select necessary columns from lv_rda_results_df
cv_rda_results_selected <- cv_rda_results_df[, c("PropVar", "Factor", "Trait")]

# Combine the selected columns
combined_results_df <- rbind(all_results_selected, cv_rda_results_selected)

combined_results_df
desired_order <- c("RDA", "SM", "SN", "LDMC", "SLA", "LA", "PDMC", "SPA", "DA")

# Convert 'Trait' to a factor with the desired order
combined_results_df$Trait <- factor(combined_results_df$Trait, levels = desired_order)
combined_results_df
# Create a stacked bar plot for all models
ggplot(combined_results_df, aes(x = Trait, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = " Traits", y = "Proportion of Variance Explained", title = "C. capillaris - Variance Explained by Climatic Variables") +
  scale_fill_manual(values = c("seagreen3", "mediumpurple", "royalblue3", "goldenrod2", "darkorange2"), 
                    labels = c("CO2", "Temperature", "Drought", "CO2:Temperature", "(CO2 + Temperature) :Drought")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Increase x-axis text size
        axis.text.y = element_text(size = 14),  # Increase y-axis text size (trait names)
        axis.title = element_text(size = 16),   # Increase axis title font size
        plot.title = element_text(size = 18),   # Increase plot title font size
        legend.position = "bottom", 
        legend.title = element_blank())



### Plotting for a single Trait ----

anova_results <- Anova(lv_DAmodel, type = "III")

# Convert ANOVA results to a data frame
anova_results_df <- as.data.frame(anova_results)

# Calculate the proportion of variance explained for each factor
anova_results_df <- anova_results_df %>%
  mutate(PropVar = `Sum Sq` / sum(`Sum Sq`),
         Factor = rownames(anova_results_df)) %>%
  filter(Factor != "(Intercept)" & Factor != "Residuals") # Exclude intercept and residuals

# Create a stacked bar plot
ggplot(anova_results_df, aes(x = 1, y = PropVar, fill = Factor)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = NULL, y = "Proportion of Variance Explained", title = "Proportion of Variance Explained by Predictors") +
  scale_fill_manual(values = c("grey80", "grey50", "black", "red", "blue"), # Adjust colors to match the example
                    labels = c("CO2", "Temperature", "Drought", "CO2:Temperature", "CO2:Drought")) +
  coord_flip() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom", 
        legend.title = element_blank())
