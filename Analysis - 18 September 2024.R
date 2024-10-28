# Installing Packages ----
install.packages("corrplot")
install.packages("ggrepel")
install.packages("Hmisc")
install.packages("VennDiagram")
install.packages("devtools")
install.github("vqv/ggbiplot") # Need to load devtools package before installing
install.packages("ggeffects")
install.packages("igraph")
install.packages("qgraph")
install.packages("vegan")
library(readr) # For easier importing of datasets 
library(lme4) # For linear models 
library(lmerTest) # For linear models
library(ggplot2) # For better plots
library(vegan)
library(dplyr) # For Variance Partitioning Barplot, for %>%
library(corrplot) # For Correlation Matrix and Corellogram
library(igraph) # For Correlation Matrix
library(qgraph) # For Correlation Matrix
library(qgraph)
library(igraph)
library(Hmisc)
library(ggrepel)
library(tidyr)
library(ggbiplot)
library(devtools)
library(ggbiplot)
library(plyr)
library(ggeffects)
library(car) # For Variance Partitioning Barplot, for Anova
library(reshape2) # For correlogram
library(vegan)
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
### Display Area ----
ggplot(AusL, aes(x = AusL$Temperature, y = AusL$`Display Area (DA)(cm2)`))+
  geom_boxplot() +
  labs(
    title = "Lotus - Effect of Temperature on Display Area",
    x = "Temperature Level",
    y = "Display Area (DA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


ggplot(AusL, aes(x = CO2, y = `Display Area (DA)(cm2)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Lotus - Effect of CO2 Level on DA by Temperature",
    x = "CO2 Level",
    y = "Display Area (DA)(cm2)",
    fill = "Temperature"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


ggplot(AusL, aes(x = CO2, y = `Display Area (DA)(cm2)`, fill = Drought)) +
  geom_boxplot() +
  facet_wrap(~ Drought) +
  labs(
    title = "Lotus - Effect of CO2 and Temperature on DA by Drought Level",
    x = "CO2 + Temperature Level",
    y = "Display Area (DA)(cm2)",
    fill = "Drought"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

### Specific Petal Area ----
ggplot(AusL, aes(x = AusL$Drought, y = AusL$`Specific Petal Area (SPA)(cm2/g)`))+
  geom_boxplot() +
  labs(
    title = "Lotus - Effect of Drought on Specific Petal Area",
    x = "Drought Level",
    y = "Specific Petal Area (SPA)(cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

### Leaf Area ----
ggplot(AusL, aes(x = AusL$CO2, y = AusL$`Leaf Area (LA)(cm2)`))+
  geom_boxplot() +
  labs(
    title = "Lotus - Effect of CO2 on Leaf Area",
    x = "CO2 Level",
    y = "Leaf Area (LA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

ggplot(AusL, aes(x = CO2, y = `Leaf Area (LA)(cm2)`, fill = Drought)) +
  geom_boxplot() +
  facet_wrap(~ Drought) +
  labs(
    title = "Lotus - Effect of CO2 and Temperature on LA by Drought Level",
    x = "CO2 + Temperature Level",
    y = "Leaf Area (LA)(cm2)",
    fill = "Drought"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


### Specific Leaf Area ----
ggplot(AusL, aes(x = CO2,y= `Specific Leaf Area (SLA)(cm2/g)`)) +
  geom_boxplot() +
  labs(
    title = "Lotus - Effect of CO2 Level on SLA",
    x = "CO2 Level",
    y = "Specific Leaf Area (SLA) (cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


ggplot(AusL, aes(x = AusL$Drought, y= AusL$`Specific Leaf Area (SLA)(cm2/g)`)) +
  geom_boxplot()+
  labs(
    title = "Lotus - Effect of Drought on SLA",
    x = "Drought Level",
    y = "Specific Leaf Area (SLA) (cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


ggplot(AusL, aes(x = CO2, y = `Specific Leaf Area (SLA)(cm2/g)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Lotus - Effect of CO2 Level on SLA by Temperature Level",
    x = "CO2 Level",
    y = "Specific Leaf Area (SLA) (cm2/g)",
    fill = "Temperature"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


### Leaf Dry Matter Content ----

ggplot(AusL, aes(x = CO2,y= `Leaf Dry Matter Content (LDMC)(g/g)`)) +
  geom_boxplot() +
  labs(
    title = "Lotus - Effect of CO2 Level on LDMC",
    x = "CO2 Level",
    y = "Leaf Dry Matter Content (LDMC) (g/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

ggplot(AusL, aes(x = Drought, y= AusL$`Leaf Dry Matter Content (LDMC)(g/g)`)) +
  geom_boxplot()+
  labs(
    title = "Lotus - Effect of Drought on LDMC",
    x = "Drought Level",
    y = "Leaf Dry Matter Content (LDMC) (g/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

ggplot(AusL, aes(x = CO2, y = `Leaf Dry Matter Content (LDMC)(g/g)`, fill = Drought)) +
  geom_boxplot() +
  facet_wrap(~ Drought) +
  labs(
    title = "Lotus - Effect of CO2 and Temperature on LDMC by Drought",
    x = "CO2 + Temperature Level",
    y = "Leaf Dry Matter Content (LDMC)(g/g)",
    fill = "Drought"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
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

### Display Area ----
ggplot(AusC, aes(x = AusC$Drought, y= `Display Area (DA)(cm2)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Drought on DA",
    x = "Drought Level",
    y = "Display Area (DA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

### Specific Petal Area ----
ggplot(AusC, aes(x = AusC$Drought, y= `Specific Petal Area (SPA)(cm2/g)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Drought on SPA",
    x = "Drought Level",
    y = "Specific Petal Area (SPA)(cm2/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

ggplot(AusC, aes(x = CO2, y = `Specific Petal Area (SPA)(cm2/g)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Crepis - Effect of CO2 Level on SPA by Temperature",
    x = "CO2 Level",
    y = "Specific Petal Area (SPA)(cm2/g)",
    fill = "Temperature"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

### Petal Dry Matter Content ----
ggplot(AusC, aes(x = AusC$Drought, y= `Petal Dry Matter Content (PDMC)(g/g)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Drought on PDMC",
    x = "Drought Level",
    y = "Petal Dry Matter Content (PDMC)(g/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

ggplot(AusC, aes(x = CO2, y = `Petal Dry Matter Content (PDMC)(g/g)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Crepis - Effect of CO2 Level on PDMC by Temperature",
    x = "CO2 Level",
    y = "Petal Dry Matter Content (PDMC)(g/g)",
    fill = "Temperature"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

View(AusL)
### Leaf Area ----
ggplot(AusC, aes(x = AusC$Temperature, y= `Leaf Area (LA)(cm2)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Temperature on LA",
    x = "Temperature Level",
    y = "Leaf Area (LA)(cm2)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

### Leaf Dry Matter Content ----
ggplot(AusC, aes(x = AusC$Drought, y= `Leaf Dry Matter Content (LDMC)(g/g)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Drought on LDMC",
    x = "Drought Level",
    y = "Leaf Dry Matter Content (g/g)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  ) 

### Seed Number ----
ggplot(AusC, aes(x = AusC$Temperature, y= `Number of Seeds (SN)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Temperature on SN",
    x = "Temperature Level",
    y = "Number of Seeds (SN)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

ggplot(AusC, aes(x = AusC$Drought, y= `Number of Seeds (SN)`)) +
  geom_boxplot()+
  labs(
    title = "Crepis - Effect of Drought on SN",
    x = "Drought Level",
    y = "Number of Seeds (SN)"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )


ggplot(AusC, aes(x = CO2, y = `Number of Seeds (SN)`, fill = Temperature)) +
  geom_boxplot() +
  facet_wrap(~ Temperature) +
  labs(
    title = "Crepis - Effect of CO2 on SN by Temperature",
    x = "CO2 Level",
    y = "Number of Seeds (SN)",
    fill = "Temperature"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    axis.text.x = element_text(size = 12),   # Increase x-axis text size
    axis.text.y = element_text(size = 12)    # Increase y-axis text size
  )

## Plotting Using GGpredict ----
ggpredict(c_sq_DAmodel)
pr=ggpredict(c_sq_DAmodel, terms = c("Drought", "Temperature"), type = "fixed")
pr
plot(pr, dodge = 1, , colors = c("red", "blue" ), ci_style = "errorbar")+
  labs(title = NULL,
       x = "Drought",
       y = "DA")
?plot.ggeffects

# Principal Component Analysis ----
## Creation of PCA models ----
View(AusPC)
lmv <- prcomp(AusPC[AusPC$Species == "Lotus", c(9,18,20:23)], scale = TRUE)
cmv <- prcomp(AusPC[AusPC$Species == "Crepis", c(9,18,20:23)], scale = TRUE)
lcav <- prcomp(AusPC[,c(16, 18,19,20)], center = TRUE, scale = TRUE)
names(AusPC)

colnames(AusPC)
## Plotting of PCA models ---- 
names(AusPC)[names(AusPC) == "Specific Petal Area (SPA)(cm2/g)"] <- "SPA"
names(AusPC)[names(AusPC) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "PDMC"
names(AusPC)[names(AusPC) == "Specific Leaf Area (SLA)(cm2/g)"] <- "SLA"
names(AusPC)[names(AusPC) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "LDMC"
names(AusPC)[names(AusPC) == "Number of Seeds (SN)"] <- "SN"
names(AusPC)[names(AusPC) == "Seed Mass (SM)(g)"] <- "SM"
names(AusPC)[names(AusPC) == "Display Area (DA)(cm2)"] <- "DA"
names(AusPC)[names(AusPC) == "Leaf Area (LA)(cm2)"] <- "LA"

biplot(lmv, scale = 0)
biplot(cmv, scale = 1, alpha = 0)
biplot(lcav, scale = 0)
scores(lmv, display = "species")
?biplot
lmv_scores <- lmv$x
lmv_scores[,3]
summary(lmv)
lmv
### PCA plots using ggbiplot ----

ggbiplot(lmv, alpha = 0, varname.size = 5) +
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(x = "PC1 (35.4%)", y = "PC2 (29.1%)")

ggbiplot(cmv, alpha = 0, varname.size = 5) + 
  xlim(-2, 2) + 
  ylim(-2, 2) +
  labs(x = "PC1 (33.6%)", y = "PC2 (21.6%)")

## Extracting PC Scores ----
summary(cmv)
lcav$x
AusPC1 <- cbind(AusPC, lcav$x[,1:4])
View(AusPC1)
View(AusC)
## Plotting with GG plot ----
ggplot(AusPC, aes(PC1,PC2, colour = Species, fill = Species)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")

## Extracting PCA Results----
### For Lotus ----

lotus_pca_trait_loadings <- lmv$rotation*lmv$sdev
lotus_pca_eigenvalues <- lmv$sdev^2
lotus_pca_variance <- lotus_pca_eigenvalues/sum(lotus_pca_eigenvalues)
lotus_pca_cum_variance <- cumsum(lotus_pca_variance)
lotus_pca_trait_loadings
lotus_pca_cum_variance
summary(cmv)

write.csv(lotus_pca_trait_loadings, "lotus_pca_trait_loadings.csv", row.names = FALSE)
write.csv(lotus_pca_eigenvalues, "lotus_pca_eigenvalues.csv", row.names = FALSE)
write.csv(lotus_pca_variance, "lotus_pca_variance.csv")
write.csv(lotus_pca_cum_variance, "lotus_pca_cum_variance.csv")

### For Crepis ----
cmv
crepis_pca_trait_loadings <- cmv$rotation*cmv$sdev
crepis_pca_eigenvalues <- cmv$sdev^2
crepis_pca_variance <- crepis_pca_eigenvalues/sum(crepis_pca_eigenvalues)
crepis_pca_cum_variance <- cumsum(crepis_pca_variance)
crepis_pca_trait_loadings
crepis_pca_cum_variance

write.csv(crepis_pca_trait_loadings, "crepis_pca_trait_loadings.csv", row.names = FALSE)
write.csv(crepis_pca_eigenvalues, "crepis_pca_eigenvalues.csv", row.names = FALSE)
write.csv(crepis_pca_variance, "crepis_pca_variance.csv")
write.csv(crepis_pca_cum_variance, "crepis_pca_cum_variance.csv")



# Proportion of variance explained
proportion_variance_explained <- eigenvalues / sum(eigenvalues)
proportion_variance_explained

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
names(AusRDL)[names(AusRDL) == "Specific Petal Area (SPA)(cm2/g)"] <- "SPA"
names(AusRDL)[names(AusRDL) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "PDMC"
names(AusRDL)[names(AusRDL) == "Specific Leaf Area (SLA)(cm2/g)"] <- "SLA"
names(AusRDL)[names(AusRDL) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "LDMC"
names(AusRDL)[names(AusRDL) == "Display Area (DA)(cm2)"] <- "DA"
names(AusRDL)[names(AusRDL) == "Leaf Area (LA)(cm2)"] <- "LA"
names(AusRDL)[names(AusRDL) == "CO2"] <- "C"
names(AusRDL)[names(AusRDL) == "Temperature"] <- "T"
names(AusRDL)[names(AusRDL) == "Drought"] <- "D"


rda.l <- rda(AusRDL[,c(9,18,20:23)] ~ CO2*Temperature*Drought, data = AusRDL, scale = TRUE)
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

# Rename the fourth row
rownames(bp_scores)[4] <- "(CT):D"

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
names(AusRDC)[names(AusRDC) == "Specific Petal Area (SPA)(cm2/g)"] <- "SPA"
names(AusRDC)[names(AusRDC) == "Petal Dry Matter Content (PDMC)(g/g)"] <- "PDMC"
names(AusRDC)[names(AusRDC) == "Specific Leaf Area (SLA)(cm2/g)"] <- "SLA"
names(AusRDC)[names(AusRDC) == "Leaf Dry Matter Content (LDMC)(g/g)"] <- "LDMC"
names(AusRDC)[names(AusRDC) == "Display Area (DA)(cm2)"] <- "DA"
names(AusRDC)[names(AusRDC) == "Leaf Area (LA)(cm2)"] <- "LA"
names(AusRDC)[names(AusRDC) == "Number of Seeds (SN)"] <- "SN"
names(AusRDC)[names(AusRDC) == "Seed Mass (SM)(g)"] <- "SM"
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



# Pearson Correlation ----
## Lotus ----
cor.test(AusL$`Specific Petal Area (SPA)(cm2/g)`, AusL$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson") # Significant
cor.test(AusL$`Petal Mass per Area (PMA)(g/cm2)`, AusL$`Petal Dry Matter Content (PDMC) (g/g)`, method = "pearson") # Significant
cor.test(AusL$`Specific Leaf Area (SLA)(cm2/g)`, AusL$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson") # Significant
cor.test(AusL$`Specific Petal Area (SPA)(cm2/g)`, AusL$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusL$`Petal Mass per Area (PMA)(g/cm2)`, AusL$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusL$`Petal Dry Matter Content (PDMC)(g/g)`, AusL$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")
cor.test(AusL$`Display Area (DA)(cm2)`, AusL$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson") #Significant
cor.test(AusL$`Display Area (DA)(cm2)`, AusL$`Specific Petal Area (SPA)(cm2/g)`, method = "pearson") # Significant
cor.test(AusL$`Display Area (DA)(cm2)`, AusL$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusL$`Display Area (DA)(cm2)`, AusL$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")
cor.test(AusL$`Specific Leaf Area (SLA)(cm2/g)`, AusL$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson")
cor.test(AusL$`Specific Petal Area (SPA)(cm2/g)`, AusL$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")
cor.test(AusL$`Petal Dry Matter Content (PDMC)(g/g)`, AusL$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")



## Crepis ----
cor.test(AusC$`Specific Petal Area (SPA)(cm2/g)`, AusC$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson") # Significant
cor.test(AusC$`Petal Mass Per Area (PMA)(g/cm2)`, AusC$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson") # Significant
cor.test(AusC$`Specific Leaf Area (SLA)(cm2/g)`, AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson") 
cor.test(AusC$`Specific Petal Area (SPA)(cm2/g)`, AusC$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusC$`Petal Mass Per Area (PMA)(g/cm2)`, AusC$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, AusC$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson")
cor.test(AusC$`Number of Seeds (SN)`, AusC$`Seed Mass (SM)(g)`, method = "pearson")
cor.test(AusC$`Number of Seeds (SN)`, AusC$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson") # Significant
cor.test(AusC$`Number of Seeds (SN)`, AusC$`Leaf Dry Matter Content (g/g)`, method = "pearson") # Significant
cor.test(AusC$`Number of Seeds (SN)`, AusC$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusC$`Display Area (DA)(cm2)`, AusC$`Petal Dry Matter Content (PDMC)(g/g)`, method = "pearson")
cor.test(AusC$`Display Area (DA)(cm2)`, AusC$`Specific Petal Area (SPA)(cm2/g)`, method = "pearson") # Significant
cor.test(AusC$`Display Area (DA)(cm2)`, AusC$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusC$`Display Area (DA)(cm2)`, AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")
cor.test(AusC$`Display Area (DA)(cm2)`, AusC$`Number of Seeds (SN)`, method = "pearson")
cor.test(AusC$`Display Area (DA)(cm2)`, AusC$`Seed Mass (SM)(g)`, method = "pearson")
cor.test(AusC$`Specific Petal Area (SPA)(cm2/g)`, AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")
cor.test(AusC$`Specific Petal Area (SPA)(cm2/g)`, AusC$`Number of Seeds (SN)`, method = "pearson")
cor.test(AusC$`Specific Petal Area (SPA)(cm2/g)`,  AusC$`Seed Mass (SM)(g)`, method = "pearson")
cor.test(AusC$`Petal Dry Matter Content (PDMC)(g/g)`, AusC$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson")
cor.test(AusC$`Petal Dry Matter Content (PDMC)(g/g)`, AusC$`Number of Seeds (SN)`, method = "pearson")
cor.test(AusC$`Petal Dry Matter Content (PDMC)(g/g)`, AusC$`Seed Mass (SM)(g)`, method = "pearson")
cor.test(AusC$`Specific Leaf Area (SLA)(cm2/g)`, AusC$`Number of Seeds (SN)`, method = "pearson") 
cor.test(AusC$`Specific Leaf Area (SLA)(cm2/g)`, AusC$`Seed Mass (SM)(g)`, method = "pearson") 
cor.test(AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, AusC$`Number of Seeds (SN)`, method = "pearson")
cor.test(AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, AusC$`Seed Mass (SM)(g)`, method = "pearson")
cor.test(AusC$`Leaf Area (LA)(cm2)`, AusC$`Specific Leaf Area (SLA)(cm2/g)`, method = "pearson") # Significant
cor.test(AusC$`Leaf Area (LA)(cm2)`, AusC$`Leaf Dry Matter Content (LDMC)(g/g)`, method = "pearson")


# Correlation Plots ----
## Lotus ----
?ggplot
SPA_PDMC_Plot <- ggplot(AusL, aes(x = AusL$`Specific Petal Area (SPA)(cm2/g)`, y = AusL$`Petal Dry Matter Content (PDMC)(g/g)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Specific Petal Area (SPA) and Petal Dry Matter Content (PDMC)",
       x = "Specfic Petal Area (SPA) (cm2/g)",
       y = "Petal Dry Matter Content (PDMC) (g/g)") +
  theme_minimal()
print(SPA_PDMC_Plot)

SPA_DA_Plot <- ggplot(AusL, aes(x = AusL$`Specific Petal Area (SPA)(cm2/g)`, y = AusL$`Display Area (DA)(cm2)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Specific Petal Area (SPA) and Display Area (DA)",
       x = "Specfic Petal Area (SPA) (cm2/g)",
       y = "Display Area (DA)(cm2)") +
  theme_minimal()
print(SPA_DA_Plot)

DA_PDMC_Plot <- ggplot(AusL, aes(x = AusL$`Display Area (DA)(cm2)`, y = AusL$`Petal Dry Matter Content (PDMC)(g/g)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Display Area (DA) and Petal Dry Matter Content (PDMC)",
       x = "Display Area (DA)(cm2)",
       y = "Petal Dry Matter Content (PDMC)(g/g)") +
  theme_minimal()
print(DA_PDMC_Plot)

SLA_LDMC_Plot <- ggplot(AusL, aes(x = AusL$`Specific Leaf Area (SLA) (cm2/g)`, y = AusL$`Leaf Dry Matter Content (LDMC) (g/g)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Lotus - SLA and LDMC",
       x = "Specific Leaf Area (SLA) (cm2/g)",
       y = "Leaf Dry Matter Content (LDMC) (g/g)") +
  theme_minimal()
print(SLA_LDMC_Plot)

## Crepis ----

SPA_DA_Plot <- ggplot(AusC, aes(x = AusC$`Specific Petal Area (SPA)(cm2/g)`, y = AusC$`Display Area (DA)(cm2)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Specific Petal Area (SPA) and Display Area (DA)",
       x = "Specfic Petal Area (SPA)(cm2/g)",
       y = "Display Area (DA)(cm2)") +
  theme_minimal()
print(SPA_DA_Plot)

SPA_PDMC_Plot <- ggplot(AusC, aes(x = AusC$`Specific Petal Area (SPA)(cm2/g)`, y = AusC$`Petal Dry Matter Content (PDMC)(g/g)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Specific Petal Area (SPA) and Petal Dry Matter Content (PDMC)",
       x = "Specfic Petal Area (SPA) (cm2/g)",
       y = "Petal Dry Matter Content (PDMC) (g/g)") +
  theme_minimal()
print(SPA_PDMC_Plot)

SLA_LA_Plot <- ggplot(AusC, aes(x = AusC$`Specific Leaf Area (SLA)(cm2/g)`, y = AusC$`Leaf Area (LA)(cm2)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Specific Leaf Area (SLA) and Leaf Area (LA)",
       x = "Specific Leaf Area (SLA)",
       y = "Leaf Area (LA)") +
  theme_minimal()
print(SLA_LA_Plot)

SN_PDMC_Plot <- ggplot(AusC, aes(x = AusC$`Number of Seeds (SN)`, y = AusC$`Petal Dry Matter Content (PDMC)(g/g)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Crepis - SN and PDMC",
       x = "Number of Seeds (SN)",
       y = "Petal Dry Matter Content (PDMC) (g/g)") +
  theme_minimal()
print(SN_PDMC_Plot)

SN_LDMC_Plot <- ggplot(AusC, aes(x = AusC$`Number of Seeds (SN)`, y = AusC$`Leaf Dry Matter Content (g/g)`)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Crepis - SN and LDMC",
       x = "Number of Seeds (SN)",
       y = "Leaf Dry Matter Content (g/g)") +
  theme_minimal()
print(SN_LDMC_Plot)

combined_data <- AusC %>%
  select(SN = `Number of Seeds (SN)`, PDMC = `Petal Dry Matter Content (PDMC)(g/g)`, LDMC = `Leaf Dry Matter Content (g/g)`) %>%
  gather(key = "Type", value = "Content", PDMC, LDMC)

# Create the combined plot
combined_plot <- ggplot(combined_data, aes(x = SN, y = Content)) +
  geom_point() +  # Scatter plot points
  geom_smooth(method = "lm", col = "blue") +  # Regression line
  labs(title = "Crepis - SN and Dry Matter Content",
       x = "Number of Seeds (SN)",
       y = "Dry Matter Content (g/g)") +
  facet_wrap(~Type, scales = "free_y") +
  theme_minimal()

print(combined_plot)



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


# Variation Partitioning ----
## For Lotus ----

l_DA_var <- varpart(AusRDL$DA, AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_DA_var$part$indfract$Adj.R.square <- l_DA_var$part$indfract$Adj.R.square *100
plot(l_DA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)

title(main = "Variance in DA of L. corniculatus")

l_SPA_var <- varpart(AusRDL$SPA, AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_SPA_var$part$indfract$Adj.R.square <- l_SPA_var$part$indfract$Adj.R.square *100
plot(l_SPA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SPA of L. corniculatus")

l_PDMC_var <- varpart(AusRDL$PDMC, AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_PDMC_var$part$indfract$Adj.R.square <- l_PDMC_var$part$indfract$Adj.R.square *100
plot(l_PDMC_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
title(main = "Variance in PDMC of L. corniculatus")

l_LA_var <- varpart(AusRDL$LA, AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_LA_var$part$indfract$Adj.R.square <- l_LA_var$part$indfract$Adj.R.square *100
plot(l_LA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
title(main = "Variance in LA of L. corniculatus")

l_SLA_var <- varpart(AusRDL$SLA, AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_SLA_var$part$indfract$Adj.R.square <- l_SLA_var$part$indfract$Adj.R.square *100
plot(l_SLA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SLA of L. corniculatus")


l_LDMC_var <- varpart(AusRDL$LDMC, AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_LDMC_var$part$indfract$Adj.R.square <- l_LDMC_var$part$indfract$Adj.R.square *100
plot(l_LDMC_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in LDMC of L. corniculatus")

l_at_var <- varpart(AusRDL[,c(9,18,20:23)], AusRDL$CO2, AusRDL$Temperature, AusRDL$Drought)
l_at_var$part$indfract$Adj.R.square <- l_at_var$part$indfract$Adj.R.square *100
plot(l_at_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)

## For Crepis ----

c_DA_var <- varpart(AusRDC$DA, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_DA_var$part$indfract$Adj.R.square <- c_DA_var$part$indfract$Adj.R.square *100
plot(c_DA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in DA of C. capillaris")

c_SPA_var <- varpart(AusRDC$SPA, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_SPA_var$part$indfract$Adj.R.square <- c_SPA_var$part$indfract$Adj.R.square *100
plot(c_SPA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SPA of C. capillaris")

c_PDMC_var <- varpart(AusRDC$PDMC, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_PDMC_var$part$indfract$Adj.R.square <- c_PDMC_var$part$indfract$Adj.R.square *100
plot(c_PDMC_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in PDMC of C. capillaris")

c_LA_var <- varpart(AusRDC$LA, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_LA_var$part$indfract$Adj.R.square <- c_LA_var$part$indfract$Adj.R.square *100
plot(c_LA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in LA of C. capillaris")

c_SLA_var <- varpart(AusRDC$SLA, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_SLA_var$part$indfract$Adj.R.square <- c_SLA_var$part$indfract$Adj.R.square *100
plot(c_SLA_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SLA of C. capillaris")

c_LDMC_var <- varpart(AusRDC$LDMC, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_LDMC_var$part$indfract$Adj.R.square <- c_LDMC_var$part$indfract$Adj.R.square *100
plot(c_LDMC_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in LDMC of C. capillaris")

c_SN_var <- varpart(AusRDC$SN, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_SN_var$part$indfract$Adj.R.square <- c_SN_var$part$indfract$Adj.R.square *100
plot(c_SN_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 1, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SN of C. capillaris")

c_SM_var <- varpart(AusRDC$SM, AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_SM_var$part$indfract$Adj.R.square <- c_SM_var$part$indfract$Adj.R.square *100
plot(c_SM_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SM of C. capillaris")

c_at_var <- varpart(AusRDC[,c(9,10,12:17)],AusRDC$CO2, AusRDC$Temperature, AusRDC$Drought)
c_at_var$part$indfract$Adj.R.square <- c_at_var$part$indfract$Adj.R.square *100
plot(c_at_var,
     Xnames = c("CO2", "Temperature", "Drought"), # name the partitions
     bg = c("seagreen3", "mediumpurple", "royalblue3"), alpha = 80, # colour the circles
     digits = 2, # only show 2 digits
     cex = 1.5)
title(main = "Variance in SM of C. capillaris")


# Correlation Networks ----
## For Lotus ----
View(Aus)
l_cor
l_corl_cor <- cor(AusRDL[,c(9,18,20:23)])
l_network <- graph_from_adjacency_matrix(l_corl_cor, weighted = TRUE, mode = "undirected", diag = FALSE)
plot(l_network)
l_network_graph <- qgraph(
  l_corl_cor, 
  graph = "cor",  # Use partial correlations
  layout = "spring", 
  sampleSize = nrow(AusRDL[,c(9,18,20:23)]),
  vsize = 6, 
  cut = 0.2,  
  maximum = 0.45, 
  border.width = 1.5,
  alpha = 0.05,
  labels = c("DA", "SPA", "PDMC","LA", "SLA", "LDMC"),
  posCol = "royalblue4",
  negCol = "brown1"
)

## For Crepis ----
c_corl_cor <- cor(AusRDC[,c(9,10,12:17)])
c_corl_cor
cor(AusRDC[,c(9,10,12:17)])
c_network <- graph_from_adjacency_matrix(c_corl_cor, weighted=T, mode="undirected", diag=F)
c_network_graph <- qgraph(c_corl_cor, 
                          graph = "cor",  # Use partial correlations
                          layout = "spring", 
                          sampleSize = nrow(AusRDC[,c(9,10,12:17)]),
                          vsize = 6, 
                          cut = 0.2,  
                          maximum = 0.45, 
                          border.width = 1.5,
                          labels = c("DA", "SPA", "PDMC", "LA","SLA", "LDMC", "SN", "SM"),
                          posCol = "royalblue4",
                          negCol = "brown1"
)

summary(c_sq_DAmodel)


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
  labs(x = " Traits", y = "Proportion of Variance Explained", title = "L. corniculatus - Variance Explained by Climatic Variables") +
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
  labs(x = " Traits", y = "Proportion of Variance Explained", title = "L. corniculatus - Variance Explained by Climatic Variables") +
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
  labs(x = " Traits", y = "Proportion of Variance Explained", title = "C. capillaris - Variance Explained by Climatic Variables") +
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
