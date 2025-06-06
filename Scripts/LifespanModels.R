# Statistical Analysis Script for the manuscript "Persistent maternal age effects on 
# male offspring fitness in wild Soay sheep"

# By Sanjana Ravindran 
set.seed(2023)
# # Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(glmmTMB, DHARMa, countreg, gratia, ggpubr, gtsummary, officer, tidyverse, cowplot, mgcv, mgcv.helper, ggeffects, MetBrewer, GGally, reshape2, patchwork, performance, mgcViz)

# Load Data 
lifetime_data <- read.table("./Datasets/OffspringLifeTimeData.csv", header=T, sep="\t")

# Convert to factors 
cols <- c("Sex", "Twin", "MumAgeClass", "DadAgeClass", "MotherID", "FatherID", "BirthYear")
lifetime_data[cols] <- lapply(lifetime_data[cols], factor)

### Maternal Age effects on lifespan

# Setting number of cores to use
ctrl <- list(nthreads=8)

# Models
m.lifespan1_split <- gam(Lifespan ~  MumAge + 
                        s(MumAge, k=4, bs = 'tp', m = c(2,0)) + 
                        DadAge + 
                        Sex  + 
                        Twin + 
                        MumLifespan +
                        DadLifespan +
                        MumAgeClass +
                        DadAgeClass +
                        s(MotherID, bs="re") + 
                        s(FatherID, bs="re") + 
                        s(BirthYear, bs="re"),
                      data=lifetime_data,
                      method="REML",
                      family = nb(link="log"), control=ctrl)
save(m.lifespan1_split, file = "./Models/Lifespan_Main_MAsplit.RData")

# Model estimating linear and non-linearity in one single term
m.lifespan1 <-  update(m.lifespan1_split, .~. -(MumAge + 
                                                   s(MumAge, k=4, bs = 'tp', m = c(2,0))) + 
                             s(MumAge, k=4, bs = 'tp'))
save(m.lifespan1, file = "./Models/Lifespan_Main_MA_NoSplit.RData")

# Testing interaction with offspring Sex 
# Run Model
m.lifespan1_ml <- update(m.lifespan1, method="ML")
m.lifespan1_ml_sex <- update(m.lifespan1, .~. -s(MumAge, k=4, bs = 'tp') +
                               DadAge:Sex +
                               s(MumAge,by=Sex, k=4, bs = 'tp'),
                             method="ML")
m.lifespan1_ml_matsex <- update(m.lifespan1, .~. -s(MumAge, k=4, bs = 'tp') +
                                  s(MumAge,by=Sex, k=4, bs = 'tp'),
                                method="ML")
m.lifespan1_ml_patsex <- update(m.lifespan1, .~. +
                                  DadAge:Sex, 
                                method="ML")

save(m.lifespan1_ml, file = "./Models/Lifespan_Main_ML.RData")
save(m.lifespan1_ml_sex, file = "./Models/Lifespan_SexInteractions_ML.RData")
save(m.lifespan1_ml_matsex, file = "./Models/Lifespan_OnlyMatSexInteraction_ML.RData")
save(m.lifespan1_ml_patsex, file = "./Models/Lifespan_OnlyPatSexInteraction_ML.RData")

# Check AICs
AIC(m.lifespan1_ml)
AIC(m.lifespan1_ml_sex)
AIC(m.lifespan1_ml_matsex)
AIC(m.lifespan1_ml_patsex)

# Final Sex model with method=REML
m.lifespan1_matsex <- update(m.lifespan1_ml_matsex, method="REML")
save(m.lifespan1_matsex, file = "./Models/Lifespan_MatSexInteraction.RData")

# Use ordered factor for Sex to check if differences are significant
lifetime_data$SexOrd <- ordered(lifetime_data$Sex)

m.lifespan1_sexord <- update(m.lifespan1, .~. -s(MumAge, k=4, bs = 'tp') +
                               s(MumAge,by=SexOrd, k=4, bs = 'tp'))
save(m.lifespan1_sexord, file = "./Models/Lifespan_OrdMatSexInteraction.RData")

################################################################################################