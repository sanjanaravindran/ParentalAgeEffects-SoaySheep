# Statistical Analysis Script for the manuscript "Persistent maternal age effects on 
# male offspring fitness in wild Soay sheep"
set.seed(2023)
# # Load Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(glmmTMB, DHARMa, countreg, gratia, ggpubr, gtsummary, officer, tidyverse, cowplot, mgcv, mgcv.helper, ggeffects, MetBrewer, GGally, reshape2, patchwork, performance, mgcViz)

# Load Data 
lifetime_data <- read.table("./Datasets/OffspringLifeTimeData.csv", header=T, sep="\t")

# Convert to factors 
cols <- c("Sex", "Twin", "MumAgeClass", "DadAgeClass", "MotherID", "FatherID", "BirthYear")
lifetime_data[cols] <- lapply(lifetime_data[cols], factor)

### Maternal Age effects on LBS

# Setting number of cores to use
ctrl <- list(nthreads=8)

# Model separating linearity and non-linearity
m.lbs1_split <- gam(LBS ~  MumAge + 
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

# Save model
save(m.lbs1_split, file = "./Models/LBS_Main_MAsplit.RData")

# The above model but also accounting for offspring longevity
m.lbs1_split_lifespan <- update(m.lbs1_split, .~. + Lifespan)
# Save model
save(m.lbs1_split_lifespan, file = "./Models/LBS_Main_MAsplit_Lifespan.RData")

# Model estimating linear and non-linearity in one single term
m.lbs1 <-  update(m.lbs1_split, .~. -(MumAge + 
                                   s(MumAge, k=4, bs = 'tp', m = c(2,0))) + 
                     s(MumAge, k=4, bs = 'tp'))
save(m.lbs1 , file = "./Models/LBS_Main_MA_NoSplit.RData")

# Testing interaction with offspring sex 
# Run Model
m.lbs1_ml <- update(m.lbs1, method="ML")
m.lbs1_ml_sex <- update(m.lbs1, .~. -s(MumAge, k=4, bs = 'tp') +
                            DadAge:Sex +
                          s(MumAge, by=Sex, k=4, bs = 'tp'),
                          method="ML")
m.lbs1_ml_matsex <- update(m.lbs1, .~. -s(MumAge, k=4, bs = 'tp') +
                                  s(MumAge,by=Sex, k=4, bs = 'tp'),
                                method="ML")
m.lbs1_ml_patsex <- update(m.lbs1, .~. +
                                  DadAge:Sex, 
                                method="ML")

# Check AICs
AIC(m.lbs1_ml)
AIC(m.lbs1_ml_sex)
AIC(m.lbs1_ml_matsex)
AIC(m.lbs1_ml_patsex)

# Final Sex model with method=REML
m.lbs1_sex_final <- update(m.lbs1_ml_sex, method="REML")
save(m.lbs1_sex_final, file = "./Models/LBS_SexInteractions.RData")

# Use ordered factor for Sex to check if differences are significant
lifetime_data$SexOrd <- ordered(lifetime_data$Sex)

m.lbs1_sexord <- update(m.lbs1, .~. + DadAge:Sex +
                          s(MumAge,by=SexOrd, k=4, bs = 'tp'))

save(m.lbs1_sexord, file = "./Models/LBS_OrdSexInteraction.RData")


# Account for offspring longevity in the model
m.lbs1_sexord_ls <- update(m.lbs1_sexord, .~. + Lifespan)
save(m.lbs1_sexord_ls, file = "./Models/LBS_OrdSexInteraction_Lifespan.RData")

##############################################################################################
