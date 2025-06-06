# Statistical Analysis Script for the manuscript "Persistent maternal age effects on 
# male offspring fitness in wild Soay sheep"
# # Load Libraries
library(tidyverse); library(mgcv); library(gratia); library(cowplot); library(MetBrewer); library(patchwork)

# Load models
load("./Models/LBS_SexInteractions.RData")

# Plotting predictoins (maternal age)
lifetime_data_plot <- lifetime_data %>%
  mutate(DadAge = replace(DadAge, DadAge >= 10, 9))  %>%
  mutate(MumAge = replace(MumAge, MumAge >= 12, 11))

# Plotting offspring sex-specific predictions
new_data_kids_mum <- with(lifetime_data, expand.grid(MumAge = seq(min(MumAge), max(MumAge), length = 1000),
                                                 DadAge = 3.237,
                                                 MumLifespan = 8.901,
                                                 DadLifespan = 5.545,
                                                 MumAgeClass = c("AdultMum"),
                                                 DadAgeClass = c("AdultDad"),
                                                 Sex = unique(Sex),
                                                 Twin = c("0"),
                                                 MotherID = MotherID[1],
                                                 FatherID=FatherID[1],
                                                 BirthYear=BirthYear[1]))

new_data_kids_dad <- with(lifetime_data, expand.grid(DadAge = seq(min(DadAge), max(DadAge), length = 1000),
                                                 MumAge = 3.879,
                                                 MumLifespan = 8.901,
                                                 DadLifespan = 5.545,
                                                 MumAgeClass = c("AdultMum"),
                                                 DadAgeClass = c("AdultDad"),
                                                 Sex = unique(Sex),
                                                 Twin = c("0"),
                                                 MotherID = MotherID[1],
                                                 FatherID=FatherID[1],
                                                 BirthYear=BirthYear[1]))

lbs_pred_kids_mum <- predict.gam(m.lbs1_sex_final, new_data_kids_mum, se.fit=T, type="link", exclude = c("s(MotherID)", "s(FatherID)", "s(BirthYear)"))
lbs_pred_kids_dad <- predict.gam(m.lbs1_sex_final, new_data_kids_dad, se.fit=T, type="link", exclude = c("s(MotherID)", "s(FatherID)", "s(BirthYear)"))


# Make predictions
ilink <- family(m.lbs1_sex_final)$linkinv
lbs_pred_kids_mum <- cbind(lbs_pred_kids_mum, new_data_kids_mum)
lbs_pred_kids_mum <- transform(lbs_pred_kids_mum, lwr_ci = ilink(fit - (1.96 * se.fit)),
                           upr_ci = ilink(fit + (1.96 * se.fit)),
                           fitted = ilink(fit))

# Make predictions
ilink <- family(m.lbs1_sex_final)$linkinv
lbs_pred_kids_dad <- cbind(lbs_pred_kids_dad, new_data_kids_dad)
lbs_pred_kids_dad <- transform(lbs_pred_kids_dad, lwr_ci = ilink(fit - (1.96 * se.fit)),
                                upr_ci = ilink(fit + (1.96 * se.fit)),
                                fitted = ilink(fit))


# Showing predictions upto mum ages 11 and dad ages 9
lbs_pred_kids_mum <- lbs_pred_kids_mum %>%
  filter(MumAge <= 11)

lbs_pred_kids_dad <- lbs_pred_kids_dad %>%
  filter(DadAge <= 9)


p1_mums_sex <- ggplot() +
  geom_jitter(data = lifetime_data_plot, aes(x = MumAge, y = LBS, group=Sex, colour=Sex), size = 1, alpha=0.3, shape=16, position = position_jitter(width = 0.4, height = 0.1, seed=123)) +
  geom_line(data= lbs_pred_kids_mum, aes(x = MumAge, y = fitted, group=Sex, colour=Sex)) +
  geom_ribbon(data= lbs_pred_kids_mum, aes(x = MumAge, y = fitted, ymin = lwr_ci, ymax = upr_ci, group=Sex, colour= Sex, fill=Sex), alpha = 0.3, colour=NA) +
  theme_cowplot() +
  theme(legend.title=element_blank()) +
  labs(y = "Offspring LBS", x="Maternal Age (in years)") +
  scale_y_continuous(breaks = seq(0,15,3), limits=c(0,17)) +
  scale_x_continuous(breaks = seq(1,11,1)) +
  scale_color_manual(values=met.brewer("Troy", 2), labels=c("Daughters", "Sons")) +
  scale_fill_manual(values=met.brewer("Troy", 2), labels=c("Daughters", "Sons"))
p1_mums_sex


p1_dads_sex <- ggplot() +
  geom_jitter(data = lifetime_data_plot, aes(x = DadAge, y = LBS, group=Sex, colour=Sex), size = 1, alpha=0.3, shape=16, position = position_jitter(width = 0.4, height = 0.1, seed=123)) +
  geom_line(data= lbs_pred_kids_dad, aes(x = DadAge, y = fitted, group=Sex, colour=Sex)) +
  geom_ribbon(data= lbs_pred_kids_dad, aes(x = DadAge, y = fitted, ymin = lwr_ci, ymax = upr_ci, group=Sex, colour= Sex, fill=Sex), alpha = 0.3, colour=NA) +
  theme_cowplot() +
  theme(legend.title=element_blank()) +
  labs(y = "Offspring LBS", x="Paternal Age (in years)") +
  scale_y_continuous(breaks = seq(0,15,3), limits=c(0,17)) +
  scale_x_continuous(breaks = seq(1,9,1)) +
  scale_color_manual(values=met.brewer("Troy", 2), labels=c("Daughters", "Sons")) +
  scale_fill_manual(values=met.brewer("Troy", 2), labels=c("Daughters", "Sons"))
p1_dads_sex

###### Getting simulataneous 95% CIs #########
## simultaneous interval for smooth of MumAge:Sex
# Final Sex model 
fd1 <- derivatives(m.lbs1_sex_final, term = "s(MumAge):SexMale", type="central", interval=c("simultaneous"))
fd1$Signif <- fd1$.derivative
fd1 <- fd1 %>%
 dplyr::mutate(Signif = case_when(is.na(.upper_ci * .lower_ci) ~ NA_real_,
                                        .upper_ci * .lower_ci >=0 ~ Signif,
                                         TRUE ~ NA_real_))
fd1 <- fd1 %>%
  mutate(group = case_when(Signif > 0 ~ "positive",
                           Signif < 0 ~ "negative",
                           TRUE ~ NA_character_))

fd1 <- fd1 %>%
  filter(MumAge <= 11)


p2_fd <- ggplot(fd1, aes(x = MumAge, y = .derivative)) +
  geom_line() +
  geom_line(mapping = aes(y = .upper_ci), lty = "dashed") +
  geom_line(mapping = aes(y = .lower_ci), lty = "dashed") +
  geom_line(mapping = aes(y = Signif, group=group,  colour=group), lwd = 1.3) +
  theme_cowplot() +
  ylab("First derivative") +
  xlab("Maternal Age (in years)") +
  scale_x_continuous(breaks=seq(1,11,2), limits=c(1,11)) +
  # ggtitle("Pointwise in males") +
  geom_hline(yintercept=0, linetype="dotted", color = "purple") +
  scale_y_continuous(breaks=seq(from=-1,to=1,by=0.4), limits = c(-1, 1)) +
  scale_colour_discrete(na.translate = F, labels=c("Negative", "Positive")) + 
  theme(legend.title=element_blank()) +
  labs(subtitle="LBS of sons")
p2_fd


## simultaneous interval for smooth of x2
fd2 <- derivatives(m.lbs1_sex_final, term = "s(MumAge):SexFemale", type="central", interval=c("simultaneous"))
fd2$Signif <- fd2$.derivative
fd2 <- fd2 %>%
  mutate(Signif = case_when(is.na(.upper_ci * .lower_ci) ~ NA_real_,
                            .upper_ci * .lower_ci >=0 ~ Signif,
                            TRUE ~ NA_real_))
fd2 <- fd2 %>%
  mutate(group = case_when(Signif > 0 ~ "positive",
                           Signif < 0 ~ "negative",
                           TRUE ~ NA_character_))
fd2 <- fd2 %>%
  mutate(MumAge <= 11)

p3_fd <- ggplot(fd2, aes(x = MumAge, y = .derivative)) +
  geom_line() +
  geom_line(mapping = aes(y = .upper_ci), lty = "dashed") +
  geom_line(mapping = aes(y = .lower_ci), lty = "dashed") +
  geom_line(mapping = aes(y = Signif), lwd = 1.3) +
  theme_cowplot() +
  ylab("First derivative") +
  xlab("Maternal Age (in years)") +
  scale_x_continuous(breaks=seq(1,11,2), limits=c(1,11)) +
  scale_y_continuous(breaks=seq(from=-1,to=1,by=0.4), limits = c(-1, 1)) +
  geom_hline(yintercept=0, linetype="dotted", color = "purple") +
  scale_colour_discrete(na.translate = F) + 
  theme(legend.title=element_blank()) +
  labs(subtitle = "LBS of daughters")
p3_fd


plots_final <- (p1_mums_sex + p1_dads_sex + plot_layout(guides = "collect", axis_titles = "collect"))/(p3_fd + p2_fd + plot_layout(axis_titles = "collect")) + plot_annotation(tag_levels = 'A')
plots_final

# ggsave("./Figures/Aug24_PAE_LBS_Fig2.tiff",  plots_final, dpi=300, height = 6, width= 8)


