#-------------------------------------------------------------------------------
# Project: ITS analysis for impact of COVID-19 on antibiotic prescribing
# Author: Olivia Bryant
# Date created: 26/09/2022
# Date updated: 13/10/2022
# Notes: ITS analysis of COVID-19 on antibiotic prescriptions
#-------------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(tsModel)
library(lubridate)
library(scales)
library(MASS)
library(Epi)
library(tseries)
library(ggtext)

# load csv
prescriptions_data <- read.csv("ITS_dataset_toJul22.csv")

startdate = min(prescriptions_data$date)

# format data with lockdown markers
prescriptions_data <- prescriptions_data %>%
  mutate(studymonth = month(date),
         date = as.Date(date, format="%d/%m/%Y"),
         Appointments = as.numeric(gsub(",","",prescriptions_data$Appointments)))

study_data <- prescriptions_data %>%
         mutate(
           time = row_number(),
           intervention1 = ifelse(prescriptions_data$date < as.Date("01/03/2020", format = "%d/%m/%Y"), 0, 1),
           intervention2 = ifelse(prescriptions_data$date < as.Date("01/07/2021", format = "%d/%m/%Y"), 0, 1))


ldn1_centre <- study_data$time[min(which(study_data$intervention1 == 1))]
ldn2_centre <- study_data$time[min(which(study_data$intervention2 == 1))]

study_data <- study_data %>%
  mutate(time_after_inter1  = (time - ldn1_centre) + 1,
         time_after_inter2 = (time - ldn2_centre) + 1)

study_data$time_after_inter1[study_data$time_after_inter1 < 0] <- 0
study_data$time_after_inter2[study_data$time_after_inter2 < 0] <- 0

write.csv(study_data, file = "ITS_processed_dataset.csv")
#-------------------------------------------------------------------------------
# Negative Binomial ITS for all antibiotics per appointment
#-------------------------------------------------------------------------------
data = study_data
outcome_measure = study_data$Appointments
drug = study_data$Item_all.antibiotics

negbinom_model1 <- glm.nb(drug ~ offset(log(outcome_measure)) + time + intervention1 +
                            time_after_inter1 + intervention2 + time_after_inter2
                          + as.factor(studymonth), 
                          data = filter(data, !is.na(intervention1)))

jpeg(file="pacf_plots/pcaf_model1")
pacf(residuals(negbinom_model1, type = "deviance"))
dev.off()

negbinom_lagres <- lag(residuals(negbinom_model1)) %>% as.numeric()
res1 <- residuals(negbinom_model1, type = "deviance")

# fit full model with lagged residuals 
negbinom_model2 <- glm.nb(drug ~ offset(log(outcome_measure)) + time + intervention1 +
                              time_after_inter1 + intervention2 + time_after_inter2
                            + as.factor(studymonth) + negbinom_lagres, 
                            data = filter(data, !is.na(intervention1)))

jpeg(file="pacf_plots/pcaf_model2")
pacf(residuals(negbinom_model2, type = "deviance"))
dev.off()

negbinom_lagres_timing <- bind_cols("time" = data$time,
                                    "negbinom_lagres" = negbinom_lagres)

outcome_pred_zeroed <- data %>%
  left_join(negbinom_lagres_timing, by = "time") %>%
  mutate_at("negbinom_lagres", ~(. = 0))

outcome_pred_nointervention <- data %>%
  mutate_at("intervention1", ~(. = 0)) %>%
  mutate_at("intervention2", ~(. = 0)) %>%
  mutate_at("time_after_inter1", ~(. = 0)) %>%
  mutate_at("time_after_inter2", ~(. = 0))

pred_noLockdown <- predict(negbinom_model2, 
                           newdata = outcome_pred_nointervention, 
                           se.fit = TRUE, interval="confidence") 

pred_twointerruption <- predict(negbinom_model2,
                                newdata = outcome_pred_zeroed,
                                se.fit = TRUE, interval="confidence") 
# save predictions from model
preds <- pred_twointerruption$fit
stbp_preds <- pred_twointerruption$se.fit

# save predictions from model
preds_0 <- pred_noLockdown$fit
stbp_preds_0 <- pred_noLockdown$se.fit

df_se <- bind_cols(stbp_preds = stbp_preds, 
                   pred = preds,
                   preds_0 = preds_0,
                   stbp_preds_0 = stbp_preds_0,
                   denom = outcome_measure) %>%
  mutate(
    #CIs
    upp = pred + (1.96 * stbp_preds),
    low = pred - (1.96 * stbp_preds),
    upp_noLdn = preds_0 + (1.96 * stbp_preds_0),
    low0_noLdn = preds_0 - (1.96 * stbp_preds_0),
    
    predicted_vals = (exp(pred) * 100) / denom,
    predicted_vals_no = (exp(preds_0) * 100) / denom,
    
    # transform CIs
    uci = (exp(upp) * 100) / denom,
    lci = (exp(low) * 100) / denom,
    uci_0 = (exp(upp_noLdn) * 100) / denom,
    lci_0 = (exp(low0_noLdn) * 100) / denom
  )

outcome_plot <- bind_cols(outcome_pred_nointervention, df_se)

parameter_estimates <- as.data.frame(ci.exp(negbinom_model2))
vals_to_print <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates))

write.csv(vals_to_print, file = "estimates/all_antibiotics_appointments.csv")

bkg_colour <- "white"
plot <- ggplot(study_data, aes(x = date, y = (Item_all.antibiotics/Appointments) * 100)) +
  geom_line(col = "gray60") +
  geom_line(data = outcome_plot, aes(y = predicted_vals, color = "ITS model"), size = 1) +
  geom_ribbon(data = outcome_plot, aes(ymin = lci, ymax = uci), fill = alpha(4, 0.4), lty = 0) +
  geom_line(data = filter(outcome_plot, date >= as.Date("01/03/2020", format="%d/%m/%Y")), 
            aes(y = predicted_vals_no, color = "Expected trend based \n on pre-COVID-19 rates"), 
            col = 2, lty = 2, size = 1) + 
  labs(y = "Items per 100 Appointments", title = "Impact of COVID-19 on Antibiotic Prescribing Rates in Primary Care in England") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), linetype = "dashed") + 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), linetype = "dashed") + 
  scale_color_manual(values = colours)

ggsave(filename = "plots/all_antibiotics_appointments",width = 8, height = 5, dpi = 150, units = "in", device='png')

plot2 <- ggplot(study_data, aes(x = date, y = (Item_all.antibiotics/Appointments) * 100)) +
  geom_line(col = "gray60") +
  geom_line(data = outcome_plot, aes(y = predicted_vals), col = 4, lty = 2) +
  geom_ribbon(data = outcome_plot, aes(ymin = lci, ymax = uci), fill = alpha(4, 0.4), lty = 0) +
  labs(y = "Items per 100 Appointments") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        legend.position = "top",
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), linetype = "dashed")+ 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), linetype = "dashed")

ggsave(filename = "plots/observed_all_antibiotics_appointments",width = 8, height = 5, dpi = 150, units = "in", device='png')

plot3 <- ggplot(study_data, aes(x = date, y = (Item_all.antibiotics/Appointments) * 100)) +
  geom_line(col = "gray60") +
  geom_line(data = outcome_plot, aes(y = predicted_vals_no), col = 2, lty = 2) +
  geom_ribbon(data = outcome_plot, aes(ymin = lci_0, ymax = uci_0), fill = alpha(2, 0.4), lty = 2) +
  labs(y = "Items per 100 Appointments") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        legend.position = "top",
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), linetype = "dashed")+ 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), linetype = "dashed")

ggsave(filename = "plots/counterfactual_all_antibiotics_appointments",width = 8, height = 5, dpi = 150, units = "in", device='png')


#-------------------------------------------------------------------------------
# Negative Binomial ITS for all antibiotics per appointment: APRIL 2020 LAG
#-------------------------------------------------------------------------------
data = study_data
outcome_measure = study_data$Appointments
drug = study_data$Item_all.antibiotics

april_index <- which(data$date == as.Date("01/04/2020", format = "%d/%m/%Y"))
data <- data[-april_index,]
outcome_measure <- outcome_measure[-april_index]
drug <- drug[-april_index]

negbinom_model1 <- glm.nb(drug ~ offset(log(outcome_measure)) + time + intervention1 +
                            time_after_inter1 + intervention2 + time_after_inter2
                          + as.factor(studymonth), 
                          data = filter(data, !is.na(intervention1)))

jpeg(file = "pacf_plots/APRIL20_EXCL_pacf_model1")
pacf(residuals(negbinom_model1, type = "deviance"))
dev.off()

negbinom_lagres <- lag(residuals(negbinom_model1)) %>% as.numeric()
res1 <- residuals(negbinom_model1, type = "deviance")

# fit full model with lagged residuals 
negbinom_model2 <- glm.nb(drug ~ offset(log(outcome_measure)) + time + intervention1 +
                            time_after_inter1 + intervention2 + time_after_inter2
                          + as.factor(studymonth) + negbinom_lagres, 
                          data = filter(data, !is.na(intervention1)))

jpeg(file = "pacf_plots/APRIL20_EXCL_pacf_model2")
pacf(residuals(negbinom_model2, type = "deviance"))
dev.off()

negbinom_lagres_timing <- bind_cols("time" = data$time,
                                    "negbinom_lagres" = negbinom_lagres)

outcome_pred_zeroed <- data %>%
  left_join(negbinom_lagres_timing, by = "time") %>%
  mutate_at("negbinom_lagres", ~(. = 0))

outcome_pred_nointervention <- data %>%
  mutate_at("intervention1", ~(. = 0)) %>%
  mutate_at("intervention2", ~(. = 0)) %>%
  mutate_at("time_after_inter1", ~(. = 0)) %>%
  mutate_at("time_after_inter2", ~(. = 0))

pred_noLockdown <- predict(negbinom_model2, 
                           newdata = outcome_pred_nointervention, 
                           se.fit = TRUE, interval="confidence") 

pred_twointerruption <- predict(negbinom_model2,
                                newdata = outcome_pred_zeroed,
                                se.fit = TRUE, interval="confidence") 
# save predictions from model
preds <- pred_twointerruption$fit
stbp_preds <- pred_twointerruption$se.fit

# save predictions from model
preds_0 <- pred_noLockdown$fit
stbp_preds_0 <- pred_noLockdown$se.fit


df_se <- bind_cols(stbp_preds = stbp_preds, 
                   pred = preds,
                   preds_0 = preds_0,
                   stbp_preds_0 = stbp_preds_0,
                   denom = outcome_measure) %>%
  mutate(
    #CIs
    upp = pred + (1.96 * stbp_preds),
    low = pred - (1.96 * stbp_preds),
    upp_noLdn = preds_0 + (1.96 * stbp_preds_0),
    low0_noLdn = preds_0 - (1.96 * stbp_preds_0),
    
    predicted_vals = (exp(pred) * 100) / denom,
    predicted_vals_no = (exp(preds_0) * 100) / denom,
    
    # transform CIs
    uci = (exp(upp) * 100) / denom,
    lci = (exp(low) * 100) / denom,
    uci_0 = (exp(upp_noLdn) * 100) / denom,
    lci_0 = (exp(low0_noLdn) * 100) / denom
  )

outcome_plot <- bind_cols(outcome_pred_nointervention, df_se)

parameter_estimates <- as.data.frame(ci.exp(negbinom_model2))
vals_to_print <- parameter_estimates %>%
  mutate(var = rownames(parameter_estimates))

write.csv(vals_to_print, file = "estimates/APRIL20_EXCL_all_antibiotics_appointments.csv")

outcome_plot <- outcome_plot %>%
  add_row(date = as.Date("01/04/2020", format = "%d/%m/%Y"))

colours <- c("Expected trend \nbased on pre-COVID-19 \nrates" = "red3", 
             "ITS model" = "cadetblue3", 
             "Data from OpenPrescribing" = "gray60")
bkg_colour <- "white"
plot <- ggplot(study_data, aes(x = date, y = (Item_all.antibiotics/Appointments) * 100)) +
  geom_line(col = "gray60") +
  geom_line(data = outcome_plot, aes(y = predicted_vals, color = "ITS model"), size = 1) +
  geom_ribbon(data = outcome_plot, aes(ymin = lci, ymax = uci), fill = alpha(4, 0.4), lty = 0) +
  geom_line(data = filter(outcome_plot, date >= as.Date("01/03/2020", format="%d/%m/%Y")), 
            aes(y = predicted_vals_no, color = "Expected trend \nbased on pre-COVID-19 \nrates"), 
            col = 2, lty = 2, size = 1) + 
  labs(y = "Items per 100 Appointments", title = "Impact of COVID-19 on Antibiotic Prescribing Rates in Primary Care in England") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), linetype = "dashed") + 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), linetype = "dashed") + 
  annotate("rect", # lag
                          fill = "gray", 
                          alpha = 0.4,
                          xmin = as.Date("2020-03-01"), 
                          xmax = as.Date("2020-05-01"),
                          ymin = -Inf, ymax = Inf) + 
  scale_color_manual(values = colours)

  
ggsave(filename = "plots/APRIL20_EXCL_appointments_all_antibiotics", width = 8, height = 5, dpi = 150, units = "in", device='png')

plot2 <- ggplot(study_data, aes(x = date, y = (Item_all.antibiotics/Appointments) * 100)) +
  geom_line(col = "gray60") +
  geom_line(data = outcome_plot, aes(y = predicted_vals), col = 4, lty = 2) +
  geom_ribbon(data = outcome_plot, aes(ymin = lci, ymax = uci), fill = alpha(4, 0.4), lty = 0) +
  labs(y = "Items per 100 Appointments") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        legend.position = "top",
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), linetype = "dashed")+ 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), linetype = "dashed") +
  annotate("rect", # lag
        fill = "gray", 
        alpha = 0.4,
        xmin = as.Date("2020-03-01"), 
        xmax = as.Date("2020-05-01"),
        ymin = -Inf, ymax = Inf)

ggsave(filename = "plots/APRIL20_EXCL_observation_model",width = 8, height = 5, dpi = 150, units = "in", device='png')

plot3 <- ggplot(study_data, aes(x = date, y = (Item_all.antibiotics/Appointments) * 100)) +
  geom_line(col = "gray60") +
  geom_line(data = outcome_plot, aes(y = predicted_vals_no), col = 2, lty = 2) +
  geom_ribbon(data = outcome_plot, aes(ymin = lci_0, ymax = uci_0), fill = alpha(2, 0.4), lty = 2) +
  labs(y = "Items per 100 Appointments") +
  theme_classic() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        legend.position = "top",
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), linetype = "dashed") + 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), linetype = "dashed") + 
  annotate("rect", # lag
                            fill = "gray", 
                            alpha = 0.4,
                            xmin = as.Date("2020-03-01"), 
                            xmax = as.Date("2020-05-01"),
                            ymin = -Inf, ymax = Inf)

ggsave(filename = "plots/APRIL20_EXCL_counterfactual", width = 8, height = 5, dpi = 150, units = "in", device='png')


#-------------------------------------------------------------------------------
# Plot type of GP appointment over time
#-------------------------------------------------------------------------------
appointments_data <- read.csv("GP_appointments_data.csv")

appointments_data <- appointments_data %>%
  mutate(total = Telephone + Face_to_face + Home_visit + video_online) %>%
  mutate(pc_telephone = (Telephone/total) * 100) %>%
  mutate(date = as.Date(Date, format = "%d/%m/%Y"))

ggplot(appointments_data, aes(x = date, y = pc_telephone)) +
  geom_line() +
  theme(axis.title = element_text(size = 12), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        legend.position = "top",
        plot.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.background = element_rect(fill = bkg_colour, colour =  NA),
        legend.background = element_rect(fill = bkg_colour, colour = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = bkg_colour, colour =  NA),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_line(size=.2, color=rgb(0,0,0,0.2)),
        panel.grid.minor.y = element_line(size=.2, color=rgb(0,0,0,0.2)) ,
        panel.grid.major.y = element_line(size=.2, color=rgb(0,0,0,0.3))) + 
  labs(y = "% of GP Appointments via Telephone", 
       caption ="Unknown appointment type excluded. Red lines represent the start/end of COVID-19 restrictions.",
       title = "Proportion of GP Appointments Conducted via Telephone in England") +
  scale_x_date(date_labels = "%b %y", breaks = breaks_pretty(10),labels = scales::label_date_short()) +
  geom_vline(xintercept = as.Date("01/03/2020", format = "%d/%m/%Y"), color="darkred", linetype = "dashed") + 
  geom_vline(xintercept = as.Date("01/07/2021", format = "%d/%m/%Y"), color="darkred", linetype = "dashed")
ggsave(filename = "plots/GP_appointment_type", width = 8, height = 5, dpi = 150, units = "in", device='png')
