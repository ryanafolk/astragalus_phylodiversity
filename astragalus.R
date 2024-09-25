########################
# Load data sets
########################

library(dplyr)
library(readr)
# Load environmental data; starting with species richness to ensure as many pixels are recovered as possible
env <- list.files(path="./environmental_data", pattern = "richness", full.names = TRUE) %>% lapply(read_csv) %>% bind_rows 
env <- data.frame(env)
env[env == -9999] <- NA
env$x <- round(env$x, digit = 1)
env$y <- round(env$y, digit = 1)
env %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> env # soil type excluded at this step (not summarized correctly in source file)
env <- as.data.frame(env)

# Add RPD randomizations
rand_RPD <- read.csv("./Astragalus_50km_csvsToShare/rand_RPD.csv")
rand_RPD$x <- round(rand_RPD$x, digit = 1)
rand_RPD$y <- round(rand_RPD$y, digit = 1)
rand_RPD %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> rand_RPD
rand_RPD <- as.data.frame(rand_RPD)
# Add significance column
rand_RPD$RPD_significance <- as.factor(ifelse(rand_RPD$value < 0.05, "Low", ifelse(rand_RPD$value > 0.95, "High", "NS")))
rand_RPD$value <- NULL

combined <- merge(rand_RPD, env, by = c("x", "y"))

# Add RPD
RPD <- read.csv("./Astragalus_50km_csvsToShare/RPD.csv")
names(RPD) <- c("x", "y", "RPD")
RPD$x <- round(RPD$x, digit = 1)
RPD$y <- round(RPD$y, digit = 1)
RPD %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> RPD
RPD <- as.data.frame(RPD)

combined <- merge(combined, RPD, by = c("x", "y"))

# Add CANAPE
CANAPE <- read.csv("./Astragalus_50km_csvsToShare/CANAPE.csv")
names(CANAPE) <- c("x", "y", "CANAPE")
CANAPE$x <- round(CANAPE$x, digit = 1)
CANAPE$y <- round(CANAPE$y, digit = 1)
CANAPE %>% group_by(x, y) %>% summarize_if(is.character, max) -> CANAPE
CANAPE <- as.data.frame(CANAPE)
CANAPE$CANAPE <- as.factor(CANAPE$CANAPE)

combined <- merge(combined, CANAPE, by = c("x", "y"))

## Add SR
SR <- read.csv("./Astragalus_50km_csvsToShare/richness.csv")
names(SR) <- c("x", "y", "SR")
SR$x <- round(SR$x, digit = 1)
SR$y <- round(SR$y, digit = 1)
SR %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> SR
SR <- as.data.frame(SR)

combined <- merge(combined, SR, by = c("x", "y"))

## Add regionalization

region <- read.csv("./Astragalus_phylogenetic_clustering/clusterOutputNexus_COLOUR.csv")
names(region) <- c("x", "y", "region")
region$x <- round(region$x, digit = 1)
region$y <- round(region$y, digit = 1)
region %>% group_by(x, y) %>% summarize_if(is.numeric, mean, na.rm = TRUE) -> region
region <- as.data.frame(region)
region[region$region == 4294967296, ] <- NA # missing data value

combined <- merge(combined, region, by = c("x", "y"))


# Normalize entire data frame
combined.temp <- combined
combined.scaled <- rapply(combined.temp, scale, c("numeric","integer"), how="replace")
combined.scaled$SR <- combined$SR
combined.scaled <- as.data.frame(combined.scaled)
combined.scaled$y <- combined.temp$y
combined.scaled$x <- combined.temp$x


########################
# Model
########################


# SR model
linear_model_complex <- lm(SR ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, data = combined.scaled)
vif(linear_model_complex)
linear_model_complex <- lm(SR ~ BIOCLIM_3 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, data = combined.scaled)
sort(abs(linear_model_complex$coefficients))
# Top 5 predictors by GLM normalized coefficient
linear_model_simple <- lm(SR ~ BIOCLIM_3 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced, data = combined.scaled)
library(lme4)
mixed_model_complex <- lmer(SR ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
library(car)
vif(mixed_model_complex)
# Reassign model to object after removing VIF > 4
mixed_model_complex <- lmer(SR ~ BIOCLIM_3 + BIOCLIM_12 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
sort(abs(fixef(mixed_model_complex)))
# Top 5 predictors by LMM normalized coefficient
mixed_model_simple <- lmer(SR ~ BIOCLIM_3 + BIOCLIM_12 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
mixed_model_noenvironment <- lmer(SR ~ (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)

AIC(linear_model_complex)
AIC(linear_model_simple)
AIC(mixed_model_complex)
AIC(mixed_model_simple)
AIC(mixed_model_noenvironment)
# Complex mixed model favored

summary(mixed_model_complex)




library(MuMIn)
r.squaredGLMM(mixed_model_complex)



# RPD model
linear_model_complex <- lm(RPD ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, data = combined.scaled)
vif(linear_model_complex)
linear_model_complex <- lm(RPD ~ BIOCLIM_3 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, data = combined.scaled)
sort(abs(linear_model_complex$coefficients))
# Top 5 predictors by GLM normalized coefficient
linear_model_simple <- lm(RPD ~ BIOCLIM_3 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced, data = combined.scaled)
mixed_model_complex <- lmer(RPD ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
vif(mixed_model_complex)
# Reassign model to object after removing VIF > 4
mixed_model_complex <- lmer(RPD ~ BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
sort(abs(fixef(mixed_model_complex)))
# Top 5 predictors by LMM normalized coefficient
mixed_model_simple <- lmer(RPD ~ BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)
mixed_model_noenvironment <- lmer(RPD ~ (1 | y) + (1 | x), na.action = na.omit, data = combined.scaled)

AIC(linear_model_complex)
AIC(linear_model_simple)
AIC(mixed_model_complex)
AIC(mixed_model_simple)
AIC(mixed_model_noenvironment)
# Complex mixed model favored

summary(mixed_model_complex)




library(MuMIn)
r.squaredGLMM(mixed_model_complex)


# CANAPE significance model
combined.scaled$CANAPE_significant <- combined.scaled$CANAPE
levels(combined.scaled$CANAPE_significant)[levels(combined.scaled$CANAPE_significant)=="Neo"] <-"Sig"
levels(combined.scaled$CANAPE_significant)[levels(combined.scaled$CANAPE_significant)=="Mixed"] <-"Sig"
levels(combined.scaled$CANAPE_significant)[levels(combined.scaled$CANAPE_significant)=="Paleo"] <-"Sig"
levels(combined.scaled$CANAPE_significant)
# Reduce factor to 2 levels; otherwise glm guesses and does this silently
# Top 5 variables in simple models chosen based on coefficients in complex logit model
logit_complex <- glm(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, family = binomial(link='logit'), data = combined.scaled)
vif(logit_complex)
logit_complex <- glm(CANAPE_significant ~ BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_12 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, family = binomial(link='logit'), data = combined.scaled)
sort(abs(logit_complex$coefficients))
logit_simple <- glm(CANAPE_significant ~ BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_12 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced, family = binomial(link='logit'), data = combined.scaled)
# Be warned, these logit models will take about 1 minute and 0.5 GB RAM.
mixed_model_complex <- glmer(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
vif(mixed_model_complex)
mixed_model_complex <- glmer(CANAPE_significant ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
sort(abs(fixef(mixed_model_complex)))
mixed_model_simple <- glmer(CANAPE_significant ~ BIOCLIM_3 + BIOCLIM_17 + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 2e5))) # Stronger likelihood search options per warnings + documentation
mixed_model_noenvironment <- glmer(CANAPE_significant ~ (1 | y) + (1 | x), family=binomial(link='logit'), na.action = na.omit, data = combined.scaled, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun = 2e5)))

AIC(logit_complex)
AIC(logit_simple)
AIC(mixed_model_complex)
AIC(mixed_model_simple)
AIC(mixed_model_noenvironment)

# Mixed model simple favored
summary(mixed_model_simple)

library(MuMIn)
r.squaredGLMM(mixed_model_simple)



########################
# Some plots
########################

library(ggplot2)
levels(combined.scaled$RPD_significance)[levels(combined.scaled$RPD_significance)=="NS"] <- NA # NS excluded by the metrics
levels(combined$RPD_significance)[levels(combined$RPD_significance)=="NS"] <- NA # NS excluded by the metrics

# Violin plots


ggplot(combined, aes(x = RPD_significance, y = BIOCLIM_3/10, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Isothermality vs. RPD significance", x="RPD significance", y = "BIO3 (isothermality)")
ggplot(combined, aes(x = RPD_significance, y = BIOCLIM_7/10, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Temperature annual range vs. RPD significance", x="RPD significance", y = "BIO7 (temperature annual range in °C)")
ggplot(combined, aes(x = RPD_significance, y = BIOCLIM_17, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Precipitation in dry season \nvs. RPD significance", x="RPD significance", y = "BIO17 (precipitation\n of driest quarter; mm)")
ggplot(combined, aes(x = RPD_significance, y = BIOCLIM_1/10, fill = RPD_significance)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Mean annual temperature vs. RPD significance", x="RPD significance", y = "BIO1 (mean annual temperature in in °C")

# Outliers trimmed if too far out
ggplot(combined, aes(x = CANAPE, y = BIOCLIM_3/10, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Isothermality vs. CANAPE significance", x="CANAPE significance", y = "BIO3 (isothermality)")
ggplot(combined, aes(x = CANAPE, y = BIOCLIM_17, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Precipitation of dry season \nvs. CANAPE significance", x="CANAPE significance", y = "BIO17 (precipitation\n of driest quarter; mm)") + ylim(quantile(combined$BIOCLIM_17, 0.025, na.rm = TRUE), quantile(combined$BIOCLIM_17, 0.975, na.rm = TRUE))
ggplot(combined, aes(x = CANAPE, y = ISRICSOILGRIDS_new_average_nitrogen_reduced, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="Nitrogen content vs. CANAPE significance", x="CANAPE significance", y = "Nitrogen content (cg/kg)") + ylim(quantile(combined$ISRICSOILGRIDS_new_average_nitrogen_reduced, 0.025, na.rm = TRUE), quantile(combined$ISRICSOILGRIDS_new_average_nitrogen_reduced, 0.975, na.rm = TRUE))
ggplot(combined, aes(x = CANAPE, y = ISRICSOILGRIDS_new_average_phx10percent_reduced/10, fill = CANAPE)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + labs(title="pH vs. CANAPE significance", x="CANAPE significance", y = "pH") + ylim(quantile(combined$ISRICSOILGRIDS_new_average_phx10percent_reduced/10, 0.025, na.rm = TRUE), quantile(combined$ISRICSOILGRIDS_new_average_phx10percent_reduced/10, 0.975, na.rm = TRUE))


# Hex plots
# Species richness response
ggplot(combined[combined$SR > 2, ], aes(x = BIOCLIM_3/10, y = SR)) + geom_hex(bins = 20) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Isothermality vs. species richness", x="BIO3 (isothermality)", y = "SR")
ggplot(combined, aes(x = BIOCLIM_1/10, y = SR)) + geom_hex(bins = 18) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Mean annual temperature (°C) \nvs. species richness", x="Aridity index", y = "SR")

# Latitude
# SR vs. latitude
ggplot(combined, aes(x = y, y = SR)) + geom_hex(bins = 18) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Species richness vs. latitude", y="Species richness", x = "Latitude")
# SR vs. latitude, north hemisphere
ggplot(combined[combined$y > 0, ], aes(x = SR, y = y)) + geom_hex(bins = 18) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Species richness vs. latitude, \nNorthern Hemisphere", y="Latitude", x = "Species richness")
# SR vs. latitude, south hemisphere
ggplot(combined[combined$y < 0, ], aes(x = SR, y = y)) + geom_hex(bins = 8) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="Species richness vs. latitude, \nSouthern Hemisphere", y="Latitude", x = "Species richness")
# SR vs. latitude, east hemisphere, logarithmic
ggplot(combined[combined$x > -27, ], aes(x = SR, y = y)) + geom_hex(bins = 18) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ log(x)) + labs(title="Species richness vs. latitude, \nEastern Hemisphere", y="Latitude", x = "Species richness")
# SR vs. latitude, north america, logarithmic
ggplot(combined[combined$x < -27 & combined$y > 0, ], aes(x = y, y = SR)) + geom_hex(bins = 15) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ log(x)) + labs(title="Species richness vs. latitude, \nNorth America", y="Species richness", x = "Latitude")


# RPD vs. latitude
# RPD vs. latitude
ggplot(combined, aes(x = y, y = RPD)) + geom_hex(bins = 30) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="RPD vs. latitude", y="RPD", x = "Latitude")
# RPD vs. latitude, north hemisphere
ggplot(combined[combined$y > 0, ], aes(x = RPD, y = y)) + geom_hex(bins = 18) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="RPD vs. latitude, \nNorthern Hemisphere", y="RPD", x = "Species richness")
# RPD vs. latitude, south hemisphere
ggplot(combined[combined$y < 0, ], aes(x = RPD, y = y)) + geom_hex(bins = 8) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ x) + labs(title="RPD vs. latitude, \nSouthern Hemisphere", y="RPD", x = "Species richness")
# RPD vs. latitude, east hemisphere, logarithmic
ggplot(combined[combined$x > -27, ], aes(x = RPD, y = y)) + geom_hex(bins = 18) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ log(x)) + labs(title="RPD vs. latitude, \nEastern Hemisphere", y="Latitude", x = "RPD") + xlim(c(0,1))
# RPD vs. latitude, north america, logarithmic
ggplot(combined[combined$x < -27 & combined$y > 0, ], aes(x = y, y = RPD)) + geom_hex(bins = 15) + scale_fill_continuous(type = "viridis") + theme_bw() + geom_smooth(method='lm', formula = y ~ log(x)) + labs(title="RPD vs. latitude, \nNorth America", y="RPD", x = "Latitude")




combined.scaled.reduced <- combined.scaled
combined.scaled.reduced[is.na(combined.scaled.reduced)] <- 0

library(MASS)
discriminant <- lda(region ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, data = combined.scaled.reduced, na.action="na.omit")


# Classification success
discriminant.jackknife <- lda(region ~ aridity_index_UNEP + BIOCLIM_1 + BIOCLIM_3 + BIOCLIM_4 + BIOCLIM_12 + BIOCLIM_7 + BIOCLIM_17 + GTOPO30_ELEVATION + ISRICSOILGRIDS_new_average_nitrogen_reduced + ISRICSOILGRIDS_new_average_phx10percent_reduced + ISRICSOILGRIDS_new_average_soilorganiccarboncontent_reduced + ISRICSOILGRIDS_new_average_sandpercent_reduced + ISRICSOILGRIDS_new_average_coarsefragmentpercent_reduced + LandCover_1_Needleleaf + LandCover_3_Deciduousbroadleaf, data = combined.scaled.reduced, CV = TRUE)
ct <- table(combined.scaled$region, discriminant.jackknife$class)
sum(diag(prop.table(ct)))



# Predict species by the discriminant function
discriminant.prediction <- predict(discriminant)

# Create dataframe for plotting
plotdata <- data.frame(type = combined$region, lda = discriminant.prediction$x)
plotdata$type <- as.factor(plotdata$type)

plotdata.americas <- plotdata[plotdata$type %in% c("3","4","5","7"),]

library(ggplot2)
ggplot(plotdata.americas) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 1) + scale_color_brewer(palette="BrBG")

plotdata.asia <- plotdata[plotdata$type %in% c("5","6","10","11","12"),]
ggplot(plotdata.asia) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 1) + scale_color_brewer(palette="BrBG")

plotdata.europe <- plotdata[plotdata$type %in% c("1","2","8","9"),]
ggplot(plotdata.europe) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 1) + scale_color_brewer(palette="BrBG")


