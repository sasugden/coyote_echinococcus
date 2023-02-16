#### PREPARE WORKSPACE ####
library(vegan)
library(ggplot2)
library(MuMIn)
library(magrittr)
library(caret)

# Create a function that mean-centers and scales only predictor variables
scale.vars <- function(data_frame){
  for(i in c("age.cement", "anth.vol", "anth.dig.vol", "anth.indig.vol", "sm.vol",
             "new.pca", "old.pca", "ungulate.vol", "food.vol.washed",
             "diet.Shannon", 
             "prey.vol", "med.vol", "bird.vol", "insect.vol", "wood.veg.vol", "herb.veg.vol",
             "veg.vol", "fruit.vol", "d13C", "d15N")){
    data_frame[,i] <- c(scale(data_frame[,i]))
  }
  return(data_frame)
}

#### DATA IMPORT AND PROCESSING ####
# Import data
sample_data <- read.csv("~/Downloads/coyote_parasite_2022_09_25.csv")

# Clean data
sample_data <- subset(sample_data, cause.of.death !="decay") # Remove the four samples that had decayed before analysis
sample_data <- subset(sample_data, !(sampleID %in% c("C098", "C100"))) # Remove two samples with no PCR or qPCR results

# Convert dummy variables (0/1) to factor variables
sample_data$echino.worms.present <- factor(as.character(sample_data$echino.worms.present),
                                           levels=c(0,1))
sample_data$bioact <- factor(as.character(sample_data$bioact),
                             levels=c(0,1))
sample_data$empty <- factor(as.character(sample_data$empty),
                            levels=c(0,1))

# Cap worm counts at 10,000 and then log-transform
sample_data$ech.cap <- sample_data$ech.count
sample_data$ech.cap[sample_data$ech.cap > 10000] <- 10000
sample_data$log1.worm <- log(sample_data$ech.cap + 1)

# Divide age classes based on the median age
sample_data$age.class <- rep("young")
sample_data$age.class[sample_data$age.cement > 1.78] <- "old"
sample_data$age.class <- factor(sample_data$age.class, 
                                levels=c("young","old"))

# Create a grouping variable based on location, age, and infection status
sample_data$group1 <- rep(0)
sample_data$group1[sample_data$location=="urban" & sample_data$age.class=="young" & sample_data$bioact==1] <- "urb_yng_1"
sample_data$group1[sample_data$location=="urban" & sample_data$age.class=="young" & sample_data$bioact==0] <- "urb_yng_0"
sample_data$group1[sample_data$location=="urban" & sample_data$age.class=="old" & sample_data$bioact==1] <- "urb_old_1"
sample_data$group1[sample_data$location=="urban" & sample_data$age.class=="old" & sample_data$bioact==0] <- "urb_old_0"
sample_data$group1[sample_data$location=="rural" & sample_data$age.class=="young" & sample_data$bioact==1] <- "rur_yng_1"
sample_data$group1[sample_data$location=="rural" & sample_data$age.class=="young" & sample_data$bioact==0] <- "rur_yng_0"
sample_data$group1[sample_data$location=="rural" & sample_data$age.class=="old" & sample_data$bioact==1] <- "rur_old_1"
sample_data$group1[sample_data$location=="rural" & sample_data$age.class=="old" & sample_data$bioact==0] <- "rur_old_0"

sample_data$group1 <- factor(sample_data$group1, levels=c("urb_yng_1", "urb_yng_0", "urb_old_1", "urb_old_0",
                                                          "rur_yng_1", "rur_yng_0", "rur_old_1", "rur_old_0"))

sample_data %>%
  dplyr::group_by(group1) %>%
  dplyr::count()

# Create a grouping variable based only on location and age
sample_data$group2 <- rep(0)
sample_data$group2[sample_data$location=="urban" & sample_data$age.class=="young"] <- "urb_yng"
sample_data$group2[sample_data$location=="urban" & sample_data$age.class=="old"] <- "urb_old"
sample_data$group2[sample_data$location=="rural" & sample_data$age.class=="young"] <- "rur_yng"
sample_data$group2[sample_data$location=="rural" & sample_data$age.class=="old"] <- "rur_old"

sample_data %>%
  dplyr::group_by(group2) %>%
  dplyr::count()

# Calculate composite condition metric by performing PCA on body condition variables
condition <- sample_data[,c("mass","length","girth","kfi")]
condition <- as.data.frame(scale(condition))
condition.pca <- capscale(condition ~ 1)
sample_data$new.pca <- as.vector(scores(condition.pca, display="sites", choices=c(1)))

# Replace NA values with zero in the diet data
sample_data$med.vol[is.na(sample_data$med.vol)=="TRUE"] <- 0
sample_data$bird.vol[is.na(sample_data$bird.vol)=="TRUE"] <- 0

# Calculate prey volume and diet diversity
sample_data$prey.vol <- sample_data$sm.vol + sample_data$med.vol + sample_data$ungulate.vol + 
  sample_data$bird.vol + sample_data$insect.vol

sample_data$diet.Shannon <- diversity(sample_data[,c("anth.dig.vol", "anth.indig.vol",
                                                     "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
                                                     "veg.vol", "fruit.vol")], index="shannon")

# Add model weights to the metadata.
sample_data$weights.l <- rep(1)
sample_data$weights.l[sample_data$location=="urban"] <- 1.731707

# Differences in age, diet, and condition as a function of location: Levene's test followed by ANOVA ####
basic.stats <- data.frame(matrix(ncol=8, nrow=0))

for(i in c("new.pca", "kfi", "splenic.index",
           "food.vol.washed", "diet.Shannon",
           "anth.vol", "anth.dig.vol", "anth.indig.vol",
           "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
           "veg.vol", "herb.veg.vol", "wood.veg.vol", 
           "fruit.vol", "d13C", "d15N")){
  levene <- car::leveneTest(sample_data[,i] ~ sample_data$location)
  
  if(levene[1,3] < 0.05){
    test <- oneway.test(sample_data[,i] ~ location, sample_data, var.equal=FALSE)
  }
  
  if(levene[1,3] >= 0.05){
    test <- oneway.test(sample_data[,i] ~ location, sample_data, var.equal=TRUE)
  }
  
  basic.stats[i,1] <- i
  basic.stats[i,2] <- levene[1,1]
  basic.stats[i,3] <- levene[1,2]
  basic.stats[i,4] <- levene[1,3]
  
  basic.stats[i,5] <- test$statistic
  basic.stats[i,6] <- test$parameter[1]
  basic.stats[i,7] <- test$parameter[2]
  basic.stats[i,8] <- test$p.value
}

colnames(basic.stats) <- c("variable", "levene F", "levene df", "levene p",
                           "ANOVA F", "ANOVA num df", "ANOVA denom df", "ANOVA p")

write.csv(basic.stats, "~/Documents/Deanna_MS_data/basic_loc_stats.csv")

#### MODELS PREDICTING ECHINOCOCCUS MULTILOCULARIS INFECTION ####
## Summary statistics: mean values of diet/health variables for infected vs. uninfected coyotes ####
grp <- sample_data %>% dplyr::group_by(bioact) %>% dplyr::summarise_if(is.numeric, .funs=c(mean, sd))
grp <- as.data.frame(t(grp))
colnames(grp) <- c("mean_uninfected", "mean_infected")
grp$id <- rownames(grp)
grp <- tidyr::separate(data = grp, col = id, into = c("variable", "calc"), sep="_")
grp$calc[grp$calc=="fn1"] <- "mean"
grp$calc[grp$calc=="fn2"] <- "sd"

grp1 <- subset(grp, calc=="mean")
grp2 <- subset(grp, calc=="sd")
colnames(grp2) <- c("sd_uninfected", "sd_infected", "variable", "calc")
grp1$calc <- NULL
grp2$calc <- NULL

grp <- merge(grp1, grp2, by="variable", all=TRUE)
rm(grp1, grp2)

sample_data_scale <- scale.vars(sample_data)

## Univariate models ####
# (a) Logistic regression models predicting infection prevalence. ####
univ.logit <- data.frame(matrix(ncol=4, nrow=0))

for(i in c("location", "age.cement", "sex", "new.pca", "kfi", "splenic.index",
           "food.vol.washed", "diet.Shannon",
           "anth.vol", "anth.dig.vol", "anth.indig.vol",
           "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
           "veg.vol", "herb.veg.vol", "wood.veg.vol", 
           "fruit.vol", "d13C", "d15N")){
  
  m1 <- glm(bioact ~ sample_data_scale[,i],
            family="binomial"(link="logit"),
            weights=weights.l,
            sample_data_scale)
  m2 <- glm(bioact ~ 1, family="binomial"(link="logit"),
            weights = weights.l,
            sample_data_scale)
  
  test <- lrtest(m1, m2)
  
  univ.logit[i, 1] <- i
  univ.logit[i, 2] <- AIC(m1)
  univ.logit[i, 3] <- test[2,4]
  univ.logit[i, 4] <- test[2,3]
  univ.logit[i, 5] <- test[2,5]
}
colnames(univ.logit) <- c("variable", "AIC", "chisq", "df", "p")

# (b) Negative binomial models predicting infection intensity. ####
univ.nbinom <- data.frame(matrix(ncol=4, nrow=0))

for(i in c("location", "age.cement", "sex", "new.pca", "kfi", "splenic.index",
           "food.vol.washed", "diet.Shannon",
           "anth.vol", "anth.dig.vol", "anth.indig.vol",
           "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
           "veg.vol", "herb.veg.vol", "wood.veg.vol", 
           "fruit.vol", "d13C", "d15N")){
  
  m1 <- MASS::glm.nb(log1.worm ~ sample_data_scale[,i],
                     weights = weights.l,
                     sample_data_scale)
  m2 <- MASS::glm.nb(log1.worm ~ 1, 
                     weights = weights.l,
                     sample_data_scale)
  
  test <- lrtest(m1, m2)
  
  univ.nbinom[i, 1] <- i
  univ.nbinom[i, 2] <- AIC(m1)
  univ.nbinom[i, 3] <- test[2,4]
  univ.nbinom[i, 4] <- test[2,3]
  univ.nbinom[i, 5] <- test[2,5]
}
colnames(univ.nbinom) <- c("variable", "nbinom.AIC", "nbinom.chisq", "nbinom.df", "nbinom.p")

univ.final.results <- merge(grp, univ.logit, by="variable", all.x=FALSE, all.y=TRUE)
univ.final.results <- merge(univ.final.results, univ.nbinom, by="variable", all=TRUE)

write.csv(univ.final.results, "~/Documents/Deanna_MS_data/univariate.regressions.weighted.csv")

## Univariate supporting statistics: chi-squared, Wilcox, Spearman correlation ####
# chi-squared tests (prevalence vs. categorical variables)
chisq.test(sample_data$bioact, sample_data$location, correct = FALSE)
chisq.test(sample_data$bioact, sample_data$sex, correct = FALSE)

# Wilcox rank-sum tests (prevalence vs. numerical variables)
numeric_predictors <- c("age.cement", "new.pca", "splenic.index", "kfi",
                        "diet.Shannon", "food.vol.washed", 
                        "anth.vol", "anth.dig.vol", "anth.indig.vol",
                        "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
                        "veg.vol", "fruit.vol", "d13C", "d15N")

wilcox_test <- data.frame(matrix(ncol=3))

for(i in c(1:length(numeric_predictors))){
  test <- wilcox.test(sample_data[,numeric_predictors[i]] ~ bioact, sample_data)
  
  wilcox_test[i,1] <- numeric_predictors[i]
  wilcox_test[i,2] <- test$statistic
  wilcox_test[i,3] <- test$p.value
  
  rm(test)
}

colnames(wilcox_test) <- c("variable", "W", "p")
write.csv(wilcox_test, "~/Documents/Deanna_MS_data/wilcox_outputs.csv")

# Wilcoxon rank-sum tests (intensity vs. categorical variables)
wilcox.test(ech.cap ~ location, sample_data)
wilcox.test(ech.cap ~ sex, sample_data)

# Spearman's rank correlation (intensity vs. continuous variables)

spearman_data <- sample_data[,c("ech.cap", "log1.worm", numeric_predictors)]

spearman_corr <- Hmisc::rcorr(as.matrix(spearman_data), type="spearman")

for(i in c(1:3)){
  spearman_corr[[i]] <- spearman_corr[[i]][,c(1:2)]
  spearman_corr[[i]] <- subset(spearman_corr[[i]], rownames(spearman_corr[[i]]) %in% numeric_predictors)
}

spearman_corr <- merge(spearman_corr$r, spearman_corr$P, by=0, all=TRUE)

write.csv(spearman_corr, "~/Documents/Deanna_MS_data/univariate.spearman.corrs.csv")

## Two and three-way interaction models ####
# (a) Models interacting with age (two-way) ####
age.interact <- data.frame(matrix(ncol=21, nrow=0))

for(i in c("location", "sex", "new.pca", "kfi", "splenic.index",
           "food.vol.washed", "diet.Shannon",
           "anth.vol", "anth.dig.vol", "anth.indig.vol",
           "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
           "veg.vol", "herb.veg.vol", "wood.veg.vol", 
           "fruit.vol", "d13C", "d15N")){
  
  # For infection status:
  # Create logistic regressions for (1) age * focal variable, (2) age + focal variable, (3) age only, and (4) null model
  logit1 <- glm(bioact ~ age.cement * sample_data_scale[,i],
            family="binomial"(link="logit"),
            sample_data_scale, weights = weights.l)
  logit2 <- glm(bioact ~ age.cement + sample_data_scale[,i], family="binomial"(link="logit"), sample_data_scale, weights = weights.l)
  logit3 <- glm(bioact ~ age.cement, family="binomial"(link="logit"), sample_data_scale, weights = weights.l)
  logitnull <- glm(bioact ~ 1, family="binomial"(link="logit"), sample_data_scale, weights = weights.l)
  
  logit.test1 <- lrtest(logit1, logit2)
  logit.test2 <- lrtest(logit1, logit3)
  logit.test3 <- lrtest(logit1, logitnull)
  
  # For infection intensity:
  # Create nbinom regressions for (1) age * focal variable, (2) age + focal variable, (3) age only, and (4) null model 
  nb1 <- MASS::glm.nb(log1.worm ~ age.cement * sample_data_scale[,i], sample_data_scale, weights = weights.l)
  nb2 <- MASS::glm.nb(log1.worm ~ age.cement + sample_data_scale[,i], sample_data_scale, weights = weights.l)
  nb3 <- MASS::glm.nb(log1.worm ~ age.cement, sample_data_scale, weights = weights.l)
  nbnull <- MASS::glm.nb(log1.worm ~ 1, sample_data_scale, weights = weights.l)
  
  nb.test1 <- lrtest(nb1, nb2)
  nb.test2 <- lrtest(nb1, nb3)
  nb.test3 <- lrtest(nb1, nbnull)
  
  # Store the outputs of the likelihood ratio tests in a data frame
  age.interact[i, 1] <- i
  age.interact[i, 2] <- AIC(logit1)
  age.interact[i, 3] <- logit.test1[2,4]
  age.interact[i, 4] <- logit.test1[2,3]
  age.interact[i, 5] <- logit.test1[2,5]
  
  age.interact[i, 6] <- logit.test2[2,4]
  age.interact[i, 7] <- logit.test2[2,3]
  age.interact[i, 8] <- logit.test2[2,5]
  
  age.interact[i, 9] <- logit.test3[2,4]
  age.interact[i, 10] <- logit.test3[2,3]
  age.interact[i, 11] <- logit.test3[2,5]
  
  age.interact[i, 12] <- AIC(nb1)
  age.interact[i, 13] <- nb.test1[2,4]
  age.interact[i, 14] <- nb.test1[2,3]
  age.interact[i, 15] <- nb.test1[2,5]
  
  age.interact[i, 16] <- nb.test2[2,4]
  age.interact[i, 17] <- nb.test2[2,3]
  age.interact[i, 18] <- nb.test2[2,5]
  
  age.interact[i, 19] <- nb.test3[2,4]
  age.interact[i, 20] <- nb.test3[2,3]
  age.interact[i, 21] <- nb.test3[2,5]
  
  rm(logit.test1, logit.test2, logit.test3, nb.test1, nb.test2, nb.test3,
     logit1, logit2, logit3, logitnull, nb1, nb2, nb3, nbull)
}

colnames(age.interact) <- c("variable", "logit.AIC",
                           "logit.chi.vs.plain", "df1", "logit.p.vs.plain",
                           "logit.chi.vs.age", "df2", "logit.p.vs.age",
                           "logit.chi.vs.null", "df3", "logit.p.vs.null",
                           "nb.AIC",
                           "nb.chi.vs.plain", "df4", "nb.p.vs.plain",
                           "nb.chi.vs.age", "df5", "nb.p.vs.age",
                           "nb.chi.vs.null", "df6", "nb.p.vs.null")

# (b) Models interacting with location (two-way) ####
loc.interact <- data.frame(matrix(ncol=19, nrow=0))

for(i in c("age.cement", "sex", "new.pca", "kfi", "splenic.index",
           "food.vol.washed", "diet.Shannon",
           "anth.vol", "anth.dig.vol", "anth.indig.vol",
           "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
           "veg.vol", "herb.veg.vol", "wood.veg.vol", 
           "fruit.vol", "d13C", "d15N")){
  
  # For infection status:
  # Create logistic regressions for (1) location * focal variable, (2) location + focal variable, (3) location only, and (4) null model
  logit1 <- glm(bioact ~ location * sample_data_scale[,i],
                family="binomial"(link="logit"),
                sample_data_scale, weights = weights.l)
  logit2 <- glm(bioact ~ location + sample_data_scale[,i], family="binomial"(link="logit"), sample_data_scale, weights = weights.l)
  logit3 <- glm(bioact ~ location, family="binomial"(link="logit"), sample_data_scale, weights = weights.l)
  logitnull <- glm(bioact ~ 1, family="binomial"(link="logit"), sample_data_scale, weights = weights.l)
  
  logit.test1 <- lrtest(logit1, logit2)
  logit.test2 <- lrtest(logit1, logit3)
  logit.test3 <- lrtest(logit1, logitnull)
  
  # For infection intensity:
  # Create nbinom regressions for (1) location * focal variable, (2) location + focal variable, (3) location only, and (4) null model
  nb1 <- MASS::glm.nb(log1.worm ~ location * sample_data_scale[,i], sample_data_scale, weights = weights.l)
  nb2 <- MASS::glm.nb(log1.worm ~ location + sample_data_scale[,i], sample_data_scale, weights = weights.l)
  nb3 <- MASS::glm.nb(log1.worm ~ location, sample_data_scale, weights = weights.l)
  nbnull <- MASS::glm.nb(log1.worm ~ 1, sample_data_scale, weights = weights.l)
  
  nb.test1 <- lrtest(nb1, nb2)
  nb.test2 <- lrtest(nb1, nb3)
  nb.test3 <- lrtest(nb1, nbnull)
  
  # Store the outputs of the likelihood ratio tests in a data frame
  loc.interact[i, 1] <- i
  loc.interact[i, 2] <- AIC(logit1)
  loc.interact[i, 3] <- logit.test1[2,4]
  loc.interact[i, 4] <- logit.test1[2,3]
  loc.interact[i, 5] <- logit.test1[2,5]
  
  loc.interact[i, 6] <- logit.test2[2,4]
  loc.interact[i, 7] <- logit.test2[2,3]
  loc.interact[i, 8] <- logit.test2[2,5]
  
  loc.interact[i, 9] <- logit.test3[2,4]
  loc.interact[i, 10] <- logit.test3[2,3]
  loc.interact[i, 11] <- logit.test3[2,5]
  
  loc.interact[i, 12] <- AIC(nb1)
  loc.interact[i, 13] <- nb.test1[2,4]
  loc.interact[i, 14] <- nb.test1[2,3]
  loc.interact[i, 15] <- nb.test1[2,5]
  
  loc.interact[i, 16] <- nb.test2[2,4]
  loc.interact[i, 17] <- nb.test2[2,3]
  loc.interact[i, 18] <- nb.test2[2,5]
  
  loc.interact[i, 19] <- nb.test3[2,4]
  loc.interact[i, 20] <- nb.test3[2,3]
  loc.interact[i, 21] <- nb.test3[2,5]
  
  rm(logit.test1, logit.test2, logit.test3, nb.test1, nb.test2, nb.test3,
     logit1, logit2, logit3, logitnull, nb1, nb2, nb3, nbull)
}

colnames(loc.interact) <- c("variable", "logit.AIC",
                            "logit.chi.vs.plain", "df1", "logit.p.vs.plain",
                            "logit.chi.vs.loc", "df2", "logit.p.vs.loc",
                            "logit.chi.vs.null", "df3", "logit.p.vs.null",
                            "nb.AIC",
                            "nb.chi.vs.plain", "df4", "nb.p.vs.plain",
                            "nb.chi.vs.loc", "df5", "nb.p.vs.loc",
                            "nb.chi.vs.null", "df6", "nb.p.vs.null")

write.csv(age.interact, "~/Documents/Deanna_MS_data/models.age.interact.weighted.csv")
write.csv(loc.interact, "~/Documents/Deanna_MS_data/models.loc.interact.weighted.csv")

# (c) Models interacting with age and location (three-way interactions) ####
tw.interact <- data.frame(matrix(ncol=15, nrow=0))

for(i in c("sex", "new.pca", "kfi", "splenic.index",
           "food.vol.washed", "diet.Shannon",
           "anth.vol", "anth.dig.vol", "anth.indig.vol",
           "prey.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "insect.vol",
           "veg.vol", "herb.veg.vol", "wood.veg.vol", 
           "fruit.vol", "d13C", "d15N")){
  
  logit1 <- glm(bioact ~ location * age.cement * sample_data_scale[,i],
                family="binomial"(link="logit"),
                weights = weights.l,
                sample_data_scale)
  logit2 <- glm(bioact ~ location * age.cement, 
                family="binomial"(link="logit"),
                weights = weights.l,
                sample_data_scale)
  logit3 <- glm(bioact ~ age.cement * sample_data_scale[,i], 
                family="binomial"(link="logit"),
                weights = weights.l,
                sample_data_scale)
  logit4 <- glm(bioact ~ age.cement,
                family="binomial"(link="logit"),
                weights = weights.l,
                sample_data_scale)
  logit5 <- glm(bioact ~ location * sample_data_scale[,i],
                family="binomial"(link="logit"),
                weights = weights.l,
                sample_data_scale)
  logit6 <- glm(bioact ~ location,
                family="binomial"(link="logit"),
                weights = weights.l,
                sample_data_scale)
  
  logit.test1 <- lrtest(logit3, logit4)
  logit.test2 <- lrtest(logit5, logit6)
  logit.test3 <- lrtest(logit1, logit2)
  
  nb1 <- MASS::glm.nb(log1.worm ~ location * age.cement * sample_data_scale[,i], 
                      weights = weights.l,
                      sample_data_scale)
  nb2 <- MASS::glm.nb(log1.worm ~ location * age.cement,
                      weights = weights.l,
                      sample_data_scale)
  nb3 <- MASS::glm.nb(log1.worm ~ age.cement * sample_data_scale[,i],
                      weights = weights.l,
                      sample_data_scale)
  nb4 <- MASS::glm.nb(log1.worm ~ age.cement,
                      weights = weights.l,
                      sample_data_scale)
  nb5 <- MASS::glm.nb(log1.worm ~ location * sample_data_scale[,i],
                      weights = weights.l,
                      sample_data_scale)
  nb6 <- MASS::glm.nb(log1.worm ~ location,
                      weights = weights.l,
                      sample_data_scale)
  
  nb.test1 <- lrtest(nb3, nb4)
  nb.test2 <- lrtest(nb5, nb6)
  nb.test3 <- lrtest(nb1, nb2)
  
  # Tests of age interactions vs. age alone
  tw.interact[i, 1] <- i
  tw.interact[i, 2] <- AIC(logit3)
  tw.interact[i, 3] <- logit.test1[2,4]
  tw.interact[i, 4] <- logit.test1[2,3]
  tw.interact[i, 5] <- logit.test1[2,5]
  
  tw.interact[i, 6] <- NA
  
  tw.interact[i, 7] <- AIC(nb3)
  tw.interact[i, 8] <- nb.test1[2,4]
  tw.interact[i, 9] <- nb.test1[2,3]
  tw.interact[i, 10] <- nb.test1[2,5]
  
  tw.interact[i, 11] <- NA
  
  # Tests of location itneractions vs. locaiton alone
  tw.interact[i, 12] <- AIC(logit5)
  tw.interact[i, 13] <- logit.test2[2,4]
  tw.interact[i, 14] <- logit.test2[2,3]
  tw.interact[i, 15] <- logit.test2[2,5]
  
  tw.interact[i, 16] <- NA
  
  tw.interact[i, 17] <- AIC(nb5)
  tw.interact[i, 18] <- nb.test2[2,4]
  tw.interact[i, 19] <- nb.test2[2,3]
  tw.interact[i, 20] <- nb.test2[2,5]
  
  tw.interact[i, 21] <- NA
  
  # Tests of three-way itneractions vs. age/location alone
  tw.interact[i, 22] <- AIC(logit1)
  tw.interact[i, 23] <- logit.test3[2,4]
  tw.interact[i, 24] <- logit.test3[2,3]
  tw.interact[i, 25] <- logit.test3[2,5]
  
  tw.interact[i, 26] <- NA

  tw.interact[i, 27] <- AIC(nb1)
  tw.interact[i, 28] <- nb.test3[2,4]
  tw.interact[i, 29] <- nb.test3[2,3]
  tw.interact[i, 30] <- nb.test3[2,5]
  
  rm(logit.test1, logit.test2, logit.test3, nb.test1, nb.test2, nb.test3,
     logit1, logit2, logit3, logit4, nb1, nb2, nb3, nb4)
}

colnames(tw.interact) <- c("variable", "logit.AIC",
                           "logit.chi.vs.ageloc", "df1", "logit.p.vs.ageloc", "blank1",
                            "logit.chi.vs.age", "df2", "logit.p.vs.age",
                            "logit.chi.vs.loc", "df3", "logit.p.vs.loc",
                           "logit.AIC.nointeract", "logit.chi.vs.nointeract", "df4", "logit.p.vs.nointeract",
                            "nb.AIC",
                           "nb.chi.vs.ageloc", "df5", "nb.p.vs.ageloc",
                            "nb.chi.vs.age", "df6", "nb.p.vs.age",
                            "nb.chi.vs.loc", "df7", "nb.p.vs.loc",
                           "nb.AIC.nointeract", "nb.chi.vs.nointeract", "df8", "nb.chi.vs.nointeract")

write.csv(tw.interact, "~/Documents/Deanna_MS_data/models.tw.interact.WEIGHTED.csv")

#### DIET CLUSTER ANALYSIS (STOMACH CONTENT ORDINATIONS) ####
# Remove samples that don't have any stomach contents.
stomach_data <- subset(sample_data, food.vol.washed > 0)

# Mean-center and scale data
stomach_data_scale <- scale.vars(sample_data)

# Define variable sets to be used in the different PCA analyses
v_sets <- list(
  v_set_4 = c("anth.vol", "sm.vol", "med.vol", "ungulate.vol", "bird.vol", "veg.vol", "fruit.vol"),
  v_set_5 = c("anth.dig.vol", "anth.indig.vol", "sm.vol", "med.vol", "ungulate.vol",
              "bird.vol", "veg.vol", "fruit.vol"),
  v_set_6 = c("anth.dig.vol", "anth.indig.vol", "sm.vol", "med.vol", "ungulate.vol",
              "bird.vol", "insect.vol", "herb.veg.vol", "wood.veg.vol", "fruit.vol"))

# Run three distance/transformation combinations...
# (1) Euclidean distance on mean-centered and scaled data
# (2) Bray-Curtis distance on relative abundances (calculated for *each variable set*)
# (3) Aitchison distance

# Create places to store data
stomach.permanova <- list()
stomach.final.df <- list()
stomach.final.ords <- list()
stomach.envfit <- list()
stomach.spp.scr <- list()

for(i in c(1:length(v_sets))){ # For each of the different variable sets being tested...
  stomach.permanova[[i]] <- list()
  stomach.final.df[[i]] <- list()
  stomach.final.ords[[i]] <- list()
  stomach.envfit[[i]] <- list()
  stomach.spp.scr[[i]] <- list()
  
  for(p in c(1:3)){ # For each of the distance/transformation combinations being tested...
    
    # Prepare data. Start with Euclidean distance.
    if(p==1){  
      temp_data <- stomach_data_scale[,c("sampleID", "location", "age.cement", "new.pca",
                                         "log1.worm", "bioact", "food.vol.washed",
                                         v_sets[[i]])]
      dist.matrix <- vegdist(temp_data[,v_sets[[i]]], method="euclidean")
    }
    
    # For Bray-Curtis data, convert to relative abundance of the total food volume
    if(p==2){
      temp_data <- stomach_data[,c("sampleID", "location", "age.cement", "new.pca",
                                   "log1.worm", "bioact", "food.vol.washed",
                                   v_sets[[i]])]
      temp_data <- subset(temp_data, rowSums(temp_data[,c(8:ncol(temp_data))]) > 0)
      
      temp_data[,c(8:ncol(temp_data))] <- decostand(temp_data[,c(8:ncol(temp_data))], 
                                                    method="total", margin=1)
      
      dist.matrix <- vegdist(temp_data[,v_sets[[i]]], method="bray")
    }
    
    # Aitchison distance
    if(p==3){
      temp_data <- stomach_data[,c("sampleID", "location", "age.cement", "new.pca", 
                                   "log1.worm", "bioact", "food.vol.washed",
                                   v_sets[[i]])]
      
      temp_data <- subset(temp_data, rowSums(temp_data[,c(8:ncol(temp_data))]) > 0)
      temp_data[,c(8:ncol(temp_data))] <- decostand(temp_data[,c(8:ncol(temp_data))], 
                                                    method="total", margin=1)
      
      for(j in c(8:ncol(temp_data))){
        temp_data[,j] <- temp_data[,j] + 0.5*min(temp_data[,c(8:ncol(temp_data))][temp_data[,c(8:ncol(temp_data))] > 0])
      }
     
      temp_data[,c(8:ncol(temp_data))] <- compositions::clr(temp_data[,c(8:ncol(temp_data))])
      
      dist.matrix <- vegdist(temp_data[,v_sets[[i]]], method="euclidean")
    }
    
    # Run PERMANOVA
    tbl <- adonis(dist.matrix ~ bioact + log1.worm + location * age.cement, temp_data)$aov.tab
    
    # Calculate ordinations
    stomach.ord <- capscale(dist.matrix ~ 1)
    temp_data <- cbind(temp_data,
                         scores(stomach.ord, choices=c(1:2), display="sites"))

    set.seed(100)
    vf <- envfit(stomach.ord, temp_data[,c("location", "age.cement", "bioact", "log1.worm",
                                           "food.vol.washed", v_sets[[i]])])
    
    # Extract envfit data into a data frame for plotting
    spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
    spp.scrs <- cbind(variable = rownames(spp.scrs), 
                      spp.scrs,
                      r2 = vf$vectors$r,
                      p = vf$vectors$pvals)
    
    # Save all data into defined tables
    stomach.permanova[[i]][[p]] <- tbl
    stomach.final.df[[i]][[p]] <- temp_data
    stomach.final.ords[[i]][[p]] <- stomach.ord
    stomach.envfit[[i]][[p]] <- vf
    stomach.spp.scr[[i]][[p]] <- spp.scrs
    
    rm(tbl, temp_data, stomach.ord, vf, spp.scrs, dist.matrix)
  }
}

# Look at the PERMANOVA results for all the different facets
tmp <- list()
for(i in c(1:3)){
  tmp[[i]] <- rbind(stomach.permanova[[i]][[1]],
                    stomach.permanova[[i]][[2]],
                    stomach.permanova[[i]][[3]])
}

permanova.master <- cbind(tmp[[1]], tmp[[2]], tmp[[3]])
rm(tmp)

openxlsx::write.xlsx(permanova.master, "stomach.permanova.master.xlsx",
                     row.names=TRUE)

#### THREE-WAY INTERACTION MODELS (FIGURES 3 and 5) ####
## Prevalence model ####
# Define the model
model <- glm(bioact ~ age.cement * location * (anth.dig.vol + sm.vol + d13C + d15N + anth.indig.vol),
             family="binomial"(link="logit"),
             sample_data_scale,
             weights = weights.l,
             na.action="na.fail")

# Evaluate all model subsets
model.dredge <- dredge(model) # three top models with weights, six top models without weights

# Average coefficients across top models
model.average <- model.avg(model.dredge, subset = delta < 2, beta="partial.sd")

# Calculate confidence intervals
temp2 <- cbind(confint(model.average, level=0.95, full=TRUE)[,1],
               confint(model.average, level=0.5, full=TRUE)[,1],
               model.average[["coefficients"]][1,],
               confint(model.average, level=0.5, full=TRUE)[,2],
               confint(model.average, level=0.95, full=TRUE)[,2])
temp2 <- as.data.frame(temp2)
colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")

# Calculate the sum of model weights of all models in which each predictor appears
model.dredge <- as.data.frame(model.dredge)
weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:(ncol(model.dredge)-5))){
  temp <- model.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(model.dredge[,c(2:(ncol(model.dredge)-5))])
colnames(weights) <- "cum_weight"

temp2$merge <- rownames(temp2)
weights$merge <- as.factor(rownames(weights))

temp2 <- data.frame(lapply(temp2, function(x) {
  gsub("locationurban", "location", x)
}))

# Combine all data together
results <- merge(temp2, weights, by="merge", all=TRUE)
results <- subset(results, merge !="(Intercept)")
results <- subset(results, merge !="food.vol.washed")
results <- results[order(-results$cum_weight),]
results$merge <- factor(results$merge)
results$merge <- forcats::fct_inorder(results$merge)

for(i in c(2:6)){
  results[,i] <- as.numeric(as.character(results[,i]))
}

results$model <- rep("prevalence")

## Intensity model ####
model <- MASS::glm.nb(log1.worm ~ age.cement * location * (anth.dig.vol + anth.indig.vol + sm.vol + d13C + d15N),
                      
                      sample_data_scale, 
                      weights = weights.l,
                      na.action="na.fail")
model.dredge <- dredge(model) # three top models with weights, six top models without weights

model.average <- model.avg(model.dredge, subset = delta < 2, beta="partial.sd")

temp2 <- cbind(confint(model.average, level=0.95, full=TRUE)[,1],
               confint(model.average, level=0.5, full=TRUE)[,1],
               model.average[["coefficients"]][1,],
               confint(model.average, level=0.5, full=TRUE)[,2],
               confint(model.average, level=0.95, full=TRUE)[,2])
temp2 <- as.data.frame(temp2)
colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")

model.dredge <- as.data.frame(model.dredge)
weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:(ncol(model.dredge)-5))){
  temp <- model.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(model.dredge[,c(2:(ncol(model.dredge)-5))])
colnames(weights) <- "cum_weight"

temp2$merge <- rownames(temp2)
weights$merge <- as.factor(rownames(weights))

temp2 <- data.frame(lapply(temp2, function(x) {
  gsub("locationurban", "location", x)
}))

results2 <- merge(temp2, weights, by="merge", all=TRUE)
results2 <- subset(results2, merge !="(Intercept)")
results2 <- results2[order(-results2$cum_weight),]
results2$merge <- factor(results2$merge)
results2$merge <- forcats::fct_inorder(results2$merge)

for(i in c(2:6)){
  results2[,i] <- as.numeric(as.character(results2[,i]))
}

results2$model <- rep("intensity")

results <- rbind(results, results2)

## Rename coefficients in the final data frame ####
levels(results$merge)[levels(results$merge)=="age.cement"] <- "age"
levels(results$merge)[levels(results$merge)=="age.cement:location"] <- "age x location"
levels(results$merge)[levels(results$merge)=="anth.dig.vol"] <- "dig. anthro"
levels(results$merge)[levels(results$merge)=="anth.dig.vol:location"] <- "x location"
levels(results$merge)[levels(results$merge)=="anth.indig.vol"] <- "indig. anthro"
levels(results$merge)[levels(results$merge)=="sm.vol"] <- "rodents"
levels(results$merge)[levels(results$merge)=="age.cement:sm.vol"] <- "__ x age"
levels(results$merge)[levels(results$merge)=="age.cement:anth.dig.vol"] <- "x age"
levels(results$merge)[levels(results$merge)=="age.cement:anth.dig.vol:location"] <- "x age x location"
levels(results$merge)[levels(results$merge)=="age.cement:anth.indig.vol"] <- "_ x age"
levels(results$merge)[levels(results$merge)=="anth.indig.vol:location"] <- "_ x location"
levels(results$merge)[levels(results$merge)=="location:sm.vol"] <- "__ x location"
levels(results$merge)[levels(results$merge)=="age.cement:location:sm.vol"] <- "__ x age x location"
levels(results$merge)[levels(results$merge)=="age.cement:anth.indig.vol:location"] <- "_ x age x location"

levels(results$merge)[levels(results$merge)=="age.cement:d13C"] <- "___ x age"
levels(results$merge)[levels(results$merge)=="age.cement:d15N"] <- "____ x age"
levels(results$merge)[levels(results$merge)=="d13C:location"] <- "___ x location"
levels(results$merge)[levels(results$merge)=="d15N:location"] <- "____ x location"
levels(results$merge)[levels(results$merge)=="age.cement:d15N:location"] <- "____ x age x location"
levels(results$merge)[levels(results$merge)=="age.cement:d13C:location"] <- "___ x age x location"

results$merge <- factor(results$merge,
                        levels = c("age", "location", "age x location", 
                                   "dig. anthro", "x age", "x location", "x age x location",
                                   "indig. anthro", "_ x age", "_ x location", "_ x age x location",
                                   "rodents", "__ x age", "__ x location", "__ x age x location",
                                   "d13C", "___ x age", "___ x location", "___ x age x location",
                                   "d15N", "____ x age", "____ x location", "____ x age x location"))

results$model <- as.character(results$model)
results$model[results$model == "prevalence"] <- "status"
results$model <- factor(results$model, levels = c("status", "intensity"))

results$facet <- rep("context")
results$facet[results$merge %in% c("dig. anthro", "x age", "x location", "x age x location")] <- "dig. anthro"
results$facet[results$merge %in% c("indig. anthro", "_ x age", "_ x location", "_ x age x location")] <- "indig. anthro"
results$facet[results$merge %in% c("rodents", "__ x age", "__ x location", "__ x age x location")] <- "rodents"
results$facet[results$merge %in% c("d13C", "___ x age", "___ x location", "___ x age x location")] <- "d13C"
results$facet[results$merge %in% c("d15N", "____ x age", "____ x location", "____ x age x location")] <- "d15N"

results$facet <- factor(results$facet, levels=c("context",
                                                "dig. anthro",
                                                "indig. anthro",
                                                "rodents",
                                                "d13C", "d15N"))

results$merge <- as.character(results$merge)
results$merge[results$merge == "____ x age"] <- "x age"
results$merge[results$merge == "___ x age"] <- "x age"
results$merge[results$merge == "__ x age"] <- "x age"
results$merge[results$merge == "_ x age"] <- "x age"
results$merge[results$merge == "____ x location"] <- "x location"
results$merge[results$merge == "___ x location"] <- "x location"
results$merge[results$merge == "__ x location"] <- "x location"
results$merge[results$merge == "_ x location"] <- "x location"
results$merge[results$merge == "____ x age x location"] <- "x age x location"
results$merge[results$merge == "___ x age x location"] <- "x age x location"
results$merge[results$merge == "__ x age x location"] <- "x age x location"
results$merge[results$merge == "_ x age x location"] <- "x age x location"

results$merge <- factor(results$merge, levels=c("age", "location", "age x location",
                                                "dig. anthro", "indig. anthro", "rodents", "d13C", "d15N",
                                                "x age", "x location", "x age x location"))

## Infection prevalence: create a data frame showing interaction results for each of the five interacting variables ####
# This data frame will be used to create the panels in Figures 3 and 5.
# For each of the five dietary predictors, the code below calculates the curves for each of the three ages (mean +/- sd) in each of the
# two sampled locations. It therefore produces six curves for each variable = the 30 curves shown in Figures 3a and 5a.
library(effects)

# Digestible anthropogenic food
model <- glm(bioact ~ age.cement * location * anth.dig.vol,
             family="binomial"(link="logit"),
             sample_data,
             weights = weights.l,
             na.action="na.fail")

handpick <- effect('age.cement*location*anth.dig.vol', mod = model,
                   xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                  anth.dig.vol=seq(0, 360, by=0.05)),
                   location = c("urban", "rural"),
                   se=TRUE, confidence.level=.95, typical=mean)
handpick <- as.data.frame(handpick)
handpick$age.cement <- as.factor(as.character(handpick$age.cement))

handpick$age.cement <- factor(handpick$age.cement,
                              levels=c(0.256, 2.47, 4.79))

# Rodents
model <- glm(bioact ~ age.cement * location * sm.vol,
             family="binomial"(link="logit"),
             sample_data,
             weights = weights.l,
             na.action="na.fail")

handpick2 <- effect('age.cement*location*sm.vol', mod = model,
                   xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                  sm.vol=seq(-0.5, 325, by=0.05)),
                   location = c("urban", "rural"),
                   se=TRUE, confidence.level=.95, typical=mean)
handpick2 <- as.data.frame(handpick2)
handpick2$age.cement <- as.factor(as.character(handpick2$age.cement))

handpick2$age.cement <- factor(handpick2$age.cement,
                              levels=c(0.256, 2.47, 4.79))

# Indigestible anthropogenic food
model <- glm(bioact ~ age.cement * location * anth.indig.vol,
             family="binomial"(link="logit"),
             sample_data,
             weights = weights.l,
             na.action="na.fail")

handpick3 <- effect('age.cement*location*anth.indig.vol', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   anth.indig.vol=seq(min(sample_data$anth.indig.vol),
                                              max(sample_data$anth.indig.vol), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick3 <- as.data.frame(handpick3)
handpick3$age.cement <- as.factor(as.character(handpick3$age.cement))

handpick3$age.cement <- factor(handpick3$age.cement,
                               levels=c(0.256, 2.47, 4.79))

# d13C isotopes
model <- glm(bioact ~ age.cement * location * d13C,
             family="binomial"(link="logit"),
             sample_data,
             weights = weights.l,
             na.action="na.fail")

handpick4 <- effect('age.cement*location*d13C', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   d13C=seq(min(sample_data$d13C),
                                                      max(sample_data$d13C), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick4 <- as.data.frame(handpick4)
handpick4$age.cement <- as.factor(as.character(handpick4$age.cement))

handpick4$age.cement <- factor(handpick4$age.cement,
                               levels=c(0.256, 2.47, 4.79))

# d15N isotopes
model <- glm(bioact ~ age.cement * location * d15N,
             family="binomial"(link="logit"),
             sample_data,
             weights = weights.l,
             na.action="na.fail")

handpick5 <- effect('age.cement*location*d15N', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   d15N=seq(min(sample_data$d15N),
                                            max(sample_data$d15N), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick5 <- as.data.frame(handpick5)
handpick5$age.cement <- as.factor(as.character(handpick5$age.cement))

handpick5$age.cement <- factor(handpick5$age.cement,
                               levels=c(0.256, 2.47, 4.79))

handpick$variable <- rep("anth.dig.vol")
handpick2$variable <- rep("sm.vol")
handpick3$variable <- rep("anth.indig.vol")
handpick4$variable <- rep("d13C")
handpick5$variable <- rep("d15N")

# Scale all the x axes so that they appear on the same scale.
# In this way, variables with different values (or ranges of values) can appear on the same axis.
scaling <- function(x){(x-min(x))/(max(x)-min(x))}

handpick$xaxis <- scaling(handpick$anth.dig.vol)
handpick2$xaxis <- scaling(handpick2$sm.vol)
handpick3$xaxis <- scaling(handpick3$anth.indig.vol)
handpick4$xaxis <- scaling(handpick4$d13C)
handpick5$xaxis <- scaling(handpick5$d15N)

handpick$type <- rep("stomach contents 1")
handpick2$type <- rep("stomach contents 2")
handpick3$type <- rep("stomach contents 3")
handpick4$type <- rep("stable isotope 1")
handpick5$type <- rep("stable isotope 2")

## Combine all interaction plots into a single data frame ####
test1 <- handpick[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test2 <- handpick2[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test3 <- handpick3[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test4 <- handpick4[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test5 <- handpick5[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test <- rbind(test1, test2, test3, test4, test5)

test$variable <- factor(test$variable, levels=c("anth.dig.vol", "anth.indig.vol", "sm.vol",
                                                "d13C", "d15N"))

levels(test$age.cement)[levels(test$age.cement)=="0.256"] <- "Younger (0.26 yr)"
levels(test$age.cement)[levels(test$age.cement)=="2.47"] <- "Mean (2.47 yr)"
levels(test$age.cement)[levels(test$age.cement)=="4.79"] <- "Older (4.79 yr)"

test$location <- factor(test$location, levels=c("urban", "rural"))

levels(test$variable)[levels(test$variable)=="anth.dig.vol"] <- "digestible anthropogenic food"
levels(test$variable)[levels(test$variable)=="anth.indig.vol"] <- "indigestible anthropogenic food"
levels(test$variable)[levels(test$variable)=="sm.vol"] <- "rodents"

## Infection intensity: repeat the code above for the models predicting infection intensity ####
# Digestible anthropogenic food
model <- glm.nb(log1.worm ~ age.cement * location * anth.dig.vol,
             sample_data,
             weights = weights.l,
             na.action="na.fail")

nb.handpick <- effect('age.cement*location*anth.dig.vol', mod = model,
                   xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                  anth.dig.vol=seq(0, 360, by=0.05)),
                   location = c("urban", "rural"),
                   se=TRUE, confidence.level=.95, typical=mean)
nb.handpick <- as.data.frame(nb.handpick)
nb.handpick$age.cement <- as.factor(as.character(nb.handpick$age.cement))

nb.handpick$age.cement <- factor(nb.handpick$age.cement,
                              levels=c(0.256, 2.47, 4.79))

# Rodents
model <- glm.nb(log1.worm ~ age.cement * location * sm.vol,
             sample_data,
             weights = weights.l,
             na.action="na.fail")

nb.handpick2 <- effect('age.cement*location*sm.vol', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   sm.vol=seq(-0.5, 325, by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
nb.handpick2 <- as.data.frame(nb.handpick2)
nb.handpick2$age.cement <- as.factor(as.character(nb.handpick2$age.cement))

nb.handpick2$age.cement <- factor(nb.handpick2$age.cement,
                               levels=c(0.256, 2.47, 4.79))

# Indigestible anthropogenic food
model <- glm.nb(log1.worm ~ age.cement * location * anth.indig.vol,
             sample_data,
             weights = weights.l,
             na.action="na.fail")

nb.handpick3 <- effect('age.cement*location*anth.indig.vol', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   anth.indig.vol=seq(min(sample_data$anth.indig.vol),
                                                      max(sample_data$anth.indig.vol), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
nb.handpick3 <- as.data.frame(nb.handpick3)
nb.handpick3$age.cement <- as.factor(as.character(nb.handpick3$age.cement))

nb.handpick3$age.cement <- factor(nb.handpick3$age.cement,
                               levels=c(0.256, 2.47, 4.79))

# d13C isotopes
model <- glm.nb(log1.worm ~ age.cement * location * d13C,
             sample_data,
             weights = weights.l,
             na.action="na.fail")

nb.handpick4 <- effect('age.cement*location*d13C', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   d13C=seq(min(sample_data$d13C),
                                            max(sample_data$d13C), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
nb.handpick4 <- as.data.frame(nb.handpick4)
nb.handpick4$age.cement <- as.factor(as.character(nb.handpick4$age.cement))

nb.handpick4$age.cement <- factor(nb.handpick4$age.cement,
                               levels=c(0.256, 2.47, 4.79))

# d15N isotopes
model <- glm.nb(log1.worm ~ age.cement * location * d15N,
              sample_data,
             weights = weights.l,
             na.action="na.fail")

nb.handpick5 <- effect('age.cement*location*d15N', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   d15N=seq(min(sample_data$d15N),
                                            max(sample_data$d15N), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
nb.handpick5 <- as.data.frame(nb.handpick5)
nb.handpick5$age.cement <- as.factor(as.character(nb.handpick5$age.cement))

nb.handpick5$age.cement <- factor(nb.handpick5$age.cement,
                               levels=c(0.256, 2.47, 4.79))

nb.handpick$variable <- rep("anth.dig.vol")
nb.handpick2$variable <- rep("sm.vol")
nb.handpick3$variable <- rep("anth.indig.vol")
nb.handpick4$variable <- rep("d13C")
nb.handpick5$variable <- rep("d15N")

scaling <- function(x){(x-min(x))/(max(x)-min(x))}

nb.handpick$xaxis <- scaling(nb.handpick$anth.dig.vol)
nb.handpick2$xaxis <- scaling(nb.handpick2$sm.vol)
nb.handpick3$xaxis <- scaling(nb.handpick3$anth.indig.vol)
nb.handpick4$xaxis <- scaling(nb.handpick4$d13C)
nb.handpick5$xaxis <- scaling(nb.handpick5$d15N)

nb.handpick$type <- rep("stomach contents 1")
nb.handpick2$type <- rep("stomach contents 2")
nb.handpick3$type <- rep("stomach contents 3")
nb.handpick4$type <- rep("stable isotope 1")
nb.handpick5$type <- rep("stable isotope 2")

nb.test1 <- nb.handpick[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
nb.test2 <- nb.handpick2[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
nb.test3 <- nb.handpick3[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
nb.test4 <- nb.handpick4[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
nb.test5 <- nb.handpick5[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
nb.test <- rbind(nb.test1, nb.test2, nb.test3, nb.test4, nb.test5)

nb.test$variable <- factor(nb.test$variable, levels=c("anth.dig.vol", "anth.indig.vol", "sm.vol",
                                                "d13C", "d15N"))

levels(nb.test$age.cement)[levels(nb.test$age.cement)=="0.256"] <- "Younger (0.26 yr)"
levels(nb.test$age.cement)[levels(nb.test$age.cement)=="2.47"] <- "Mean (2.47 yr)"
levels(nb.test$age.cement)[levels(nb.test$age.cement)=="4.79"] <- "Older (4.79 yr)"

nb.test$location <- factor(nb.test$location, levels=c("urban", "rural"))

levels(nb.test$variable)[levels(nb.test$variable)=="anth.dig.vol"] <- "digestible anthropogenic food"
levels(nb.test$variable)[levels(nb.test$variable)=="anth.indig.vol"] <- "indigestible anthropogenic food"
levels(nb.test$variable)[levels(nb.test$variable)=="sm.vol"] <- "rodents"

nb.test$fit[nb.test$fit > 21] <- NA

#### K-FOLD CROSS VALIDATION: LOGISTIC REGRESSION #####
# Recreate the logistic regression model
model <- glm(bioact ~ age.cement * location * (anth.dig.vol + anth.indig.vol + sm.vol + d13C + d15N),
             family="binomial"(link="logit"),
             sample_data_scale,
             weights = sample_data_scale$weights.l,
             na.action="na.fail")
model.dredge <- dredge(model)

# Extract all the top-ranked models 
list.of.models <- get.models(model.dredge, subset = delta < 2)

# Calculate coefficients and confidence intervals for each model (standardized by partial standard deviation)
std.coefs <- list()

for(i in c(1:length(list.of.models))){
  std.coefs[[i]] <- as.data.frame(std.coef(list.of.models[[i]], partial.sd = TRUE))
  colnames(std.coefs[[i]]) <- c("estimate", "se", "df")
  std.coefs[[i]]$conf <- std.coefs[[i]]$se * qnorm(0.975)
  std.coefs[[i]]$low <- std.coefs[[i]]$estimate - std.coefs[[i]]$conf
  std.coefs[[i]]$high <- std.coefs[[i]]$estimate + std.coefs[[i]]$conf
  std.coefs[[i]][,c("se","df")] <- NULL
  
  colnames(std.coefs[[i]]) <- paste0(colnames(std.coefs[[i]]), i)
  
  std.coefs[[i]]$merge <- rownames(std.coefs[[i]])
}

std.coefs.table <- merge(std.coefs[[1]], std.coefs[[2]], by="merge", all=TRUE)

for(i in c(3:length(std.coefs))){
  std.coefs.table <- merge(std.coefs.table, std.coefs[[i]], by="merge", all=TRUE)
}

rownames(std.coefs.table) <- std.coefs.table$merge
std.coefs.table$merge <- NULL

write.csv(std.coefs.table, "~/Documents/Deanna_MS_data/logit_coefs.csv")

# Perform cross-validation of each model to assess model accuracy #
model.eval <- data.frame(matrix(ncol=6, nrow=0))
ctrl <- trainControl(method = "cv", number = 5)

for(i in c(1:length(list.of.models))){
  set.seed(100)
  modelo1 <- train(list.of.models[[i]]$formula,
                   data = sample_data_scale, 
                   weights = sample_data_scale$weights.l,
                   method = "glm", family = binomial,
                   trControl = ctrl)
  
  vector <- c(modelo1$results, pR2(list.of.models[[i]])[[4]])
  vector <- as.vector(vector)
  names(vector) <- c("parameter", "accuracy", "kappa", "accuracySD", "kappaSD", "mcfadden")
  
  model.eval <- rbind(model.eval, vector)
}

write.csv(model.eval, "logit_model_fit.csv")

## Check the null, global, and age/location only models for AIC and accuracy estimates ####
model.table <- as.data.frame(model.dredge)

# Look at information for the null model
model.1 <- model.table
model.1 <- subset(model.1, is.na(age.cement)=="TRUE" & is.na(location)=="TRUE" & 
                    is.na(d13C)=="TRUE" & is.na(d15N)=="TRUE" & is.na(anth.dig.vol)=="TRUE" &
                    is.na(anth.indig.vol)=="TRUE" & is.na(sm.vol)=="TRUE")
model.1 <- get.models(model.dredge, subset = "1")
model.1 <- model.1[[1]]
pR2(model.1)

# Get information for the age/location only model
model.2 <- model.table
model.2  <- subset(model.2, is.na(age.cement)=="FALSE" & is.na(location)=="FALSE" &
                     is.na(d13C)=="TRUE" & is.na(d15N)=="TRUE" &
                     is.na(anth.dig.vol)=="TRUE" & is.na(anth.indig.vol)=="TRUE" &
                     is.na(sm.vol)=="TRUE" & is.na(`age.cement:location`)=="FALSE")
model.2 <- get.models(model.dredge, subset = rownames(model.2))
model.2 <- model.2[[1]]
pR2(model.2)

# Get information for the global model
model.3 <- model.table
model.3  <- subset(model.3, is.na(age.cement)=="FALSE" & is.na(location)=="FALSE" &
                     is.na(d13C)=="FALSE" & is.na(d15N)=="FALSE" &
                     is.na(anth.dig.vol)=="FALSE" & is.na(anth.indig.vol)=="FALSE" &
                     is.na(sm.vol)=="FALSE" & is.na(`age.cement:d15N:location`)=="FALSE" &
                     is.na(`age.cement:anth.indig.vol:location`)=="FALSE" &
                     is.na(`age.cement:location:sm.vol`)=="FALSE" &
                     is.na(`age.cement:d13C:location`)=="FALSE" &
                     is.na(`age.cement:anth.dig.vol:location`)=="FALSE" &
                     is.na(`age.cement:d15N:location`)=="FALSE" &
                     is.na(`age.cement:anth.dig.vol`)=="FALSE")
model.3 <- get.models(model.dredge, subset = rownames(model.3))
model.3 <- model.3[[1]]
pR2(model.3)

# Calculate accuracy and kappa for the second two models (age/location and global)
model.list <- list(model.2, model.3)

model.eval <- data.frame(matrix(ncol=6, nrow=0))

for(i in c(1:length(model.list))){
  set.seed(100)
  modelo1 <- train(model.list[[i]]$formula,
                   data = sample_data_scale, 
                   weights = sample_data_scale$weights.l,
                   method = "glm", family = binomial,
                   trControl = ctrl)
  
  vector <- c(modelo1$results, pR2(model.list[[i]])[[4]])
  vector <- as.vector(vector)
  names(vector) <- c("parameter", "accuracy", "kappa", "accuracySD", "kappaSD", "mcfadden")
  
  model.eval <- rbind(model.eval, vector)
}

write.csv(model.eval, "logit_model_fit_ctrl.csv")


#### K-FOLD CROSS VALIDATION: NEGATIVE BINOMIAL REGRESSION ####
# Recreate the negative binomial regression model 
model <- glm.nb(log1.worm ~ age.cement * location * (anth.dig.vol + anth.indig.vol + sm.vol + d13C + d15N),
                sample_data_scale,
                weights = sample_data_scale$weights.l,
                na.action="na.fail")
model.dredge <- dredge(model) # three top models with weights, six top models without weights

list.of.models <- get.models(model.dredge, subset = delta < 2)

# Retrieve coefficients
std.coefs <- list()

for(i in c(1:length(list.of.models))){
  std.coefs[[i]] <- as.data.frame(std.coef(list.of.models[[i]], partial.sd = TRUE))
  colnames(std.coefs[[i]]) <- c("estimate", "se", "df")
  std.coefs[[i]]$conf <- std.coefs[[i]]$se * qnorm(0.975)
  std.coefs[[i]]$low <- std.coefs[[i]]$estimate - std.coefs[[i]]$conf
  std.coefs[[i]]$high <- std.coefs[[i]]$estimate + std.coefs[[i]]$conf
  std.coefs[[i]][,c("se","df")] <- NULL
  
  colnames(std.coefs[[i]]) <- paste0(colnames(std.coefs[[i]]), i)
  
  std.coefs[[i]]$merge <- rownames(std.coefs[[i]])
}

std.coefs.table <- merge(std.coefs[[1]], std.coefs[[2]], by="merge", all=TRUE)

for(i in c(3:length(std.coefs))){
  std.coefs.table <- merge(std.coefs.table, std.coefs[[i]], by="merge", all=TRUE)
}

rownames(std.coefs.table) <- std.coefs.table$merge
std.coefs.table$merge <- NULL

write.csv(std.coefs.table, "~/Documents/Deanna_MS_data/nbinom_coefs.csv")

# Cross validation
ctrl <- trainControl(method = "cv", number = 5)

model.eval <- data.frame(matrix(ncol=8, nrow=0))

for(i in c(1:length(list.of.models))){
  set.seed(100)
  modelo1 <- train(list.of.models[[i]]$call$formula,
                   data = df, 
                   weights = sample_data_scale$weights.l,
                   method = "glm.nb",
                   trControl = ctrl)
  
  vector <- c(modelo1$results[2,], pR2(list.of.models[[i]])[[4]])
  vector <- as.vector(vector)
  names(vector) <- c("link", "RMSE", "RSq", "MAE", "RMSE.SD", "RSq.SD", "MAE.SD", "mcfadden")
  
  model.eval <- rbind(model.eval, vector)
}

write.csv(model.eval, "~/Documents/Deanna_MS_data/nbinom_model_eval.csv")

## Check the null, global, and age/location only models for AIC and accuracy estimates ####
# Look at information for the null model
model.table <- as.data.frame(model.dredge)

model.1 <- model.table
model.1 <- subset(model.1, is.na(age.cement)=="TRUE" & is.na(location)=="TRUE" & 
                    is.na(d13C)=="TRUE" & is.na(d15N)=="TRUE" & is.na(anth.dig.vol)=="TRUE" &
                    is.na(anth.indig.vol)=="TRUE" & is.na(sm.vol)=="TRUE")
model.1 <- get.models(model.dredge, subset = "1")
model.1 <- model.1[[1]]
pR2(model.1)

# Get information for the age/location only model
model.2 <- model.table
model.2  <- subset(model.2, is.na(age.cement)=="FALSE" & is.na(location)=="FALSE" &
                     is.na(d13C)=="TRUE" & is.na(d15N)=="TRUE" &
                     is.na(anth.dig.vol)=="TRUE" & is.na(anth.indig.vol)=="TRUE" &
                     is.na(sm.vol)=="TRUE" & is.na(`age.cement:location`)=="FALSE")
model.2 <- get.models(model.dredge, subset = rownames(model.2))
model.2 <- model.2[[1]]
pR2(model.2)

# Get information for the global model
model.3 <- model.table
model.3  <- subset(model.3, is.na(age.cement)=="FALSE" & is.na(location)=="FALSE" &
                     is.na(d13C)=="FALSE" & is.na(d15N)=="FALSE" &
                     is.na(anth.dig.vol)=="FALSE" & is.na(anth.indig.vol)=="FALSE" &
                     is.na(sm.vol)=="FALSE" & is.na(`age.cement:d15N:location`)=="FALSE" &
                     is.na(`age.cement:anth.indig.vol:location`)=="FALSE" &
                     is.na(`age.cement:location:sm.vol`)=="FALSE" &
                     is.na(`age.cement:d13C:location`)=="FALSE" &
                     is.na(`age.cement:anth.dig.vol:location`)=="FALSE" &
                     is.na(`age.cement:d15N:location`)=="FALSE" &
                     is.na(`age.cement:anth.dig.vol`)=="FALSE")
model.3 <- get.models(model.dredge, subset = rownames(model.3))
model.3 <- model.3[[1]]
pR2(model.3)

# Calculate accuracy and kappa for the second two models (age/location and global)
model.list <- list(model.2, model.3)

model.eval <- data.frame(matrix(ncol=8, nrow=0))

for(i in c(1:length(model.list))){
  set.seed(100)
  modelo1 <- train(model.list[[i]]$call$formula,
                   data = df, 
                   weights = sample_data_scale$weights.l,
                   method = "glm.nb",
                   trControl = ctrl)
  
  vector <- c(modelo1$results[2,], pR2(model.list[[i]])[[4]])
  vector <- as.vector(vector)
  names(vector) <- c("link", "RMSE", "RSq", "MAE", "RMSE.SD", "RSq.SD", "MAE.SD", "mcfadden")
  
  model.eval <- rbind(model.eval, vector)
}


write.csv(model.eval, "logit_model_fit_ctrl.csv")





#### EFFECTS OF CONDITION ON INFECTION STATUS (FIGURE 2) ####
# Prevalence: age/location interaction #
model <- glm(bioact ~ age.cement*location,
                         family="binomial"(link="logit"),
                         sample_data,
                         weights = weights.l,
                         na.action="na.fail")

# Generate interaction data frame
cd.logit1 <- effect('age.cement*location', mod = model,
                    xlevels = list(age.cement=seq(min(sample_data$age.cement),
                                               max(sample_data$age.cement),
                                               by=0.1)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
cd.logit1 <- as.data.frame(cd.logit1)

# Intensity: age/location interaction #
model <- glm.nb(log1.worm ~ age.cement*location,
             sample_data,
             weights = weights.l,
             na.action="na.fail")

# Generate interaction data frame
cd.nb1 <- effect('age.cement*location', mod = model,
                    xlevels = list(age.cement=seq(min(sample_data$age.cement),
                                                  max(sample_data$age.cement),
                                                  by=0.1)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
cd.nb1 <- as.data.frame(cd.nb1)

# Prevalence: spleen/location/age interaction
model <- glm(bioact ~ age.cement * location * splenic.index,
             family="binomial"(link="logit"),
             sample_data,
             weights = weights.l,
             na.action="na.fail")

cd.logit2 <- effect('age.cement*location*splenic.index', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   splenic.index=seq(min(sample_data$splenic.index),
                                               max(sample_data$splenic.index),
                                               by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
cd.logit2 <- as.data.frame(cd.logit2)
cd.logit2$age.cement <- as.factor(as.character(cd.logit2$age.cement))

cd.logit2$age.cement <- factor(cd.logit2$age.cement,
                               levels=c(0.256, 2.47, 4.79))

# Intensity: spleen/location/age interaction
model <- glm.nb(log1.worm ~ age.cement * location * splenic.index,
             sample_data,
             weights = weights.l,
             na.action="na.fail")

cd.nb2 <- effect('age.cement*location*splenic.index', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   splenic.index=seq(min(sample_data$splenic.index),
                                                     max(sample_data$splenic.index),
                                                     by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
cd.nb2 <- as.data.frame(cd.nb2)
cd.nb2$age.cement <- as.factor(as.character(cd.nb2$age.cement))

cd.nb2$age.cement <- factor(cd.nb2$age.cement,
                               levels=c(0.256, 2.47, 4.79))

#### Differences in infection status based on proximity to the Leduc dump ####
sample_data_geo <- sample_data

sample_data_geo$site.tag[sample_data_geo$site.tag=="F"] <- "E"
sample_data_geo$site.tag[sample_data_geo$site.tag=="D"] <- "E"
sample_data_geo$site.tag[sample_data_geo$site.tag=="G"] <- "E"

ggplot(sample_data_geo, aes(x=site.tag, y=log(anth.indig.vol + 0.1))) +
  geom_boxplot()

ggplot(sample_data_geo, aes(x=site.tag, y=age.cement)) +
  geom_boxplot()

prev.geography <- sample_data_geo %>%
  dplyr::group_by(site.tag, bioact) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

ggplot(subset(prev.geography, bioact==1),
       aes(x=site.tag, y=freq)) + geom_bar(stat="identity")

temp <- subset(sample_data_geo, is.na(site.tag) == "FALSE")

chisq.test(temp$site.tag, temp$bioact) # No significant differences, p = 0.2077

temp.rural <- subset(temp, location=="rural")

chisq.test(temp.rural$site.tag, temp.rural$bioact) # No significant differences, p = 0.1672

summary(aov(anth.indig.vol ~ site.tag, temp)) # p = 0.0506
summary(aov(anth.indig.vol ~ site.tag, temp.rural)) # p = 0.472
summary(aov(anth.dig.vol ~ site.tag, temp)) # p = 0.134
summary(aov(age.cement ~ site.tag, temp)) # p = 0.227
summary(aov(age.cement ~ site.tag, temp.rural)) # p = 0.188

summary(aov(anth.indig.vol ~ site.tag, sample_data_geo, na.action="na.omit")) # p = 0.035
summary(aov(anth.indig.vol ~ site.tag, sample_data_geo.rural, na.action="na.omit")) # p = 0.438


sample_data_geo <- subset(sample_data_geo, is.na(site.tag)=="FALSE")
sample_data_geo <- subset(sample_data_geo, site.tag !="Z")

sample_data_geo$site.tag <- factor(sample_data_geo$site.tag)
levels(sample_data_geo$site.tag)[levels(sample_data_geo$site.tag)=="A"] <- "Edmonton\n(urban)"
levels(sample_data_geo$site.tag)[levels(sample_data_geo$site.tag)=="B"] <- "Other rural"
levels(sample_data_geo$site.tag)[levels(sample_data_geo$site.tag)=="C"] <- "Other rural"
levels(sample_data_geo$site.tag)[levels(sample_data_geo$site.tag)=="E"] <- "Leduc dump\nsite"

sample_data_geo$site.tag <- factor(sample_data_geo$site.tag,
                                   levels = c("Edmonton\n(urban)", "Leduc dump\nsite", "Other rural"))

sample_data_geo_melt <- sample_data_geo[,c("site.tag","bioact","anth.indig.vol","ech.cap","log1.worm",
                                           "age.cement")]

sample_data_geo_melt <- reshape2::melt(sample_data_geo_melt)