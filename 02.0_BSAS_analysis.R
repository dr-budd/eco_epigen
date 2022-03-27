## set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## LOAD LIBRARIES ----

library(betareg) 
library(car) 
require(devtools) 
library(cowplot) 
library(grid) 
library(RColorBrewer) 
library(ggpubr) 
library(AICcmodavg) 
library(data.table)
library(multcomp)
library(biostat)
library(tidyverse)

## IMPORT DATA ----

## import data from geneious 
data_unmod<-read.csv("02.1_geneious_data_2021.csv")

## import gene and cpg information
gene_cpg <- read.csv("02.2_gene_CpG_info.csv")

## import fish ID numbers 
fish_data<-read.csv("02.3_Barra_IDs.csv")

## import sampling information (location etc) from LTMP data
barra_data <- read.csv("02.4_barra_data.csv")
 
## FORMAT DATA ---- 

data <- data_unmod %>% 
  ## drop gene annotation information from 'type' column
  filter(Type!="gene") %>% 
  ## remove annotations that are not polymorphisms
  filter(Type=="Polymorphism") %>% 
  ## remove SNPs that didn't occur at CpG sites
  na.omit(CpG_site) %>%
  ## remove non-C-T changes
  filter(Change=="C -> T") %>%
  ## remove anything less than 500 reads 
  filter(., Coverage>499) %>%
  ## add gene and cpg information
  mutate(Minimum = as.numeric(Minimum)) %>%
  left_join(., gene_cpg, by = "Minimum") %>%
  ## add fish IDs 
  left_join(., fish_data, by = "Track.Name") %>%
  rename(Fish = LTMPno) %>%
  ## remove % signs and calculate prop meth
  mutate(Reference.Frequency = sub("\\%.*","\\", Reference.Frequency)) %>%
  mutate(Reference.Frequency = as.numeric(as.character(Reference.Frequency))) %>%
  mutate(prop_meth = Reference.Frequency/100) %>%
  ## remove barra from low frequency catch locations 
  filter(., Fish!="#251") %>%
  filter(., Fish!="#254") %>%
  filter(., Fish!="#110") %>%
  ## add sampling location info
  left_join(., barra_data %>% 
              select(JulieRegion, SexCode, Sex, tl, 
                     yr, AgeClass, YearClass, MonthCaught, LTMPno) %>%
              rename(Fish = LTMPno), 
            by = "Fish") %>%
  ## remove intersex? fish
  filter(., Sex!="Transitional") %>% 
  droplevels(.) %>%
  ## delete unwanted columns
  select(-c(Variant.Raw.Frequency, Variant.Frequency, 
            Name, Change, Type, Track.Name)) %>%
  ## rename some others
  rename(Total_length = tl, 
         Percent_methylation = Reference.Frequency) %>%
  ## create size bins
  mutate(Length_binned = cut(Total_length, 
                             breaks = c(0,60,70,80,90,100,120), 
                             labels = c("50-60","60-70","70-80","80-90","90-100","v100"))) %>%
  ## make factors factors
  mutate(Sex = as.factor(Sex), 
         Gene = as.factor(Gene), 
         JulieRegion = as.factor(JulieRegion), 
         CpG_site = Position) %>%
  ## edit region and sex names so that east coast males come first alphabetically
  mutate(JulieRegionEC = gsub("Mid-northern GoC", "MidNorthernGulf", JulieRegion) %>%
           gsub("Southern GoC", "SouthernGulf", .) %>%
           gsub("Wet-tropics East Coast", "EastCoastWetTropics", .)) %>%
  mutate(SexM = gsub("Female", "xFemale", Sex)) %>%
  ## alternatively, just reorder the factor here (although df edits can cause headaches)
  mutate(JulieRegionEC = as.factor(JulieRegionEC)) %>%
  mutate(SexM = as.factor(SexM))

## STATS ----

## RUN BETA REGRESSION MODELS ----

for (g in levels(data$Gene)) {
  data_gene <- filter(data, Gene==g)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionEC+SexM+Total_length+CpG_site:SexM+
                                                       JulieRegionEC:SexM+JulieRegionEC:Total_length+
                                                       SexM:Total_length+CpG_site:SexM:Total_length+JulieRegionEC:SexM:Total_length), 
                              link="logit", data=data_gene)) ## same for all
  ## CpG automatically gets set as #1
  SumBetaRegMod_OR_Coef_temp<-as.data.frame(summary(betareg_mod_temp)$coefficients) ## just coefficients
  frmla <- formula(betareg_mod_temp$formula)
  model_formula <- c(paste(frmla[2], frmla[1], frmla[3]))
  model_name <- paste("SumBetaRegMod_OR_Coef", g, sep=("_"))
  print(model_name)
  print("pseudo.r.squared")
  print(betareg_mod_temp$pseudo.r.squared)
  data_subset_name <- c(paste(g))
  assign(model_name, cbind(SumBetaRegMod_OR_Coef_temp, data_subset_name, model_formula))
  ## Anova
  AnovaBetaRegMod_OR_temp<-Anova(betareg_mod_temp)
  Anova_model_name <-paste("AnovaBetaRegMod_OR", g, sep=("_"))
  assign(Anova_model_name, cbind(AnovaBetaRegMod_OR_temp, data_subset_name, model_formula))
  
  ## by sex (within gene)
  for (s in levels(data_gene$Sex)) {
    data_sex <- filter(data_gene, Sex==s) 
    data_sex$CpG_site <- factor(data_sex$CpG_site)
    data_sex$Sex <- factor(data_sex$Sex)
    ## no significant 3 way interactions for cyp19a1a sex model, so 2-way
    betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionEC+Total_length+JulieRegionEC:Total_length), 
                                link="logit", data=data_sex))
    SumBetaRegMod_OR_Coef_temp<-as.data.frame(summary(betareg_mod_temp)$coefficients)
    frmla <- formula(betareg_mod_temp$formula)
    model_formula <- c(paste(frmla[2], frmla[1], frmla[3]))
    model_name <- paste("SumBetaRegMod_OR_Coef", g, s, sep=("_"))
    print(model_name)
    print(betareg_mod_temp$pseudo.r.squared)
    data_subset_name <- c(paste(g, s, sep=" "))
    assign(model_name, cbind(SumBetaRegMod_OR_Coef_temp, data_subset_name, model_formula))
    ## Anova
    AnovaBetaRegMod_OR_temp<-Anova(betareg_mod_temp)
    Anova_model_name <-paste("AnovaBetaRegMod_OR", g, s, sep=("_"))
    assign(Anova_model_name, cbind(AnovaBetaRegMod_OR_temp, data_subset_name, model_formula))
  }
}

## CREATE SUMMARY TABLE ----

##  If i had time, I would do this with greater elegance... **

## cyp19a1a
SumBetaRegMod_OR_Coefs_cyp19a1a <- rbind(SumBetaRegMod_OR_Coef_cyp19a1a,
                                       SumBetaRegMod_OR_Coef_cyp19a1a_Female,
                                       SumBetaRegMod_OR_Coef_cyp19a1a_Male)
## correct p-values (by gene)
SumBetaRegMod_OR_Coefs_cyp19a1a$BH <- p.adjust(SumBetaRegMod_OR_Coefs_cyp19a1a$"mean.Pr...z..", 
                                               method = "BH")

## esr1
SumBetaRegMod_OR_Coefs_esr1 <- rbind(SumBetaRegMod_OR_Coef_esr1,
                                     SumBetaRegMod_OR_Coef_esr1_Female,
                                     SumBetaRegMod_OR_Coef_esr1_Male)

## correct p-values (by gene)
SumBetaRegMod_OR_Coefs_esr1$BH <- p.adjust(SumBetaRegMod_OR_Coefs_esr1$"mean.Pr...z..", 
                                           method = "BH")

## nr5a2
SumBetaRegMod_OR_Coefs_nr5a2 <- rbind(SumBetaRegMod_OR_Coef_nr5a2,
                                      SumBetaRegMod_OR_Coef_nr5a2_Female,
                                      SumBetaRegMod_OR_Coef_nr5a2_Male)
## correct p-values (by gene)
SumBetaRegMod_OR_Coefs_nr5a2$BH <- p.adjust(SumBetaRegMod_OR_Coefs_nr5a2$"mean.Pr...z..", 
                                            method = "BH")

## dmrt1
SumBetaRegMod_OR_Coefs_dmrt1 <- rbind(SumBetaRegMod_OR_Coef_dmrt1,
                                      SumBetaRegMod_OR_Coef_dmrt1_Female,
                                      SumBetaRegMod_OR_Coef_dmrt1_Male)

## correct p-values (by gene)
SumBetaRegMod_OR_Coefs_dmrt1$BH <- p.adjust(SumBetaRegMod_OR_Coefs_dmrt1$"mean.Pr...z..", 
                                            method = "BH")

SumBetaRegMod_OR_Coefs<-rbind(SumBetaRegMod_OR_Coefs_cyp19a1a, 
                              SumBetaRegMod_OR_Coefs_esr1, 
                              SumBetaRegMod_OR_Coefs_nr5a2, 
                              SumBetaRegMod_OR_Coefs_dmrt1)

##ADD SIGNIFICANCE CODES
SumBetaRegMod_OR_Coefs$sig_code <- ifelse(SumBetaRegMod_OR_Coefs$"mean.Pr...z.." <=0.001, "***",
                                          ifelse(SumBetaRegMod_OR_Coefs$"mean.Pr...z.." <=0.01, "**",
                                                 ifelse(SumBetaRegMod_OR_Coefs$"mean.Pr...z.." <=0.05, "*",
                                                        ifelse(SumBetaRegMod_OR_Coefs$"mean.Pr...z.." <=0.1, ".",
                                                               " "))))

## ADD BH CORRECTED SIGNIFICANCE CODES 
SumBetaRegMod_OR_Coefs$sig_code2 <- ifelse(SumBetaRegMod_OR_Coefs$BH <=0.001, "***",
                                           ifelse(SumBetaRegMod_OR_Coefs$BH <=0.01, "**",
                                                  ifelse(SumBetaRegMod_OR_Coefs$BH <=0.05, "*",
                                                         ifelse(SumBetaRegMod_OR_Coefs$BH <=0.1, ".",
                                                                " "))))


## EDIT OUTPUT

## delete phi estimate values
DropColumns <- c("precision.Estimate", "precision.Std..Error", "precision.z.value",
                 "precision.Pr...z..")
SumBetaRegMod_OR_Coefs <- SumBetaRegMod_OR_Coefs[ , !(names(SumBetaRegMod_OR_Coefs) %in% DropColumns)]

## rename remaining columns 
ColNames <- c("Estimate", "Std. Error", "Z value", "Pr(>|z|)", "Data subset", "Beta regression model", 
              "BH adjusted p-values", "Significance code", "Adjusted significance code")
SumBetaRegMod_OR_Coefs <- setNames(SumBetaRegMod_OR_Coefs, ColNames)

## reorder columns
SumBetaRegMod_OR_Coefs <- SumBetaRegMod_OR_Coefs[c("Data subset", "Beta regression model", "Estimate", 
                                                   "Std. Error", "Z value", "Pr(>|z|)", "Significance code",  
                                                   "BH adjusted p-values", "Adjusted significance code")]

##  EXPORT CSV
write.csv(SumBetaRegMod_OR_Coefs, file="table_summary_betareg_mods.csv")

## CREATE ANOVA TABLE ----

## cyp19a1a
AnovaBetaRegMod_ORs_cyp19a1a <- rbind(AnovaBetaRegMod_OR_cyp19a1a,
                                    AnovaBetaRegMod_OR_cyp19a1a_Female,
                                    AnovaBetaRegMod_OR_cyp19a1a_Male)
## correct p-values (by gene)
AnovaBetaRegMod_ORs_cyp19a1a$BH <- p.adjust(AnovaBetaRegMod_ORs_cyp19a1a$"Pr(>Chisq)", method = "BH")

## esr1
AnovaBetaRegMod_ORs_esr1 <- rbind(AnovaBetaRegMod_OR_esr1,
                                  AnovaBetaRegMod_OR_esr1_Female,
                                  AnovaBetaRegMod_OR_esr1_Male)

## correct p-values (by gene)
AnovaBetaRegMod_ORs_esr1$BH <- p.adjust(AnovaBetaRegMod_ORs_esr1$"Pr(>Chisq)", method = "BH")

## nr5a2
AnovaBetaRegMod_ORs_nr5a2 <- rbind(AnovaBetaRegMod_OR_nr5a2,
                                   AnovaBetaRegMod_OR_nr5a2_Female,
                                   AnovaBetaRegMod_OR_nr5a2_Male)
## correct p-values (by gene)
AnovaBetaRegMod_ORs_nr5a2$BH <- p.adjust(AnovaBetaRegMod_ORs_nr5a2$"Pr(>Chisq)", method = "BH")

## dmrt1
AnovaBetaRegMod_ORs_dmrt1 <- rbind(AnovaBetaRegMod_OR_dmrt1,
                                   AnovaBetaRegMod_OR_dmrt1_Female,
                                   AnovaBetaRegMod_OR_dmrt1_Male)

## correct p-values (by gene)
AnovaBetaRegMod_ORs_dmrt1$BH <- p.adjust(AnovaBetaRegMod_ORs_dmrt1$"Pr(>Chisq)", method = "BH")

AnovaBetaRegMod_ORs<-rbind(AnovaBetaRegMod_ORs_cyp19a1a, 
                           AnovaBetaRegMod_ORs_esr1, 
                           AnovaBetaRegMod_ORs_nr5a2, 
                           AnovaBetaRegMod_ORs_dmrt1)



## ADD SIGNIFICANCE CODES
AnovaBetaRegMod_ORs$sig_code <- ifelse(AnovaBetaRegMod_ORs$"Pr(>Chisq)" <=0.001, "***",
                                       ifelse(AnovaBetaRegMod_ORs$"Pr(>Chisq)" <=0.01, "**",
                                              ifelse(AnovaBetaRegMod_ORs$"Pr(>Chisq)" <=0.05, "*",
                                                     ifelse(AnovaBetaRegMod_ORs$"Pr(>Chisq)" <=0.1, ".",
                                                            " "))))


## ADD CORRECTED SIGNIFICANCE CODES
AnovaBetaRegMod_ORs$sig_code_BH <- ifelse(AnovaBetaRegMod_ORs$BH <=0.001, "***",
                                          ifelse(AnovaBetaRegMod_ORs$BH <=0.01, "**",
                                                 ifelse(AnovaBetaRegMod_ORs$BH <=0.05, "*",
                                                        ifelse(AnovaBetaRegMod_ORs$BH <=0.1, ".",
                                                               " "))))


## EXPORT CSV
write.csv(AnovaBetaRegMod_ORs, file="table_anova_betareg_mods.csv")

## INDIVIDUAL LINEAR HYPOTHESIS TESTING ----

## SEPARATE BY SEX

## east coast as baseline

temp_lh_table_EC=NULL

for (g in levels(data$Gene)) {
  data_gene <- filter(data, Gene==g)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  for (s in levels(data_gene$Sex)) {
    data_sex <- filter(data_gene, Sex==s) 
    data_sex$CpG_site <- factor(data_sex$CpG_site)
    data_sex$Sex <- factor(data_sex$Sex)
    ## no significant 3 way interactions for cyp19a1a sex model, so 2-way
    betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionEC+Total_length+JulieRegionEC:Total_length), 
                                link="logit", data=data_sex))
    coefcompare<-paste(round(summary(betareg_mod_temp)$coefficients$mean["JulieRegionECSouthernGulf", "Estimate"], digits=4), 
                       round(summary(betareg_mod_temp)$coefficients$mean["JulieRegionECMidNorthernGulf", "Estimate"], digits=4), 
                       sep=" = ")
    compare<-"JulieRegionECSouthernGulf = JulieRegionECMidNorthernGulf" # create character vector giving the hypothesis in symbolic form for linearHypothesis
    model_name <- paste(g, s, compare, sep=("_"))
    ## test individual linear hypothesis
    ## NB cant test linear hypothesis for the baseline (e.g. anything with EC, -83 (cyp19a1a), male)
    temp_lh<-(linearHypothesis(betareg_mod_temp, compare))
    temp_lh2<-cbind(s, g, compare, coefcompare, temp_lh)
    temp_lh_table_EC<-rbind(temp_lh_table_EC, temp_lh2)
  }
}


## MNG as baseline 
data$JulieRegionMNG <- data$JulieRegion
data$JulieRegionMNG <- gsub("Mid-northern GoC", "MidNorthernGulf", data$JulieRegionMNG)
data$JulieRegionMNG <- gsub("Southern GoC", "SouthernGulf", data$JulieRegionMNG)
data$JulieRegionMNG <- gsub("Wet-tropics East Coast", "WetTropicsEastCoast", data$JulieRegionMNG)
data$JulieRegionMNG=as.factor(data$JulieRegionMNG)
levels(data$JulieRegionMNG)

temp_lh_table_MNG=NULL

for (g in levels(data$Gene)) {
  data_gene <- filter(data, Gene==g)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  for (s in levels(data_gene$Sex)) {
    data_sex <- filter(data_gene, Sex==s) 
    data_sex$CpG_site <- factor(data_sex$CpG_site)
    data_sex$Sex <- factor(data_sex$Sex)
    ## no significant 3 way interactions for cyp19a1a sex model, so 2-way
    betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionMNG+Total_length+JulieRegionMNG:Total_length), 
                                link="logit", data=data_sex))
    coefcompare<-paste(round(summary(betareg_mod_temp)$coefficients$mean["JulieRegionMNGSouthernGulf", "Estimate"], digits=4), 
                       round(summary(betareg_mod_temp)$coefficients$mean["JulieRegionMNGWetTropicsEastCoast", "Estimate"], digits=4),  
                       sep=" = ")
    compare<-"JulieRegionMNGSouthernGulf = JulieRegionMNGWetTropicsEastCoast"
    model_name <- paste(s, g, compare, sep=("_"))
    ## test individual linear hypothesis
    ## NB cant test linear hypothesis for the baseline (e.g. anything with EC, -83 (cyp19a1a), male)
    temp_lh<-(linearHypothesis(betareg_mod_temp, compare))
    temp_lh2<-cbind(g, s, compare, coefcompare, temp_lh)
    temp_lh_table_MNG<-rbind(temp_lh_table_MNG, temp_lh2)
  }
}

## SG as baseline
data$JulieRegionSG <- data$JulieRegionEC
data$JulieRegionSG <- gsub("SouthernGulf", "ASouthernGulf", data$JulieRegionSG)
data$JulieRegionSG <- as.factor(data$JulieRegionSG)

levels(data$JulieRegionSG)

temp_lh_table_SG=NULL

for (g in levels(data$Gene)) {
  data_gene <- filter(data, Gene==g)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  for (s in levels(data_gene$Sex)) {
    data_sex <- filter(data_gene, Sex==s) 
    data_sex$CpG_site <- factor(data_sex$CpG_site)
    data_sex$Sex <- factor(data_sex$Sex)
    ## no significant 3 way interactions for cyp19a1a sex model, so 2-way
    betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionSG+Total_length+JulieRegionSG:Total_length), 
                                link="logit", data=data_sex))
    coefcompare<-paste(round(summary(betareg_mod_temp)$coefficients$mean["JulieRegionSGMidNorthernGulf", "Estimate"], digits=4), 
                       round(summary(betareg_mod_temp)$coefficients$mean["JulieRegionSGEastCoastWetTropics", "Estimate"], digits=4), 
                       sep=" = ")
    compare<-"JulieRegionSGMidNorthernGulf = JulieRegionSGEastCoastWetTropics"
    model_name <- paste(g, s, compare, sep=("_"))
    ## test individual linear hypothesis
    ## NB cant test linear hypothesis for the baseline (e.g. anything with EC, -83 (cyp19a1a), male)
    temp_lh<-(linearHypothesis(betareg_mod_temp, compare))
    temp_lh2<-cbind(s, g, compare, coefcompare, temp_lh)
    temp_lh_table_SG<-rbind(temp_lh_table_SG, temp_lh2)
  }
}

lh_table_all<-rbind(temp_lh_table_EC, temp_lh_table_MNG, temp_lh_table_SG)
lh_table_all<-na.omit(lh_table_all)
lh_table_all$PrChisqBH <- p.adjust(lh_table_all$"Pr(>Chisq)", method = "BH")


## ADD SIGNIFICANCE CODES
lh_table_all$BHsig_code <- ifelse(lh_table_all$PrChisqBH <=0.001, "***",
                                          ifelse(lh_table_all$PrChisqBH <=0.01, "**",
                                                 ifelse(lh_table_all$PrChisqBH <=0.05, "*",
                                                        ifelse(lh_table_all$PrChisqBH <=0.1, ".",
                                                               " "))))

lh_table_all$"FDR adjusted p-value"<-paste(round(lh_table_all$PrChisqBH, digits = 4), lh_table_all$BHsig_code, sep="")
lh_table_all<-lh_table_all[,-c(8:10)]
lh_table_all$Chisq<-round(lh_table_all$Chisq, digits=4)

names(lh_table_all)[(names(lh_table_all) == "s")] <- "Sex" 
names(lh_table_all)[(names(lh_table_all) == "g")] <- "Amplicon" 
names(lh_table_all)[(names(lh_table_all) == "compare")] <- "Comparison" 
names(lh_table_all)[(names(lh_table_all) == "coefcompare")] <- "Comparison coefficients"
names(lh_table_all)[(names(lh_table_all) == "Res.Df")] <- "Residual degrees of freedom"
names(lh_table_all)[(names(lh_table_all) == "Df")] <- "Degrees of freedom"
names(lh_table_all)[(names(lh_table_all) == "Chisq")] <- "Chi-square statistic"

lh_table_all$Comparison<-gsub("Julie.*SouthernGulf", "southern GoC", lh_table_all$Comparison) ##'.' for character '*' for any number of character
lh_table_all$Comparison<-gsub("Julie.*MidNorthernGulf", "mid-northern GoC", lh_table_all$Comparison)
lh_table_all$Comparison<-gsub("Julie.*EastCoastWetTropics", "north Qld east coast", lh_table_all$Comparison)
lh_table_all$Comparison<-gsub("Julie.*WetTropicsEastCoast", "north Qld east coast", lh_table_all$Comparison)

lh_table_all<-
  lh_table_all %>% 
  arrange(desc(Sex),
          Amplicon,
          factor(lh_table_all$Comparison, levels = c("southern GoC = mid-northern GoC", 
                                        "southern GoC = north Qld east coast", 
                                        "mid-northern GoC = north Qld east coast")))

names(lh_table_all)[(names(lh_table_all) == "Comparison")] <- "Comparison (null hypothesis)"

write.csv(lh_table_all, "table_linear_hypoth.csv", row.names=FALSE) 

## BETAREG PLOTS (length) ----

## create small legend function
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)))+#,
    #color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = 10), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}


## define CpG 
cpg<-"CpG4" #"CpG#" you want to look at
n<-4 #as above 
# run loop
for (g in levels(data$Gene)) { #subset data by gene
  data_gene <- filter(data, Gene==g)
  data_gene <- droplevels(data_gene)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  #by sex (within gene)
  for (s in levels(data_gene$Sex)) { # subset data by sex ####
    data_sex <- filter(data_gene, Sex==s) 
    data_sex <- droplevels(data_sex)
    data_sex$CpG_site <- factor(data_sex$CpG_site)
    data_sex$Sex <- factor(data_sex$Sex)
    # model data using betareg package ####
    betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionEC+Total_length+JulieRegionEC:Total_length), 
                                link="logit", data=data_sex))
    # create titles ####
    p1_ttle<-paste(g, s, sep="_")#, "CpG4")  
    CpG_legendttle<-paste(g, "CpG_site", sep="_")
    # predict betareg ####
    # CpG bit ####
    c<-data_sex$CpG_site[n] #for selecting CpG 
    # make fake data #
    fakedataEC = data.frame(Total_length = 55:105, 
                            CpG_site = c,
                            JulieRegionEC = 'EastCoastWetTropics') 
    fakedataMNG = data.frame(Total_length = 55:105, 
                             CpG_site = c,
                             JulieRegionEC = 'MidNorthernGulf')
    fakedataSG = data.frame(Total_length = 55:105, 
                            CpG_site = c,
                            JulieRegionEC = 'SouthernGulf')
    # make predictions based on fake data #
    predicted_MNG <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                   newdata = fakedataMNG), 
                                TL.pred=55:105) #MNG
    predicted_SG <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                  newdata = fakedataSG), 
                               TL.pred=55:105) #SG
    predicted_EC <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                  newdata = fakedataEC), 
                               TL.pred=55:105) #EC  
    # create "p1" plot with CpG legend only and assign p1_ttle ("gene_sex") ####
    p<-ggplot(data = data_sex, 
              aes(Total_length, prop_meth, color=JulieRegionEC))+
      geom_point(aes(shape=CpG_site), size=2)+
      scale_shape_manual(values=1:nlevels(data_sex$CpG_site))+
      scale_color_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"))+ # colour brewer palette Set2
      guides(color="none", shape=guide_legend(nrow=1))+ 
      labs(shape = "CpG\nsite")+ 
      geom_line(color='#66c2a5', data = predicted_EC, aes(x=TL.pred, y=RFP.pred))+
      geom_line(color='#fc8d62', data = predicted_MNG, aes(x=TL.pred, y=RFP.pred))+
      geom_line(color='#8da0cb', data = predicted_SG, aes(x=TL.pred, y=RFP.pred))+
      ggtitle(paste(s, "fish", sep=" "))+
      scale_x_continuous(limits=c(50, 110), "Total length (cm)",
                         breaks=c(50, 60, 70, 80, 90, 100, 110))+
      scale_y_continuous("Proportion methylated", 
                         sec.axis = sec_axis(~., name = g))+
      theme_bw()+
      theme(axis.text.y.right = element_blank(), 
            axis.ticks.y.right = element_blank(),
            axis.title.y.right = element_text(face = "italic"), 
            plot.title = element_text(hjust = 0.5, size=12),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
            legend.position="right")
    p1<-addSmallLegend(p)
    assign(p1_ttle, p1) ## assign gene_sex_nolegends title to plot
    ## create "pl2" plot with region legend only and extract "l2" legend (once per everything...) ####
    pl2<-ggplot(data = data_sex, 
                aes(Total_length, prop_meth, color=JulieRegionEC))+
      geom_point()+
      scale_shape_manual(values=1:nlevels(data_sex$CpG_site))+
      geom_line(color='#66c2a5', data = predicted_EC, aes(x=TL.pred, y=RFP.pred))+
      geom_line(color='#fc8d62', data = predicted_MNG, aes(x=TL.pred, y=RFP.pred))+
      geom_line(color='#8da0cb', data = predicted_SG, aes(x=TL.pred, y=RFP.pred))+
      xlab("Total length (cm)")+
      scale_x_continuous(limits=c(50, 110), 
                         breaks=c(50, 60, 70, 80, 90, 100, 110))+
      scale_y_continuous(limits=c(0, 1))+
      scale_colour_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"), ## colour brewer palette Set2
                          labels = c("North Qld east coast", "Mid-northern GoC", "Southern GoC"))+ ## change legend names
      theme_bw()+
      labs(color = "Region")+ ## change legend title
      theme(legend.position = "bottom", 
            legend.margin=margin(c(0,0,0,0)))
    l2<-get_legend(pl2) ## extract legend 2 
  }
}

## remove redundant axis information

nr5a2_Males<-nr5a2_Male+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(), 
                              axis.ticks.x=element_blank(), 
                              axis.title.y.right=element_blank())

nr5a2_Females<-nr5a2_Female+theme(axis.title.x=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(),
                                  axis.title.y=element_blank())

dmrt1_Males<-dmrt1_Male+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(), 
                              axis.ticks.x=element_blank(), 
                              axis.title.y.right=element_blank(), 
                              plot.title=element_blank())

dmrt1_Females<-dmrt1_Female+theme(axis.title.x=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(), 
                                  axis.title.y=element_blank(),
                                  plot.title=element_blank())

cyp19a1a_Males<-cyp19a1a_Male+theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(), 
                                    axis.ticks.x=element_blank(), 
                                    axis.title.y.right=element_blank(), 
                                    plot.title=element_blank())

cyp19a1a_Females<-cyp19a1a_Female+theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank(), 
                                        axis.title.y=element_blank(),
                                        plot.title=element_blank())

esr1_Males<-esr1_Male+theme(axis.title.y.right=element_blank(), 
                            plot.title=element_blank())

esr1_Females<-esr1_Female+theme(axis.title.y=element_blank(),
                                plot.title=element_blank())

all_genes <- ggarrange(ggarrange(nr5a2_Males, nr5a2_Females, common.legend = TRUE, legend="right"),
                       ggarrange(dmrt1_Males, dmrt1_Females, common.legend = TRUE, legend="right"),
                       ggarrange(cyp19a1a_Males, cyp19a1a_Females, common.legend = TRUE, legend="right"),
                       ggarrange(esr1_Males, esr1_Females, common.legend = TRUE, legend="right"),
                       widths = c(19.4, 18),
                       heights = c(10.5, 9, 9, 11),
                       nrow=4, ncol=1)

## add legend 
all_genes_leg <- plot_grid(all_genes, l2, ## add region legend 
                           nrow = 2,
                           rel_heights = c(14.5, 0.5)) ## change me
all_genes_leg

## save plot
save_plot("figure_all_nolims.pdf", 
          all_genes_leg,
          base_height = 29/3,
          base_width = 21/2.5)

## BETAREG PLOTS (age) ----

## define CpG   
cpg<-"CpG4"
n<-4 ## as above 
## run loop 
for (g in levels(data$Gene)) { ## subset data by gene
  data_gene <- filter(data, Gene==g)
  data_gene <- droplevels(data_gene)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  ## by sex (within gene)
  for (s in levels(data_gene$Sex)) { ## subset data by sex 
    data_sex <- filter(data_gene, Sex==s) 
    data_sex <- droplevels(data_sex)
    data_sex$CpG_site <- factor(data_sex$CpG_site)
    data_sex$Sex <- factor(data_sex$Sex)
    betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionEC+AgeClass+JulieRegionEC:AgeClass), 
                                link="logit", data=data_sex))
    p1_ttle<-paste(g, s, sep="_")
    CpG_legendttle<-paste(g, "CpG_site", sep="_")
    c<-data_sex$CpG_site[n]
    fakedataEC = data.frame(AgeClass = 2:15, 
                            CpG_site = c,
                            JulieRegionEC = 'EastCoastWetTropics') 
    fakedataMNG = data.frame(AgeClass = 2:15, 
                             CpG_site = c,
                             JulieRegionEC = 'MidNorthernGulf')
    fakedataSG = data.frame(AgeClass = 2:15, 
                            CpG_site = c,
                            JulieRegionEC = 'SouthernGulf')
    predicted_MNG <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                   newdata = fakedataMNG), 
                                A.pred=2:15) #MNG
    predicted_SG <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                  newdata = fakedataSG), 
                               A.pred=2:15) #SG
    predicted_EC <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                  newdata = fakedataEC), 
                               A.pred=2:15) #EC  
    p<-ggplot(data = data_sex, 
              aes(AgeClass, prop_meth, color=JulieRegionEC))+
      geom_point(aes(shape=CpG_site), size=2)+
      scale_shape_manual(values=1:nlevels(data_sex$CpG_site))+
      scale_color_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"))+ 
      guides(color="none", shape=guide_legend(nrow=1))+ 
      labs(shape = "CpG\nsite")+ 
      geom_line(color='#66c2a5', data = predicted_EC, aes(x=A.pred, y=RFP.pred))+
      geom_line(color='#fc8d62', data = predicted_MNG, aes(x=A.pred, y=RFP.pred))+
      geom_line(color='#8da0cb', data = predicted_SG, aes(x=A.pred, y=RFP.pred))+
      ggtitle(paste(s, "fish", sep=" "))+
      scale_x_continuous(limits=c(2, 15), "Age (years)"
                         )+
      scale_y_continuous("Proportion methylated", 
                         sec.axis = sec_axis(~., name = g))+
      theme_bw()+
      theme(axis.text.y.right = element_blank(), 
            axis.ticks.y.right = element_blank(),
            axis.title.y.right = element_text(face = "italic"), 
            plot.title = element_text(hjust = 0.5, size=12),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
            legend.position="right")
    p1<-addSmallLegend(p)
    assign(p1_ttle, p1)
    pl2<-ggplot(data = data_sex, 
                aes(AgeClass, prop_meth, color=JulieRegionEC))+
      geom_point()+
      scale_shape_manual(values=1:nlevels(data_sex$CpG_site))+
      geom_line(color='#66c2a5', data = predicted_EC, aes(x=A.pred, y=RFP.pred))+
      geom_line(color='#fc8d62', data = predicted_MNG, aes(x=A.pred, y=RFP.pred))+
      geom_line(color='#8da0cb', data = predicted_SG, aes(x=A.pred, y=RFP.pred))+
      xlab("Age (years)")+
      scale_x_continuous(limits=c(2, 15)
                         )+
      scale_y_continuous(limits=c(0, 1))+
      scale_colour_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"),
                          labels = c("North Qld east coast", "Mid-northern GoC", "Southern GoC"))+
      theme_bw()+
      labs(color = "Region")+
      theme(legend.position = "bottom", 
            legend.margin=margin(c(0,0,0,0)))
    l2<-get_legend(pl2)
  }
}

## remove redundant axis information

nr5a2_Males<-nr5a2_Male+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(), 
                              axis.ticks.x=element_blank(), 
                              axis.title.y.right=element_blank())

nr5a2_Females<-nr5a2_Female+theme(axis.title.x=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(),
                                  axis.title.y=element_blank())

dmrt1_Males<-dmrt1_Male+theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(), 
                              axis.ticks.x=element_blank(), 
                              axis.title.y.right=element_blank(), 
                              plot.title=element_blank())

dmrt1_Females<-dmrt1_Female+theme(axis.title.x=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(), 
                                  axis.title.y=element_blank(),
                                  plot.title=element_blank())

cyp19a1a_Males<-cyp19a1a_Male+theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(), 
                                    axis.ticks.x=element_blank(), 
                                    axis.title.y.right=element_blank(), 
                                    plot.title=element_blank())

cyp19a1a_Females<-cyp19a1a_Female+theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank(), 
                                        axis.title.y=element_blank(),
                                        plot.title=element_blank())

esr1_Males<-esr1_Male+theme(axis.title.y.right=element_blank(), 
                            plot.title=element_blank())

esr1_Females<-esr1_Female+theme(axis.title.y=element_blank(),
                                plot.title=element_blank())

all_genes <- ggarrange(ggarrange(nr5a2_Males, nr5a2_Females, common.legend = TRUE, legend="right"),
                       ggarrange(dmrt1_Males, dmrt1_Females, common.legend = TRUE, legend="right"),
                       ggarrange(cyp19a1a_Males, cyp19a1a_Females, common.legend = TRUE, legend="right"),
                       ggarrange(esr1_Males, esr1_Females, common.legend = TRUE, legend="right"),
                       widths = c(19.4, 18),
                       heights = c(10.5, 9, 9, 11),
                       nrow=4, ncol=1)

all_genes_leg <- plot_grid(all_genes, l2,
                           nrow = 2,
                           rel_heights = c(14.5, 0.5))
all_genes_leg

## save plot
save_plot("figure_all_nolims_age.pdf", 
          all_genes_leg,
          base_height = 29/3,
          base_width = 21/2.5)

## BETAREG PLOTS (plot by gene, colour by region) ----

for (g in levels(data$Gene)) {
  data_gene <- filter(data, Gene==g)
  data_gene <- droplevels(data_gene)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  betareg_mod_temp <-(betareg(prop_meth~(CpG_site+JulieRegionEC+Total_length+JulieRegionEC:Total_length), 
                              link="logit", data=data_gene))
  p1_ttle<-paste(g, sep="_")
  CpG_legendttle<-paste(g, "CpG_site", sep="_")
  c<-data_gene$CpG_site[n]
  fakedataEC = data.frame(Total_length = 55:105, 
                          CpG_site = c,
                          JulieRegionEC = 'EastCoastWetTropics')
  fakedataMNG = data.frame(Total_length = 55:105, 
                           CpG_site = c,
                           JulieRegionEC = 'MidNorthernGulf')
  fakedataSG = data.frame(Total_length = 55:105, 
                          CpG_site = c,
                          JulieRegionEC = 'SouthernGulf')
  predicted_MNG <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                 newdata = fakedataMNG), 
                              TL.pred=55:105) #MNG
  predicted_SG <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                newdata = fakedataSG), 
                             TL.pred=55:105) #SG
  predicted_EC <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                newdata = fakedataEC), 
                             TL.pred=55:105) #EC  
  p<-ggplot(data = data_gene, 
            aes(Total_length, prop_meth, color=JulieRegionEC))+
    geom_point(aes(shape=CpG_site), size=2)+
    labs(shape = "CpG site")+ 
    scale_shape_manual(values=1:nlevels(data_gene$CpG_site))+
    scale_color_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"))+ 
    guides(color="none", shape=guide_legend(nrow=1))+ 
    geom_line(color='#66c2a5', data = predicted_EC, aes(x=TL.pred, y=RFP.pred))+
    geom_line(color='#fc8d62', data = predicted_MNG, aes(x=TL.pred, y=RFP.pred))+
    geom_line(color='#8da0cb', data = predicted_SG, aes(x=TL.pred, y=RFP.pred))+
    xlab("Total length (cm)")+
    scale_x_continuous(limits=c(50, 110), "Total length (cm)",
                       breaks=c(50, 60, 70, 80, 90, 100, 110))+
    scale_y_continuous("Proportion methylated", limits=c(0, 1),
                       sec.axis = sec_axis(~., name = g))+
    theme_bw()+
    theme(axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          axis.title.y.right = element_text(face = "italic"), 
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          legend.position="bottom", 
          legend.margin=margin(c(0,0,0,0)))
  p1<-addSmallLegend(p)
  assign(p1_ttle, p1)
  pl2<-ggplot(data = data_gene, 
              aes(Total_length, prop_meth, color=JulieRegionEC))+
    geom_point()+
    scale_shape_manual(values=1:nlevels(data_gene$CpG_site))+
    geom_line(color='#66c2a5', data = predicted_EC, aes(x=TL.pred, y=RFP.pred))+
    geom_line(color='#fc8d62', data = predicted_MNG, aes(x=TL.pred, y=RFP.pred))+
    geom_line(color='#8da0cb', data = predicted_SG, aes(x=TL.pred, y=RFP.pred))+
    xlab("Total length (cm)")+
    ylab("Proportion Methlyation")+
    scale_x_continuous(limits=c(50, 110), 
                       breaks=c(50, 60, 70, 80, 90, 100, 110))+
    scale_y_continuous(limits=c(0, 1))+
    scale_colour_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"), 
                        labels = c("North Qld east coast", "Mid-northern GoC", "Southern GoC"))+ 
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position="bottom", 
          legend.margin=margin(c(0,0,0,0)))
  l2<-get_legend(pl2)
}


all_genes_both_sexes <- ggarrange(dmrt1+theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank()), 
                                  nr5a2+theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.title.y=element_blank()), 
                                  cyp19a1a, 
                                  esr1+theme(axis.title.y=element_blank()), 
                  widths = c(19.4, 18),
                  heights = c(9, 10),
                  nrow=2, ncol=2)

all_genes_both_sexes

all_genes_both_sexes_leg <- plot_grid(all_genes_both_sexes, l2,
                           nrow = 2,
                           rel_heights = c(14.5, 0.5)) 
all_genes_both_sexes_leg

save_plot("facet_gene_colour_region.pdf", 
          all_genes_both_sexes_leg,
          base_height = 29/5,
          base_width = 21/2.5)

## BETAREG PLOTS (plot by gene, colour by sex) ----

for (g in levels(data$Gene)) {
  data_gene <- filter(data, Gene==g)
  data_gene <- droplevels(data_gene)
  data_gene$CpG_site <- factor(data_gene$CpG_site)
  betareg_mod_temp <-(betareg(prop_meth~(CpG_site+Total_length), 
                              link="logit", data=data_gene))
  p1_ttle<-paste(g, sep="_")
  CpG_legendttle<-paste(g, "CpG_site", sep="_")
  c<-data_gene$CpG_site[n] 
  fakedata = data.frame(Total_length = 55:105, 
                          CpG_site = c)
  predicted <- data.frame(RFP.pred = predict(betareg_mod_temp, 
                                                 newdata = fakedata), 
                              TL.pred=55:105)
  p<-ggplot(data = data_gene, 
            aes(Total_length, prop_meth, color=Sex))+
    geom_point(aes(shape=CpG_site), size=2)+
    labs(shape = "CpG site")+
    scale_shape_manual(values=1:nlevels(data_gene$CpG_site))+
    scale_color_manual(values=c("#fc8d62", "#8da0cb"))+
    guides(color="none", shape=guide_legend(nrow=1))+ 
    geom_line(color='black', data = predicted, aes(x=TL.pred, y=RFP.pred))+
    xlab("Total length (cm)")+
    scale_x_continuous(limits=c(50, 110), "Total length (cm)",
                       breaks=c(50, 60, 70, 80, 90, 100, 110))+
    scale_y_continuous("Proportion methylated",
                       sec.axis = sec_axis(~., name = g))+
    theme_bw()+
    theme(axis.text.y.right = element_blank(), 
          axis.ticks.y.right = element_blank(),
          axis.title.y.right = element_text(face = "italic"), 
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          legend.position="right", 
          legend.margin=margin(c(0,0,0,0)))
  p1<-addSmallLegend(p)
  assign(p1_ttle, p1)
  pl2<-ggplot(data = data_gene, 
              aes(Total_length, prop_meth, colour=Sex))+
    scale_colour_manual(values=c("#fc8d62", "#8da0cb"), 
                        labels = c("Male", "Female"))+ 
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position="bottom", 
          legend.margin=margin(c(0,0,0,0)))
  l2<-get_legend(pl2) 
}

all_genes_no_region_or_sex <- ggarrange(nr5a2+theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(), 
                                                    axis.ticks.x=element_blank()),
                                        dmrt1+theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(), 
                                                    axis.ticks.x=element_blank()),
                                        cyp19a1a+theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(), 
                                        axis.ticks.x=element_blank()),
                                  esr1, 
                                  widths = c(19.4, 18),
                                  heights = c(9, 10),
                                  nrow=4, ncol=1)

all_genes_no_region_or_sex

all_genes_no_region_or_sex_leg <- plot_grid(all_genes_no_region_or_sex, 
                                            l2,
                                      nrow = 2,
                                      rel_heights = c(14.5, 0.5))
all_genes_no_region_or_sex_leg

## save plot
save_plot("facet_gene_colour_sex.pdf",
          all_genes_no_region_or_sex_leg,
          base_height = 29/5, ## more than half
          base_width = 21/2.5)

## Mann-Whitney-Wilcoxon ----
## FOR SEX BY GENE

for (i in levels(data$Gene)) {
  data_gene <- filter(data, Gene==i)
  name_temp<-paste(i)
  print(name_temp)
  temp_MWW<-with(data_gene, wilcox.test(prop_meth ~ Sex, alpha=0.05, p.adj="BH")) #Wilcoxon rank sum test (aka Mann-Whiteney U test)
  print(temp_MWW)
  print('p-value')
  print(temp_MWW$p.value)
  Zscore<-qnorm(temp_MWW$p.value)
  print('Zscore')
  print(Zscore)
  rscore<- abs(Zscore)/sqrt(length(data_gene)) ## sometimes calculated as an effect size for M-W test
  print('rscore')
  print(rscore)
  data_gene_F<-filter(data_gene, Sex=='Female') ## female data
  print('Female Median')
  print(median(data_gene_F$prop_meth))
  print(IQR(data_gene_F$prop_meth))
  data_gene_M<-filter(data_gene, Sex=='Male') ## male data
  print('Male Median')
  print(median(data_gene_M$prop_meth))
  print(IQR(data_gene_M$prop_meth))
}

## FOR SEX AND REGION BY GENE

data$RegionSex <- factor(paste0(data$JulieRegion, "-", data$Sex))
data$SexRegion <- factor(paste0(data$Sex, "_", data$JulieRegion))
data$SexRegionND <- gsub("-", ".", data$SexRegion)
data$SexRegionND <- factor(data$SexRegionND,
                           levels=c("Male_Southern GoC", 
                                    "Male_Mid.northern GoC", 
                                    "Male_Wet.tropics East Coast", 
                                    "Female_Southern GoC", 
                                    "Female_Mid.northern GoC", 
                                    "Female_Wet.tropics East Coast"))


PWT_cld <- NULL
for (i in levels(data$Gene)) {
  data_gene <- filter(data, Gene==i)
  name_temp<-paste(i)
  print(name_temp)
  temp_MWW<-with(data_gene, pairwise.wilcox.test(prop_meth, SexRegionND, alpha=0.05, p.adjust.method="BH")) ## Wilcoxon rank sum test (aka Mann-Whiteney U test)
  print(summary(temp_MWW))
  temp_df <- make_cld(temp_MWW)
  temp_df$Gene <- paste(i)
  PWT_cld <- rbind(PWT_cld, temp_df)
}

PWT_cld_tab <- PWT_cld %>%
  mutate(Sex = gsub("_.*", "", group)) %>%
  mutate(Region = gsub(".*_", "", group) %>%
           gsub("\\.", "-", .)) %>%
  mutate(JulieRegions = gsub("Mid-northernGoC", "Mid-northern GoC", Region) %>%
           gsub("SouthernGoC", "Southern GoC", .) %>%
           gsub("Wet-tropicsEastCoast", "North Qld east coast", .)) %>%
  mutate(x = gsub("Male_SouthernGoC", "0.75", group) %>%
           gsub("Male_Mid.northernGoC", "1.0", .) %>%
           gsub("Male_Wet.tropicsEastCoast", "1.25", .) %>%
           gsub("Female_SouthernGoC", "1.75", .) %>%
           gsub("Female_Mid.northernGoC", "2.0", .) %>%
           gsub("Female_Wet.tropicsEastCoast", "2.25", .)) %>%
  mutate(x = as.numeric(x)) %>%
  select(Gene, JulieRegions, Sex, Region, cld, x) %>%
  mutate(JulieRegions_FL = factor(JulieRegions, 
                                 levels = c("Southern GoC", 
                                            "Mid-northern GoC", 
                                            "North Qld east coast")))
    
PWT_cld_tab    

for (i in levels(data$Gene)) {
  print(i)
  temp_aov<-aov(prop_meth ~ SexRegion, data = filter(data, Gene==i))
  print(summary(temp_aov))
  my_glht<-glht(temp_aov, linfct = mcp(SexRegion = "Tukey"), test = adjusted("bonferroni"))
  my_cld<-cld(my_glht)
  my_cld_df<-as.data.frame(my_cld$mcletters$Letters, rownames=NULL)
}

## HYPOTHESIS FIGURE ----
## male/female methylation plot

an_text <- data.frame(label=c("a", "b", "a", "b", "a", "b", "a", "b"), 
                      x = c(1, 2, 1, 2, 1, 2, 1, 2),
                      y = c(1.05, 1.05, 0.35, 0.35, 1.05, 1.05, 0.9, 0.9),
                      Sex = c("Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female"), 
                      Gene = c("cyp19a1a", "cyp19a1a", "dmrt1", "dmrt1", "esr1", "esr1", "nr5a2", "nr5a2"))

## order factors
data$gene_fac<-factor(data$Gene, levels=c("nr5a2", "dmrt1", "cyp19a1a", "esr1"))
data$sex_fac<-factor(data$Sex, levels=c("Male", "Female"))
an_text$gene_fac<-factor(an_text$Gene, levels=c("nr5a2", "dmrt1", "cyp19a1a", "esr1"))
an_text$sex_fac<-factor(an_text$Sex, levels=c("Male", "Female"))

## make long for combine figure
bps3<-ggplot(data, aes(sex_fac, prop_meth, fill=sex_fac))+
  geom_boxplot()+
  geom_jitter(pch=21, alpha=0.3, size=1, width=0.35)+
  facet_wrap(~gene_fac, scales="free", ncol=1, strip.position="right")+
  ylab("Proportion methylated")+
  xlab("Sex")+
  geom_text(an_text, mapping=aes(x=x, y=y, label=label), 
            fontface = 3, size = 3, 
            show.legend=FALSE)+
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text = element_text(face="italic", size=10), 
        plot.title = element_text(hjust = 0.5, size=11))+
  scale_fill_manual(values=c("#8da0cb", "#fc8d62"))+
  ggtitle("A. Initial observation")+
  guides(fill="none")

# bps3

## fake data (hypothesis) plot

## for fake data ranges
for (g in levels(data$gene_fac)) {
  temp_name <- paste0(g, "_nums")
  temp_data <- filter(data, gene_fac == g)
  assign(temp_name, seq(from = min(temp_data$prop_meth), 
                        to = max(temp_data$prop_meth), 
                        by = (max(temp_data$prop_meth-min(temp_data$prop_meth))/10)))
  
}

## create df
fake_data_exp<-data.frame(gene=c("nr5a2", "nr5a2", "nr5a2", "nr5a2", "nr5a2",
                                 "nr5a2", "nr5a2", "nr5a2", "nr5a2", "nr5a2", "nr5a2", 
                                 "dmrt1", "dmrt1", "dmrt1", "dmrt1", "dmrt1",
                                 "dmrt1", "dmrt1", "dmrt1", "dmrt1", "dmrt1", "dmrt1", 
                                 "cyp19a1a", "cyp19a1a", "cyp19a1a", "cyp19a1a", "cyp19a1a",
                                 "cyp19a1a", "cyp19a1a", "cyp19a1a", "cyp19a1a", "cyp19a1a", "cyp19a1a", 
                                 "esr1", "esr1", "esr1", "esr1", "esr1",
                                 "esr1", "esr1", "esr1", "esr1", "esr1", "esr1"), 
                          total_length=c(55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, ## 11
                                         55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 
                                         55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 
                                         55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105), 
                          sex=c("Male", "Male", "Male", "Male", "Male", "Male",
                                "Female", "Female", "Female","Female", "Female",
                                "Male", "Male", "Male", "Male", "Male", "Male",
                                "Female", "Female", "Female","Female", "Female",
                                "Male", "Male", "Male", "Male", "Male", "Male",
                                "Female", "Female", "Female","Female", "Female",
                                "Male", "Male", "Male", "Male", "Male", "Male",
                                "Female", "Female", "Female","Female", "Female"),
                          meth=c(nr5a2_nums, dmrt1_nums, rev(cyp19a1a_nums), rev(esr1_nums)))

## make factors
fake_data_exp$gene_fac<-factor(fake_data_exp$gene, levels=c("nr5a2", "dmrt1", "cyp19a1a", "esr1"))
fake_data_exp$sex_fac<-factor(fake_data_exp$sex, levels=c("Male", "Female"))

## fake data expected points
fake_data_exp_points<-rbind.data.frame(fake_data_exp, fake_data_exp, 
                                       fake_data_exp, fake_data_exp,
                                       fake_data_exp, fake_data_exp,
                                       fake_data_exp, fake_data_exp)

fdep<-ggplot(fake_data_exp_points, aes(total_length, meth))+
  geom_line()+
  geom_jitter(aes(fill=sex_fac), width=2, height=0.05, pch=21, alpha=0.4, size=1.5)+
  scale_fill_manual(values=c("#8da0cb", "#fc8d62"))+
  facet_wrap(~gene_fac, scales="free", ncol=1, strip.position="right")+
  ylab("Proportion methylated")+
  theme_bw()+
  xlab("Total length (cm)")+
  theme(strip.background = element_blank(), 
        strip.text = element_text(face="italic", size=10), 
        plot.title = element_text(hjust = 0.5, size=11))+
  ggtitle("B. Simplified hypothesis")+
  guides(fill="none")

# fdep

## regression plot

simple_lm_both<-ggplot(data, aes(x=Total_length, y=prop_meth))+
  geom_point(aes(fill=sex_fac), pch=21, alpha=0.4, size=1.5)+
  geom_smooth(method="glm",
              colour = "black", 
              se=FALSE, size=0.5)+
  stat_regline_equation(size = 2, 
                        vjust = 0.75, hjust = 0.4, label.x = 60, label.y.npc = "bottom",
                        aes(label = paste0(..rr.label..)))+
  geom_smooth(method="glm",
              aes(colour=sex_fac),
              lty=2,
              se=FALSE, size=0.5)+
  facet_wrap(~gene_fac, scales="free", ncol=1, strip.position="right")+
  ylab("Proportion methylated")+
  xlab("Total length (cm)")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face="italic", size=10),
        plot.title = element_text(hjust = 0.5, size=11))+
  scale_fill_manual(values=c("#8da0cb", "#fc8d62"))+
  scale_colour_manual(values=c("#8da0cb", "#fc8d62"))+
  ggtitle("C. Preliminary result")+
  guides(fill="none", colour="none")

simple_lm_both

legend_plot <- ggplot(data, aes(x=Total_length, y=prop_meth))+
  geom_point(aes(fill=sex_fac), pch=21, size=1.5)+
  labs(fill="Sex")+
  scale_fill_manual(values=c("#8da0cb", "#fc8d62"))+
  theme_bw()+
  theme(legend.position="bottom", 
        legend.margin=margin(0,0,0,0))

legend<-get_legend(legend_plot)

pohpr<-ggarrange(bps3+theme(strip.text = element_blank(), 
                     plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")), 
          fdep+theme(axis.title.y = element_blank(), 
                     strip.text = element_blank(),
                     plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")), 
          simple_lm_both+theme(axis.title.y = element_blank(),
                          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")),
          nrow=1, ncol=3, 
          widths = c(1.12, 1, 1.12))

pohpr_leg <- plot_grid(pohpr, legend, nrow = 2, rel_heights = c(14.5, 0.5))
pohpr_leg

## save plot
save_plot("figure_hypothesis_more.pdf", 
          pohpr_leg,
          base_height = 29/3,
          base_width = 21/2.5)

## AIC/BIC/VIF ----


AIC_BIC<-NULL
for (Amplicon in levels(data$Gene)) {
  data_gene <- filter(data, Gene==Amplicon) 
  betareg_mod_temp_age <- (betareg(prop_meth ~ 
                                     (CpG_site+JulieRegionEC+SexM+AgeClass)^3, 
                                   link="logit", data=data_gene))
  betareg_mod_temp_size <- (betareg(prop_meth ~ 
                                      (CpG_site+JulieRegionEC+SexM+Total_length)^3, 
                                    link="logit", data=data_gene))
  betareg_mod_temp_size_age <- (betareg(prop_meth ~ 
                                          (CpG_site+JulieRegionEC+SexM+Total_length+AgeClass)^4,
                                        link="logit", data=data_gene))
  Cand.models <- list(betareg_mod_temp_age, betareg_mod_temp_size, betareg_mod_temp_size_age)
  Modnames <- c("Age", "Total length", "Total length and age")
  AIC_temp<-aictab(cand.set = Cand.models, modnames = Modnames,
                   second.ord = TRUE, nobs = NULL, 
                   sort = FALSE) 
  BIC_temp<-bictab(cand.set = Cand.models, modnames = Modnames,
                   nobs = NULL, 
                   sort = FALSE)
  AIC_BIC_temp<-merge.data.frame(AIC_temp, BIC_temp, 
                                 by=c("Modnames", "K", "LL"), 
                                 suffixes=c(" AIC"," BIC"))
  AIC_BIC_temp<-cbind(Amplicon, AIC_BIC_temp)
  AIC_BIC<-rbind(AIC_BIC_temp, AIC_BIC)
}

names(AIC_BIC)[names(AIC_BIC)=="Modnames"]<-"Covariate(s)"
names(AIC_BIC)[names(AIC_BIC)=="LL"]<-"Log-likelihood"
names(AIC_BIC)[names(AIC_BIC)=="LL"]<-"Log-likelihood"
names(AIC_BIC)[names(AIC_BIC)=="Delta_AICc"]<-"Delta AICc"
names(AIC_BIC)[names(AIC_BIC)=="Delta_BIC"]<-"Delta BIC"
names(AIC_BIC)[names(AIC_BIC)=="ModelLik AIC"]<-"Relative likelihood (AICc)"
names(AIC_BIC)[names(AIC_BIC)=="ModelLik BIC"]<-"Relative likelihood (BIC)"
names(AIC_BIC)[names(AIC_BIC)=="AICcWt"]<-"AIC model weight"
names(AIC_BIC)[names(AIC_BIC)=="BICWt"]<-"BIC model weight"

write.csv(AIC_BIC, "AICBIC_weights.csv")

## AIC BIC TEST - MALES/FEMALEs ----

AIC_BIC<-NULL 

for (Sex in levels(data$Sex)) {
  data_sex <- filter(data, Sex==Sex)
  data_sex$Gene <- factor(data_sex$Gene)
  
  for (Amplicon in levels(data_sex$Gene)) {
    data_gene <- filter(data_sex, Gene==Amplicon)
    betareg_mod_temp_age <- (betareg(prop_meth ~ 
                                       (CpG_site+JulieRegionEC+AgeClass)^3, 
                                     link="logit", data=data_gene))
    
    betareg_mod_temp_size <- (betareg(prop_meth ~ 
                                        (CpG_site+JulieRegionEC+Total_length)^3, 
                                      link="logit", data=data_gene))
    
    betareg_mod_temp_size_age <- (betareg(prop_meth ~ 
                                            (CpG_site+JulieRegionEC+SexM+Total_length+AgeClass)^4,
                                          link="logit", data=data_gene))
    Cand.models <- list(betareg_mod_temp_age, betareg_mod_temp_size, betareg_mod_temp_size_age)
    Modnames <- c("Age", "Total length", "Total length and age")
    AIC_temp<-aictab(cand.set = Cand.models, modnames = Modnames,
                 second.ord = TRUE, nobs = NULL, 
                 sort = FALSE)
    BIC_temp<-bictab(cand.set = Cand.models, modnames = Modnames,
                 nobs = NULL, 
                 sort = FALSE)
    AIC_BIC_temp<-merge.data.frame(AIC_temp, BIC_temp, 
                                   by=c("Modnames", "K", "LL"), 
                                   suffixes=c(" AIC"," BIC"))
    AIC_BIC_temp<-cbind(Sex, Amplicon, AIC_BIC_temp)
    AIC_BIC<-rbind(AIC_BIC_temp, AIC_BIC)
  }
}
names(AIC_BIC)[names(AIC_BIC)=="Modnames"]<-"Covariate(s)"
names(AIC_BIC)[names(AIC_BIC)=="LL"]<-"Log-likelihood"
names(AIC_BIC)[names(AIC_BIC)=="LL"]<-"Log-likelihood"
names(AIC_BIC)[names(AIC_BIC)=="Delta_AICc"]<-"Delta AICc"
names(AIC_BIC)[names(AIC_BIC)=="Delta_BIC"]<-"Delta BIC"
names(AIC_BIC)[names(AIC_BIC)=="ModelLik AIC"]<-"Relative likelihood (AICc)"
names(AIC_BIC)[names(AIC_BIC)=="ModelLik BIC"]<-"Relative likelihood (BIC)"
names(AIC_BIC)[names(AIC_BIC)=="AICcWt"]<-"AIC model weight"
names(AIC_BIC)[names(AIC_BIC)=="BICWt"]<-"BIC model weight"

write.csv(AIC_BIC, "AICBIC_weights_mf.csv")

## age/size correlation ----

data_distinct <- data %>% distinct(Fish, .keep_all = TRUE)

# create function 
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    theme_minimal()+
    stat_smooth(method = "lm", col = "red") +
    labs(caption = paste("Adjusted R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P-value =",signif(summary(fit)$coef[2,4], 5)))
}

model <- lm(AgeClass ~ Total_length + JulieRegion, data = data_distinct)
g<-ggplotRegression(model)
g
g+xlab("Total length (cm)")+ylab("Age (years)")

## calculate VIF values for model (no interaction terms) ----

#VIFs start at 1 and have no upper limit. 
#A value of 1 indicates that there is no correlation between this independent variable and any others. 
#VIFs between 1 and 5 suggest that there is a moderate correlation, but it is not severe enough to warrant corrective measures. 
#VIFs greater than 5 represent critical levels of multicollinearity where the coefficients are poorly estimated, and the p-values are questionable.

# If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, 
# then that predictor is more related to the other predictors than it is to the response.
VIF_table<-NULL
for (Amplicon in levels(data$Gene)) {
  data_gene <- filter(data, Gene==Amplicon)
  betareg_mod_temp_age <- (betareg(prop_meth ~ 
                                     (CpG_site+JulieRegionEC+SexM+AgeClass), 
                                   link="logit", data=data_gene)) 
  betareg_mod_temp_size <- (betareg(prop_meth ~ 
                                      (CpG_site+JulieRegionEC+SexM+Total_length), 
                                    link="logit", data=data_gene)) 
  betareg_mod_temp_size_age <- (betareg(prop_meth ~ 
                                          (CpG_site+JulieRegionEC+SexM+Total_length+AgeClass),
                                        link="logit", data=data_gene))
  VIF_temp<-as.data.frame(vif(betareg_mod_temp_size_age))
  VIF_temp<-setDT(VIF_temp, keep.rownames = TRUE)[]
  VIF_temp<-cbind(Amplicon, VIF_temp) 
  VIF_table<-rbind(VIF_temp, VIF_table)
}

VIF_table$rn<-gsub("CpG_site", "CpG site", VIF_table$rn)
VIF_table$rn<-gsub("JulieRegionEC", "Region", VIF_table$rn)
VIF_table$rn<-gsub("SexM", "Sex", VIF_table$rn)
VIF_table$rn<-gsub("Total_length", "Total length", VIF_table$rn)
VIF_table$rn<-gsub("AgeClass", "Age", VIF_table$rn)

names(VIF_table)[names(VIF_table)=="rn"]<-"Covariate"
names(VIF_table)[names(VIF_table)=="GVIF"]<-"Generalised VIF (GVIF)"
names(VIF_table)[names(VIF_table)=="Df"]<-"Degrees of freedom (Df)"

write.csv(VIF_table, "VIF_table.csv", row.names=FALSE)


