library(tidyverse)
library(viridis)
library(MASS)
setwd("C:/Users/fc809/Downloads/")
getwd()
# Load data
all_patients_genotype_collapsed_annotated_wide <- read.csv("all_patients_genotype_collapsed_annotated_wide_v2.csv")

# Keep only rows with tel_length > 0
a <- all_patients_genotype_collapsed_annotated_wide %>% filter(tel_length > 0)

# Replace NAs with FALSE
a[is.na(a)] <- F

# Fit negative binomial regression including 'multihit'
model <- glm.nb(tel_length ~ ASXL1 + TET2 + DNMT3A + TP53 + JAK2 + Multihit + patient, data = a)
parameters<-
summary(glm.nb(tel_length~ASXL1+TET2+DNMT3A+TP53+JAK2+Multihit+patient,data=a))$coefficient[2:7,1:2]
Genes<-substr(rownames(parameters),1,str_length(rownames(parameters))-4)
parameters<-as.tibble(parameters)
parameters$Gene<-Genes
parameters$Estimate_l<-parameters$Estimate-1.96*parameters$`Std. Error`
parameters$Estimate_h<-parameters$Estimate+1.96*parameters$`Std. Error`

ggplot(parameters,aes(x=reorder(Gene,Estimate),y=Estimate,ymin=Estimate_l,ymax=Estimate_h))+geom_point()+geom_errorbar(width=0)+theme_classic()+geom_hline(yintercept=0,linetype="dashed")+labs(x="",y="effect on telomere length")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

