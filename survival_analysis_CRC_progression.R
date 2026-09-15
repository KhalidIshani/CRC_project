### Survival Modeling for "Discrete Transcriptional States Define Immune Remodeling and 
### Dynamic CMS Transitions in Colorectal Cancer"
###
### Authors: Katie Storey, Zejie Yan
### Date: 2026-09-14

library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

wd = "~/Documents/Research/CRC/"
setwd(wd)

# load patient data for survival analysis
data <- read_excel(paste0(wd,"survivalData_TUMOR.xlsx"))

### clean up data
data <- data %>%
  mutate(
    os_status_num  = as.numeric(sub("^([01]).*", "\\1", OS_STATUS)),
    pfs_status_num = as.numeric(sub("^([01]).*", "\\1", PFS_STATUS))
  )

data <- data %>%
  mutate(
    OS_MONTHS  = as.numeric(OS_MONTHS),
    PFS_MONTHS = as.numeric(PFS_MONTHS)
  )

data <- data %>%
  mutate(
    diff_os_pfs = OS_MONTHS - PFS_MONTHS,
    group_var   = factor(Predicted_Class)  
  )

## functions to recode T,N,M,AJCC stages

recode_T <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("TIS")               ~ 0,
    x_chr %in% c("T1")                ~ 1,
    x_chr %in% c("T2")                ~ 2,
    x_chr %in% c("T3")                ~ 3,
    x_chr %in% c("T4", "T4A", "T4B")  ~ 4,
    TRUE                              ~ NA_real_
  )
}

recode_N <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("N0")                          ~ 0,
    x_chr %in% c("N1", "N1A", "N1B", "N1C")     ~ 1,
    x_chr %in% c("N2", "N2A", "N2B")            ~ 2,
    TRUE                                        ~ NA_real_
  )
}

recode_M <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("M0")                      ~ 0,
    x_chr %in% c("M1", "M1A", "M1B")        ~ 1,
    x_chr %in% c("MX")                      ~ NA_real_,
    TRUE                                    ~ NA_real_
  )
}

recode_ajcc_stage <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("STAGE I", "STAGE IA") ~ 1,
    x_chr %in% c("STAGE II", "STAGE IIA", "STAGE IIB", "STAGE IIC") ~ 2,
    x_chr %in% c("STAGE III", "STAGE IIIA", "STAGE IIIB", "STAGE IIIC") ~ 3,
    x_chr %in% c("STAGE IV", "STAGE IVA", "STAGE IVB") ~ 4,
    TRUE ~ NA_real_
  )
}

data <- data %>%
  mutate(
    AJCC_STAGE_NUM = recode_ajcc_stage(AJCC_PATHOLOGIC_TUMOR_STAGE),
    PATH_T_NUM     = recode_T(PATH_T_STAGE),
    PATH_N_NUM     = recode_N(PATH_N_STAGE),
    PATH_M_NUM     = recode_M(PATH_M_STAGE)
  )


# remove the NAs from data
data_clean <- data[data$group_var != "None", ]

# drop the unused level so it doesn't show up in the plot
data_clean$group_var <- droplevels(as.factor(data_clean$group_var))


### K-M curve for overall survival
fit_os <- survfit(
  Surv(OS_MONTHS, os_status_num) ~ group_var,
  data = data_clean
)

plot_os = ggsurvplot(
  fit_os,
  data         = data_clean,
  #pval         = TRUE,        
  risk.table   = TRUE,
  legend.title = "Group",
  xlab = "Overall survival (months)",
  ylab = "Survival probability"
)

plot_os

# arrange ggsurvplot to save
arranged_plot <- arrange_ggsurvplots(list(plot_os), print = FALSE)

# save plot
ggsave(filename = paste0(wd, "immune_exposed_tumor_classifier_OS.png"),
       plot = arranged_plot,
       width = 15, # potentially adjust width/height to desired size
       height = 7,
       units = "in",
       dpi = 300) 


### K-M curve for progression-free survival
fit_pfs <- survfit(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var,
  data = data_clean
)

plot_pfs = ggsurvplot(
  fit_pfs,
  data         = data_clean,
  #pval         = TRUE,
  risk.table   = TRUE,
  legend.title = "Group",
  xlab = "Progression-free survival (months)",
  ylab = "Progression-free survival probability"
)

plot_pfs

# arrange ggsurvplot to save
arranged_plot <- arrange_ggsurvplots(list(plot_pfs), print = FALSE)

ggsave(filename = paste0(wd, "immune_exposed_tumor_classifier_PFS.png"),
       plot = arranged_plot,
       width = 15, # potentially adjust width/height to desired size
       height = 7,
       units = "in",
       dpi = 300) 


### Cox Proportional Hazards model

## make "Polyps" the reference group
data_clean$group_var <- relevel(factor(data_clean$group_var), ref = "Polyps")
## use below to make "AKPS-Tu" the reference group instead:
#data_clean$group_var <- relevel(factor(data_clean$group_var), ref = "AKPS_Tu")

## unadjusted
# OS
cox_os <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var,
  data = data_clean
)
summary(cox_os)

#PFS
cox_pfs <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var,
  data = data_clean
)
summary(cox_pfs)


## adjusting for AJCC stage

## set "Polyps" is the reference group:
data_clean$group_var <- relevel(factor(data_clean$group_var), ref = "Polyps")
## or use this to compare to AKPS:
#data_clean$group_var <- relevel(factor(data_clean$group_var), ref = "AKPS_Tu")

# OS
cox_os_ajcc <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var +
    AJCC_STAGE_NUM,
  data = data_clean
)
summary(cox_os_ajcc)

# PFS
cox_pfs_ajcc <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var +
    AJCC_STAGE_NUM,
  data = data_clean
)
summary(cox_pfs_ajcc)


## adjusting for clinical factors 

# clean up clinical factors and check for number of NAs
table(data$AGE, useNA = "ifany")
data_clean$AGE = as.numeric(data_clean$AGE)

data_clean$SEX <- factor(data_clean$SEX)
table(data_clean$SEX)
# make "Female" the reference group
data_clean$SEX <- relevel(data_clean$SEX, ref = "Female")

data_clean$RACE <- factor(data_clean$RACE)
table(data_clean$RACE)
# make "White" the reference group
data_clean$RACE <- relevel(data_clean$RACE, ref = "White")

## all clinical factors + AJCC stage
# OS
cox_os_ajcc_clin <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var +
    AJCC_STAGE_NUM + AGE + SEX + RACE,
  data = data_clean
)
summary(cox_os_ajcc_clin)

# PFS
cox_pfs_ajcc_clin <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var +
    AJCC_STAGE_NUM + AGE + SEX + RACE,
  data = data_clean
)
summary(cox_pfs_ajcc_clin)


## all clinical factors without AJCC
## OS
cox_os_clin <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var +
    AGE + SEX + RACE,
  data = data_clean
)
summary(cox_os_clin)


#PFS
cox_pfs_clin <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var +
    AGE + SEX + RACE,
  data = data_clean
)
summary(cox_pfs_clin)


## individual clinical faxtors
## AGE
# OS
cox_os_clin <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var +
    AGE,
  data = data_clean
)
summary(cox_os_clin)

# PFS
cox_pfs_clin <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var +
    AGE,
  data = data_clean
)
summary(cox_pfs_clin)

## SEX
# OS
cox_os_clin <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var +
    SEX,
  data = data_clean
)
summary(cox_os_clin)

# PFS
cox_pfs_clin <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var +
    SEX,
  data = data_clean
)
summary(cox_pfs_clin)

## RACE
# OS
cox_os_clin <- coxph(
  Surv(OS_MONTHS, os_status_num) ~ group_var +
    RACE,
  data = data_clean
)
summary(cox_os_clin)

#PFS
cox_pfs_clin <- coxph(
  Surv(PFS_MONTHS, pfs_status_num) ~ group_var +
    RACE,
  data = data_clean
)
summary(cox_pfs_clin)



