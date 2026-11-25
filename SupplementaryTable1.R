######################################################################
# Supplementary Table S1 -- survival analyses of classifier-defined
# colorectal cancer states in TCGA-COAD (updated classifiers)
#
# In vivo whole-tumor classifier only. For each endpoint (OS, PFS):
#   * global log-rank test across the three classes ("None" excluded)
#   * Cox models fit on all three classes (reference releveled per
#     comparison; zero follow-up times retained), unadjusted and adjusted
#     (age, sex, race, AJCC stage, and combinations)
# Output: one tidy table (CSV + console) to paste back for formatting.
######################################################################

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(readr); library(purrr)
})

base_dir <- "/Users/khalidishani/Desktop/Summer Research - Dr. Hwang/Colon Cancer Project/Genome Medicine Revision"
out_csv  <- file.path(base_dir, "Supplementary Tables/SuppTableS1_survival_results.csv")
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

######################################################################
# CONFIG 1 -- input files (merged predictions + clinical data)
######################################################################
inputs <- list(
  "In vivo whole-tumor" = list(
    file        = file.path(base_dir, "DataSets/Figure 8/Figure 8E-H/survivalData_TUMOR.csv"),
    class_col   = "Predicted_Class_mouseLoadings",
    labels      = c(Polyps = "Polyp-like", AKP_Tu = "AKP_Tu-like", AKPS_Tu = "AKPS_Tu-like"),
    comparisons = list(c("AKP_Tu", "Polyps"), c("AKPS_Tu", "Polyps"), c("AKP_Tu", "AKPS_Tu"))
  )
)
# Each comparison is c(test, reference): HR = test vs reference.

######################################################################
# CONFIG 2 -- clinical column names (cBioPortal TCGA defaults; run
# names(read_csv(<file>)) once and edit if yours differ)
######################################################################
cols <- list(
  os_time  = "OS_MONTHS",   os_status  = "OS_STATUS",
  pfs_time = "PFS_MONTHS",  pfs_status = "PFS_STATUS",
  age      = "AGE",         sex        = "SEX",
  race     = "RACE",        stage      = "AJCC_PATHOLOGIC_TUMOR_STAGE"
)

adjustments <- list(
  "Unadjusted"                    = character(0),
  "Age"                           = "age",
  "Sex"                           = "sex",
  "Race"                          = "race",
  "Age + sex + race"              = c("age", "sex", "race"),
  "AJCC stage"                    = "stage",
  "AJCC stage + age + sex"        = c("stage", "age", "sex"),
  "AJCC stage + age + sex + race" = c("stage", "age", "sex", "race")
)

######################################################################
# Helpers
######################################################################
# "1:DECEASED" / "0:LIVING" / numeric -> 0/1
parse_status <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  as.integer(sub(":.*$", "", as.character(x)))
}

na_tokens <- c("", "NA", "[Not Available]", "[Not Evaluated]", "[Unknown]",
               "[Discrepancy]", "Not Reported", "Unknown", "[Not Applicable]")

clean_stage <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- sub("^STAGE\\s+", "", s)
  s <- sub("[ABC]$", "", s)
  factor(ifelse(s %in% c("I", "II", "III", "IV"), s, NA), levels = c("I", "II", "III", "IV"))
}

clean_cat <- function(x) {
  x <- as.character(x); x[x %in% na_tokens] <- NA
  factor(x)
}

prep_data <- function(cfg) {
  raw <- read_csv(cfg$file, show_col_types = FALSE)
  miss <- setdiff(c(cfg$class_col, unlist(cols)), names(raw))
  if (length(miss)) stop(cfg$file, " lacks column(s): ", paste(miss, collapse = ", "),
                         "\nAvailable: ", paste(names(raw), collapse = ", "))
  d <- tibble(
    class      = as.character(raw[[cfg$class_col]]),
    os_time    = suppressWarnings(as.numeric(raw[[cols$os_time]])),
    os_event   = parse_status(raw[[cols$os_status]]),
    pfs_time   = suppressWarnings(as.numeric(raw[[cols$pfs_time]])),
    pfs_event  = parse_status(raw[[cols$pfs_status]]),
    age        = suppressWarnings(as.numeric(raw[[cols$age]])),
    sex        = clean_cat(raw[[cols$sex]]),
    race       = clean_cat(raw[[cols$race]]),
    stage      = clean_stage(raw[[cols$stage]])
  ) %>% filter(class %in% names(cfg$labels))
  # collapse rare race categories (<10) into "Other" to keep models estimable
  rc <- table(d$race); rare <- names(rc)[rc < 10]
  if (length(rare)) d$race <- factor(ifelse(d$race %in% rare, "Other", as.character(d$race)))
  d$class <- factor(d$class, levels = names(cfg$labels))
  d
}

# Matches the user's process: all classified tumors with non-missing time and
# status are kept (zero follow-up times included), and the global test is the
# log-rank / Cox score test across all three classes.
global_logrank <- function(d, time, event) {
  dd <- d[!is.na(d[[time]]) & !is.na(d[[event]]), ]
  sd <- survdiff(Surv(dd[[time]], dd[[event]]) ~ class, data = dd)
  tibble(n = nrow(dd), events = sum(dd[[event]]),
         HR = NA_real_, lower = NA_real_, upper = NA_real_,
         p = pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE), ph_p = NA_real_)
}

# One Cox model with all three classes; `ref` is set as the reference level and
# the HR for `test` vs `ref` is extracted (AKP vs AKPS uses AKPS as reference).
pairwise_cox <- function(d, time, event, test, ref, covars) {
  dd <- d
  keep <- !is.na(dd[[time]]) & !is.na(dd[[event]])
  for (v in covars) keep <- keep & !is.na(dd[[v]])
  dd <- dd[keep, ]
  dd$grp <- relevel(droplevels(factor(as.character(dd$class))), ref = ref)
  dd <- droplevels(dd)
  empty <- tibble(n = nrow(dd), events = sum(dd[[event]]), HR = NA_real_,
                  lower = NA_real_, upper = NA_real_, p = NA_real_, ph_p = NA_real_)
  if (nlevels(dd$grp) < 2 || sum(dd[[event]]) < 5) return(empty)
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ grp",
                         if (length(covars)) paste0(" + ", paste(covars, collapse = " + ")) else ""))
  fit <- tryCatch(coxph(f, data = dd), error = function(e) NULL)
  if (is.null(fit)) return(empty)
  s  <- summary(fit)
  rn <- paste0("grp", test)
  if (!rn %in% rownames(s$coefficients)) return(empty)
  ph <- tryCatch(cox.zph(fit)$table["grp", "p"], error = function(e) NA_real_)
  tibble(n = s$n, events = s$nevent,
         HR = s$conf.int[rn, "exp(coef)"],
         lower = s$conf.int[rn, "lower .95"], upper = s$conf.int[rn, "upper .95"],
         p = s$coefficients[rn, "Pr(>|z|)"], ph_p = ph)
}

######################################################################
# Run
######################################################################
endpoints <- list(OS = c("os_time", "os_event"), PFS = c("pfs_time", "pfs_event"))

results <- imap_dfr(inputs, function(cfg, classifier) {
  d <- prep_data(cfg)
  message(classifier, ": ", nrow(d), " classified tumors; class counts: ",
          paste(names(table(d$class)), table(d$class), sep = "=", collapse = ", "))
  imap_dfr(endpoints, function(ep, ep_name) {
    glob <- global_logrank(d, ep[1], ep[2]) %>%
      mutate(Comparison = "Global comparison", Model = "Log-rank")
    pairs <- map_dfr(cfg$comparisons, function(cmp) {
      imap_dfr(adjustments, function(covars, adj_name) {
        pairwise_cox(d, ep[1], ep[2], cmp[1], cmp[2], covars) %>%
          mutate(Comparison = paste(cfg$labels[cmp[1]], "vs", cfg$labels[cmp[2]]),
                 Model = adj_name)
      })
    })
    bind_rows(glob, pairs) %>% mutate(Classifier = classifier, Endpoint = ep_name)
  })
})

results <- results %>%
  mutate(`95% CI` = ifelse(is.na(lower), "—", sprintf("%.3f–%.3f", lower, upper)),
         HR_fmt   = ifelse(is.na(HR), "—", sprintf("%.3f", HR)),
         P_fmt    = ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                                                  sprintf("%.3f", p)))) %>%
  select(Classifier, Endpoint, Comparison, Model, n, events,
         HR = HR_fmt, `95% CI`, `P value` = P_fmt, PH_test_p = ph_p, HR_raw = HR, p_raw = p)

write_csv(results, out_csv)
message("Saved: ", out_csv)

options(width = 200)
print(as.data.frame(results %>% select(-HR_raw, -p_raw)), row.names = FALSE)