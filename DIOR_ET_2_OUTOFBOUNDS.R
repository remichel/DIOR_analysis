# OUT OF BOUND ANALSIS --------------------------------------------------------
#
# Very rarely, the online fixation control misses to reject a trial in which
# fixation was not kept within the allowed boundaries. Here, we identify these
# trials and save them in a separate list to load them in our behavioral analysis
# to exclude them post-hoc.

# LOAD PACKAGES ----------------------------------------------------------------
library(rmTools)
libraries("dplyr", "tidyr")


# CONFIG -----------------------------------------------------------------------
STUDYpath   <- "insert_path_here"
RDATpath    <- "/Out/DIO_"
maxDeg      <- c(-2,2)
fixWait     <- .3 # in s
subjectlist <- seq(1,44,1)
outofbounds <- sapply(as.character(subjectlist),function(x) NULL)
maxSOA      <- 1.042 # in s

# LOAD ET DATA AND IDENTIFY TRIALS ---------------------------------------------
for(i in subjectlist){

  load(paste0(STUDYpath, i, RDATpath, i, ".Rdat"))

  # prepare dataset for analysis
  clean <- et %>%
    # filter to test trials, exclude practice
    filter(section == 1) %>%
    # filter to non aborted trials
    filter(aborted == FALSE) %>%
    # filter to fix control period
    mutate(fixtimeoff = case_when(solotrigger == ">TTargetOff"~ trialtime, # target trials
                                  solotrigger == ">TCueOff"~ trialtime,  # catch trials
                            TRUE~ 0  )) %>%
    mutate(fixtimeoff = case_when(soa == "0" ~  fixtimeoff + maxSOA,
                                  TRUE ~ fixtimeoff  )) %>%
    group_by(trial) %>%
    filter(trialtime > min(trialtime) + fixWait) %>%
    filter(trialtime < max(fixtimeoff)) %>%
    ungroup()


  # if any, which trials were out of bounds?
  if(any((clean$out_of_bounds == 1))){
    outofbounds[[i]] <-  unique(clean$trial[which(clean$out_of_bounds == 1)])
  }
}

# SAVE IDENTIFIED TRIALS -------------------------------------------------------

# convert out of bounds list into df
outofbounds <- outofbounds %>%
  tibble::enframe() %>%
  tidyr::unnest(value) %>%
  dplyr::rename(subject = name, trial = value)

# save to RDAT
save(file = paste0(STUDYpath, "outofbounds.Rdat"), outofbounds)
