# ET PREPROCESSING -------------------------------------------------------------
#
# This script preprocesses raw .edf files and prepares the ET data to analyze
# out-of-boundary fixations which were not detected by the online fixation control.

# INSTALL PACKAGES ------------------------------------------------------------
# If you run this script for the first time, you might need to install packages

# install.packages(c('devtools','doParallel','data.table','ggplot2','lme4','dplyr','afex'))
# devtools::install_github('hedderik/pR')
# devtools::install_github('wanjam/wmR')
# devtools::install_github('hedderik/hvR')
# devtools::install_github('jashubbard/edfR')
# devtools::install_github("jwdink/eyetrackingR")

# LOAD PACKAGES ----------------------------------------------------------------

library(rmTools)
libraries(c('data.table', 'dplyr', 'hvR', 'pR', 'doParallel', 'edfR', 'stringr', "tidyr"))

# CONFIG -----------------------------------------------------------------------

# folders
STUDYpath   <- "insert_path_here"
RDATpath    <- "/Temp/"
OUTpath     <- "/Out/"
EDFpattern  <- "DIO_\\d.edf" 
# pattern used to match your filenames (see https://github.com/rstudio/cheatsheets/raw/master/strings.pdf)

# config
subjectlist <- seq(1,44,1)
max_cores   <- 35
sfreq       <- 500
maxdev      <- 2 # maximum allowed distance from screen scenter in degrees
x           <- 1920/2 # screen center x
y           <- 1080/2 # screen center y
pixperdeg   <- 55.1734

# PREPROCESSING ----------------------------------------------------------------

# prepare parallel processing
n_cores   <- min(length(subjectlist), max_cores) # start a core for each subject, or a max of max_cores
cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)

# create list of all edf files
edfs <- paste0(STUDYpath, subjectlist, '/DIO_', subjectlist, '.edf')

# LOOP OVER SUBJECTS -----------------------------------------------------------
foreach(i = 1:length(subjectlist), 
        .packages = c('data.table','edfR', 'dplyr', "tidyr")) %dopar% {

  # get current edf file name
  fname  <- edfs[i]
  # get current subject number
  subject <- subjectlist[i]

  # create output folders
  if (!dir.exists(paste0(STUDYpath, subject, RDATpath))) {dir.create(paste0(STUDYpath, subject, RDATpath))}
  if (!dir.exists(paste0(STUDYpath, subject, OUTpath))) {dir.create(paste0(STUDYpath, subject, OUTpath))}

  # edf.all parses 'all' possible data.types in the file. You can also just
  # parse parts of it. Usually that's not gonna be necessary, though.
  curdat               <- edf.all(fname, samples = TRUE, eventmask = TRUE)

  # The file header contains some meta-info (tracked eye, Hz, etc). get that:
  curdat$recordingInfo <- edf.recordings(fname)
  curdat$recordingInfo <- curdat$recordingInfo[1,]  # we only need the info once

  message("Save raw dataset as .Rdat ...")
  # Save as Rdat
  save(file = paste0(STUDYpath, subject, RDATpath, "DIO_", subject, ".Rdat"), curdat)

  message("Extracting relevant variables, incl. triggers...")
  # Extract relevant variables
  et <- as.data.frame(curdat$samples)
  # extract triggers and get rid of calibration messages
  msg <- as.data.frame(curdat$messages)
  msg$msgstart <- lapply(strsplit(msg$msg, ''), `[[`, 1)
  msg <- msg %>%
    filter(msgstart == ">") %>%
    select(-c(msgstart)) %>%
    separate(col = msg,
             into = c("trigger", paste("triggerinfo", seq(1,11, 1), sep = "_")),
             sep = "\\\\t",
             fill = "right" 
             # only TData has 11 rows of information, all other triggers have less. 
             # if trigger has less then 11 infos, e.g. only 6, fill rows 7-11 with NA
             ) %>%
    # extract session info
    mutate(section = as.numeric(case_when(trigger == ">TStart" ~ triggerinfo_1,
                                        TRUE ~ "999"))) %>%
    mutate(section = na_if(section, "999")) %>%
    # extract trial info
    mutate(trial = as.numeric(case_when(trigger == ">TStart" ~ triggerinfo_3,
                             TRUE ~ "0"))) %>%
    mutate(trial = na_if(trial, "0")) %>%
    # extract SOA info (position 7 in ">TData")
    mutate(soa = triggerinfo_7) %>%
    # extract RT info from ">TResp" trigger (position 5)
    mutate(rt = case_when(trigger == ">TResp" ~ as.numeric(triggerinfo_5),
                          TRUE ~ 0)) %>%
    mutate(rt = na_if(rt, 0)) %>%
    # if aborted, soa is stored in triggerinfo 6
    mutate(soa = case_when(trigger == ">TAbort" ~ triggerinfo_6,
                           TRUE ~ soa))

  # merge data frames by time
  et <- merge(et,
             msg,
             by = c("time"),
             all = T) 
  # for binocular, we need all.x and all.y = T, otherwise messages might be skipped.


  et <- et %>%
    # create a trigger variable, which will have an entry every row, 
    # and a solo trigger variable which codes the time when the trigger was sent
    mutate(solotrigger = trigger) %>%
    # fill all relevant variables so that they have entries in every row, 
    # trigger and session are  coded in TStart trigger, therefore filling downwards
    fill(trigger, .direction = "down") %>%
    fill(section, .direction = "down") %>%
    fill(trial, .direction = "down") %>%
    # soa and rt are coded in the TData trigger which is sent in the end, 
    # so we need to fill upwards here
    fill(soa, .direction = "up") %>%
    fill(rt, .direction = "up") %>%
    # exclude all triggers, and therefore datapoints, that arent necessary 
    # for our analyyes
    filter(trigger %in% c(">TStart", ">TCueOn", ">TCueOff", ">TTargetOn", ">TTargetOff", ">TResp", ">TAbort")) %>%
    # as some trials were aborted and repeated, trial numbers can occur more 
    # than once, thus we assign a unique trial id for each sent TStart
    group_by(solotrigger) %>%
    # for each group of solo triggers, incl TStart, we use row numbers as trial_ids
    mutate(trial_id = row_number()) %>%
    ungroup() %>%
    # if solotrigger is trial start, keep the newly assigned trial_id, if not, 
    # it is meaningless and therefore is set to 0 (and then to NA 
    # in a follwing step)
    mutate(trial_id = case_when(solotrigger == ">TStart"~ trial_id,
                                TRUE~ as.integer(0)   )) %>%
    # set all other entries to NA so that we can use fill to add the trial id 
    # to every row
    mutate(trial_id = na_if(trial_id, 0)) %>%
    fill(trial_id, .direction = "down") %>%
    # for each trial, code whether a trial was aborted in the new variable 
    # validTrial, code not-aborted trials as NA, so that we can use fill 
    # afterwards
    group_by(trial_id) %>%
    mutate(validTrial = case_when(solotrigger == ">TAbort" ~ FALSE,
                               TRUE ~ NA )) %>%
    fill(validTrial, .direction = "updown") %>%
    ungroup() %>%
    # recode the validTrial variable  into an aborted variable to be more intuitive
    mutate(aborted = case_when(validTrial == FALSE~  TRUE,
                               TRUE~ FALSE))

  
  
  # compute trial time vector (relative to cue onset OR trial onset if no cue was presented)
  et <- et %>%
    # apply for each trial separately
    group_by(trial_id) %>%
    mutate(trialtime = (time - min(time)) / 1000) %>%
    ungroup() %>%
    # set cue time to 0 if no cue was presentd
    mutate(cuetime = case_when(solotrigger == ">TCueOn" ~ trialtime,
                               TRUE ~ 0)) %>%
    group_by(trial_id) %>%
    mutate(trialtime = trialtime - max(cuetime)) %>%
    ungroup()

  
  # compute deviation from screen center in degrees vis angle 
  et <- et %>%
    mutate(xdevLdeg = (gxL - x) / pixperdeg,
           xdevRdeg = (gxR - x) / pixperdeg,
           ydevLdeg = (gyL - y) / pixperdeg,
           ydevRdeg = (gyR - y) / pixperdeg) 
  
  
  # compute fixations' euclidean distance from screen center 
  et <- et %>% 
    rowwise() %>%
    mutate(devL = sqrt(xdevLdeg^2 + ydevLdeg^2),
           devR = sqrt(xdevRdeg^2 + ydevRdeg^2)) %>%
    ungroup()
  
  
  # bring deviations into long format
  et <- et %>%
    pivot_longer(c("devL", "devR"),
                 names_to = "coordinates", 
                 values_to = "deviation") 


  # create an out-of-bound variable
  et <- et %>%
    mutate(out_of_bounds = case_when(abs(deviation) > maxdev ~ 1,
                                     TRUE ~ 0 ))


  message("Save preprocessed data...")
  save(file = paste0(STUDYpath, subject, OUTpath, "DIO_", subject, ".Rdat"), et)



}

parallel::stopCluster(cl = cl)
