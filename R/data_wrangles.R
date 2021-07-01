get_participant_data <- function(subjects, sessions, data_path, folstr) {
  # this function loads each participant's data
  # and concatenates them into a longform dataset
  fn <- 'sub-0%d_ses-0%d_task-learnAttExp_events.tsv'
  get_subject_strings <- function(i) {
    get_session_strings <- function(j) dir(sprintf(paste(data_path, "sub-0%d", "ses-0%d", folstr, sep = "/"), i, j), pattern=sprintf(fn, i, j), full.names = TRUE)
    do.call(cbind, lapply(sessions, get_session_strings))
  }
  files <- do.call(rbind, lapply(subjects, get_subject_strings))
  rownames(files) <- subjects
  colnames(files) <- sessions
  resplog <- function(i, j) read.table(files[as.character(i),as.character(j)], sep = "\t", header = TRUE)
  d <- do.call(rbind, lapply(subjects, function(i) do.call(rbind, lapply(sessions, function (j) resplog(i, j)))))
  
  # get trials to getthe reward condition 
  fn <- 'sub-0%d_ses-0%d_task-learnAttExp_trls.csv'
  files <- do.call(rbind, lapply(subjects, get_subject_strings))      
  rownames(files) <- subjects
  colnames(files) <- sessions
  eventlog <- function(i, j) {
    e = read.table(files[as.character(i),as.character(j)], sep = ",", header = TRUE)
    e$sub = i
    e$sess= j
    e
  }
  e <- do.call(rbind, lapply(subjects, function(i) do.call(rbind, lapply(sessions, function (j) eventlog(i, j)))))
  names(e)[1] = "t"
  e <- e[, -c(2, 4:8 )]  
  d <- inner_join(d, e, by=c("sub", "sess", "t"))
  # ALLOCATE CUES
  # 1 = left arrow, 2 = right arrow, 3 = bidirectional arrow
  d$cert <- NA
  d$cert[ d$loc == 1 & d$cue == 1 ] = ".8"
  d$cert[ d$loc == 2 & d$cue == 1 ] = ".2"
  d$cert[ d$loc == 2 & d$cue == 2 ] = ".8"
  d$cert[ d$loc == 1 & d$cue == 2 ] = ".2"  
  d$cert[ d$cue == 3 ] = ".5"
  d$cert <- as.factor(d$cert)
  
  d$loc <- as.factor(d$loc)
  levels(d$loc) <- c("left", "right")
  d$sub <- as.factor(d$sub)
  d$rew <- as.factor(d$rew)
  levels(d$rew) <- c("5", "50")
  d$reward_type <- as.factor(d$reward_type)
  levels(d$reward_type) <- c("htgt/hdst", "htgt/ldst", "ltgt/ldst", "ltgt/hdst")
  
  d <- d %>% filter(t > 120) # get rid of initial trials
  
  d
}


get_mri_data <- function(subjects, sessions, data_path, TRs) {
  # this function loads each participant's data
  # and concatenates them into a longform dataset
  
  get_subject_strings <- function(i, k, fn) { # i = sub, j = session, k = TR
    get_session_strings <- function(j, k) {
      t <- try(dir(sprintf(paste(data_path, "sub-0%d", "ses-0%d", "func", sep = "/"), i, j), pattern=sprintf(fn, i, j, k), full.names = TRUE))
      if (!any(length(t))){
        print("no-file")
      } else {
        dir(sprintf(paste(data_path, "sub-0%d", "ses-0%d", "func", sep = "/"), i, j), pattern=sprintf(fn, i, j, k), full.names = TRUE)
      }
    }
    do.call(cbind, lapply(sessions, get_session_strings, k=k))
  }
  
  get_trls_by_TR <- function(k, subjects, sessions, data_path){
    fn <- 'sub-0%d_ses-0%d_task-learnAtt_acq-TR%d_trls.txt'
    files <- do.call(rbind, lapply(subjects, get_subject_strings, k=k, fn=fn))
    rownames(files) <- subjects
    colnames(files) <- sessions
    
    resplog <- function(i, j) {
      if (files[as.character(i),as.character(j)]!="no-file"){
        tmp = read.table(files[as.character(i),as.character(j)], sep = ",", header = TRUE)
        tmp$sub = i
        tmp$sess = j
        tmp
      }
    }
    
    d <- do.call(rbind, lapply(subjects, function(i) do.call(rbind, lapply(sessions, function (j) resplog(i, j)))))
    d$TR = k
    d
  }
  
  d = do.call(rbind, lapply(TRs, get_trls_by_TR, subjects=subjects, sessions=sessions, data_path=data_path))
  
  # get trials to getthe reward condition 
  get_rew_by_TR <- function(k, subjects, sessions, data_path){
    fn <- 'sub-0%d_ses-0%d_task-learnAtt_acq-TR%d_trls.tsv'
    files <- do.call(rbind, lapply(subjects, get_subject_strings, k=k, fn=fn)) 
    rownames(files) <- subjects
    colnames(files) <- sessions
    eventlog <- function(i, j, k) {
      if (files[as.character(i),as.character(j)]!="no-file"){
        e = read.table(files[as.character(i),as.character(j)], sep = "\t", header = TRUE)
        e$sub = i
        e$sess= j
        e$TR = k
        e
      }
    }
    get.event.logs <- function(TR, subjects, sessions){
      k = TR
      do.call(rbind, lapply(subjects, function(i, k) do.call(rbind, lapply(sessions, function (j, k) eventlog(i, j, k), k=TR)), k=TR))
    }
    e = get.event.logs(k, subjects=subjects, sessions=sessions)
    e
  } 
  e = do.call(rbind, lapply(TRs, get_rew_by_TR, subjects=subjects, sessions=sessions, data_path=data_path))
  e <- e %>% select(-c("rew", "loc", "co1", "co2", "or")) 
  names(d)[1] = "t"
  d <- inner_join(d, e, by=c("sub", "sess", "t", "TR"))
  # ALLOCATE CUES
  # 1 = left arrow, 2 = right arrow, 3 = bidirectional arrow
  d$cert <- NA
  d$cert[ d$position == 0 & d$cue == 1 ] = ".8"
  d$cert[ d$position == 1 & d$cue == 1 ] = ".2"
  d$cert[ d$position == 1 & d$cue == 2 ] = ".8"
  d$cert[ d$position == 0 & d$cue == 2 ] = ".2"  
  d$cert[ d$cue == 3 ] = ".5"
  d$cert <- as.factor(d$cert)
  
  d$position <- as.factor(d$position)
  levels(d$position) <- c("left", "right")
  d$sub <- as.factor(d$sub)
  d$reward_trial <- as.factor(d$reward_trial)
  levels(d$reward_trial) <- c("5", "50")
  d$reward_type <- as.factor(d$reward_type)
  levels(d$reward_type) <- c("ltgt/ldst", "ltgt/hdst", "htgt/hdst", "htgt/ldst")
  d
}



clean.sub.data <- function(subject, data, sd_reject=2.5, RT_min=0.1){
  # this function takes the subject index (defined in the variable subject)
  # and all the raw data (data)
  # it filters the data to get the individual subject's data, then trims to
  # get the correct RTs that are > .1 s, and are < 3 * sd from the median
  # for each certainty and reward condition
  sub.data <- data %>% filter(sub == subject) %>%
    filter(rt > RT_min) %>%
    filter(resp == 1) %>%
    group_by(cert, reward_type) %>%
    filter(rt < median(rt) + sd_reject*sd(rt)) 
}

plot.behav <- function(data, dv, iv, grp, ylims, cols){
  # compute mean of data across the variables of interest
  sum.dat <- data %>% group_by(eval(parse(text=grp)), eval(parse(text=iv))) %>%
                  summarise(mu=mean(eval(parse(text=dv))))
  names(sum.dat) <- c(grp, iv, "mu")
  sum.dat$reward_type <- factor(sum.dat$reward_type, c("htgt/ldst", "htgt/hdst", "ltgt/ldst", "ltgt/hdst"))
  sum.dat$cert <- factor(sum.dat$cert, c(".8", ".5", ".2"))
  data$reward_type <- factor(data$reward_type, c("htgt/ldst", "htgt/hdst", "ltgt/ldst", "ltgt/hdst"))
  data$cert <- factor(data$cert, c(".8", ".5", ".2"))
  
  alpha = 1/3
  ggplot(sum.dat, aes_string(x=iv, y="mu", col=grp)) +
    geom_line(aes_string(group=grp), size=1.1) + 
    geom_line(data[data$reward_type == "htgt/ldst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[1]) +
    geom_line(data[data$reward_type == "htgt/hdst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[2]) +
    geom_line(data[data$reward_type == "ltgt/ldst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[3]) +
    geom_line(data[data$reward_type == "ltgt/hdst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[4]) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) + 
    ylab(dv) + xlab(iv) + ylim(ylims) +
    theme(panel.border = element_blank(), 
          panel.grid.major =   element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.margin=margin(6,0,6,0),
          legend.title=element_blank()) 
}

plot.mri <- function(data, dv, iv, grp, facet, ylims, cols){
  # compute mean of data across the variables of interest
  sum.dat <- data %>% group_by(eval(parse(text=grp)), eval(parse(text=iv)), eval(parse(text=facet))) %>%
                      summarise(mu=mean(eval(parse(text=dv))))
  names(sum.dat) <- c(grp, iv, facet, "mu")
  sum.dat$reward_type <- factor(sum.dat$reward_type, c("htgt/ldst", "htgt/hdst", "ltgt/ldst", "ltgt/hdst"))
  sum.dat$cert <- factor(sum.dat$cert, c(".8", ".5", ".2"))
  # now reorder the un-summarised data
  data$reward_type <- factor(data$reward_type, c("htgt/ldst", "htgt/hdst", "ltgt/ldst", "ltgt/hdst"))
  data$cert <- factor(data$cert, c(".8", ".5", ".2"))
  
  # splitting the data into each reward condition for overlaying on plot
  alpha = 1/4
  ggplot(sum.dat, aes_string(x=iv, y="mu", col=grp)) +
    geom_line(aes_string(group=grp), size=1.1) +
    geom_line(data[data$reward_type == "htgt/ldst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[1]) +
    geom_line(data[data$reward_type == "htgt/hdst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[2]) +
    geom_line(data[data$reward_type == "ltgt/ldst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[3]) +
    geom_line(data[data$reward_type == "ltgt/hdst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[4]) +
    facet_wrap(.~eval(parse(text=facet)), nrow=1) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) + 
    ylab(dv) + xlab(iv) + ylim(ylims) +
    theme(panel.border = element_blank(), 
          panel.grid.major =   element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.title=element_blank()) +
    theme_cowplot()
}

plot.mri.simp <- function(data, dv, iv, grp, ylims, cols){
  # compute mean of data across the variables of interest
  sum.dat <- data %>% group_by(eval(parse(text=grp)), eval(parse(text=iv))) %>%
              summarise(mu=mean(eval(parse(text=dv))))
  names(sum.dat) <- c(grp, iv, "mu")
  sum.dat$cert <- factor(sum.dat$cert, c(".8", ".5", ".2"))
  # now reorder the un-summarised data
  data$cert <- factor(data$cert, c(".8", ".5", ".2"))
  
  # splitting the data into each reward condition for overlaying on plot
  alpha = 1/4
  ggplot(sum.dat, aes_string(x=iv, y="mu", col=grp)) +
    geom_line(aes_string(group=grp), size=1.1) +
    geom_line(data[data$reward_type == "htgt/ldst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[1]) +
    geom_line(data[data$reward_type == "htgt/hdst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[2]) +
    geom_line(data[data$reward_type == "ltgt/ldst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[3]) +
    geom_line(data[data$reward_type == "ltgt/hdst", ], mapping=aes_string(x=iv, y=dv, group="sub"), alpha=alpha, colour=cols[4]) +
    facet_wrap(.~eval(parse(text=facet)), nrow=1) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) + 
    ylab(dv) + xlab(iv) + ylim(ylims) +
    theme(panel.border = element_blank(), 
          panel.grid.major =   element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.title=element_blank()) +
    theme_cowplot()
}