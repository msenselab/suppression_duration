circleFun <- function(center = c(0,0), r = 2, npoints = 100){
  t <- seq(0,2*pi,length.out = npoints)
  x <- center[1] + r * cos(t)
  y <- center[2] + r * sin(t)
  return(data.frame(x = x, y = y))
}

readData <- function(file) {
  dat <- read.csv(file)
  
  addmissing <- function(dat, list) {
    for(n in list) {
      if(!(n %in% colnames(dat))) {
        dat[n] <- NA
      }
    }
    return(dat) 
  }
  dat <- addmissing(dat, c("fix_check_3.started", "fix_check_3.stopped", "text_feedback.started",
                           "text_feedback.stopped", "text_feedback_2.stopped"))
  return(dat)
}

fix_pattern <- function(tar_fix) {
  if (1 %in% tar_fix && 2 %in% tar_fix) {
    tar1_fix <- which(tar_fix==1)
    tar2_fix <- which(tar_fix==2)
    first_tar_fixed <- ifelse(min(min(tar1_fix),min(tar2_fix)) %in% tar1_fix,1,2)
    if(first_tar_fixed == 1) {
      if(max(tar1_fix)>min(tar2_fix)) {
        return("Refixation")
      } else {
        return("Both fix.")
      }
    } else {
      if(max(tar2_fix)>min(tar1_fix)) {
        return("Refixation")
      } else {
        return("Both fix.")
      }
    }
  } else if (1 %in% tar_fix || 2 %in% tar_fix) {
    return("Single tar. fix.")
  } else {
    return("No fix.")
  }
}

compare <- function(filelist) {
  if (length(filelist) == 0) {
    return;
  }
  largest <- c()
  smallest <- names(read.csv(filelist[1]))
  for(f in filelist) {
    n <- names(read.csv(f))
    largest <- union(largest,n)
    smallest <- intersect(smallest,n)
  }
  return(setdiff(largest,smallest))
}


within_error <- function(data, var=rt, name="rt") {
  var = enquo(var)
  mdat <- data %>% group_by(participant) %>% summarize(mrt = mean(!!var))
  sdata <- full_join(data, mdat, by = c('participant')) %>% 
    mutate(drt = !!var - mrt) %>% 
    summarize(mrt = mean(!!var), se_rt = sd(drt)/sqrt(n()-1), ci_rt = se_rt * qt(0.975,n()-1)) 
  names <- names(sdata)
  names[(length(names)-2):length(names)] <-  c(paste0("m_",name),paste0("se_",name),paste0("ci_",name))
  sdata <- set_names(sdata, names)
  return(sdata)
}
