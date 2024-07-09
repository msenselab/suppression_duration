library(tidyverse)
library(cowplot)
source("functions.R")

# ---- Load data ----

b_data_files <- c('behavioral/ds_eye_track_exp_2024-02-06_13h32.24.549.csv',
                  'behavioral/2_eye_track_exp_2024-02-06_14h57.08.086.csv',
                  'behavioral/AI_eye_track_exp_2024-02-06_16h35.17.298.csv',
                  'behavioral/ES_eye_track_exp_2024-02-08_13h01.40.124.csv',
                  'behavioral/cml_eye_track_exp_2024-02-08_14h57.30.496.csv',
                  'behavioral/AK_eye_track_exp_2024-02-08_16h06.12.429.csv',
                  'behavioral/JM_eye_track_exp_2024-02-08_17h32.52.232.csv',
                  'behavioral/TS_eye_track_exp_2024-02-09_10h36.30.013.csv',
                  'behavioral/CAN_eye_track_exp_2024-02-09_11h54.07.842.csv',
                  'behavioral/NGS_eye_track_exp_2024-02-09_14h55.39.801.csv',
                  'behavioral/di_eye_track_exp_2024-02-09_16h20.31.132.csv',
                  'behavioral/by_eye_track_exp_2024-02-12_12h17.15.003.csv',
                  'behavioral/DA_eye_track_exp_2024-02-12_13h39.16.268.csv',
                  'behavioral/CIW_eye_track_exp_2024-02-15_16h49.45.708.csv',
                  'behavioral/MSR_eye_track_exp_2024-02-16_10h38.39.135.csv',
                  'behavioral/s1_eye_track_exp_2024-01-19_13h07.43.352.csv')

bdata <- reduce(lapply(b_data_files, readData), rbind)
bdata <- bdata %>% mutate(blocktype = ifelse(is.na(trials.thisN), "Practice", "Experiment"),
                          tno = ifelse(is.na(trials.thisN), practice_trials.thisRepN,
                                       trials.thisN)) %>% filter(blocktype=="Experiment")

# Add target locations and distractor location
bdata <- bdata %>% 
  pivot_longer(cols = c(ori10, ori11, ori12, ori13, ori14, ori16, ori17, ori18, ori19, ori20),
               names_to = "loc", values_to = "ori") %>% group_by(tno, participant) %>% 
  mutate(tar1_pos = which(ori %in% c(12, -12, 168, -168))[1],
         tar2_pos = which(ori %in% c(12, -12, 168, -168))[2], 
         dist_pos = ifelse(is_empty(which(ori %in% c(270, 90))),
                           0 ,which(ori %in% c(270, 90)))) %>% spread(loc, ori) %>%
  group_by(participant) %>%
  mutate(dist_region = ifelse(dist_pos == 0, "Absent", 
                              ifelse(dist_freq=="freq_bottom", ifelse(dist_pos > 5,
                                                                      "Freq.", "Rare"),
                                     ifelse(dist_pos > 5, "Rare", "Freq."))),
         responses = factor(responses, levels=c("resp", "inhibit"), labels=c("Go", "No-go")),
         distractor = ifelse(dist_pos > 0, "Present", "Absent"),
         dist_rep = ifelse(dist_pos == 0, "Absent", ifelse(dist_pos == lag(dist_pos), "Repeat", "Change")))

subs <- unique(bdata$participant)
bdata$distractor <- factor(bdata$distractor, levels = c("Absent", "Present"), ordered = TRUE)

# ---- behavioral data analysis ----

rt_data <- bdata %>% group_by(dist_region, distractor, participant) %>% 
  filter(response.corr == 1, responses=="Go", response.rt > 0.2) %>% 
  summarize(rt=mean(response.rt)*1000) 

we_rt_fig <- within_error(rt_data) %>% 
  ggplot(aes(x=dist_region, y = m_rt, linetype=distractor)) + 
  geom_point() + geom_errorbar(aes(ymin=m_rt-ci_rt, ymax=m_rt+ci_rt), width = 0.2) + 
  labs(x="Distractor region", y="Mean RT (ms)") + coord_cartesian(ylim = c(950,1350)) + 
  theme_classic() + theme(legend.position = "None") 

pc_data <-  bdata %>% group_by(dist_region, distractor, responses, participant)  %>% 
  summarize(pc=mean(response.corr)) 

we_pc_fig_go <- within_error(pc_data, pc, "pc") %>% filter(responses=="Go")  %>% 
  ggplot(aes(x=dist_region, y = m_pc, linetype=distractor)) +
  geom_point(color="black", aes(fill=dist_region)) + 
  geom_errorbar(aes(ymin=m_pc-ci_pc, ymax=m_pc+ci_pc), width = 0.2) + theme_classic() + 
  labs(x="Distractor region", y="Accuracy") + coord_cartesian(ylim = c(0.85, 1)) +
  theme(legend.position = "None") 

nrep_rt_data <- bdata %>% filter(dist_rep %in% c("Absent", "Change")) %>% group_by(dist_region, participant) %>% 
  filter(response.corr == 1, responses=="Go", response.rt > 0.2) %>% 
  summarize(rt=mean(response.rt)*1000) %>% mutate(dist_present = ifelse(dist_region == "Absent", "Absent", "Present")) %>%
  group_by(dist_region, dist_present)

nrep_rt_fig_e2 <- within_error(nrep_rt_data) %>%
  ggplot(aes(x=dist_region, y = m_rt, linetype = dist_present)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=m_rt-ci_rt, ymax=m_rt+ci_rt), width = 0.2) + theme_classic() + 
  labs(x="Distractor region", y="Reaction time (ms)", linetype = "Distractor") + 
  theme(legend.position = "bottom")

snrep_rt_data_e1 <- read.csv('../data exp1/behavioral//nrep_rt_e1.csv')

nrep_rt_fig_e1 <- snrep_rt_data_e1 %>%
  ggplot(aes(x=dist_region, y = m_rt, linetype = dist_present)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=m_rt-ci_rt, ymax=m_rt+ci_rt), width = 0.2) + theme_classic() + 
  labs(x="Distractor region", y="Reaction time (ms)", linetype = "Distractor") + 
  theme(legend.position = "bottom")

nrep_rt_fig <- plot_grid(nrep_rt_fig_e1, nrep_rt_fig_e2, nrow = 1, ncol=2, labels=c("A","B"))

ggsave('figures/nrep_rt.png', nrep_rt_fig, width = 7, height = 3)

# ---- Load processed data ----

all_fixations <- readRDS('eye-tracking summary/fixation_summary.rds')
all_fixations <- full_join(all_fixations, select(bdata, participant, tno, dist_rep), by = c('participant', 'tno')) 
all_saccades <- readRDS('eye-tracking summary/saccade_summary.rds')
all_saccades <- full_join(all_saccades, select(bdata, participant, tno, dist_rep), by = c('participant', 'tno')) 

# ---- Classify fixations ----
item_theta <- c(30, 60, 90, 120, 150, 210, 240, 270, 300, 330)*pi/180
item_x <- 6*cos(item_theta)
item_y <- 6*sin(item_theta)
all_fixations <- all_fixations %>% group_by(tno, participant) %>%
  mutate(tar1_X = item_x[tar1_pos], tar1_Y = item_y[tar1_pos], 
         tar1_dist = sqrt((tar1_X - average_gaze_x)^2+ 
                            (tar1_Y - average_gaze_y)^2),
         tar2_X = item_x[tar2_pos], tar2_Y = item_y[tar2_pos], 
         dist_X = ifelse(dist_pos > 0, item_x[dist_pos], NA), 
         dist_Y = ifelse(dist_pos > 0, item_y[dist_pos], NA),
         tar2_dist = sqrt((tar2_X - average_gaze_x)^2+ 
                            (tar2_Y - average_gaze_y)^2),
         fix_dist = sqrt(average_gaze_x^2+average_gaze_y^2),
         dist_dist = sqrt((dist_X - average_gaze_x)^2+ 
                            (dist_Y - average_gaze_y)^2),
         fix_type = ifelse((tar1_dist < 3 & !(!is.na(dist_dist) & tar1_dist > dist_dist)) |
                             (tar2_dist < 3 &!(!is.na(dist_dist) & tar2_dist > dist_dist)), "Target", 
                           ifelse(fix_dist < 3, "Fixation", 
                                  ifelse(!is.na(dist_dist) & dist_dist < 3, "Distractor", 
                                         "Other"))),
         tar_fixed = ifelse(fix_type != "Target", NA, ifelse(tar1_dist < tar2_dist, 1, 2)),
         fix_pattern = fix_pattern(tar_fixed[cat=="Search"]),
         responses = factor(responses, levels=c("Go", "No-go")))

all_fixations$fix_pattern <-
  factor(all_fixations$fix_pattern, levels = 
           c("No fix." , "Single tar. fix.", "Both fix." , "Refixation"), 
         ordered = TRUE)

# ---- First saccade to target ----

# Find first target fixations
first_tar_fix <- filter(all_fixations, cat=="Search", fix_type=="Target") %>% 
  select(tno, fix_count, participant) %>% 
  group_by(tno, participant) %>% summarize(first_tar_fix = min(fix_count))

all_fixations <- full_join(all_fixations, first_tar_fix, 
                           by=c("tno", "participant"))

first_tar_fixed <- all_fixations %>% group_by(participant, tno) %>% 
  filter(fix_count == first_tar_fix) %>% summarize(first_tar_fixed = tar_fixed)

all_fixations <- full_join(all_fixations, first_tar_fixed, 
                           by=c("tno", "participant"))

both_tar_fixed <- all_fixations %>% group_by(tno, participant) %>% 
  filter(cat=="Search", fix_type=="Target", fix_count > first_tar_fix, 
         tar_fixed!=first_tar_fixed) %>% summarize(both_tar_fixed = min(fix_count))

all_fixations <- full_join(all_fixations, both_tar_fixed, 
                           by=c("tno", "participant"))

# Find first target fixation time
first_fix_time <- all_fixations %>% filter(fix_count==first_tar_fix) %>% 
  select(tno, participant, first_tar_fix, first_fix_time = start_time)

all_fixations <- full_join(all_fixations, first_fix_time %>% select(-first_tar_fix),
                          by=c("tno", "participant"))

all_saccades <- full_join(all_saccades, first_fix_time, 
                          by=c("tno", "participant"))

# Find both target fixation time

both_fix_time <- all_fixations %>% filter(fix_count==both_tar_fixed) %>% 
  select(tno, participant, both_tar_fixed, both_fix_time = start_time)

all_fixations <- full_join(all_fixations, both_fix_time %>% select(-both_tar_fixed),
                           by=c("tno", "participant"))

all_saccades <- full_join(all_saccades, both_fix_time, 
                          by=c("tno", "participant"))

# Find first saccade to target
first_tar_sacc <- all_saccades %>% filter(start_time < first_fix_time) %>% 
  group_by(tno, participant) %>% filter(start_time==max(start_time), cat=="Search") %>% 
  mutate(latency=(start_time - search_start_time)*1000)

pj = position_dodge(width = 0.4)
tar_sacc_lat <- first_tar_sacc %>% group_by(dist_region, responses, participant, tno) %>%
  filter(response.corr==1) %>% summarize(latency=mean(latency)) %>% summarize(latency=mean(latency))

# ---- First saccade to distractor ----

# Find first distractor fixations
first_dist_fix <- filter(all_fixations, cat=="Search", fix_type=="Distractor") %>% 
  select(tno, fix_count, participant) %>% 
  group_by(tno, participant) %>% summarize(first_dist_fix = min(fix_count))

all_fixations <- full_join(all_fixations, first_dist_fix, 
                           by=c("tno", "participant"))

# Find first distractor fixation time
first_dist_fix_time <- all_fixations %>% filter(fix_count==first_dist_fix) %>% 
  select(tno, participant, first_dist_fix, first_dist_fix_time = start_time)

all_fixations <- full_join(all_fixations, first_dist_fix_time %>% select(-first_dist_fix),
                           by=c("tno", "participant"))

all_saccades <- full_join(all_saccades, first_dist_fix_time, 
                          by=c("tno", "participant"))

# Find first saccade to distractor

first_dist_sacc <- all_saccades %>% filter(start_time < first_dist_fix_time) %>% 
  group_by(tno, participant) %>% filter(start_time==max(start_time), cat=="Search") %>% 
  mutate(latency=(start_time - search_start_time)*1000)

pj = position_dodge(width = 0.4)
dist_sacc_lat <- first_dist_sacc %>% group_by(dist_region, responses, participant, tno) %>%
  filter(response.corr==1) %>% summarize(latency=mean(latency)) %>% summarize(latency=mean(latency))

# ---- Bar plot of proportion of first target fixations ----

first_fix_dat_all <- all_fixations %>% group_by(participant, tno) %>% 
  filter(cat=="Search", fix_type!="Fixation", 
         response.corr == 1) %>% filter(fix_count==min(fix_count)) %>% 
  group_by(dist_region, responses, participant) 

first_fix_dat_all$fix_type = factor(first_fix_dat_all$fix_type, 
                                    levels = c("Distractor", "Target", "Other"),
                                    ordered = TRUE)

pj = position_dodge(width = 0.6)
fig_first_fix_bar <- first_fix_dat_all %>% filter(dist_region!="Absent") %>% 
  group_by(dist_region, responses, participant, tno) %>% summarize(fix_type=unique(fix_type)) %>%
  ggplot(aes(x=1, fill=fix_type, 
             y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  geom_histogram(stat="count", color="black") + theme_classic() + 
  facet_grid(dist_region~responses) + theme(strip.background = element_blank()) + 
  labs(x="Distractor region", y="Proportion first fixation", fill="Fixated object") + 
  coord_flip() + scale_y_reverse(breaks=c(0,0.25,0.5,0.75,1), 
                                 labels=c("1","0.75","0.5","0.25","0")) + 
  scale_fill_grey() +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(), 
        legend.position = "bottom")  #remove y axis ticks

ggsave("./figures/first_fix_bar.png", fig_first_fix_bar, width = 7.5, height = 3)

# ---- Plot target location effect ----

tar_loc_dat <- all_fixations %>% 
  mutate(tar1_region = ifelse(dist_freq=="freq_bottom", 
                              ifelse(tar1_pos > 5, "Frequent", "Rare"), 
                              ifelse(tar1_pos > 5, "Rare", "Frequent")), 
         tar2_region = ifelse(dist_freq=="freq_bottom", 
                              ifelse(tar2_pos > 5, "Frequent", "Rare"),
                              ifelse(tar2_pos > 5, "Rare", "Frequent")), 
         tar_regions = ifelse(tar1_region==tar2_region, tar1_region, 
                              "Freq./Rare"), 
         first_region_fixed = ifelse(first_tar_fixed==1, tar1_region, tar2_region)) %>% 
  group_by(participant, tno) 

tar_loc_dat$tar_regions <- factor(tar_loc_dat$tar_regions, levels=c("Frequent", 
                                   "Freq./Rare", "Rare"),
                                  labels=c("Freq. dist.", 
                                           "Freq./Rare dist.", "Rare dist."),
                                  ordered =TRUE)

tar_loc_rt_dat <- tar_loc_dat %>% 
  filter(response.corr == 1, 
         responses=="Go", dist_region=="Absent", !is.na(tar_regions)) %>% 
  group_by(tar_regions, responses, participant, tno) %>%
  summarize(rt=mean(response.rt)*1000) %>%
  summarize(rt=mean(rt))

star_loc_rt_dat <- tar_loc_rt_dat %>% 
  summarize(mrt=mean(rt), se_rt=sd(rt)/sqrt(n()-1))

fig_tar_loc_we <- within_error(tar_loc_rt_dat) %>% ggplot(aes(x=tar_regions, y=m_rt)) + geom_point() + 
  geom_errorbar(aes(ymin=m_rt-ci_rt, ymax=m_rt+ci_rt),width=0.2) + theme_classic() +
  labs(x="Target regions", y="Mean RT (ms)")  + scale_y_continuous(limits = c(850,1250))

we_rt_pc_fig <- plot_grid(we_rt_fig, we_pc_fig_go, fig_tar_loc_we, nrow = 1, ncol=3, 
                       labels=c("A","B","C"))

ggsave('./figures/we_rt-pc.png', we_rt_pc_fig, width = 9.5, height = 3)

# ---- Classify distractor fixations ----
  
fig_fix_pattern <- all_fixations %>% group_by(responses, participant, tno) %>%
    filter(!is.na(responses)) %>% summarize(fix_pattern=unique(fix_pattern)) %>%
    ggplot(aes(x=1, fill=fix_pattern, y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
    geom_histogram(stat="count", color="black") + facet_grid(responses~.) + theme_classic() +
    labs(x="Response", y="Proportion of trials", fill = "Fix. pattern") + 
    coord_flip() + scale_y_reverse(breaks=c(0,0.25,0.5,0.75,1), 
                                   labels=c("1","0.75","0.5","0.25","0")) + 
    scale_fill_grey() + 
    theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank())  #remove y axis ticks
  
ggsave("./figures/fix_pattern.png", fig_fix_pattern, width = 6, height = 2)

dist_fixations <- all_fixations %>% mutate(dfix_type = 
                           ifelse(cat=="Search" & fix_type=="Distractor",
                                  ifelse(is.na(first_tar_fix) | fix_count < first_tar_fix, "Before", 
                                                  ifelse(is.na(both_tar_fixed) | fix_count < both_tar_fixed, 
                                                               "Between", "After")), NA)) 

# Distractor fixation proportions

sdist_fix_prop <- filter(dist_fixations, cat=="Search", response.corr==1, !is.na(both_tar_fixed)) %>% 
  group_by(dist_region, responses, participant, tno) %>% filter(dist_region!="Absent") %>% 
  summarize(before_count = sum(dfix_type=="Before", na.rm=T)>0, 
            between_count = sum(dfix_type=="Between", na.rm=T)>0, 
            after_count = sum(dfix_type=="After", na.rm=T)>0) %>%
  summarize(before_prop = mean(before_count), after_prop = mean(after_count), 
            between_prop = mean(between_count))

we_ssdist_fix_prop <- cbind(within_error(sdist_fix_prop, before_prop, 'before'), 
                            within_error(sdist_fix_prop, between_prop, 'between') %>%
                              ungroup() %>% select(-c(responses, dist_region)),
                            within_error(sdist_fix_prop, after_prop, 'after') %>%
                              ungroup() %>% select(-c(responses, dist_region)))

we_fix_prop_spread <- we_ssdist_fix_prop %>% 
  pivot_longer(cols = c(m_before, m_between, m_after, 
                        ci_before, ci_between, ci_after),
               names_to=c("stat_type", "fix_type"), names_sep='_') %>%
  spread(stat_type, value)

we_fix_prop_spread$fix_type <- factor(we_fix_prop_spread$fix_type, 
                                   levels = c("before", "between", "after"),
                                   labels = c("Before", "Between", "After"),
                                   ordered = TRUE)

we_fig_fix_prop <- we_fix_prop_spread %>%
  ggplot(aes(x=dist_region, y=m, color=responses,
             group=responses)) + 
  geom_point() + geom_errorbar(aes(ymin=m-ci,ymax=m+ci), width=0.3) + 
  geom_line() + theme_classic() + facet_wrap(~fix_type) + 
  theme(strip.background = element_blank()) + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Prop. of distractor fixations", color="Response")

ggsave('./figures/we_fix_prop.png', we_fig_fix_prop, width = 6.5, height=2.8)

# ---- Repeated distractor fixations -----

srep_dist_fix_prop <- filter(dist_fixations, cat=="Search", response.corr==1, !is.na(both_tar_fixed)) %>% 
  group_by(dist_region, responses, participant, tno) %>% filter(dist_region!="Absent") %>% 
  summarize(before_fix = sum(dfix_type=="Before", na.rm=T)>0, 
            between_fix = sum(dfix_type=="Between", na.rm=T)>0, 
            after_fix = sum(dfix_type=="After", na.rm=T)>0) %>% 
  filter(between_fix==0) %>%
  group_by(before_fix, dist_region, responses, participant) %>%
  summarize(between_prop=mean(between_fix), after_prop=mean(after_fix)) 

good_subs <- filter(summarize(group_by(srep_dist_fix_prop, participant), N=n()), N==8)$participant

srep_dist_fix_prop <- filter(srep_dist_fix_prop, participant %in% good_subs)

we_ssrep_dist_fix_prop <- within_error(srep_dist_fix_prop, after_prop, 'after')

we_ssrep_dist_fix_prop$before_fix <- factor(we_ssrep_dist_fix_prop$before_fix,
                                            levels=c(T, F), labels=c("Fix.", "No fix."))

fig_dist_fix_after_we <- we_ssrep_dist_fix_prop %>%
  ggplot(aes(x=dist_region, y=m_after, color=responses, shape=before_fix, linetype=before_fix,
             group=interaction(responses, before_fix))) + 
  geom_point() + geom_errorbar(aes(ymin=m_after-ci_after,ymax=m_after+ci_after), width=0.3) + 
  geom_line() + theme_classic()  + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Prop. of post-target dist. fix.", color="Resp.", 
       shape="Pre-tar.", linetype="Pre-tar.") + theme(legend.position = "bottom")

ggsave('./figures/we_dist_fix_after.png', fig_dist_fix_after_we, width=4, height=3)

# ---- Distractor fixation prop. without repetitions ----

sdist_fix_prop_norep <- filter(dist_fixations, cat=="Search", response.corr==1) %>% 
  group_by(dist_region, responses, participant, tno) %>% filter(dist_rep=="Change") %>% 
  summarize(before_count = sum(dfix_type=="Before", na.rm=T)>0, 
            between_count = sum(dfix_type=="Between", na.rm=T)>0, 
            after_count = sum(dfix_type=="After", na.rm=T)>0) %>%
  summarize(before_prop = mean(before_count), after_prop = mean(after_count), 
            between_prop = mean(between_count))

we_ssdist_fix_prop_norep <- cbind(within_error(sdist_fix_prop_norep, before_prop, 'before'), 
                                  within_error(sdist_fix_prop_norep, between_prop, 'between') %>%
                                    ungroup() %>% select(-c(responses, dist_region)),
                                  within_error(sdist_fix_prop_norep, after_prop, 'after') %>%
                                    ungroup() %>% select(-c(responses, dist_region)))

we_fix_prop_spread_norep <- we_ssdist_fix_prop_norep %>% 
  pivot_longer(cols = c(m_before, m_between, m_after, 
                        ci_before, ci_between, ci_after),
               names_to=c("stat_type", "fix_type"), names_sep='_') %>%
  spread(stat_type, value)

we_fix_prop_spread_norep$fix_type <- factor(we_fix_prop_spread_norep$fix_type, 
                                            levels = c("before", "between", "after"),
                                            labels = c("Before", "Between", "After"),
                                            ordered = TRUE)

we_fig_fix_prop_norep_e2 <- we_fix_prop_spread_norep %>%
  ggplot(aes(x=dist_region, y=m, color=responses,
             group=responses)) + 
  geom_point() + geom_errorbar(aes(ymin=m-ci,ymax=m+ci), width=0.3) + 
  geom_line() + theme_classic() + facet_wrap(~fix_type) + 
  theme(strip.background = element_blank()) + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Prop. of distractor fixations", color="Response")

ggsave('./figures/we_fix_prop_norep.png', we_fig_fix_prop_norep_e2, width = 6.5, height=2.8)
