library(tidyverse)
library(cowplot)
source("functions.R")

# ---- Load data ----

b_data_files <- c('behavioral/001_eye_tracking_exp_2021_Nov_08_0830.csv',
                  'behavioral/002_eye_tracking_exp_2021_Nov_06_1955.csv',
                  'behavioral/003_eye_tracking_exp_2021_Nov_08_1002.csv',
                  'behavioral/004_eye_tracking_exp_2021_Nov_08_1130.csv',
                  'behavioral/005_eye_tracking_exp_2021_Nov_08_1304.csv',
                  'behavioral/006_eye_tracking_exp_2021_Nov_08_1505.csv',
                  'behavioral/007_eye_tracking_exp_2021_Nov_08_1616.csv',
                  'behavioral/008_eye_tracking_exp_2021_Nov_08_1735.csv',
                  'behavioral/009_eye_tracking_exp_2021_Nov_08_1913.csv',
                  'behavioral/010_eye_tracking_exp_2021_Nov_09_1137.csv',
                  'behavioral/011_eye_tracking_exp_2021_Nov_09_1400.csv',
                  'behavioral/012_eye_tracking_exp_2021_Nov_10_1011.csv',
                  'behavioral/013_eye_tracking_exp_2021_Nov_10_1137.csv',
                  'behavioral/014_eye_tracking_exp_2021_Nov_10_1340.csv',
                  'behavioral/015_eye_tracking_exp_2021_Nov_10_1439.csv')

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
                           0 ,which(ori %in% c(270, 90)))) %>% spread(loc, ori) %>% group_by(participant) %>%
  mutate(dist_region = ifelse(dist_pos == 0, "Absent", 
                              ifelse(dist_freq=="freq_bottom", ifelse(dist_pos > 5,
                                                                      "Freq.", "Rare"),
                                     ifelse(dist_pos > 5, "Rare", "Freq."))), 
         respon = factor(respon, levels=c("resp", "inhibit"), labels=c("Go", "No-go")),
         dist_rep = ifelse(dist_pos == 0, "Absent", ifelse(dist_pos == lag(dist_pos), "Repeat", "Change")))

subs <- unique(bdata$participant)

# ---- behavioral data analysis ----

rt_data <- bdata %>% group_by(dist_region, participant) %>% 
  filter(response.corr == 1, respon=="Go", response.rt > 0.2) %>% 
  summarize(rt=mean(response.rt)*1000) 

err_data <-  bdata %>% group_by(dist_region, respon, participant)  %>% 
  summarize(err=1-mean(response.corr)) 

we_rt_fig <- within_error(mutate(rt_data, distractor = ifelse(dist_region == "Absent", "Absent", "Present")) %>% 
                            group_by(dist_region, distractor)) %>% 
  ggplot(aes(x=dist_region, y = m_rt, linetype=distractor)) + 
  geom_point() + geom_errorbar(aes(ymin=m_rt-ci_rt, ymax=m_rt+ci_rt), width = 0.2) + 
  labs(x="Distractor region", y="Mean RT (ms)") + coord_cartesian(ylim = c(950,1400)) + 
  theme_classic() + theme(legend.position = "None") 

we_err_fig <- within_error(mutate(err_data, distractor = ifelse(dist_region == "Absent", "Absent", "Present"),
                                  acc = 1 - err) %>% 
                             group_by(dist_region, distractor), acc, 'acc')  %>% 
  ggplot(aes(x=dist_region, y = m_acc, linetype=distractor)) + 
  geom_point() + geom_errorbar(aes(ymin=m_acc-ci_acc, ymax=m_acc+ci_acc), width = 0.2) + 
  labs(x="Distractor region", y="Accuracy") + coord_cartesian(ylim = c(0.85,1)) + 
  theme_classic() + theme(legend.position = "None") 

nrep_rt_data <- bdata %>% filter(dist_rep %in% c("Absent", "Change")) %>% group_by(dist_region, participant) %>% 
  filter(response.corr == 1, respon=="Go", response.rt > 0.2) %>% 
  summarize(rt=mean(response.rt)*1000) %>% mutate(dist_present = ifelse(dist_region == "Absent", "Absent", "Present")) %>% 
  group_by(dist_region, dist_present) 

write.csv(within_error(nrep_rt_data), './behavioral/nrep_rt_e1.csv')

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
         respon = factor(respon, levels=c("resp", "inhibit"), labels=c("Go", "No-go")))

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
tar_sacc_lat <- first_tar_sacc %>% group_by(dist_region, respon, participant, tno) %>%
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
dist_sacc_lat <- first_dist_sacc %>% group_by(dist_region, respon, participant, tno) %>%
  filter(response.corr==1) %>% summarize(latency=mean(latency)) %>% summarize(latency=mean(latency))

# ---- Bar plot of proportion of first target fixations ----

first_fix_dat_all <- all_fixations %>% group_by(participant, tno) %>% 
  filter(cat=="Search", fix_type!="Fixation", 
         response.corr == 1) %>% filter(fix_count==min(fix_count)) %>% 
  group_by(dist_region, respon, participant) 

first_fix_dat_all$fix_type = factor(first_fix_dat_all$fix_type, 
                                    levels = c("Distractor", "Target", "Other"),
                                    ordered = TRUE)

pj = position_dodge(width = 0.6)
fig_first_fix_bar <- first_fix_dat_all %>% filter(dist_region!="Absent") %>% 
  group_by(dist_region, respon, participant, tno) %>% summarize(fix_type=unique(fix_type)) %>%
  ggplot(aes(x=1, fill=fix_type, 
             y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
  geom_histogram(stat="count", color="black") + theme_classic() + 
  facet_grid(dist_region~respon) + theme(strip.background = element_blank()) + 
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
                                  ordered = TRUE)

tar_loc_rt_dat <- tar_loc_dat %>% 
  filter(response.corr == 1, 
         respon=="Go", dist_region=="Absent", !is.na(tar_regions)) %>% 
  group_by(tar_regions, respon, participant, tno) %>%
  summarize(rt=mean(response.rt)) %>%
  summarize(rt=mean(rt))

star_loc_rt_dat <- within_error(tar_loc_rt_dat %>% mutate(rt=rt*1000))

fig_tar_loc_we <- star_loc_rt_dat  %>% ggplot(aes(x=tar_regions, y=m_rt)) + geom_point() + 
  geom_errorbar(aes(ymin=m_rt-ci_rt, ymax=m_rt+ci_rt),width=0.2) + theme_classic() +
  labs(x="Target regions", y="Mean RT (ms)") + coord_cartesian(ylim = c(950,1200))

we_rt_pc_fig <- plot_grid(we_rt_fig, we_err_fig, fig_tar_loc_we, labels = c("A", "B", "C"), ncol = 3)

ggsave('./figures/we_rt-pc.png', we_rt_pc_fig, width = 9, height = 3)

# ---- Plot first target saccade latency on only trials with first fixation- on target ----

first_search_fix <- all_fixations %>% group_by(participant, tno) %>% 
  filter(cat=="Search") %>% summarize(first_search_fix = min(fix_count))

first_tar_sacc_lat <- full_join(first_tar_sacc, first_search_fix) %>%
  group_by(dist_region, participant, tno) %>% filter(response.corr == 1, latency > 50, 
                                                     latency < 500,
                                                     first_tar_fix==first_search_fix) %>%
  summarize(latency=mean(latency)) %>% summarize(latency=mean(latency))

first_dist_sacc_lat <- full_join(first_dist_sacc, first_search_fix) %>%
  group_by(dist_region, participant, tno) %>% filter(response.corr == 1, latency > 50, 
                                                     latency < 500,
                                                     first_dist_fix==first_search_fix) %>%
  summarize(latency=mean(latency)) %>% summarize(latency=mean(latency))

goodsubs <- intersect(intersect(filter(first_tar_sacc_lat, dist_region=="Absent")$participant,
                               filter(first_tar_sacc_lat, dist_region=="Freq.")$participant),
                      filter(first_tar_sacc_lat, dist_region=="Rare")$participant)

first_tar_sacc_lat$type <- "Target"
first_dist_sacc_lat$type <- "Distractor"
first_sacc_lat <- rbind(filter(first_tar_sacc_lat, participant %in% goodsubs),
                        filter(first_dist_sacc_lat, participant %in% goodsubs))

we_sacc_lat_fig <- within_error(group_by(first_sacc_lat, dist_region, type), latency, 'lat') %>%
  ggplot(aes(x=dist_region, y=m_lat, color=type, group=type)) + geom_point(position=pj, size=2) + 
  geom_errorbar(aes(ymin=m_lat-ci_lat, ymax=m_lat+ci_lat),width=0.4, position=pj) +
  geom_line(position=pj) +  theme_classic() + scale_color_grey(start=0.4, end=0.6) +
  labs(x="Distractor region", y="Saccadic latency (ms)")

ggsave("figures/we_sacc_lat.png", we_sacc_lat_fig, width = 4, height = 2.5)

# ---- Classify distractor fixations ----

fig_fix_pattern <- all_fixations %>% group_by(respon, participant, tno) %>%
    filter(!is.na(respon)) %>% summarize(fix_pattern=unique(fix_pattern)) %>%
    ggplot(aes(x=1, fill=fix_pattern, y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..])) +
    geom_histogram(stat="count", color="black") + facet_grid(respon~.) + theme_classic() +
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
  group_by(dist_region, respon, participant, tno) %>% filter(dist_region!="Absent") %>% 
  summarize(before_count = sum(dfix_type=="Before", na.rm=T)>0, 
            between_count = sum(dfix_type=="Between", na.rm=T)>0, 
            after_count = sum(dfix_type=="After", na.rm=T)>0) %>%
  summarize(before_prop = mean(before_count), after_prop = mean(after_count), 
            between_prop = mean(between_count))

we_ssdist_fix_prop <- cbind(within_error(sdist_fix_prop, before_prop, 'before'), 
                            within_error(sdist_fix_prop, between_prop, 'between') %>%
                              ungroup() %>% select(-c(respon, dist_region)),
                            within_error(sdist_fix_prop, after_prop, 'after') %>%
                              ungroup() %>% select(-c(respon, dist_region)))

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
  ggplot(aes(x=dist_region, y=m, color=respon,
             group=respon)) + 
  geom_point() + geom_errorbar(aes(ymin=m-ci,ymax=m+ci), width=0.3) + 
  geom_line() + theme_classic() + facet_wrap(~fix_type) + 
  theme(strip.background = element_blank()) + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Prop. of distractor fixations", color="Response")

ggsave('./figures/we_fix_prop.png', we_fig_fix_prop, width = 6.5, height=2.8)

# ---- Distractor fixation durations first ----

sdist_fix_dur_first <- filter(dist_fixations, cat=="Search",  fix_type=="Distractor", response.corr==1, !is.na(both_tar_fixed)) %>% 
  group_by(dist_region, dfix_type, respon, participant, tno) %>% filter(fix_count==min(fix_count))%>% 
  summarize(dur = sum(duration)*1000) %>% summarize(dur = mean(dur)) %>% 
  filter(dfix_type!="Between") %>% spread(dfix_type, dur)

missing_data_subs <- union(unique(filter(sdist_fix_dur_first, is.na(After))$participant), 
                           unique(filter(sdist_fix_dur_first, is.na(Before))$participant))

sdist_fix_dur_first <- filter(sdist_fix_dur_first, !(participant %in% missing_data_subs)) %>%
  gather(type, dur, -c(dist_region, respon, participant))

we_ssdist_fix_dur_first <- within_error(group_by(sdist_fix_dur_first, dist_region, respon, type), dur, 'dur')

we_ssdist_fix_dur_first$type <- factor(we_ssdist_fix_dur_first$type,
                                    levels=c("Before", "After"),
                                    ordered=TRUE)

we_dist_fix_dur_fig <- we_ssdist_fix_dur_first %>%
  ggplot(aes(x=dist_region, y=m_dur, color=respon,
             group=respon)) + 
  geom_point() + geom_errorbar(aes(ymin=m_dur-ci_dur,ymax=m_dur+ci_dur), width=0.3) + 
  geom_line() + theme_classic() + facet_wrap(~type) + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Distractor fix. dur. (ms)", color="Response") +
  theme(strip.background = element_blank()) + lims(y=c(140,330))

ggsave('./figures/we_fix_dur_first.png', we_dist_fix_dur_fig, width = 5, height=2.8)

# ---- Target fixation durations first vs. second ----

first_tar_fix_dur <- all_fixations %>% filter(fix_type=="Target") %>%
  filter(fix_count==first_tar_fix) %>% 
  group_by(dist_region, respon, participant) %>% 
  filter(response.corr == 1) %>%
  summarize(dur=mean(duration)*1000) 

second_tar_fix_dur <- all_fixations %>% filter(fix_type=="Target") %>%
  filter(fix_count==both_tar_fixed) %>% 
  group_by(dist_region, respon, participant) %>% 
  filter(response.corr == 1) %>%
  summarize(dur=mean(duration)*1000) 

first_tar_fix_dur$tar_fixed <- "First"
second_tar_fix_dur$tar_fixed <- "Second"
both_tar_fix_dur <- rbind(first_tar_fix_dur, second_tar_fix_dur)

pj = position_dodge(width = 0.4)
we_tar_fix_dur_fig <- within_error(group_by(both_tar_fix_dur,respon,dist_region,tar_fixed), dur, 'dur') %>%
  ggplot(aes(x = dist_region, y = m_dur, color=respon, group=respon)) + geom_point(position=pj) + 
  geom_line(position=pj) + theme_classic() + facet_grid(~tar_fixed) +
  geom_errorbar(aes(ymin=m_dur - ci_dur, ymax=m_dur + ci_dur), width=0.2, position=pj) +
  labs(x="Distractor region", y="Target fix. dur. (ms)") +
  theme(strip.background = element_blank()) + scale_color_grey(start=0.2, end=0.6) + 
  lims(y=c(160,350))

we_fix_dur_fig <- plot_grid(we_dist_fix_dur_fig, we_tar_fix_dur_fig, nrow = 2, ncol=1, 
                         labels=c("A","B"))

ggsave("figures/we_fix_dur.png", we_fix_dur_fig, width=5,height=5)

# ---- Repeated distractor fixations -----

srep_dist_fix_prop <- filter(dist_fixations, cat=="Search", response.corr==1, !is.na(both_tar_fixed)) %>% 
  group_by(dist_region, respon, participant, tno) %>% filter(dist_region!="Absent") %>% 
  summarize(before_fix = sum(dfix_type=="Before", na.rm=T)>0, 
            between_fix = sum(dfix_type=="Between", na.rm=T)>0, 
            after_fix = sum(dfix_type=="After", na.rm=T)>0) %>% 
  filter(between_fix==0) %>%
  group_by(before_fix, dist_region, respon, participant) %>%
  summarize(between_prop=mean(between_fix), after_prop=mean(after_fix)) 

good_subs <- intersect(unique(filter(srep_dist_fix_prop, before_fix==TRUE, dist_region=="Rare",respon=="No-go"))$participant,
                       unique(filter(srep_dist_fix_prop, before_fix==TRUE, dist_region=="Rare",respon=="Go"))$participant)

srep_dist_fix_prop <- filter(srep_dist_fix_prop, participant %in% good_subs)

we_ssrep_dist_fix_prop <- within_error(srep_dist_fix_prop, after_prop, 'after')

we_ssrep_dist_fix_prop$before_fix <- factor(we_ssrep_dist_fix_prop$before_fix,
                                         levels=c(T, F), labels=c("Fix.", "No fix."))

fig_dist_fix_after_we <- we_ssrep_dist_fix_prop %>%
  ggplot(aes(x=dist_region, y=m_after, color=respon, shape=before_fix, linetype=before_fix,
             group=interaction(respon, before_fix))) + 
  geom_point() + geom_errorbar(aes(ymin=m_after-ci_after,ymax=m_after+ci_after), width=0.3) + 
  geom_line() + theme_classic()  + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Prop. of post-target distractor fixations", color="Response", 
       shape="Pre-target", linetype="Pre-target")

ggsave('./figures/we_dist_fix_after.png', fig_dist_fix_after_we, width=4, height=3)

# ---- Distractor fixation prop. without repetitions ----

sdist_fix_prop_norep <- filter(dist_fixations, cat=="Search", response.corr==1, !is.na(both_tar_fixed)) %>% 
  group_by(dist_region, respon, participant, tno) %>% filter(dist_rep=="Change") %>% 
  summarize(before_count = sum(dfix_type=="Before", na.rm=T)>0, 
            between_count = sum(dfix_type=="Between", na.rm=T)>0, 
            after_count = sum(dfix_type=="After", na.rm=T)>0) %>%
  summarize(before_prop = mean(before_count), after_prop = mean(after_count), 
            between_prop = mean(between_count))

we_ssdist_fix_prop_norep <- cbind(within_error(sdist_fix_prop_norep, before_prop, 'before'), 
                            within_error(sdist_fix_prop_norep, between_prop, 'between') %>%
                              ungroup() %>% select(-c(respon, dist_region)),
                            within_error(sdist_fix_prop_norep, after_prop, 'after') %>%
                              ungroup() %>% select(-c(respon, dist_region)))

we_fix_prop_spread_norep <- we_ssdist_fix_prop_norep %>% 
  pivot_longer(cols = c(m_before, m_between, m_after, 
                        ci_before, ci_between, ci_after),
               names_to=c("stat_type", "fix_type"), names_sep='_') %>%
  spread(stat_type, value)

we_fix_prop_spread_norep$fix_type <- factor(we_fix_prop_spread_norep$fix_type, 
                                      levels = c("before", "between", "after"),
                                      labels = c("Before", "Between", "After"),
                                      ordered = TRUE)

we_fig_fix_prop_norep <- we_fix_prop_spread_norep %>%
  ggplot(aes(x=dist_region, y=m, color=respon,
             group=respon)) + 
  geom_point() + geom_errorbar(aes(ymin=m-ci,ymax=m+ci), width=0.3) + 
  geom_line() + theme_classic() + facet_wrap(~fix_type) + 
  theme(strip.background = element_blank()) + scale_color_grey(start=0.2, end=0.6) +
  labs(x="Distractor region", y = "Prop. of distractor fixations", color="Response")

ggsave('./figures/we_fix_prop_norep.png', we_fig_fix_prop_norep, width = 6.5, height=2.8)

