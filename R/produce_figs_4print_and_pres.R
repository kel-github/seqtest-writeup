## Written by K. Garner, 2020
#############################################################################################
# Generate plots for the sequence comparison write up that require further post-processing in
# lucid press
library(magrittr)
library(tidyverse)
library(cowplot)
library(readr)
library(wesanderson)
library(RJSONIO)
library(cowplot)
library(rlang)
library(kableExtra)
library(scales)
library(rstatix)
library(lme4) 
library(emmeans)
source("R_rainclouds.R") # for the raincloud plot
source("data_wrangles.R")



#############################################################################################
### Plot behavioural data by cue probability collapsed over value condition
### to go on plot with paradigm etc

subjects = c(1, 2, 3, 4, 5)
sessions = c(2, 3, 4)
TRs = c(700, 1510, 1920)
data_path = '../data/'
mri.beh.raw <- get_mri_data(subjects, sessions, data_path, TRs)

RT_min = .2
sd_reject = 2.5

mri.beh.clean <- mri.beh.raw %>% group_by(sub, sess, TR, cert) %>%
  filter(rt > RT_min) %>%
  filter(resp == 1) %>%
  ungroup()
mri.beh.crit <- mri.beh.clean %>% group_by(sub, sess, TR, cert) %>%
  summarise(crit = median(rt) + (sd_reject*sd(rt))) %>% ungroup()

mri.beh.clean <- inner_join(mri.beh.clean, mri.beh.crit, by=c("sub", "sess", "TR", "cert")) %>%
                 filter(rt < crit)

### calculate inv eff
mri.acc = mri.beh.raw %>% group_by(sub, TR, cert) %>%
  summarise(acc = mean(resp))
mri.inv.eff = mri.beh.clean %>% group_by(sub, TR, cert) %>%
  summarise(RT = median(rt)) %>%
  inner_join(mri.acc, sum.inv.eff, by=c("sub", "TR", "cert")) %>%
  transform(inv_eff = RT/acc)

cert.inv.eff.grp <- mri.inv.eff %>% group_by(TR, cert) %>%
  summarise(mu=mean(inv_eff),
            N=length(inv_eff),
            se=sd(inv_eff)/sqrt(N))
fig.cols = wes_palette("IsleofDogs1")
iv = "cert"
grp = "TR"
cols = fig.cols[c(1,3,4)]
ylb = "IE"
ylims=c(0,5)

plt.crt.sum <- function(data, iv, grp, cols, ylb, ylims){
  data$TR <- as.factor(data$TR)
  ggplot(data, aes_string(x=iv, y="mu", col=grp, group="TR")) +
    geom_line(lwd=1.5) +
    geom_errorbar(aes(ymin=mu-se, ymax=mu+se), lwd=1.5, width=.2) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) + 
    ylab(ylb) + xlab(iv) + ylim(ylims) +
    theme_cowplot() +
    theme(panel.border = element_blank(), 
          panel.grid.major =   element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) 
}
plt.crt.sum(cert.inv.eff.grp, iv, grp, cols, ylb, c(0.5,.85))


#############################################################################################




plt.crt.sum(cert.inv.eff.grp, iv, grp, cols, ylb, c(0.5,.85))

cert.acc = mri.beh.raw %>% group_by(sub, TR, cert) %>%
  summarise(acc = mean(resp))
cert.inv.eff = mri.beh.clean %>% group_by(sub, TR, cert) %>%
  summarise(RT = median(rt)) %>%
  inner_join(mri.acc, sum.inv.eff, by=c("sub", "TR", "cert")) %>%
  transform(inv_eff = RT/acc)



#############################################################################################
### Plot behavioural data by cue probability collapsed over value condition
tmp <- mri.inv.eff %>% group_by(sub, cert, TR) %>% 
                       summarise(RT = mean(RT),
                                 acc = mean(acc),
                                 inv_eff = mean(inv_eff))

tmp.sum <- tmp %>% group_by(cert, TR) %>%
           summarise(mu=mean(acc))
names(tmp.sum ) <- c("cert", "TR", "mu")
tmp.sum$cert <- factor(tmp.sum$cert, c(".8", ".5", ".2"))
# now reorder the un-summarised data
tmp$cert <- factor(tmp$cert, c(".8", ".5", ".2"))

# splitting the data into each reward condition for overlaying on plot
yval = "acc"
alpha = 1/3
cols = wes_palette("IsleofDogs1")
tp <- ggplot(tmp.sum, aes_string(x="cert", y="mu", group=1)) +
  geom_line(size=1.1, color=cols[4]) +
  geom_line(tmp[tmp$sub == "1", ], mapping=aes_string(x="cert", y=yval, group=1), alpha=alpha, colour=cols[1]) +
  geom_line(tmp[tmp$sub == "2", ], mapping=aes_string(x="cert", y=yval, group=1), alpha=alpha, colour=cols[2]) +
  geom_line(tmp[tmp$sub == "3", ], mapping=aes_string(x="cert", y=yval, group=1), alpha=alpha, colour=cols[3]) + 
  geom_line(tmp[tmp$sub == "4", ], mapping=aes_string(x="cert", y=yval, group=1), alpha=alpha, colour=cols[5]) + 
  geom_line(tmp[tmp$sub == "5", ], mapping=aes_string(x="cert", y=yval, group=1), alpha=alpha, colour=cols[6]) + 
  facet_wrap(.~TR, nrow=1) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) + 
  ylab("Accuracy") + xlab("Cue Certainty") + ylim(0.6,1) +
  theme(panel.border = element_blank(), 
        panel.grid.major =   element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank()) + theme_cowplot()

ggsave("../images/Fig_ACC.png", plot=tp, width = 10, height = 5, units="cm")
rm(tmp)
rm(tmp.sum)
rm(tp)

############################ CNR values
h <- draw.sub.plts(CNR, C="hand") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../images/Fig_CNR_hand.png", plot=h, width = 10, height = 12, units="cm")
ggsave("../images/Fig_CNR_hand.pdf", plot=h, width = 10, height = 12, units="cm")


cp <- draw.sub.plts(CNR, C="tgtLoc") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("../images/Fig_CNR_tgtLoc.png", plot=cp, width = 10, height = 12, units="cm")
ggsave("../images/Fig_CNR_tgtLoc.pdf", plot=cp, width = 10, height = 12, units="cm")

########################### TR Table
kable(t(TRs), caption="Table 1. Compared sequences") %>% save_kable(file="../images/test.png")
