## Written by K. Garner, 2020
#############################################################################################
# Generate plots for the sequence comparison write up that require further post-processing in
# lucid press
# this code needs to be run along side the document in order to produce outputs


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
mri.dat.4.filt <- mri.beh.clean %>% group_by(sub, sess, TR, cert) %>%
  summarise(m = median(rt),
            s = sd_reject*sd(rt),
            f = m + s) %>%
  ungroup()
mri.beh.clean <- inner_join(mri.beh.clean, mri.dat.4.filt, by=c("sub", "sess", "TR", "cert")) %>%
  filter(rt < f)

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
ggsave(filename="../images/behav_by_Prot.png", plot=last_plot(), 
                 width = 2.45, height = 2.95, units = "in", dpi=300)

#############################################################################################
# plot distributions of CNR data 
CNR_dens <- ggplot(CNR, aes(x=CNR, fill=sub)) +
  geom_density(alpha=0.4, trim=TRUE) +
  facet_grid(rows=vars(TR)) +
  xlim(c(5.44, 25)) + xlab("CNR") + 
  theme_cowplot() +
  scale_colour_manual(values=wes_palette("FantasticFox1")[1:5]) +
  scale_fill_manual(values=wes_palette("FantasticFox1")[1:5]) +
  theme( axis.text.y = element_blank(),
         strip.background = element_blank(),
         strip.text.x = element_blank())
ggsave(filename='~/Dropbox/documents/MC-Docs/seqtest-writeup/images/thrsh_comp_densities.png', 
       plot=last_plot(), width= 2.45, height= 2.95, units="in", dpi=300)

CNR.p <- roi_cnr %>% filter(roi == "IPS" | roi == "CN" | roi == "Put") %>%
  group_by(sub, TR, roi) %>% summarise(R=mean(CNR)) %>%
  ggplot(aes(x=TR, y=R, group=sub)) +
  geom_line(aes(color=sub), lwd=1.1, alpha=.75) +
  xlab(expression(italic("P"))) + 
  facet_wrap(~roi, nrow=1) +
  scale_colour_manual(values=c(wes_palette("FantasticFox1")[c(2:5)])) +
  ylab("CNR") + 
  theme(strip.text.x = element_text(size=8)) + 
  theme_cowplot() +
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         legend.position = "none",
         text = element_text(size = 10),
         axis.text.x = element_text(size = 8),
         axis.text.y = element_text(size = 8))
ggsave(filename='~/Dropbox/documents/MC-Docs/seqtest-writeup/images/CNR_by_key_reg.png', 
       plot=CNR.p, width = 3.5, height = 2.45, units="in", dpi=300)

#############################################################################################
### Plot FIR by selected regions for each sequence
mu_firs %>% filter(roi == "IPS" | roi == "CN" | roi == "Put") %>%
  ggplot(aes(x=order, y=beta, group=sub)) +
  geom_line(aes(color=sub), lwd=1.1, alpha=.75) +
  xlab("t") +
  scale_x_continuous(breaks=seq(2,18,by=2),
                     labels = c("2","","","","10","","","","18")) +
  facet_grid(vars(TR), vars(roi), scales="free_y") +
  scale_colour_manual(values=c(wes_palette("FantasticFox1")[c(2:5)])) +
  ylab(expression(beta)) + theme_cowplot() +
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         strip.text.y = element_blank(),
         legend.position = "none",
         text = element_text(size = 10),
         axis.text.x = element_text(size = 8),
         axis.text.y = element_text(size = 8))
ggsave(filename='~/Dropbox/documents/MC-Docs/seqtest-writeup/images/FIR_by_reg.png', 
       plot=last_plot(), width=3.5, height=2.95, units="in", dpi=300)

#############################################################################################
### Plot tSNR data by key regions
tSNR.p <- mutSNR %>% filter(roi == "IPS" | roi == "CN" | roi == "Put") %>%
                    ggplot(aes(x=TR, y=tSNR, group=sub)) +
                    geom_line(aes(color=sub), lwd=1.1, alpha=.75) +
                    xlab(expression(italic("P"))) + 
                    facet_wrap(~roi, nrow=1) +
                    scale_colour_manual(values=c(wes_palette("FantasticFox1"))) +
                    ylab("tSNR") + 
                    theme(strip.text.x = element_text(size=8)) + 
                    theme_cowplot() +
                    theme( strip.background = element_blank(),
                           strip.text.x = element_blank(),
                           legend.position = "none",
                           text = element_text(size = 10),
                           axis.text.x = element_text(size = 8),
                           axis.text.y = element_text(size = 8))
ggsave(filename='~/Dropbox/documents/MC-Docs/seqtest-writeup/images/tSNR_by_key_reg.png', 
       plot=tSNR.p, width = 3.5, height = 2.45, units="in", dpi=300)


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
