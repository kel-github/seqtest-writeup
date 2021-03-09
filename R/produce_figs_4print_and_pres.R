## Written by K. Garner, 2020
#############################################################################################
# This code has copies of the plotting code in seq-test-writeup.Rmd that can be run alongside
# the document. Is useful for making small plots for talks etc that won't appear in the main
# writeup.
# Note: a lot of it is hardcoded, as its just a scratch pad to make bespoke images from
# the data


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
