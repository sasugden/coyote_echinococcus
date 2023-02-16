#### Figure 1: Infection status/intensity #####
intensityrow <- sample_data
intensityrow <- intensityrow[,c("ech.cap", "log1.worm", "location", "age.class", "new.pca")]

intensityrow$health.class <- rep("below\nmedian")
intensityrow$health.class[intensityrow$new.pca > median(intensityrow$new.pca)] <- "above\nmedian"
intensityrow$health.class <- factor(intensityrow$health.class)

intensityrow$log1.worm <- NULL
intensityrow$new.pca <- NULL

intensityrow <- reshape2::melt(intensityrow)

intensityrow1 <- intensityrow
intensityrow1$facet <- rep("location")
intensityrow1$axis <- intensityrow1$location
intensityrow1[,c("location", "age.class", "health.class")] <- NULL

intensityrow2 <- intensityrow
intensityrow2$facet <- rep("age class")
intensityrow2$axis <- intensityrow2$age.class
intensityrow2[,c("location", "age.class", "health.class")] <- NULL

intensityrow3 <- intensityrow
intensityrow3$facet <- rep("health class")
intensityrow3$axis <- intensityrow3$health.class
intensityrow3[,c("location", "age.class", "health.class")] <- NULL

intensityrow <- rbind(intensityrow1, intensityrow2, intensityrow3)

intensityrow$facet <- factor(intensityrow$facet, levels=c("location", "age class", "health class"))

prevrow <- sample_data
prevrow <- prevrow[,c("bioact", "location", "age.class", "new.pca")]

prevrow$health.class <- rep("below average")
prevrow$health.class[prevrow$new.pca > median(prevrow$new.pca)] <- "above average"
prevrow$health.class <- factor(prevrow$health.class)

prevrow$new.pca <- NULL

prevrow1 <- prevrow %>%
  dplyr::group_by(location) %>%
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(location, bioact) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(location, bioact, perc) %>%
  unique() %>%
  as.data.frame()

prevrow2 <- prevrow %>%
  dplyr::group_by(age.class) %>%
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(age.class, bioact) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(age.class, bioact, perc) %>%
  unique() %>%
  as.data.frame()

prevrow3 <- prevrow %>%
  dplyr::group_by(health.class) %>%
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(health.class, bioact) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(health.class, bioact, perc) %>%
  unique() %>%
  as.data.frame()

prevrow1$facet <- rep("location")
prevrow2$facet <- rep("age class")
prevrow3$facet <- rep("health class")
colnames(prevrow1) <- c("axis", "bioact", "perc", "facet")
colnames(prevrow2) <- c("axis", "bioact", "perc", "facet")
colnames(prevrow3) <- c("axis", "bioact", "perc", "facet")

prevrow <- rbind(prevrow1, prevrow2, prevrow3)
prevrow$facet <- factor(prevrow$facet, levels=c("location", "age class", "health class"))

prevrow$axis <- factor(prevrow$axis, levels=c("urban", "rural", "young", "old", "below average", "above average"))
intensityrow$axis <- factor(intensityrow$axis, levels=c("urban", "rural", "young", "old", "below\nmedian", "above\nmedian"))

p1 <- ggplot(subset(prevrow, bioact==1), aes(x=axis, y=perc)) + geom_bar(stat="identity", fill="darkgrey") + 
  facet_grid(~facet, scales="free_x") +
  theme_bw() + plot_theme + scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0))) +
  labs(y="Infection prevalence (%)") + theme(axis.text.x=element_blank(),
                                             axis.title.x=element_blank())

p2 <- ggplot(intensityrow, aes(x=axis, y=value+1)) + geom_boxplot(fill="darkgrey") +
  facet_grid(~facet, scales="free_x") + scale_y_log10() +
  theme_bw() + plot_theme +
  labs(y="Infection intensity (# of worms)", x="Group") +
  theme(strip.text=element_blank(), strip.background=element_blank())

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

p1$widths <- p2$widths

grid.arrange(p1, p2, nrow=2, heights=c(0.475, 0.525))

comp <- arrangeGrob(p1, p2, nrow=2, heights=c(0.475,0.525))

ggsave("~/Documents/Deanna_MS_data/intro.bars.jpg", comp, width=4.5, height=4, units="in", dpi=300)

#### Figure 2: Age/location and spleen/age/location interactions as predictors of infection status ####
plot1 <- ggplot(cd.logit1, aes(x=age.cement, y=fit, color=location)) +
  geom_line(size=1) +
  scale_color_manual(values=c("forestgreen","purple")) +
  guides(color = guide_legend(ncol=2, title.position="top", title.hjust=0.5, title = "Dietary predictor")) +
  theme_bw() + plot_theme +
  labs(x="Age (yr)", y="Probability of infection", tag="a") +
  theme(legend.position="none")

plot2 <- ggplot(cd.nb1, aes(x=age.cement, y=fit, color=location)) +
  geom_line(size=1) +
  scale_color_manual(values=c("forestgreen","purple")) +
  guides(color = guide_legend(ncol=1, title.position="top", title.hjust=0.5, title = "Location")) +
  theme_bw() + plot_theme +
  labs(x="Age (yr)", y="Natural log worm count", tag="b") +
  theme(legend.position="bottom")

cd.logit2$response <- rep("infection status")
cd.nb2$response <- rep("infection intensity")
cd.nb2$fit <- log(cd.nb2$fit, base=10)

scaling <- function(x){(x-min(x))/(max(x)-min(x))}

cd.logit2$fit <- scaling(cd.logit2$fit)
cd.nb2$fit <- scaling(cd.nb2$fit)

cd.test <- rbind(cd.logit2, cd.nb2)

levels(cd.test$age.cement)[levels(cd.test$age.cement)=="0.256"] <- "Younger (0.26 yr)"
levels(cd.test$age.cement)[levels(cd.test$age.cement)=="2.47"] <- "Mean (2.47 yr)"
levels(cd.test$age.cement)[levels(cd.test$age.cement)=="4.79"] <- "Older (4.79 yr)"

cd.test$location <- factor(cd.test$location, levels=c("urban", "rural"))

cd.test$response <- factor(cd.test$response, levels=c("infection status", "infection intensity"))

plot3 <- ggplot(cd.test, aes(x=splenic.index, y=fit, color=response)) +
  geom_line(size=1) + facet_grid(location ~ age.cement, scales="free_y") +
  guides(color = guide_legend(ncol=1, title.position="top", title.hjust=0.5, title = "Response variable")) +
  theme_bw() + plot_theme +
  labs(x="Spleen mass (normalized by body mass)", y="Infection status (pink) or intensity (blue)", tag="c") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom")

grid.arrange(plot1, plot2, plot3, 
             layout_matrix=rbind(c(1,3),c(2,3)), widths=c(0.25,0.75), heights=c(0.4,0.55))


plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)

plot2$widths <- plot1$widths


comp <- arrangeGrob(plot1, plot2, plot3, 
                    layout_matrix=rbind(c(1,3),c(2,3)), widths=c(0.25,0.75), heights=c(0.4,0.55))

ggsave("~/Documents/Deanna_MS_data/age_and_spleens.jpg", comp, width=6.5, height=5, units="in", dpi=300)

#### Figure 3: Predictors of infection status ####
plot2 <- ggplot(test, aes(x=xaxis, y=fit, color=variable, linetype=type)) +
  geom_line(size=1) + facet_grid(location ~ age.cement) +
  scale_linetype_manual(values=c("dashed", "dashed", "solid","solid","solid"),
                        guide = "none") +
  guides(color = guide_legend(ncol=2, title.position="top", title.hjust=0.5, title = "Dietary predictor")) +
  theme_bw() + plot_theme +
  labs(x="Predictor value (min to max)", y="Probability of infection", tag="a") +
  theme(axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        legend.position="bottom")

plot1 <- ggplot(subset(results, model=="status")) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  geom_linerange(aes(x=merge, y=psd_Coef, ymin=psd_25, ymax=psd_75),
                 size=2, position = position_dodge(width=0.5)) + 
  geom_linerange(aes(x=merge, y=psd_Coef, ymin=psd_2.5, ymax=psd_97.5),
                 size=0.5, position = position_dodge(width=0.5)) + 
  geom_point(aes(x=merge, y=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE,
             position = position_dodge(width=0.5)) +
  geom_vline(xintercept=c(4.5, 8.5, 12.5, 16.5, 20.5), color="black", size=1.5) +
  scale_x_discrete(limits=rev) +
  #scale_y_continuous(limits=c(-1.75, 1.75)) +
  scale_color_manual(values = c("red", "blue"), guide = guide_legend(title.position="top", title = "Dependent variable",
                                                                     title.hjust = 0.5)) +
  theme_bw() + coord_flip() +
  labs(x="Predictor", y="Coefficient", tag="b") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="bottom",
        axis.text.x=element_text(size=9, color="black"),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

grid.arrange(plot2, plot1, ncol=2, widths=c(0.6,0.4))

comp <- arrangeGrob(plot2, plot1, ncol=2, widths=c(0.6,0.4))

ggsave("~/Documents/Deanna_MS_data/status.model.figure.jpg", comp,
       width=8.5, height=6, units="in", dpi=300)
#### Figure 4: Bar graphs showing predictor values based on age, location, and infection status ####
sample_data_melt <- reshape2::melt(sample_data)
sample_data_melt <- subset(sample_data_melt, variable %in% c("sm.vol", "anth.dig.vol", "anth.indig.vol", "d13C", "d15N"))

new_plots <- sample_data_melt %>%
  dplyr::group_by(location, age.class, variable, bioact) %>%
  dplyr::summarise(avg=mean(value), sd=sd(value), se=plotrix::std.error(value))

new_plots$variable <- factor(new_plots$variable,
                             levels=c("sm.vol", "anth.dig.vol", "anth.indig.vol", "d13C", "d15N"))

levels(new_plots$variable)[levels(new_plots$variable)=="food.vol.washed"] <- "Total food"
levels(new_plots$variable)[levels(new_plots$variable)=="sm.vol"] <- "Rodents"
levels(new_plots$variable)[levels(new_plots$variable)=="anth.dig.vol"] <- "Dig. anthro."
levels(new_plots$variable)[levels(new_plots$variable)=="anth.indig.vol"] <- "Indig. anthro."

new_plots$fill <- paste0(new_plots$age.class, new_plots$bioact)
new_plots$fill <- factor(new_plots$fill, levels=c("young0", "young1", "old0", "old1"))

sample_data_melt$fill <- paste0(sample_data_melt$age.class, sample_data_melt$bioact)
sample_data_melt$fill <- factor(sample_data_melt$fill, levels=c("young0", "young1", "old0", "old1"))

p1 <- ggplot(subset(new_plots, location=="urban" & variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
             aes(x=age.class, y=avg, fill=fill)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + plot_theme +
  geom_errorbar(aes(x=age.class, ymin=avg-se, ymax=avg+se),
                position=position_dodge(width=0.9), color="black", width=0.3) +
  facet_wrap(~variable, scales="free_y", ncol=3) +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title="Urban", y="Mean volume (ml)") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p2 <- ggplot(subset(sample_data_melt, location=="urban" & variable %in% c("d13C", "d15N")),
             aes(x=age.class, y=value, fill=fill)) +
  geom_boxplot() +
  theme_bw() + plot_theme +
  facet_wrap(~variable, scales="free_y", ncol=2) +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y="Isotope value") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p3 <- ggplot(subset(new_plots, location=="rural" & variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
             aes(x=age.class, y=avg, fill=fill)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  geom_errorbar(aes(x=age.class, ymin=avg-se, ymax=avg+se),
                position=position_dodge(width=0.9), color="black", width=0.3) +
  facet_wrap(~variable, scales="free_y", ncol=3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title="Rural", y="Mean volume (ml)") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p4 <- ggplot(subset(sample_data_melt, location=="rural" & variable %in% c("d13C", "d15N")),
             aes(x=age.class, y=value, fill=fill)) +
  geom_boxplot() +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  facet_wrap(~variable, scales="free_y", ncol=2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y="Isotope value") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)

p2$heights <- p1$heights
p4$heights <- p3$heights

grid.arrange(p1, p2, p3, p4,
             ncol=2, nrow=2, widths=c(0.6, 0.4), heights=c(0.475,0.525))

legend_data <- subset(new_plots, age.class=="young")
levels(legend_data$bioact)[levels(legend_data$bioact)==0] <- "uninfected"
levels(legend_data$bioact)[levels(legend_data$bioact)==1] <- "infected"

legend_data2 <- subset(new_plots, age.class=="old")
levels(legend_data$bioact)[levels(legend_data$bioact)==0] <- "uninfected"
levels(legend_data$bioact)[levels(legend_data$bioact)==1] <- "infected"

legend1 <- cowplot::get_legend(
  ggplot(subset(legend_data, location=="rural" & 
                  variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
         aes(x=age.class, y=avg, fill=bioact)) +
    geom_bar(stat="identity", position="dodge") +
    theme_bw() + plot_theme +
    scale_fill_manual(values=c("cadetblue2", "salmon")) +
    guides(fill = guide_legend(title="Juveniles", title.hjust=0.5, title.position="top", ncol=1)) +
    theme(legend.position="bottom")
)
legend2 <- cowplot::get_legend(
  ggplot(subset(legend_data, location=="rural" & 
                  variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
         aes(x=age.class, y=avg, fill=bioact)) +
    geom_bar(stat="identity", position="dodge") +
    theme_bw() + plot_theme +
    scale_fill_manual(values=c("deepskyblue3", "orangered1")) +
    guides(fill = guide_legend(title="Adults", title.hjust=0.5, title.position="top", ncol=1)) +
    theme(legend.position="bottom")
)

grid.arrange(p1, p2, p3, p4, legend1, legend2,
             layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2),
                                 c(3,3,3,3,3,3,4,4,4,4),
                                 c(NA,NA,NA,5,5,6,6,NA,NA,NA)),
             heights=c(0.4,0.4,0.2))

comp <- arrangeGrob(p1, p2, p3, p4, legend1, legend2,
                    layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2),
                                        c(3,3,3,3,3,3,4,4,4,4),
                                        c(NA,NA,NA,5,5,6,6,NA,NA,NA)),
                    heights=c(0.4,0.4,0.2))

ggsave("~/Documents/Deanna_MS_data/diet_bars.jpg", comp, width=6.5, height=5, units="in", dpi=300)

#### Figure 5: Predictors of infection intensity ####
plot2 <- ggplot(nb.test, aes(x=xaxis, y=fit, color=variable, linetype=type)) +
  geom_line(size=1) + facet_grid(location ~ age.cement, scales="free_y") +
  scale_linetype_manual(values=c("dashed", "dashed", "solid","solid","solid"),
                        guide = "none") +
  #scale_y_continuous(limits=c(0, 20)) +
  guides(color = guide_legend(ncol=2, title.position="top", title.hjust=0.5, title = "Dietary predictor")) +
  theme_bw() + plot_theme +
  labs(x="Predictor value (min to max)", y="Infection intensity (log worm counts)", tag="a") +
  theme(axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        legend.position="bottom")

plot1 <- ggplot(subset(results, model=="intensity")) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  geom_linerange(aes(x=merge, y=psd_Coef, ymin=psd_25, ymax=psd_75),
                 size=2, position = position_dodge(width=0.5)) + 
  geom_linerange(aes(x=merge, y=psd_Coef, ymin=psd_2.5, ymax=psd_97.5),
                 size=0.5, position = position_dodge(width=0.5)) + 
  geom_point(aes(x=merge, y=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE,
             position = position_dodge(width=0.5)) +
  geom_vline(xintercept=c(4.5, 8.5, 12.5, 16.5, 20.5), color="black", size=1.5) +
  scale_x_discrete(limits=rev) +
  #scale_y_continuous(limits=c(-1.75, 1.75)) +
  scale_color_manual(values = c("red", "blue"), guide = guide_legend(title.position="top", title = "Dependent variable",
                                                                     title.hjust = 0.5)) +
  theme_bw() + coord_flip() +
  labs(x="Predictor", y="Coefficient") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.y = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="bottom",
        axis.text.x=element_text(size=9, color="black"),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))


grid.arrange(plot2, plot1, ncol=2, widths=c(0.6,0.4))

comp <- arrangeGrob(plot2, plot1, ncol=2, widths=c(0.6,0.4))

ggsave("~/Documents/Deanna_MS_data/intensity.model.figure.jpg", comp,
       width=8.5, height=6, units="in", dpi=300)

#### Figure A2, S1: Infection histograms ####
sample_data$group2 <- factor(sample_data$group2, levels=c("rur_yng", "rur_old", "urb_yng", "urb_old"))

p1 <- ggplot(subset(sample_data, ech.count==0), aes(x=ech.count, fill=group2)) + geom_histogram() +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks=c(0),
                     expand = expansion(mult = c(0, 0))) +
  labs(x=NULL, y="Number of coyotes") +
  scale_fill_manual(values=c("darkolivegreen3", "forestgreen", "mediumpurple1", "purple")) +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11), legend.position="none")

p2 <- ggplot(sample_data, aes(x=ech.count, fill=group2)) + geom_histogram() +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_log10(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x="Worm counts (number of worms)", y = NULL) +
  scale_fill_manual(values=c("darkolivegreen3", "forestgreen", "mediumpurple1", "purple")) +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11), legend.position="none",
        plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"))

grid.arrange(p1, p2, ncol=2, widths=c(0.1, 0.8))

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p1$heights <- p2$heights

ggpubr::ggarrange(p1, p2, ncol=2, widths=c(0.08, 0.8), labels=c("a"))

p3 <- ggplot(subset(sample_data, ech.count==0), aes(x=log1.worm, fill=group2)) + geom_histogram() +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(breaks=c(0),
                     expand = expansion(mult = c(0, 0))) +
  labs(x=NULL, y="Number of coyotes") +
  scale_fill_manual(values=c("darkolivegreen3", "forestgreen", "mediumpurple1", "purple")) +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11), legend.position="none")

p4 <- ggplot(subset(sample_data, ech.count > 0), aes(x=log1.worm, fill=group2)) + geom_histogram() +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  #scale_x_log10(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x="Worm counts (natural log-transformed)", y = NULL) +
  scale_fill_manual(values=c("darkolivegreen3", "forestgreen", "mediumpurple1", "purple")) +
  theme(axis.text=element_text(size=11), axis.title=element_text(size=11), legend.position="none",
        plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"))

p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)
p3$heights <- p4$heights

ggpubr::ggarrange(p3, p4, ncol=2, widths=c(0.08, 0.8), labels=c("b"))

top <- ggpubr::ggarrange(
  ggpubr::ggarrange(p1, p2, ncol=2, widths=c(0.1, 0.8), labels=c("a")),
  ggpubr::ggarrange(p3, p4, ncol=2, widths=c(0.1, 0.8), labels=c("b")), nrow=2)

levels(sample_data$group2)[levels(sample_data$group2)=="urb_yng"] <- "urban, juvenile"
levels(sample_data$group2)[levels(sample_data$group2)=="urb_old"] <- "urban, adult"
levels(sample_data$group2)[levels(sample_data$group2)=="rur_yng"] <- "rural, juvenile"
levels(sample_data$group2)[levels(sample_data$group2)=="rur_old"] <- "rural, adult"

legend <- cowplot::get_legend(ggplot(subset(sample_data, ech.count > 0), aes(x=log1.worm, fill=group2)) + geom_histogram() +
                                theme_bw() + plot_theme +
                                scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                                #scale_x_log10(expand = expansion(mult = c(0.01, 0.01))) +
                                labs(x="Worm counts (natural log-transformed)", y = NULL) +
                                scale_fill_manual(values=c("darkolivegreen3", "forestgreen", "mediumpurple1", "purple")) +
                                guides(fill=guide_legend(ncol=2, title.position="top", title.hjust=0.5, title=NULL)) +
                                theme(axis.text=element_text(size=11), axis.title=element_text(size=11), legend.position="bottom",
                                      plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt")))

top <- ggpubr::ggarrange(
  ggpubr::ggarrange(p1, p2, ncol=2, widths=c(0.11, 0.8), labels=c("a")),
  ggpubr::ggarrange(p3, p4, ncol=2, widths=c(0.11, 0.8), labels=c("b")), nrow=2)

grid.arrange(top, legend, nrow=2, heights=c(0.85, 0.1))

histograms <- arrangeGrob(top, legend, nrow=2, heights=c(0.85, 0.1))

ggsave("~/Documents/Deanna_MS_data/histograms.jpg", histograms, width=6.5, height=6, units="in", dpi=300)
#### Figure A2, S2: Composite condition metric ####
condition.ord <- as.data.frame(scores(condition.pca, choices=c(1:2), display = "sites"))
condition.ord$sampleID <- sample_data$sampleID

condition.ord <- merge(condition.ord, sample_data[,c("sampleID", "location", "age.cement", "age.class", 
                                                     "bioact", "mass", "girth", "length", "kfi")],
                       by="sampleID", all=TRUE)

cors <- condition.ord[,c("MDS1", "MDS2", "mass", "girth", "length", "kfi")]
cors <- Hmisc::rcorr(as.matrix(cors), type="pearson")

write.csv(cors$r, "~/Documents/Deanna_MS_data/cor.r.csv")
write.csv(cors$P, "~/Documents/Deanna_MS_data/cor.p.csv")

levels(condition.ord$bioact)[levels(condition.ord$bioact)=="0"] <- "uninfected"
levels(condition.ord$bioact)[levels(condition.ord$bioact)=="1"] <- "infected"

fit <- envfit(condition.pca, condition.ord[,c("mass","girth","length","kfi")])

spp.scrs.BC <- as.data.frame(scores(fit, display="vectors"))
spp.scrs.BC <- cbind(spp.scrs.BC, Rsq = fit[["vectors"]]$r,
                     pvals = fit[["vectors"]]$pvals, Variable = rownames(spp.scrs.BC))


p1 <- ggplot(condition.ord, aes(x=MDS1, y=MDS2)) + geom_point(aes(shape=location, color=bioact), size=2) +
  theme_bw() + plot_theme +
  geom_segment(data = spp.scrs.BC,
               aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="darkgrey") +
  scale_color_manual(values=c("deepskyblue3", "orangered1")) +
  geom_text(data = subset(spp.scrs.BC, Variable %in% c("length", "girth")),
            aes(x=MDS1+0.075, y=MDS2+0.075, label=Variable), size = 3.5) +
  geom_text(data = subset(spp.scrs.BC, Variable %in% c("mass")),
            aes(x=MDS1+0.075, y=MDS2-0.075, label=Variable), size = 3.5) +
  geom_text(data = subset(spp.scrs.BC, Variable %in% c("kfi")),
            aes(x=MDS1-0.05, y=MDS2-0.075, label=Variable), size = 3.5) +
  guides(color = guide_legend(title=expression(paste(italic("E. multilocularis"))), hjust=0.5, ncol=1, title.position="top"),
         shape = guide_legend(title="location", hjust=0.5, ncol=1, title.position="top")) +
  labs(x = "PC1 (67.6%)", y = "PC2 (20.6%)", tag = "b") +
  theme(legend.position="bottom")

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("",""
                   ,"Spearman's R = 0.630\np < 0.001 ",""),
  hjustvar = c(0,0,1.1,1) ,
  vjustvar = c(0,1,-0.5,1)) #<- adjust


p2 <- ggplot(condition.ord, aes(x=age.cement, y=MDS1)) + geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  # scale_y_continuous(limits=c(0,12)) +
  theme_bw() + plot_theme +
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=3.5) +
  labs(x="Cementum age (yr)", y="Condition measure (PC1)", tag="c")

grid.arrange(p1, p2, ncol=2)

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p2$heights <- p1$heights

pdata <- sample_data[,c("mass", "length", "girth", "kfi", "splenic.index")]
colnames(pdata) <- c("mass (kg)", "length (cm)", "girth (cm)", "kfi", "splenic.index (g/kg)")

GGally::ggpairs(pdata, 
                diag = NULL, 
                upper = list(continuous = wrap(ggally_cor, digits=2))) +
  theme_bw() + plot_theme +
  labs(tag="a")
theme(panel.grid=element_blank())

ptop <- ggplotGrob(ptop)


comp <- arrangeGrob(p1, p2, ncol=2)
ggsave("~/Documents/Deanna_MS_data/supp__condition.metrics.jpg", comp,
       width=6.5, height=4, units="in", dpi=300)

rm(p1, p2, comp, fit)

# Spleen mass and all the other items ###
GGally::ggpairs(sample_data[,c("mass","length","girth","kfi","splenic.index")], diag = NULL) + theme_bw()

melted <- reshape2::melt(sample_data[,c("mass","length","girth","kfi","splenic.index")])

melted <- cbind(melted, melted)
colnames(melted) <- c("xaxis", "valuex", "yaxis", "valuey")

ggplot(melted, aes(x=valuex, y=valuey)) + geom_point() + facet_grid(yaxis ~ xaxis, scales="free")



#### Figure A2, S3: Coyote age and sex distribution #####
age_plot <- sample_data
age_plot$bin <- rep(0.5)
age_plot$bin[age_plot$age.cement >= 1] <- 1.5
age_plot$bin[age_plot$age.cement >= 2] <- 2.5
age_plot$bin[age_plot$age.cement >= 3] <- 3.5
age_plot$bin[age_plot$age.cement >= 4] <- 4.5
age_plot$bin[age_plot$age.cement >= 5] <- 5.5
age_plot$bin[age_plot$age.cement >= 6] <- 6.5
age_plot$bin[age_plot$age.cement >= 7] <- 7.5
age_plot$bin[age_plot$age.cement >= 8] <- 8.5
age_plot$bin[age_plot$age.cement >= 9] <- 9.5
age_plot$bin[age_plot$age.cement >= 10] <- 10.5
age_plot$bin[age_plot$age.cement >= 11] <- 11.5

age_plot_count <- age_plot %>%
  dplyr::group_by(location, bin) %>%
  dplyr::count()

plot1 <- ggplot(age_plot_count, aes(x=bin, y=n)) + geom_bar(stat="identity") + 
  facet_wrap(~location, scales="free_y", nrow=2) +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12),
                     limits=c(0,12),
                     expand = expansion(mult=c(0, 0))) +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x="Cementum age (yr)", y="Number of individuals", tag="a")

sex_plot_count <- sample_data %>%
  dplyr::group_by(location, sex) %>%
  dplyr::count()
sex_plot_count$sex <- factor(sex_plot_count$sex, levels=c("M","F"))
levels(sex_plot_count$sex)[levels(sex_plot_count$sex)=="M"] <- "male"
levels(sex_plot_count$sex)[levels(sex_plot_count$sex)=="F"] <- "female"

plot2 <- ggplot(sex_plot_count, aes(x=sex, y=n)) + geom_bar(stat="identity") + 
  facet_wrap(~location, scales="free_y", nrow=2) +
  theme_bw() + plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x="Sex", y="Number of individuals", tag="b") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)

plot1$heights <- plot2$heights

grid.arrange(plot1, plot2, ncol=2, widths=c(0.75, 0.2))

comp <- arrangeGrob(plot1, plot2, ncol=2, widths=c(0.75, 0.2))

ggsave("~/Documents/Deanna_MS_data/age_sex_distribution.jpg", comp, width=6.5, height=4.5, units="in", dpi=300)


#### Figure A3, S1: Stomach content ordination ####
stomach_plot <- as.data.frame(scores(stomach.final.ords[[2]][[2]], choices=c(1:2), display="site"))

stomach_plot <- cbind(stomach.final.df[[2]][[2]], stomach_plot)

levels(stomach.final.df[[2]][[2]]$bioact)[levels(stomach.final.df[[2]][[2]]$bioact)=="0"] <- "uninfected"
levels(stomach.final.df[[2]][[2]]$bioact)[levels(stomach.final.df[[2]][[2]]$bioact)=="1"] <- "infected"

levels(stomach.spp.scr[[2]][[2]]$variable)[levels(stomach.spp.scr[[2]][[2]]$variable)=="ungulate.vol"] <- "ungulates"
levels(stomach.spp.scr[[2]][[2]]$variable)[levels(stomach.spp.scr[[2]][[2]]$variable)=="food.vol.washed"] <- "stomach volume"
levels(stomach.spp.scr[[2]][[2]]$variable)[levels(stomach.spp.scr[[2]][[2]]$variable)=="med.vol"] <- "meso-mammals"
levels(stomach.spp.scr[[2]][[2]]$variable)[levels(stomach.spp.scr[[2]][[2]]$variable)=="sm.vol"] <- "rodents"
levels(stomach.spp.scr[[2]][[2]]$variable)[levels(stomach.spp.scr[[2]][[2]]$variable)=="anth.dig.vol"] <- "dig. anthro"

stomach.final.df[[2]][[2]]$loc_bioact <- paste0(stomach.final.df[[2]][[2]]$location, ", ", stomach.final.df[[2]][[2]]$bioact)

stomach.plot <- ggplot(stomach.final.df[[2]][[2]], aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color=loc_bioact), size=2) +
  #stat_ellipse(aes(linetype=location), level = 0.75) +
  stat_ellipse(aes(color=loc_bioact), level = 0.75) +
  theme_bw() + plot_theme +
  labs(x = "PCoA1 (28.4%)", y="PCoA2 (20.1%)") +
  #scale_color_manual(values=c("deepskyblue3", "orangered1")) +
  #scale_linetype_manual(values=c("dotted","dashed")) +
  guides(shape = guide_legend(title.hjust=0.5),
         linetype = guide_legend(title.hjust=0.5),
         color = guide_legend(title=expression(paste(italic("E. multilocularis"))), title.hjust=0.5)) +
  
  geom_segment(data = subset(stomach.spp.scr[[2]][[2]], p < 0.05),
               aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow = arrow(length=unit(0.25, "cm")), color="darkgrey") +
  
  geom_text(data = subset(stomach.spp.scr[[2]][[2]], p < 0.05 & variable=="ungulates"),
            aes(x=MDS1+0.075, y=MDS2-0.075, label=variable), size = 3.5) +
  
  geom_text(data = subset(stomach.spp.scr[[2]][[2]], p < 0.05 & variable=="stomach volume"),
            aes(x=MDS1+0.075, y=MDS2-0.075, label=variable), size = 3.5) +
  
  geom_text(data = subset(stomach.spp.scr[[2]][[2]], p < 0.05 & variable=="meso-mammals"),
            aes(x=MDS1, y=MDS2-0.075, label=variable), size = 3.5) +  
  
  geom_text(data = subset(stomach.spp.scr[[2]][[2]], p < 0.05 & variable=="rodents"),
            aes(x=MDS1-0.075, y=MDS2-0.075, label=variable), size = 3.5) +  
  
  geom_text(data = subset(stomach.spp.scr[[2]][[2]], p < 0.05 & variable=="dig. anthro"),
            aes(x=MDS1+0.075, y=MDS2+0.075, label=variable), size = 3.5) 


ggsave("~/Documents/Deanna_MS_data/stomach_ordination.jpg", stomach.plot, width=6.5, height=5, units="in", dpi=300)


#### Figure A3, S2: Superspreaders ####
suprespread <- subset(sample_data, ech.count > 20000)
#suprespread$d13C <- -1*suprespread$d13C
suprespread <- reshape2::melt(suprespread)
suprespread <- subset(suprespread, variable %in% c("anth.dig.vol", "anth.indig.vol", "sm.vol", "d13C", "d15N", "age.cement",
                                                   "food.vol.washed", "splenic.index"))

ggplot(suprespread, aes(x=sampleID, y=value)) + geom_bar(stat="identity") + 
  facet_wrap(~variable, scales="free_y", ncol=4) +
  theme_bw() + plot_theme +
  scale_y_continuous(expand=expansion(mult = c(0, 0.1)))

#### Figure A3, S3: Infection status based on proximity to Leduc dump ####
temp <- sample_data_geo %>% 
  dplyr::group_by(site.tag) %>% 
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(site.tag, bioact) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(site.tag, bioact, perc) %>%
  unique() %>%
  as.data.frame()

temp <- subset(temp, bioact=="1")
temp$facet <- rep("Infection prevalence")

p1 <- ggplot(temp, aes(x=site.tag, y=perc, fill=site.tag)) + 
  geom_bar(stat="identity") +
  theme_bw() + plot_theme +
  facet_wrap(~facet) +
  scale_y_continuous(limits=c(0, 80), expand=expansion(mult=c(0,0))) +
  labs(x="Sample location", y="Infection prevalence (%)", tag="a") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

levels(sample_data_geo_melt$variable)[levels(sample_data_geo_melt$variable)=="ech.cap"] <- "E. multi counts"

p2 <- ggplot(subset(sample_data_geo_melt, variable=="E. multi counts"),
             aes(x=site.tag, y=value+1, fill=site.tag)) +
  geom_boxplot() + scale_y_log10() +
  theme_bw() + plot_theme +
  facet_wrap(~variable) +
  labs(x="Sample location", y="Worm count (n)", tag="b") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

levels(sample_data_geo_melt$variable)[levels(sample_data_geo_melt$variable)=="anth.indig.vol"] <- "Indig. anthro"

p3 <- ggplot(subset(sample_data_geo_melt, variable=="Indig. anthro"),
             aes(x=site.tag, y=value+1, fill=site.tag)) +
  geom_boxplot() + scale_y_log10() +
  theme_bw() + plot_theme +
  facet_wrap(~variable) +
  labs(x="Sample location", y="Volume (ml)", tag="c") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

levels(sample_data_geo_melt$variable)[levels(sample_data_geo_melt$variable)=="age.cement"] <- "Cementum age"

p4 <- ggplot(subset(sample_data_geo_melt, variable=="Cementum age"),
             aes(x=site.tag, y=value, fill=site.tag)) +
  geom_boxplot() +
  theme_bw() + plot_theme +
  facet_wrap(~variable) +
  labs(x="Sample location", y="Age (yr)", tag="d") +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))

grid.arrange(p1, p2, p3, p4, ncol=2)

comp <- arrangeGrob(p1, p2, p3, p4, ncol=2)

ggsave("~/Documents/Deanna_MS_data/worms_and_location.jpg", comp, width=5, height=6.5, units="in", dpi=300)

#### Figure A4, S1: Infection status/intensity for qPCR-positive infections #####
intensityrow <- sample_data
intensityrow <- intensityrow[,c("ech.cap", "log1.worm", "location", "age.class", "new.pca")]

intensityrow$health.class <- rep("below\nmedian")
intensityrow$health.class[intensityrow$new.pca > median(intensityrow$new.pca)] <- "above\nmedian"
intensityrow$health.class <- factor(intensityrow$health.class)

intensityrow$log1.worm <- NULL
intensityrow$new.pca <- NULL

intensityrow <- reshape2::melt(intensityrow)

intensityrow1 <- intensityrow
intensityrow1$facet <- rep("location")
intensityrow1$axis <- intensityrow1$location
intensityrow1[,c("location", "age.class", "health.class")] <- NULL

intensityrow2 <- intensityrow
intensityrow2$facet <- rep("age class")
intensityrow2$axis <- intensityrow2$age.class
intensityrow2[,c("location", "age.class", "health.class")] <- NULL

intensityrow3 <- intensityrow
intensityrow3$facet <- rep("health class")
intensityrow3$axis <- intensityrow3$health.class
intensityrow3[,c("location", "age.class", "health.class")] <- NULL

intensityrow <- rbind(intensityrow1, intensityrow2, intensityrow3)

intensityrow$facet <- factor(intensityrow$facet, levels=c("location", "age class", "health class"))

prevrow <- sample_data
prevrow <- subset(prevrow, is.na(qPCR.binary)=="FALSE")
prevrow <- prevrow[,c("qPCR.binary", "location", "age.class", "new.pca")]

prevrow$health.class <- rep("below average")
prevrow$health.class[prevrow$new.pca > median(prevrow$new.pca)] <- "above average"
prevrow$health.class <- factor(prevrow$health.class)

prevrow$new.pca <- NULL

prevrow1 <- prevrow %>%
  dplyr::group_by(location) %>%
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(location, qPCR.binary) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(location, qPCR.binary, perc) %>%
  unique() %>%
  as.data.frame()

prevrow2 <- prevrow %>%
  dplyr::group_by(age.class) %>%
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(age.class, qPCR.binary) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(age.class, qPCR.binary, perc) %>%
  unique() %>%
  as.data.frame()

prevrow3 <- prevrow %>%
  dplyr::group_by(health.class) %>%
  dplyr::mutate(countL = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(health.class, qPCR.binary) %>%
  dplyr::mutate(perc = 100*n()/countL) %>%
  dplyr::select(health.class, qPCR.binary, perc) %>%
  unique() %>%
  as.data.frame()

prevrow1$facet <- rep("location")
prevrow2$facet <- rep("age class")
prevrow3$facet <- rep("health class")
colnames(prevrow1) <- c("axis", "qPCR.binary", "perc", "facet")
colnames(prevrow2) <- c("axis", "qPCR.binary", "perc", "facet")
colnames(prevrow3) <- c("axis", "qPCR.binary", "perc", "facet")

prevrow <- rbind(prevrow1, prevrow2, prevrow3)
prevrow$facet <- factor(prevrow$facet, levels=c("location", "age class", "health class"))


prevrow$axis <- factor(prevrow$axis, levels=c("urban", "rural", "young", "old", "below average", "above average"))
intensityrow$axis <- factor(intensityrow$axis, levels=c("urban", "rural", "young", "old", "below\nmedian", "above\nmedian"))


p1 <- ggplot(subset(prevrow, qPCR.binary==1), aes(x=axis, y=perc)) + geom_bar(stat="identity", fill="darkgrey") + 
  facet_grid(~facet, scales="free_x") +
  theme_bw() + plot_theme + scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  labs(y="Infection prevalence (%)") + theme(axis.text.x=element_blank(),
                                             axis.title.x=element_blank())

p2 <- ggplot(intensityrow, aes(x=axis, y=value+1)) + geom_boxplot(fill="darkgrey") +
  facet_grid(~facet, scales="free_x") + scale_y_log10() +
  theme_bw() + plot_theme +
  labs(y="Infection intensity (# of worms)", x="Group") +
  theme(strip.text=element_blank(), strip.background=element_blank())

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)

p1$widths <- p2$widths

grid.arrange(p1, p2, nrow=2, heights=c(0.475, 0.525))

comp <- arrangeGrob(p1, p2, nrow=2, heights=c(0.475,0.525))

ggsave("~/Documents/Deanna_MS_data/qPCR__intro.bars.jpg", comp, width=4.5, height=4, units="in", dpi=300)


#### Figure A4, S2: qPCR vs. bioactive infection... as a function of age and location ####
qPCR_plot_data <- sample_data
qPCR_plot_data <- subset(qPCR_plot_data, is.na(qPCR.binary)=="FALSE")

qPCR_plot_data$group <- rep("uninfected")
qPCR_plot_data$group[qPCR_plot_data$bioact==1 & qPCR_plot_data$qPCR.binary==1] <- "biologically active"
qPCR_plot_data$group[qPCR_plot_data$bioact==0 & qPCR_plot_data$qPCR.binary==1] <- "qPCR-positive,\nno worms"

qPCR_plot_data <- qPCR_plot_data[order(-qPCR_plot_data$log1.worm, -qPCR_plot_data$qPCR.binary), ]
qPCR_plot_data$rank <- seq(1, 103, by=1)

plot1 <- ggplot(qPCR_plot_data,
                aes(x=rank, y=log1.worm, color=group)) +
  scale_color_manual(values=c("orangered1", "limegreen", "black")) +
  guides(color=guide_legend(title="Infection group", title.hjust=0.5, ncol=1, title.position="top")) +
  geom_point() +
  theme_bw() + plot_theme +
  theme(legend.position="bottom") +
  labs(x="Infection rank", y="Log-transformed worm count", tag="a")

qPCR_plot_data$bioact_qPCR <- paste0(qPCR_plot_data$bioact, "_", qPCR_plot_data$qPCR.binary)

sample_data2 %>% dplyr::group_by(bioact_qPCR) %>% dplyr::summarise(mean=mean(age.cement))

summary(aov(age.cement ~ bioact_qPCR, sample_data2))

TukeyHSD(aov(age.cement ~ bioact_qPCR, sample_data2))

qPCR_plot_data %>% dplyr::group_by(group) %>% dplyr::summarise(mean.age = mean(age.cement))

plot2 <- ggplot(qPCR_plot_data, aes(x=group, y=age.cement, fill=group)) + geom_boxplot() +
  facet_wrap(~location) +
  scale_fill_manual(values=c("orangered1", "limegreen", "darkgrey")) +
  theme_bw() + plot_theme +
  guides(fill=guide_legend(title="Infection group", title.hjust=0.5, ncol=1, title.position="top")) +
  theme(legend.position="bottom", axis.text.x=element_blank(), axis.title.x=element_blank()) +
  labs(y="Cementum age (yr)", tag="b")

plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)

grid.arrange(plot1, plot2, ncol=2)

comp <- arrangeGrob(plot1, plot2, ncol=2)

ggsave("~/Documents/Deanna_MS_data/qPCR.vs.bioact.age.jpg", comp, width=6.5, height=5, units="in", dpi=300)


#### Figure A4, S3: Figure 2 reproduced for qPCR-positive coyotes ####
# Prevalence: age/location interaction #
model <- glm(qPCR.binary ~ age.cement*location,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

cd.logit1 <- effect('age.cement*location', mod = model,
                    xlevels = list(age.cement=seq(min(sample_data$age.cement),
                                                  max(sample_data$age.cement),
                                                  by=0.1)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
cd.logit1 <- as.data.frame(cd.logit1)

ggplot(cd.logit1, aes(x=age.cement, y=fit)) +
  geom_point() + geom_line() + facet_grid(~location)

# Intensity: age/location interaction #
model <- glm.nb(log1.worm ~ age.cement*location,
                sample_data_qPCR,
                weights = weights.l,
                na.action="na.fail")

cd.nb1 <- effect('age.cement*location', mod = model,
                 xlevels = list(age.cement=seq(min(sample_data$age.cement),
                                               max(sample_data$age.cement),
                                               by=0.1)),
                 location = c("urban", "rural"),
                 se=TRUE, confidence.level=.95, typical=mean)
cd.nb1 <- as.data.frame(cd.nb1)

ggplot(cd.nb1, aes(x=age.cement, y=fit)) +
  geom_point() + geom_line() + facet_grid(~location)

# Prevalence: spleen/location/age interaction
model <- glm(qPCR.binary ~ age.cement * location * splenic.index,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

cd.logit2 <- effect('age.cement*location*splenic.index', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   splenic.index=seq(min(sample_data$splenic.index),
                                                     max(sample_data$splenic.index),
                                                     by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
cd.logit2 <- as.data.frame(cd.logit2)
cd.logit2$age.cement <- as.factor(as.character(cd.logit2$age.cement))

cd.logit2$age.cement <- factor(cd.logit2$age.cement,
                               levels=c(0.256, 2.47, 4.79))

ggplot(cd.logit2, aes(x=splenic.index, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement)

# Intensity: spleen/location/age interaction
model <- glm.nb(log1.worm ~ age.cement * location * splenic.index,
                sample_data_qPCR,
                weights = weights.l,
                na.action="na.fail")

cd.nb2 <- effect('age.cement*location*splenic.index', mod = model,
                 xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                splenic.index=seq(min(sample_data$splenic.index),
                                                  max(sample_data$splenic.index),
                                                  by=0.05)),
                 location = c("urban", "rural"),
                 se=TRUE, confidence.level=.95, typical=mean)
cd.nb2 <- as.data.frame(cd.nb2)
cd.nb2$age.cement <- as.factor(as.character(cd.nb2$age.cement))

cd.nb2$age.cement <- factor(cd.nb2$age.cement,
                            levels=c(0.256, 2.47, 4.79))

ggplot(cd.nb2, aes(x=splenic.index, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement) +
  scale_y_log10()


plot1 <- ggplot(cd.logit1, aes(x=age.cement, y=fit, color=location)) +
  geom_line(size=1) +
  scale_color_manual(values=c("forestgreen","purple")) +
  guides(color = guide_legend(ncol=2, title.position="top", title.hjust=0.5, title = "Dietary predictor")) +
  theme_bw() + plot_theme +
  labs(x="Age (yr)", y="Probability of infection", tag="a") +
  theme(legend.position="none")

plot2 <- ggplot(cd.nb1, aes(x=age.cement, y=fit, color=location)) +
  geom_line(size=1) +
  scale_color_manual(values=c("forestgreen","purple")) +
  guides(color = guide_legend(ncol=1, title.position="top", title.hjust=0.5, title = "Location")) +
  theme_bw() + plot_theme +
  labs(x="Age (yr)", y="Natural log worm count", tag="b") +
  theme(legend.position="bottom")

cd.logit2$response <- rep("infection status")
cd.nb2$response <- rep("infection intensity")
cd.nb2$fit <- log(cd.nb2$fit, base=10)

scaling <- function(x){(x-min(x))/(max(x)-min(x))}

cd.logit2$fit <- scaling(cd.logit2$fit)
cd.nb2$fit <- scaling(cd.nb2$fit)

cd.test <- rbind(cd.logit2, cd.nb2)

levels(cd.test$age.cement)[levels(cd.test$age.cement)=="0.256"] <- "Younger (0.26 yr)"
levels(cd.test$age.cement)[levels(cd.test$age.cement)=="2.47"] <- "Mean (2.47 yr)"
levels(cd.test$age.cement)[levels(cd.test$age.cement)=="4.79"] <- "Older (4.79 yr)"

cd.test$location <- factor(cd.test$location, levels=c("urban", "rural"))

cd.test$response <- factor(cd.test$response, levels=c("infection status", "infection intensity"))

plot3 <- ggplot(cd.test, aes(x=splenic.index, y=fit, color=response)) +
  geom_line(size=1) + facet_grid(location ~ age.cement, scales="free_y") +
  guides(color = guide_legend(ncol=1, title.position="top", title.hjust=0.5, title = "Response variable")) +
  theme_bw() + plot_theme +
  labs(x="Spleen nass (normalized by body mass)", y="Infection status (red) or intensity (blue)", tag="c") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom")

grid.arrange(plot1, plot2, plot3, 
             layout_matrix=rbind(c(1,3),c(2,3)), widths=c(0.25,0.75), heights=c(0.4,0.55))

plot1 <- ggplotGrob(plot1)
plot2 <- ggplotGrob(plot2)

plot2$widths <- plot1$widths

comp <- arrangeGrob(plot1, plot2, plot3, 
                    layout_matrix=rbind(c(1,3),c(2,3)), widths=c(0.25,0.75), heights=c(0.4,0.55))

ggsave("~/Documents/Deanna_MS_data/age_and_spleens_qPCR.jpg", comp, width=6.5, height=5, units="in", dpi=300)

#### Figure A4, S4: Figure 3 reproduced for qPCR-positive coyotes ####
# Run the model and calculate the weights
sample_data_qPCR <- subset(sample_data, is.na(qPCR.binary)=="FALSE")
sample_data_qPCR_scale <- scale.vars(sample_data_qPCR)

model <- glm(qPCR.binary ~ age.cement * location * (anth.dig.vol + sm.vol + d13C + d15N + anth.indig.vol),
             family="binomial"(link="logit"),
             sample_data_qPCR_scale,
             weights = weights.l,
             na.action="na.fail")
model.dredge <- dredge(model) # three top models with weights, six top models without weights

model.average <- model.avg(model.dredge, subset = delta < 2, beta="partial.sd")

temp2 <- cbind(confint(model.average, level=0.95, full=TRUE)[,1],
               confint(model.average, level=0.5, full=TRUE)[,1],
               model.average[["coefficients"]][1,],
               confint(model.average, level=0.5, full=TRUE)[,2],
               confint(model.average, level=0.95, full=TRUE)[,2])
temp2 <- as.data.frame(temp2)
colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")

model.dredge <- as.data.frame(model.dredge)
weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:(ncol(model.dredge)-5))){
  temp <- model.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(model.dredge[,c(2:(ncol(model.dredge)-5))])
colnames(weights) <- "cum_weight"

temp2$merge <- rownames(temp2)
weights$merge <- as.factor(rownames(weights))

temp2 <- data.frame(lapply(temp2, function(x) {
  gsub("locationurban", "location", x)
}))

results <- merge(temp2, weights, by="merge", all=TRUE)
results <- subset(results, merge !="(Intercept)")
results <- subset(results, merge !="food.vol.washed")
results <- results[order(-results$cum_weight),]
results$merge <- factor(results$merge)
results$merge <- forcats::fct_inorder(results$merge)

for(i in c(2:6)){
  results[,i] <- as.numeric(as.character(results[,i]))
}

results$model <- rep("prevalence")

levels(results$merge)[levels(results$merge)=="age.cement"] <- "age"
levels(results$merge)[levels(results$merge)=="age.cement:location"] <- "age * location"
levels(results$merge)[levels(results$merge)=="anth.dig.vol"] <- "dig. anthro"
levels(results$merge)[levels(results$merge)=="anth.dig.vol:location"] <- "* location"
levels(results$merge)[levels(results$merge)=="anth.indig.vol"] <- "indig. anthro"
levels(results$merge)[levels(results$merge)=="sm.vol"] <- "rodents"
levels(results$merge)[levels(results$merge)=="age.cement:sm.vol"] <- "__ * age"
levels(results$merge)[levels(results$merge)=="age.cement:anth.dig.vol"] <- "* age"
levels(results$merge)[levels(results$merge)=="age.cement:anth.dig.vol:location"] <- "* age * location"
levels(results$merge)[levels(results$merge)=="age.cement:anth.indig.vol"] <- "_ * age"
levels(results$merge)[levels(results$merge)=="anth.indig.vol:location"] <- "_ * location"
levels(results$merge)[levels(results$merge)=="location:sm.vol"] <- "__ * location"
levels(results$merge)[levels(results$merge)=="age.cement:location:sm.vol"] <- "__ * age * location"
levels(results$merge)[levels(results$merge)=="age.cement:anth.indig.vol:location"] <- "_ * age * location"

levels(results$merge)[levels(results$merge)=="age.cement:d13C"] <- "___ * age"
levels(results$merge)[levels(results$merge)=="age.cement:d15N"] <- "____ * age"
levels(results$merge)[levels(results$merge)=="d13C:location"] <- "___ * location"
levels(results$merge)[levels(results$merge)=="d15N:location"] <- "____ * location"
levels(results$merge)[levels(results$merge)=="age.cement:d15N:location"] <- "____ * age * location"
levels(results$merge)[levels(results$merge)=="age.cement:d13C:location"] <- "___ * age * location"

results$merge <- factor(results$merge,
                        levels = c("age", "location", "age * location", 
                                   "dig. anthro", "* age", "* location", "* age * location",
                                   "indig. anthro", "_ * age", "_ * location", "_ * age * location",
                                   "rodents", "__ * age", "__ * location", "__ * age * location",
                                   "d13C", "___ * age", "___ * location", "___ * age * location",
                                   "d15N", "____ * age", "____ * location", "____ * age * location"))

results$merge
results$facet <- rep("context")
results$facet[results$merge %in% c("dig. anthro", "* age", "* location", "* age * location")] <- "dig. anthro"
results$facet[results$merge %in% c("indig. anthro", "_ * age", "_ * location", "_ * age * location")] <- "indig. anthro"
results$facet[results$merge %in% c("rodents", "__ * age", "__ * location", "__ * age * location")] <- "rodents"
results$facet[results$merge %in% c("d13C", "___ * age", "___ * location", "___ * age * location")] <- "d13C"
results$facet[results$merge %in% c("d15N", "____ * age", "____ * location", "____ * age * location")] <- "d15N"

results$facet <- factor(results$facet, levels=c("context",
                                                "dig. anthro",
                                                "indig. anthro",
                                                "rodents",
                                                "d13C", "d15N"))

results$merge <- as.character(results$merge)
results$merge[results$merge == "____ * age"] <- "* age"
results$merge[results$merge == "___ * age"] <- "* age"
results$merge[results$merge == "__ * age"] <- "* age"
results$merge[results$merge == "_ * age"] <- "* age"
results$merge[results$merge == "____ * location"] <- "* location"
results$merge[results$merge == "___ * location"] <- "* location"
results$merge[results$merge == "__ * location"] <- "* location"
results$merge[results$merge == "_ * location"] <- "* location"
results$merge[results$merge == "____ * age * location"] <- "* age * location"
results$merge[results$merge == "___ * age * location"] <- "* age * location"
results$merge[results$merge == "__ * age * location"] <- "* age * location"
results$merge[results$merge == "_ * age * location"] <- "* age * location"

results$merge <- factor(results$merge, levels=c("age", "location", "age * location",
                                                "dig. anthro", "indig. anthro", "rodents", "d13C", "d15N",
                                                "* age", "* location", "* age * location"))

# Put together all the interaction plots
library(effects)

model <- glm(qPCR.binary ~ age.cement * location * anth.dig.vol,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

## Diet 1
handpick <- effect('age.cement*location*anth.dig.vol', mod = model,
                   xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                  anth.dig.vol=seq(0, 360, by=0.05)),
                   location = c("urban", "rural"),
                   se=TRUE, confidence.level=.95, typical=mean)
handpick <- as.data.frame(handpick)
handpick$age.cement <- as.factor(as.character(handpick$age.cement))

handpick$age.cement <- factor(handpick$age.cement,
                              levels=c(0.256, 2.47, 4.79))

ggplot(handpick, aes(x=anth.dig.vol, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement)

## Diet 2
model <- glm(qPCR.binary ~ age.cement * location * sm.vol,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

handpick2 <- effect('age.cement*location*sm.vol', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   sm.vol=seq(-0.5, 325, by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick2 <- as.data.frame(handpick2)
handpick2$age.cement <- as.factor(as.character(handpick2$age.cement))

handpick2$age.cement <- factor(handpick2$age.cement,
                               levels=c(0.256, 2.47, 4.79))

ggplot(handpick2, aes(x=sm.vol, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement)

## Diet 3
model <- glm(qPCR.binary ~ age.cement * location * anth.indig.vol,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

handpick3 <- effect('age.cement*location*anth.indig.vol', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   anth.indig.vol=seq(min(sample_data_qPCR$anth.indig.vol),
                                                      max(sample_data_qPCR$anth.indig.vol), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick3 <- as.data.frame(handpick3)
handpick3$age.cement <- as.factor(as.character(handpick3$age.cement))

handpick3$age.cement <- factor(handpick3$age.cement,
                               levels=c(0.256, 2.47, 4.79))

ggplot(handpick3, aes(x=anth.indig.vol, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement)

jtools::interact_plot(model, pred=anth.indig.vol, modx=age.cement, mod2=location,
                      modx.values = c(0.256, 2.47, 4.79)) 



## Diet 4
model <- glm(qPCR.binary ~ age.cement * location * d13C,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

handpick4 <- effect('age.cement*location*d13C', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   d13C=seq(min(sample_data_qPCR$d13C),
                                            max(sample_data_qPCR$d13C), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick4 <- as.data.frame(handpick4)
handpick4$age.cement <- as.factor(as.character(handpick4$age.cement))

handpick4$age.cement <- factor(handpick4$age.cement,
                               levels=c(0.256, 2.47, 4.79))

ggplot(handpick4, aes(x=d13C, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement)

jtools::interact_plot(model, pred=d13C, modx=age.cement, mod2=location,
                      modx.values = c(0.256, 2.47, 4.79)) 

## Diet 4
model <- glm(qPCR.binary ~ age.cement * location * d15N,
             family="binomial"(link="logit"),
             sample_data_qPCR,
             weights = weights.l,
             na.action="na.fail")

handpick5 <- effect('age.cement*location*d15N', mod = model,
                    xlevels = list(age.cement=c(0.256, 2.47, 4.79),
                                   d15N=seq(min(sample_data_qPCR$d15N),
                                            max(sample_data_qPCR$d15N), by=0.05)),
                    location = c("urban", "rural"),
                    se=TRUE, confidence.level=.95, typical=mean)
handpick5 <- as.data.frame(handpick5)
handpick5$age.cement <- as.factor(as.character(handpick5$age.cement))

handpick5$age.cement <- factor(handpick5$age.cement,
                               levels=c(0.256, 2.47, 4.79))

ggplot(handpick5, aes(x=d15N, y=fit)) +
  geom_point() + geom_line() + facet_grid(location ~ age.cement)

jtools::interact_plot(model, pred=d15N, modx=age.cement, mod2=location,
                      modx.values = c(0.256, 2.47, 4.79)) 


handpick$variable <- rep("anth.dig.vol")
handpick2$variable <- rep("sm.vol")
handpick3$variable <- rep("anth.indig.vol")
handpick4$variable <- rep("d13C")
handpick5$variable <- rep("d15N")

scaling <- function(x){(x-min(x))/(max(x)-min(x))}

handpick$xaxis <- scaling(handpick$anth.dig.vol)
handpick2$xaxis <- scaling(handpick2$sm.vol)
handpick3$xaxis <- scaling(handpick3$anth.indig.vol)
handpick4$xaxis <- scaling(handpick4$d13C)
handpick5$xaxis <- scaling(handpick5$d15N)

handpick$type <- rep("stomach contents 1")
handpick2$type <- rep("stomach contents 2")
handpick3$type <- rep("stomach contents 3")
handpick4$type <- rep("stable isotope 1")
handpick5$type <- rep("stable isotope 2")

test1 <- handpick[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test2 <- handpick2[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test3 <- handpick3[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test4 <- handpick4[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test5 <- handpick5[,c("age.cement", "location", "xaxis", "fit", "variable", "type")]
test <- rbind(test1, test2, test3, test4, test5)

test$variable <- factor(test$variable, levels=c("anth.dig.vol", "anth.indig.vol", "sm.vol",
                                                "d13C", "d15N"))

levels(test$age.cement)[levels(test$age.cement)=="0.256"] <- "Younger (0.26 yr)"
levels(test$age.cement)[levels(test$age.cement)=="2.47"] <- "Mean (2.47 yr)"
levels(test$age.cement)[levels(test$age.cement)=="4.79"] <- "Older (4.79 yr)"

test$location <- factor(test$location, levels=c("urban", "rural"))

levels(test$variable)[levels(test$variable)=="anth.dig.vol"] <- "digestible anthropogenic food"
levels(test$variable)[levels(test$variable)=="anth.indig.vol"] <- "indigestible anthropogenic food"
levels(test$variable)[levels(test$variable)=="sm.vol"] <- "rodents"


plot2 <- ggplot(test, aes(x=xaxis, y=fit, color=variable, linetype=type)) +
  geom_line(size=1) + facet_grid(location ~ age.cement) +
  scale_linetype_manual(values=c("dashed", "dashed", "solid","solid","solid"),
                        guide = "none") +
  guides(color = guide_legend(ncol=2, title.position="top", title.hjust=0.5, title = "Dietary predictor")) +
  theme_bw() + plot_theme +
  labs(x="Predictor value (min to max)", y="Probability of infection", tag="a") +
  theme(axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        legend.position="bottom")


plot1 <- ggplot(subset(results, model=="prevalence")) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  geom_linerange(aes(x=merge, y=psd_Coef, ymin=psd_25, ymax=psd_75),
                 size=2, position = position_dodge(width=0.5)) + 
  geom_linerange(aes(x=merge, y=psd_Coef, ymin=psd_2.5, ymax=psd_97.5),
                 size=0.5, position = position_dodge(width=0.5)) + 
  geom_point(aes(x=merge, y=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE,
             position = position_dodge(width=0.5)) +
  #geom_vline(xintercept=c(4.5, 8.5, 12.5, 16.5, 20.5), color="black", size=1.5) +
  scale_x_discrete(limits=rev, position="top") +
  scale_color_manual(values = c("red", "blue"), guide = guide_legend(title.position="top", title = "Dependent variable",
                                                                     title.hjust = 0.5)) +
  theme_bw() + coord_flip() +
  labs(x="Predictor", y="Coefficient", tag="b") +
  facet_grid(facet ~ ., scales="free_y", switch="both") +
  theme(plot.tag=element_text(face="bold"),
        panel.grid.major.x = element_line(linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position="bottom",
        panel.spacing.y=unit(0, "lines"),
        panel.border=element_rect(fill=NA, size=2, color="darkgrey"),
        axis.text.x=element_text(size=9, color="black"),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))


grid.arrange(plot2, plot1, ncol=2, widths=c(0.6,0.4))

comp <- arrangeGrob(plot2, plot1, ncol=2, widths=c(0.6,0.4))

ggsave("~/Documents/Deanna_MS_data/qPCR_status.model.figure.jpg", comp,
       width=8.5, height=6, units="in", dpi=300)


#### Figure A4, S5: Figure 5 reproduced for qPCR-positive coyotes ####
sample_data_qPCR$qPCR.binary <- factor(sample_data_qPCR$qPCR.binary,
                                       levels=c("0","1"))

sample_data_melt <- reshape2::melt(sample_data_qPCR)
sample_data_melt <- subset(sample_data_melt, variable %in% c("sm.vol", "anth.dig.vol", "anth.indig.vol", "d13C", "d15N"))

new_plots <- sample_data_melt %>%
  dplyr::group_by(location, age.class, variable, qPCR.binary) %>%
  dplyr::summarise(avg=mean(value), sd=sd(value), se=plotrix::std.error(value))

new_plots$variable <- factor(new_plots$variable,
                             levels=c("sm.vol", "anth.dig.vol", "anth.indig.vol", "d13C", "d15N"))

levels(new_plots$variable)[levels(new_plots$variable)=="food.vol.washed"] <- "Total food"
levels(new_plots$variable)[levels(new_plots$variable)=="sm.vol"] <- "Rodents"
levels(new_plots$variable)[levels(new_plots$variable)=="anth.dig.vol"] <- "Dig. anthro."
levels(new_plots$variable)[levels(new_plots$variable)=="anth.indig.vol"] <- "Indig. anthro."

new_plots$fill <- paste0(new_plots$age.class, new_plots$qPCR.binary)
new_plots$fill <- factor(new_plots$fill, levels=c("young0", "young1", "old0", "old1"))

sample_data_melt$fill <- paste0(sample_data_melt$age.class, sample_data_melt$qPCR.binary)
sample_data_melt$fill <- factor(sample_data_melt$fill, levels=c("young0", "young1", "old0", "old1"))

p1 <- ggplot(subset(new_plots, location=="urban" & variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
             aes(x=age.class, y=avg, fill=fill)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + plot_theme +
  geom_errorbar(aes(x=age.class, ymin=avg-se, ymax=avg+se),
                position=position_dodge(width=0.9), color="black", width=0.3) +
  facet_wrap(~variable, scales="free_y", ncol=3) +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title="Urban", y="Mean volume (ml)") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p2 <- ggplot(subset(sample_data_melt, location=="urban" & variable %in% c("d13C", "d15N")),
             aes(x=age.class, y=value, fill=fill)) +
  geom_boxplot() +
  theme_bw() + plot_theme +
  facet_wrap(~variable, scales="free_y", ncol=2) +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y="Isotope value") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p3 <- ggplot(subset(new_plots, location=="rural" & variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
             aes(x=age.class, y=avg, fill=fill)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  geom_errorbar(aes(x=age.class, ymin=avg-se, ymax=avg+se),
                position=position_dodge(width=0.9), color="black", width=0.3) +
  facet_wrap(~variable, scales="free_y", ncol=3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title="Rural", y="Mean volume (ml)") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p4 <- ggplot(subset(sample_data_melt, location=="rural" & variable %in% c("d13C", "d15N")),
             aes(x=age.class, y=value, fill=fill)) +
  geom_boxplot() +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("cadetblue2", "salmon", "deepskyblue3", "orangered1")) +
  facet_wrap(~variable, scales="free_y", ncol=2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y="Isotope value") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        legend.position="none")

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p3 <- ggplotGrob(p3)
p4 <- ggplotGrob(p4)

p2$heights <- p1$heights
p4$heights <- p3$heights

grid.arrange(p1, p2, p3, p4,
             ncol=2, nrow=2, widths=c(0.6, 0.4), heights=c(0.475,0.525))

legend_data <- subset(new_plots, age.class=="young")
levels(legend_data$qPCR.binary)[levels(legend_data$qPCR.binary)==0] <- "uninfected"
levels(legend_data$qPCR.binary)[levels(legend_data$qPCR.binary)==1] <- "infected"

legend_data2 <- subset(new_plots, age.class=="old")
levels(legend_data$qPCR.binary)[levels(legend_data$qPCR.binary)==0] <- "uninfected"
levels(legend_data$qPCR.binary)[levels(legend_data$qPCR.binary)==1] <- "infected"

legend1 <- cowplot::get_legend(
  ggplot(subset(legend_data, location=="rural" & 
                  variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
         aes(x=age.class, y=avg, fill=qPCR.binary)) +
    geom_bar(stat="identity", position="dodge") +
    theme_bw() + plot_theme +
    scale_fill_manual(values=c("cadetblue2", "salmon")) +
    guides(fill = guide_legend(title="Juveniles", title.hjust=0.5, title.position="top", ncol=1)) +
    theme(legend.position="bottom")
)
legend2 <- cowplot::get_legend(
  ggplot(subset(legend_data, location=="rural" & 
                  variable %in% c("Rodents", "Dig. anthro.", "Indig. anthro.")),
         aes(x=age.class, y=avg, fill=qPCR.binary)) +
    geom_bar(stat="identity", position="dodge") +
    theme_bw() + plot_theme +
    scale_fill_manual(values=c("deepskyblue3", "orangered1")) +
    guides(fill = guide_legend(title="Adults", title.hjust=0.5, title.position="top", ncol=1)) +
    theme(legend.position="bottom")
)

grid.arrange(p1, p2, p3, p4, legend1, legend2,
             layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2),
                                 c(3,3,3,3,3,3,4,4,4,4),
                                 c(NA,NA,NA,5,5,6,6,NA,NA,NA)),
             heights=c(0.4,0.4,0.2))

comp <- arrangeGrob(p1, p2, p3, p4, legend1, legend2,
                    layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2),
                                        c(3,3,3,3,3,3,4,4,4,4),
                                        c(NA,NA,NA,5,5,6,6,NA,NA,NA)),
                    heights=c(0.4,0.4,0.2))

ggsave("~/Documents/Deanna_MS_data/diet_bars_for_qPCR_infections.jpg", comp, width=6.5, height=5, units="in", dpi=300)

#### qPCR positive 


