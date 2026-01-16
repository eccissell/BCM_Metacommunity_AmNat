#CODE AUTHOR: ETHAN C. CISSELL
#DATE LAST CODE UPDATED: January 15, 2026
#####################################################
#Code for: Cissell EC, McCoy SJ. In Review @ The American Naturalist. Heterogeneous predation, community asynchrony, and metacommunity stability in cyanobacterial mats (Preprint Available Here: https://www.biorxiv.org/content/10.1101/2022.10.07.511315v2)
#####################################################
## ============================
## USER SETUP
## ============================

## Set this to the folder containing the data file
data_dir <- "/LOCAL/PATH/TO/DATA"
# NOTE: Users must update `data_dir` to point to the location of the data file on their system

#Load libraries
# Vector of required packages
packages <- c(
  "plyr", "ggbreak", "ggpubr", "readxl", "gdata", "gtools",
  "ggplot2", "reshape2", "colorspace", "gtable", "labeling",
  "Rcpp", "digest", "magrittr", "stringi", "munsell", "stringr",
  "scales", "grid", "ggthemes", "lme4", "lmerTest", "gridExtra",
  "nlme", "tidyverse", "glmmTMB", "DHARMa", "MuMIn", "Rmisc",
  "fishualize", "plotly", "mgcv", "mgcViz", "rgl", "plot3D",
  "plotrix", "dplyr", "matrixStats", "car", "survival"
)

# Install any packages that are not already installed
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))


########
##Metacommunity benthic cover
########
#Read in data
benthos=read.csv(file.path(data_dir,'metacommunity_BCM_benthic_df.csv'))

#Cleaning
benthos <- benthos %>%
  mutate(
    Date = as.factor(Date),      # keep as factor for ANOVA-style contrasts
    Depth = as.numeric(Depth),
    BCM = as.numeric(BCM)
  )

##Focus on spread of depth over time
#Calc CI
DEPTH<- summarySE(benthos, groupvars = c("Date"), measurevar = c("Depth"))

#Now Pass that to plot of means
DEPTH %>% ggplot()+
  geom_bar(aes(x = Date, y = Depth), stat= "summary", fun="mean", fill="#393E42") +
  theme_classic() +
  ylab("Depth") +
  scale_fill_brewer(palette="Paired", name = "Substrate") +
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black")) +
  theme(legend.text=element_text(size=16), legend.title = element_text(size = 18)) +
  geom_errorbar(aes(x= Date, y= Depth, ymin=Depth-ci, ymax=Depth+ci),
                width=0.09,
                color="#fe5d26",
                size=1)
#Boxplot for quick look at spread
boxplot(Depth~Date,las=1, data=benthos)

#Model to look at BCM across date and depth
lmbcmmet=lm(BCM~Date*Depth,data=benthos)
summary(lmbcmmet)
#Check Assumptions
library("DHARMa")
check_bcmmet_model1 <- simulateResiduals(fittedModel = lmbcmmet, n = 500)
plot(check_bcmmet_model1)
#Good
#Significance
library("car")
Anova(lmbcmmet,type="II",test.statistic = "F")

#Run model without interaction because interactions isn't significant, then compare with LRT
lmbcmmet2=lm(BCM~Date+Depth,data=benthos)
summary(lmbcmmet2)
#Check Assumptions
check_bcmmet_model2 <- simulateResiduals(fittedModel = lmbcmmet2, n = 500)
plot(check_bcmmet_model2)
#Good
#Now test if adding complexity is good using LRT
anova(lmbcmmet2,lmbcmmet,test="LRT") #p=0.269
#Use reduced model
#Significance of reduced model (no interaction)
Anova(lmbcmmet2,type="II",test.statistic = "F")
#Now use Tukey HSD
#Now use pairwise Wilcox with corrections to look at specific differences
pairwise.wilcox.test(benthos$BCM, benthos$Date, p.adjust.method="bonferroni")
#Only 6/11 vs. 6/13, 6/20 and 6/26 differ
#Summarize and Organize Benthic Data by Date
BCM<- summarySE(benthos, groupvars = c("Date"), measurevar = c("BCM"))
colnames(BCM)[3] <- "Mean"
BCM$Group = "BCM"
summary(BCM)

#Plot cover of BCM against date
bcm_date=BCM %>%
  ggplot()+
  geom_bar(aes(x = Date, y = Mean), stat= "identity", fill="#393E42",size=1.5,color="black") +
  theme_classic() +
  ylab("Proportional Cover of BCM") +
  scale_fill_brewer(palette="Paired", name = "Substrate") +
  theme(axis.ticks.x=element_blank()) +
  theme(axis.title = element_text(size=16),
        axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=1.1)) +
  theme(legend.text=element_text(size=16), legend.title = element_text(size = 18)) +
  geom_errorbar(aes(x= Date, y= Mean, ymin=Mean-ci, ymax=Mean+ci),
                width=0.09,
                color="black",
                size=1) +
  geom_hline(yintercept=0.1956,
             linetype="dashed",
             color="black",
             size=1.07)+
  xlab("")
bcm_date

#Get mean and SD of BCM cover across dates
mean(BCM$Mean) #0.196
sd(BCM$Mean) #0.0494

#Plot BCM cover vs. depth
benthos%>%
  ggplot()+
  geom_smooth(aes(x = Depth, y = BCM),method=lm,se=T,color="black")+
  geom_point(aes(x = Depth, y = BCM), stat= "identity") +
  theme_classic() +
  ylab("Proportional Cover of BCM") +
  xlab("Depth (m)") +
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=13, color = "black"),
        axis.text.y = element_text(size=13, color = "black"),
        axis.line.y = element_line(size=1.1),
        axis.line.x = element_line(size=1.1)) +
  theme(legend.text=element_text(size=13), legend.title = element_text(size = 18))

#Model BCM cover against depth
depthlm = lm (BCM~Depth, data = benthos)
summary(depthlm)

########
##Community persistence and coring experiment
########
#Looking at overall mat deaths
df=read.csv(file.path(data_dir,'coring_experiment_df.csv'))

#Sort by death date
df=df[order(df$Last_Date),]

#Lets plot the mean + SE (Fig. 2a)
p=ggplot(df,aes(x=Treatment,y=Binom)) +
  geom_point(alpha=0.8) +
  geom_jitter(width=0.2,height=0.07) +
  theme_classic()
p= p + stat_summary(fun.data="mean_se",fun.args = list(mult=1),
                    geom="pointrange",color="black",size=1) +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1)) +
  scale_x_discrete(labels=c("Wounded","Unmanipulated"))+
  labs(y="Probability of Survival")
p

#Cumulative sum plot of proportion dead by date
#Need to change 0 and 1 assignment for this plot that counts cumulative dead
df_zerooneswap=df %>%
  subset(Treatment=="W")%>%
  mutate(binom_swap = ifelse(Binom==1,0,1))
#Order
#Sort by death date
df_zerooneswap=df_zerooneswap[order(df_zerooneswap$Last_Date),]
#Cumulative sum plot of proportion dead by date (Fig 1c)
death_cumulsum=ggplot(data=df_zerooneswap, aes(x=Last_Date)) +
  geom_col(aes(y=binom_swap))+
  geom_line(aes(x=Last_Date,y=(cumsum(binom_swap)/26*100),group=1), inherit.aes=FALSE, size=1.1)+
  theme(panel.background = element_blank()) +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt"))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  scale_y_continuous(name="Cumulative % dead")+
  scale_x_discrete(name="Date")

#Call plot
death_cumulsum

#14 died in unmanipulated treatment across 10 different death dates
#53.8% mortality (14/26)

#Combine this plot with metacommunity persistence plot from above for Fig 1b,c
bcm_ar=ggarrange(bcm_date, NULL, death_cumulsum,
                 ncol=1, nrow=3, align="h",
                 widths=c(1,-0.1,1),
                 heights=c(1,-0.18,1.1),
                 labels=c("a","","b"))

#Use Kruskal-Wallis test to initial test for differences between treatments in start date
kruskal.test(First_Date~Treatment,data=df) #p=0.8821 ; no significant difference
#Use GLM for differences between treatments
trtglm=glm(Binom~Treatment,data=df,family=binomial)
summary(trtglm)
#Check Assumptions
check_trt_model1 <- simulateResiduals(fittedModel = trtglm, n = 500)
plot(check_trt_model1)
#Looks good

#Look for effect of start date
trtglm2=glm(Binom~Treatment*First_Date,data=df,family=binomial)
summary(trtglm2)

#Check Assumptions
check_trt_model2 <- simulateResiduals(fittedModel = trtglm2, n = 500)
plot(check_trt_model2)
#Looks good, and no significant effects

#Finally test if adding complexity is good using LRT
anova(trtglm,trtglm2,test="LRT") #p=0.02 ; Marginally better use full model

#Now test for significance
Anova(trtglm2,type="II",test.statistic = "F")

#Compute Kaplan-Meier survivorship curves to compare persistence across treatments
#Read in df
df2 <- read.csv(file.path(data_dir,'coring_experiment_df.csv'))

#Clean up df and prepare 'status' for KM analysis
df2$Treatment <- factor(df2$Treatment, levels = c("H","W"))
df2$status    <- as.integer(toupper(df2$Died) == "YES")      # 1 = died, 0 = censored

# Run Kaplan–Meier by treatment
library("survival")
kmfit <- survfit(Surv(Duration, status) ~ Treatment, data = df2)
summary(kmfit)

#Log-rank p-value 
lr   <- survdiff(Surv(Duration, status) ~ Treatment, data = df2)
pval <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
print(lr) #p = 0.07

# Extract medians for median persistence
km_medians <- quantile(kmfit, probs = 0.5)
print(km_medians)
#Cored median = 47
#Uncored median = 36

# Build a df with SE for plotting
strata_labels <- rep(names(kmfit$strata), kmfit$strata)
km_df <- tibble(
  time  = kmfit$time,
  surv  = kmfit$surv,
  se    = kmfit$std.err,                            
  Treatment = factor(sub("^Treatment=", "", strata_labels), levels = c("H","W"))
) |>
  mutate(lower = pmax(0, surv - se),
         upper = pmin(1, surv + se)) |>
  group_by(Treatment) |>
  arrange(time, .by_group = TRUE) |>
  reframe(
    time  = c(0, time),
    surv  = c(1, surv),
    lower = c(1, lower),
    upper = c(1, upper)
  ) |>
  ungroup()

# Plot KM survivor curves (Figure 2b)
decay_km_plot <-
  ggplot() +
  geom_ribbon(
    data = filter(km_df, Treatment == "W"),
    aes(x = time, ymin = lower, ymax = upper),
    fill = "#455765", alpha = 0.20
  ) +
  geom_ribbon(
    data = filter(km_df, Treatment == "H"),
    aes(x = time, ymin = lower, ymax = upper),
    fill = "#cd9fb2", alpha = 0.20
  ) +
  geom_step(
    data = km_df,
    aes(x = time, y = surv, color = Treatment),
    linewidth = 1.75
  ) +
  scale_color_manual(
    values = c(H = "#cd9fb2", W = "#455765"),
    labels = c(H = "Coring (H)", W = "Unmanipulated (W)"),
    name   = "Treatment"
  ) +
  theme(panel.background = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        axis.line.y = element_line(size = 1.1),
        axis.line.x = element_line(size = 1.1)) +
  ylab("Survival Probability") +
  xlab("Time (Days)") +
  theme(legend.position = "none")

print(decay_km_plot)

#Test for difference in survivability at median survival for unmanipulated treatment (i.e., natural expected survivability)
#Extract median
kmW  <- survfit(Surv(Duration, status) ~ 1, data = subset(df2, Treatment == "W"))
tau  <- as.numeric(summary(kmW)$table["median"])
if (is.na(tau)) stop("W median not reached; pick a prespecified tau within follow-up.")

#KM by treatment & survival at tau
km   <- survfit(Surv(Duration, status) ~ Treatment, data = df2)
sm   <- summary(km, times = tau)
labs <- sm$strata
iH   <- grep("(^H$)|(^Treatment=H$)", labs)
iW   <- grep("(^W$)|(^Treatment=W$)", labs)
SH <- sm$surv[iH]; SEH <- sm$std.err[iH]
SW <- sm$surv[iW]; SEW <- sm$std.err[iW]

#Calculate difference, SE, z, and p-values
diff   <- SH - SW
SEdiff <- sqrt(SEH^2 + SEW^2)
z      <- diff / SEdiff
p_two  <- 2 * pnorm(-abs(z))    # two-sided
p_one_HgtW <- 1 - pnorm(z)      # one-sided, alternative: H has higher survival at tau

cat(sprintf("Difference S_H - S_W = %.3f (SE=%.3f)\n", diff, SEdiff)) 
#Difference of 0.296 SE=0.142
cat(sprintf("Two-sided p = %.3f; One-sided p (H > W) = %.3f\n", p_two, p_one_HgtW))
#Two-sided p = 0.037; one-sided p = 0.018

#Run a permutational p-check to test rigor of above (relaxes assumptions)
set.seed(10000)
Sdiff_at_tau <- function(dat, tau) {
  sm <- summary(survfit(Surv(Duration, status) ~ Treatment, data = dat),
                times = tau, extend = TRUE)
  tab <- data.frame(strata = sm$strata, S = sm$surv, check.names = FALSE)
  SH <- tab$S[tab$strata %in% c("Treatment=H","H")]
  SW <- tab$S[tab$strata %in% c("Treatment=W","W")]
  SH - SW
}

obs_diff <- Sdiff_at_tau(df2, tau)

B <- 10000
perm_diffs <- replicate(B, {
  df2$Treatment <- sample(df2$Treatment)
  Sdiff_at_tau(df2, tau)
})

p_perm_two  <- mean(abs(perm_diffs) >= abs(obs_diff))
p_perm_oneH <- mean(perm_diffs >= obs_diff)  # one-sided H > W

cat(sprintf("Permutation p (two-sided) at tau=%.0f d: %.3f\n", tau, p_perm_two)) #p = 0.045
cat(sprintf("Permutation p (one-sided H>W): %.3f\n", p_perm_oneH)) #p = 0.023
#Permutational p's confirms confidence in statistical significance


########
##Fish predation pressure
########
#Read in data
matfollow=read.csv(file.path(data_dir,'mat_video_follow_counts.csv'))

## Species columns used in gather
species_cols <- c(
  "bicolor", "yellow_goat", "tp_stripe", "ip_princess",
  "ip_stripe", "queen", "surgeon", "tang", "spotted_goat",
  "ip_redband", "rock", "french", "tp_princess", "tp_redband"
)

## Add total.bites per mat (row-wise sum)
matfollow <- matfollow %>%
  mutate(
    total.bites = rowSums(across(all_of(species_cols)), na.rm = TRUE)
  )

## Convert to long format
gathered <- matfollow %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "species",
    values_to = "freq"
  )

gathered$species <- as.factor(gathered$species)
str(gathered)

#Fit model
m1=glmmTMB(data=gathered, freq ~ site+species + (1|site/mat) + offset(log(duration)),family=nbinom2,REML=TRUE)
summary(m1)
Anova(m1,type="II")
#Check model assumptions
simulationOutput <- simulateResiduals(m1)
testDispersion(simulationOutput)
plot(simulationOutput)

#Obtain SD of total bites among mats to better understand among mat variability suggested in RE
#All mats
sd(gathered$total.bites) #108.8 bites SD
mean(gathered$total.bites) #72.8 bites Mean

#Plot of bite count by species for Fig 3a
g=gathered %>%
  mutate(species=species%>%
           fct_relevel("bicolor",
                       "tp_stripe",
                       "ip_stripe",
                       "tp_princess",
                       "ip_princess",
                       "tp_redband",
                       "ip_redband",
                       "queen",
                       "surgeon",
                       "tang",
                       "rock",
                       "french",
                       "yellow_goat",
                       "spotted_goat"))
fish_species_bites=ggplot(data=g,aes(x=freq, y=as.factor(species)))+
  labs(y="Species") +
  labs(x="Bite Count") +
  geom_boxplot(notch=FALSE, fill='#A4A4A4',outlier.shape=NA)+
  geom_jitter(data = g, height=0.1, width=0, aes(alpha=8/10),size=3)+
  scale_y_discrete(labels=c("St. partitus", "TP Sc. iseri","IP Sc. iseri","TP Sc. taeniopterus","IP Sc. taeniopterus","TP Sp. aurofrenatum","IP Sp. aurofrenatum","Sc. vetula","A. bahianus","A. coeruleus","H. tricolor","Po. paru","M. martinicus","Ps. maculatus")) +
  theme(panel.grid.major = element_blank (), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.line = element_line(size=1, colour = "black"),
        axis.line.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  guides (size=FALSE) + guides(alpha=FALSE) +
  theme(axis.text.y=element_text(face="italic")) +
  scale_x_break(c(150,220))
fish_species_bites

#Summarize to obtain rate by species
gathered2=gather(matfollow,rates,rate,r_bi,r_yel,r_tstripe,r_iprince,r_qu,r_surgeon,r_tang,r_spot,r_ired,r_rock,r_french,r_tprince,r_tred)

#Pull Mean and SD
mean(gathered2$rate) #0.19
sd(gathered2$rate) #0.70

#Invert axes for consistent plotting 
# Explicitly set site order
matfollow$site <- factor(
  matfollow$site,
  levels = rev(c(
    "Bachelor",
    "The Lake",
    "Angel City",
    "Aquarius",
    "Invisibles"
  ))
)

#Figure 3b
site_bites_with=ggplot(data=matfollow, aes(x=total.bites, y=site)) +
  geom_boxplot(aes(x=total.bites,y=site, group=site),
               fill="#A4A4A4",
               size=0.75,
               color="black") +
  theme_classic() +
  theme(axis.text = element_text(size=15,color="black"),
        axis.title=element_text(size=15,color="black"),
        axis.line = element_line(size=1.1),
        axis.ticks.length=unit(4,"pt")) +
  ylab("Study Site") +
  xlab("Total Bites (#)")
site_bites_with

##END OF SCRIPT
