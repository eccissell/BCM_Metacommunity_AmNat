# BCM_Metacommunity_AmNat

# CITATION
AUTHOR 1 &amp; AUTHOR 2. (IN REVIEW) Heterogeneous predation, community asynchrony, and metacommunity stability in cyanobacterial mats

# AUTHOR NAMES AND CONTACT DETAILS
Redacted for ensuring double-blind peer review.

# SUMMARY OF STUDY
Benthic cyanobacterial mats are increasing in abundance on coral reefs, yet the processes governing their persistence across scales of ecological organization remain poorly resolved. Here, we integrated observational and experimental approaches to characterize the ecological processes driving the dynamics of these mats. Community- and metacommunity-scale dynamics of mats off the island of Bonaire, Caribbean Netherlands, were tracked for 49 days alongside quantification of predation pressure from fishes. We also tested the hypothesis that elevated predation would decrease mat persistence using an in-situ coring experiment. Cyanobacterial mat metacommunity-scale cover was temporally stable across the study despite high community-level mortality, likely stabilized by asynchrony in community extinction. Natural mat communities had a median persistence time of 36 days, with experimentally disturbed mats exhibiting significantly longer persistence. Predation rates from diverse reef fishes significantly varied among spatially-distinct cyanobacterial mat communities within a reef site, but not among reef sites.

# RESPONSIBILITY FOR DATA COLLECTION AND CODE
Author 1

# INCLUDED FILES &amp; FOLDERS  
## metacommunity_BCM_benthic_df.csv
**Description:**
Data from photoquadrat-based estimates of benthic cyanobacterial mat (BCM) cover used to quantify metacommunity-scale dynamics.

**Variables:**
- **Name** — Unique identifier for each photoquadrat
- **Date** — Sampling date of the photoquadrat
- **Depth** — Depth of the photoquadrat (in m)
- **BCM** — Proportional benthic cover of cyanobacterial mats in that photoquadrat (ranging from a value of 0–1), estimated from 50-point CoralNet annotations

## coring_experiment_df.csv
**Description:**
Community-scale (e.g., Individual cyanobacterial mats) tracking data, including from the in situ coring (simulated disturbance) experiment, used to analyze community-scale (mat) survival and persistence.

**Variables:**
- **Treatment** — Experimental treatment (W = Unmanipulated; H = Wounded)
- **Set** — Unique Identifier for each mat
- **First_Date** — Date monitoring of the mat began
- **Last_Date** — Date monitoring of the mat ended
- **Duration** — Number of days the mat was observed
- **Died** — Binary indicator of mat death during the study 
- **Binom** — Binary response variable used in GLM analyses (0 = died, 1 = survived)

## mat_video_follow_counts.csv
**Description:**
Raw data from video-based observations of reef fish interactions with cyanobacterial mats, used to quantify predation rates from reef fishes on cyanobacterial mats.

**Variables:**
- **mat** — Unique identifier for the focal mat community observed in the video
- **site** — Reef site name where the mat was observed
- **date** — Observation date
- **depth** — Depth of the focal mat (m)
- **min** — Minutes component of observation duration
- **sec** — Seconds component of observation duration
- **duration** — Total observation duration (in seconds)
- **bicolor_ind** — Number of *Stegastes partitus* individuals observed near/associated with the mat during the observation
- **bicolor** — Total bite count for bicolor damselfish [*Stegastes partitus*] (# of bites)
- **stripe** — Total bite count for striped parrotfish [*Scarus iseri*], all phase combined (# of bites)
- **tp_stripe** — Total bite count for terminal-phase striped parrotfish [*Scarus iseri*] (# of bites)
- **ip_stripe** — Total bite count for initial-phase striped parrotfish [*Scarus iseri*] (# of bites)
- **princess** — Total bite count for princess parrotfish [*Scarus taeniopterus*], all phases combined (# of bites)
- **tp_princess** — Total bite count for terminal-phase princess parrotfish [*Scarus taeniopterus*] (# of bites)
- **ip_princess** — Total bite count for initial-phase princess parrotfish [*Scarus taeniopterus*] (# of bites)
- **redband** — Total bite count for redband parrotfish [*Sparisoma aurofrenatum*], all phases combined (# of bites)
- **tp_redband** — Total bite count for terminal-phase redband parrotfish [*Sparisoma aurofrenatum*] (# of bites)
- **ip_redband** — Total bite count for intial-phase redband parrotfish [*Sparisoma aurofrenatum*] (# of bites)
- **queen** — Total bite count for queen parrotfish [*Scarus vetula*] (# of bites)
- **surgeon** — Total bite count for ocean surgeonfish [*Acanthurus bahianus*] (# of bites)
- **tang** — Total bite count for Atlantic blue tang [*Acanthurus coeruleus*] (# of bites)
- **rock** — Total bite count for rock beauty angelfish [*Holacanthus tricolor*] (# of bites)
- **french** — Total bite count for French angelfish [*Pomacanthus paru*] (# of bites)
- **yellow_goat** — Total foraging disturbance count for yellow goatfish [*Mulloidichthys martinicus*], treated as bites in analyses (# of disturbances)
- **spotted_goat** — Total foraging disturbance count for spotted goatfish [*Pseudupeneus maculatus*], treated as bites in analyses (# of disturbances)
- **r_bi** — Bite rate for bicolor damselfish [*Stegastes partitus*] (bites·min⁻¹)
- **r_yel** — Foraging disturbance rate for yellow goatfish [*Mulloidichthys martinicus] (disturbances·min⁻¹)
- **r_tstripe** — Bite rate for terminal-phase striped parrotfish [*Scarus iseri*] (bites·min⁻¹)
- **r_iprince** — Bite rate for initial-phase Scarus taeniopterus [Princess parrotfish] (bites·min⁻¹)
- **r_qu** — Bite rate for queen parrotfish [*Scarus vetula*] (bites·min⁻¹)
- **r_surgeon** — Bite rate for ocean surgeonfish [*Acanthurus bahianus*] (bites·min⁻¹)
- **r_tang** — Bite rate for Atlantic blue tang [*Acanthurus coeruleus*] (bites·min⁻¹)
- **r_spot** — Foraging disturbance rate for spotted goatfish [*Pseudupeneus maculatus*] (disturbances·min⁻¹)
- **r_ired** — Bite rate for initial-phase redband parrotfish [*Sparisoma aurofrenatum*] (bites·min⁻¹)
- **r_rock** — Bite rate for rock beauty angelfish [*Holacanthus tricolor*] (bites·min⁻¹)
- **r_french** — Bite rate for French angelfish [*Pomacanthus paru*] (bites·min⁻¹)
- **r_tprince** — Bite rate for terminal-phase princess parrotfish [*Scarus taeniopterus*] (bites·min⁻¹)
- **r_tred** — Bite rate for terminal-phase redband parrotfish [*Sparisoma aurofrenatum*] (bites·min⁻¹)

## Combined_Master_Script.R
**Description:**
R script used to conduct all analyses included in the manuscript.

- **R version:** 4.5.2

**Attached base packages:**
- grid
- stats
- graphics
- grDevices
- utils
- datasets
- base
- methods

**Included packages:**
- car 3.1-3
- carData 3.0-5
- mgcViz 0.2.1
- qgam 2.0.0
- mgcv 1.9-3
- plotly 4.11.0
- fishualize 0.2.3
- Rmisc 1.5.1
- lattice 0.22-7
- MuMIn 1.48.11
- DHARMa 0.4.7
- glmmTMB 1.1.14
- lubridate 1.9.4
- forcats 1.0.1
- dplyr 1.1.4
- purrr 1.2.0
- readr 2.1.6
- tidyr 1.3.1
- tibble 3.3.0
- tidyverse 2.0.0
- nlme 3.1-168
- gridExtra 2.3
- lmerTest 3.1-3
- lme4 1.1-38
- Matrix 1.7-4
- ggthemes 5.2.0
- scales 1.4.0
- stringr 1.6.0
- munsell 0.5.1
- stringi 1.8.7
- magrittr 2.0.4
- digest 0.6.39
- Rcpp 1.1.1
- labeling 0.4.3
- gtable 0.3.6
- colorspace 2.1-2
- reshape2 1.4.5
- gtools 3.9.5
- gdata 3.0.1
- readxl 1.4.5
- ggpubr 0.6.2
- ggplot2 4.0.1
- ggbreak 0.1.6
- plyr 1.8.9

## Representative Videos and Photos
- **Description:** Included in this directory are raw representative photographs (as .jpg files; file naming convention: "GOPROXXX.JPG") of photoquadrats used for creating "metacommunity_BCM_benthic_df.csv" and raw representative videos (as .mp4 files; file naming convention: "SITE_DATE_MATID.mp4"). Note that video quality was reduced using HandBrake v 1.10.2 to allow for upload. Full datasets are available upon request (contact ethan.c.cissell@gmail.com)

