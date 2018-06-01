library(dplyr)

### Rast & Hofer 2013, Table 5
data_all <- read.csv("Hertzog-Table5.csv", sep = ";", dec = ",")
names(data_all) <- c("Study", "y", "x", "N", "sigma2_Iy", "sigma2_Sy",
                     "sigma2_Ix", "sigma2_Sx", "sigma_IySy", "sigma_IyIx",
                     "sigma_IySx", "sigma_SyIx", "sigma_SySx", "sigma_IxSx",
                     "sigma2_Ey", "sigma2_Ex", "sigma_EyEx", "Waves", "Length")

### One variant of one study: Victoria Longitudinal 
### with x = Word Recall, y = Reaction Time
vls_RT_WRC <- data_all %>% filter(Study == "VLS", y == "RT", x == "WRC")
vls_RT_WRC_timepoints <- c(3.06, 6.08, 9.5)
vls_RT_WRC_age_at_start <- 55:85
vls_RT_WRC_x_axis <- 55:(max(vls_RT_WRC_age_at_start) + 
    ceiling(max(vls_RT_WRC_timepoints)))
vls_RT_WRC_n <- 50