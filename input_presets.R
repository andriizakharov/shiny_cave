library(dplyr)

### Rast & Hofer 2013, Table 5
data_all <- read.csv("data/Hertzog-Table5.csv", sep = ";", dec = ",")
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

acad_A_CD <- data_all %>% filter(Study == "ACAD", y == "A", x == "CD")

elsa_AF_PM <- data_all %>% filter(Study == "ELSA", y == "AF", x == "PM")


### Add custom parameter settings from Lindenberger, von Oertzen, Ghisletta, Hertzog 2011
LOGH_2011_A <- list(sigma2_Iy = 100, sigma2_Ix = 100, sigma2_Sy = 200, sigma2_Sx = 200,
                    sigma2_Ex = 25, sigma2_Ey = 25, sigma_IyIx = -3, sigma_IySy = 0,
                    sigma_IySx = 0, sigma_SyIx = 0, sigma_IxSx = 0, mean_slope_x = -14,
                    mean_slope_y = -5)
LOGH_2011_B <- list(sigma2_Iy = 100, sigma2_Ix = 100, sigma2_Sy = 200, sigma2_Sx = 200,
                    sigma2_Ex = 25, sigma2_Ey = 25, sigma_IyIx = -3, sigma_IySy = 0,
                    sigma_IySx = 0, sigma_SyIx = 0, sigma_IxSx = 0, mean_slope_x = -14,
                    mean_slope_y = -2.5)
LOGH_2011_C <- list(sigma2_Iy = 100, sigma2_Ix = 100, sigma2_Sy = 200, sigma2_Sx = 200,
                    sigma2_Ex = 25, sigma2_Ey = 25, sigma_IyIx = 60, sigma_IySy = 0,
                    sigma_IySx = 0, sigma_SyIx = 0, sigma_IxSx = 0, mean_slope_x = -14,
                    mean_slope_y = -5)