library(semper)

bivariateLgcm<-function(lgcm.x, lgcm.y, slope.x.slope.y.covariance=0, icept.x.slope.y.covariance=0,
                        icept.y.slope.x.covariance =0, icept.x.icept.y.covariance=0,
                        latent.correlations=FALSE)
{
    result <- list()
    result$lgcm.x <- lgcm.x
    result$lgcm.y <- lgcm.y
    result$slope.x.slope.y.covariance <- slope.x.slope.y.covariance
    result$icept.x.slope.y.covariance <- icept.x.slope.y.covariance
    result$icept.y.slope.x.covariance <- icept.y.slope.x.covariance
    result$icept.x.icept.y.covariance <- icept.x.icept.y.covariance
    
    result$latent.correlations <- latent.correlations
    class(result) <- "bivariateLgcm"
    
    return(result)
}

toOpenMx <- function(generic.model, name=NULL)
{
    manifest.varname <- "X"
    if (inherits(generic.model,"lgcm")) {
        
        #if (latent.correlations) stop("Not implemented yet.")
        
        lgcm <- generic.model
        if (is.null(name)) {
            name <- "Latent Growth Curve Model"
        }
        manifests <- paste(manifest.varname, 1:length(lgcm$timepoints),sep="")
        latents <- c("intercept","slope")
        
        p1 <- mxPath(from=latents[1], to=manifests, free=FALSE, values=1)
        p2 <- mxPath(from=latents[2], to=manifests, free=FALSE, values= lgcm$timepoints)
        
        p3 <- mxPath(from=manifests, to=manifests, free = TRUE, values=lgcm$residual.variance,
                     labels="residualvariance", arrows=2)
        
        p4 <- mxPath(from=c(latents[1],latents[1],latents[2]),
                     to = c(latents[1],latents[2],latents[2]), free = TRUE, arrows=2,values=
                         c(lgcm$intercept.variance, lgcm$intercept.slope.covariance,
                           lgcm$slope.variance),connect = "unique.pairs",
                     labels=c("interceptvariance","interceptslopecovariance","slopevariance"))
        
        p5 <- mxPath(from="one",to=manifests,values=0,free=FALSE)
        
        p6 <- mxPath(from="one",to=latents,values=c(lgcm$intercept.mean, lgcm$slope.mean),free=TRUE,
                     labels=c("interceptmean","slopemean"))
        
        lgcmModel <- mxModel(name,
                             type="RAM",
                             manifestVars=manifests,
                             latentVars = latents,
                             p1,p2,p3,p4,p5,p6)
        
        return(lgcmModel)
        
        
    } else if (inherits(generic.model,"qgcm")) 
    {
        
        #if (latent.correlations) stop("Not implemented yet.")
        
        qgcm <- generic.model
        if (is.null(name)) {
            name <- "Quadratic Latent Growth Curve Model"
        }
        manifests <- paste(manifest.varname, 1:length(qgcm$timepoints),sep="")
        latents <- c("intercept","slope","quadratic")
        
        p1 <- mxPath(from=latents[1], to=manifests, free=FALSE, values=1)
        p2 <- mxPath(from=latents[2], to=manifests, free=FALSE, values= qgcm$timepoints)
        
        p7 <- mxPath(from=latents[3], to=manifests, free=FALSE, values= qgcm$timepoints^2)
        
        p3 <- mxPath(from=manifests, to=manifests, free = TRUE, values=qgcm$residual.variance,
                     labels="residualvariance", arrows=2)
        
        p4 <- mxPath(from=c(latents[1],latents[1],latents[2],latents[1],latents[2],latents[3]),
                     to = c(latents[1],latents[2],latents[2],latents[3],latents[3],latents[3]), free = TRUE, arrows=2,values=
                         c(qgcm$intercept.variance, qgcm$intercept.slope.covariance,
                           qgcm$slope.variance, qgcm$intercept.quadratic.covariance,
                           qgcm$slope.quadratic.covariance, qgcm$quadratic.variance),connect = "single",
                     labels=c("interceptvariance","interceptslopecovariance","slopevariance",
                              "interceptquadraticcovariance","slopequadraticcovariance","quadraticvariance"))
        
        p5 <- mxPath(from="one",to=manifests,values=0,free=FALSE)
        
        p6 <- mxPath(from="one",to=latents,values=
                         c(qgcm$intercept.mean, qgcm$slope.mean,qgcm$quadratic.mean),free=TRUE,
                     labels=c("interceptmean","slopemean","quadraticmean"))
        
        qgcmModel <- mxModel(name,
                             type="RAM",
                             manifestVars=manifests,
                             latentVars = latents,
                             p1,p2,p3,p4,p5,p6,p7)
        
        return(qgcmModel)
        
        
    } else if (inherits(generic.model,"bivariateLgcm")) 
    {
        
        if (generic.model$latent.correlations) {
            
            lgcm.x <- generic.model$lgcm.x
            lgcm.y <- generic.model$lgcm.y
            
            manifestsx <- paste("X", 1:length(lgcm.x$timepoints),sep="")
            manifestsy <- paste("Y", 1:length(lgcm.y$timepoints),sep="")
            manifests <- c(manifestsx, manifestsy)
            
            latentsx <- c("interceptx","slopex")
            latentsy <- c("intercepty","slopey")
            latents <- c(latentsx, latentsy)
            
            latentsx.std <- c("interceptx.std","slopex.std")
            latentsy.std <- c("intercepty.std","slopey.std")
            latents <- c(latents, latentsx.std, latentsy.std)
            
            
            p1x <- mxPath(from=latentsx[1], to=manifestsx, free=FALSE, values=1)
            p2x <- mxPath(from=latentsx[2], to=manifestsx, free=FALSE, values= lgcm.x$timepoints)
            
            p3x <- mxPath(from=manifestsx, to=manifestsx, free = TRUE, values=lgcm.x$residual.variance,
                          labels="residualvariancex", arrows=2)
            
            # fixed variances
            p4x <- mxPath(from=c(latentsx.std[1],latentsx.std[2]), to=c(latentsx.std[1],latentsx.std[2]), free=FALSE,
                          arrows=2,values=1,connect="single")
            # projections to std latents
            p4xx <- mxPath(from=latentsx.std,to=latentsx, free=TRUE,arrows=1,
                           values=c(sqrt(lgcm.x$intercept.variance),
                                    sqrt(lgcm.x$slope.variance)),
                           connect="single")
            
            #   p4x <- mxPath(from=c(latentsx[1],latentsx[1],latentsx[2]),
            #                  to = c(latentsx[1],latentsx[2],latentsx[2]), free = TRUE, arrows=2,values=
            #                    c(lgcm.x$intercept.variance, lgcm.x$intercept.slope.covariance,
            #                      lgcm.x$slope.variance),connect = "unique.pairs",
            #                  labels=c("interceptvariancex","interceptslopecovariancex","slopevariancex"))
            
            p5x <- mxPath(from="one",to=manifestsx,values=0,free=FALSE)
            
            # --
            
            p1y <- mxPath(from=latentsy[1], to=manifestsy, free=FALSE, values=1)
            p2y <- mxPath(from=latentsy[2], to=manifestsy, free=FALSE, values= lgcm.y$timepoints)
            
            p3y <- mxPath(from=manifestsy, to=manifestsy, free = TRUE, values=lgcm.y$residual.variance,
                          labels="residualvariancey", arrows=2)
            
            #      p4y <- mxPath(from=c(latentsy[1],latentsy[1],latentsy[2]),
            #                    to = c(latentsy[1],latentsy[2],latentsy[2]), free = TRUE, arrows=2,values=
            #                      c(lgcm.y$intercept.variance, lgcm.y$intercept.slope.covariance,
            #                        lgcm.y$slope.variance),connect = "unique.pairs",
            #                    labels=c("interceptvariancey","interceptslopecovariancey","slopevariancey"))
            # fixed variances
            p4y <- mxPath(from=c(latentsy.std[1],latentsy.std[2]), to=c(latentsy.std[1],latentsy.std[2]), free=FALSE,
                          arrows=2,values=1,connect="single")
            # projections to std latents
            p4yy <- mxPath(from=latentsy.std,to=latentsy, free=TRUE,arrows=1,
                           values=c(sqrt(lgcm.y$intercept.variance),
                                    sqrt(lgcm.y$slope.variance)),
                           connect="single")      
            
            p5y <- mxPath(from="one",to=manifestsy,values=0,free=FALSE)
            
            ## --
            
            p1xy <- mxPath(from=latentsx.std, to=latentsy.std, connect="unique.pairs",free=TRUE,arrows=2,
                           labels=c("iceptxicepty","iceptxslopey","slopexicepty","slopexslopey"),
                           values=c(generic.model$icept.x.icept.y.covariance /
                                        (sqrt(lgcm.x$intercept.variance)*sqrt(lgcm.y$intercept.variance)),
                                    generic.model$icept.x.slope.y.covariance /
                                        (sqrt(lgcm.x$intercept.variance)*sqrt(lgcm.y$slope.variance)),
                                    generic.model$icept.y.slope.x.covariance /
                                        (sqrt(lgcm.x$slope.variance)*sqrt(lgcm.y$intercept.variance)),
                                    generic.model$slope.x.slope.y.covariance /
                                        (sqrt(lgcm.x$slope.variance)*sqrt(lgcm.y$slope.variance))
                           )
            )
            
            p2xy <- mxPath(from=latentsx.std[1],to=latentsx.std[2],arrows=2, labels=c("interceptslopecorrelationx"),
                           values=lgcm.x$intercept.slope.covariance /
                               (sqrt(lgcm.x$slope.variance)*sqrt(lgcm.x$intercept.variance))
            )
            
            p3xy <- mxPath(from=latentsy.std[1],to=latentsy.std[2],arrows=2, labels=c("interceptslopecorrelationy"),
                           values=lgcm.y$intercept.slope.covariance /
                               (sqrt(lgcm.y$slope.variance)*sqrt(lgcm.y$intercept.variance))
            )
            
            lgcmModel <- mxModel("Latent Growth Curve Model",
                                 type="RAM",
                                 manifestVars=manifests,
                                 latentVars = latents,
                                 p1x,p2x,p3x,p4x,p4xx,p5x,
                                 p1y,p2y,p3y,p4y,p4yy,p5y,
                                 p1xy, p2xy, p3xy
            )
            
            
        } 
        else 
        {
            
            lgcm.x <- generic.model$lgcm.x
            lgcm.y <- generic.model$lgcm.y
            
            manifestsx <- paste("X", 1:length(lgcm.x$timepoints),sep="")
            manifestsy <- paste("Y", 1:length(lgcm.y$timepoints),sep="")
            manifests <- c(manifestsx, manifestsy)
            
            latentsx <- c("interceptx","slopex")
            latentsy <- c("intercepty","slopey")
            latents <- c(latentsx, latentsy)
            
            
            
            
            p1x <- mxPath(from=latentsx[1], to=manifestsx, free=FALSE, values=1)
            p2x <- mxPath(from=latentsx[2], to=manifestsx, free=FALSE, values= lgcm.x$timepoints)
            
            p3x <- mxPath(from=manifestsx, to=manifestsx, free = TRUE, values=lgcm.x$residual.variance,
                          labels="residualerrorx", arrows=2)
            
            p4x <- mxPath(from=c(latentsx[1],latentsx[1],latentsx[2]),
                          to = c(latentsx[1],latentsx[2],latentsx[2]), free = TRUE, arrows=2,values=
                              c(lgcm.x$intercept.variance, lgcm.x$intercept.slope.covariance,
                                lgcm.x$slope.variance),connect = "unique.pairs",
                          labels=c("interceptvariancex","interceptslopecovariancex","slopevariancex"))
            
            p5x <- mxPath(from="one",to=manifestsx,values=0,free=FALSE)
            
            # --
            
            p1y <- mxPath(from=latentsy[1], to=manifestsy, free=FALSE, values=1)
            p2y <- mxPath(from=latentsy[2], to=manifestsy, free=FALSE, values= lgcm.y$timepoints)
            
            p3y <- mxPath(from=manifestsy, to=manifestsy, free = TRUE, values=lgcm.y$residual.variance,
                          labels="residualerrory", arrows=2)
            
            p4y <- mxPath(from=c(latentsy[1],latentsy[1],latentsy[2]),
                          to = c(latentsy[1],latentsy[2],latentsy[2]), free = TRUE, arrows=2,values=
                              c(lgcm.y$intercept.variance, lgcm.y$intercept.slope.covariance,
                                lgcm.y$slope.variance),connect = "unique.pairs",
                          labels=c("interceptvariancey","interceptslopecovariancey","slopevariancey"))
            
            
            p5y <- mxPath(from="one",to=manifestsy,values=0,free=FALSE)
            
            ## --
            
            p1xy <- mxPath(from=latentsx, to=latentsy, connect="unique.pairs",free=TRUE,arrows=2,
                           labels=c("iceptxicepty","iceptxslopey","slopexicepty","slopexslopey"),
                           values=c(generic.model$icept.x.icept.y.covariance,
                                    generic.model$icept.x.slope.y.covariance,
                                    generic.model$icept.y.slope.x.covariance, 
                                    generic.model$slope.x.slope.y.covariance
                                    
                           ))
            
            lgcmModel <- mxModel("Latent Growth Curve Model",
                                 type="RAM",
                                 manifestVars=manifests,
                                 latentVars = latents,
                                 p1x,p2x,p3x,p4x,p5x,
                                 p1y,p2y,p3y,p4y,p5y,
                                 p1xy
            )
            
        } # end if (latent.correlation)
        
    } else {
        warning("Unknown model class.")
        return(NA)
    }
}



### Define the parameters of the models
mean_icept_x <- 0
mean_slope_x <- -14
var_icept_x <- 100
var_slope_x <- 200
error_var_x <- 25
cov_icept_x_slope_x <- 0

mean_icept_y <- 0
mean_slope_y <- -5
var_icept_y <- 100
var_slope_y <- 200
error_var_y <- 25
cov_icept_y_slope_y <- 0

cov_icept_x_icept_y <- 60
cov_icept_x_slope_y <- 0
cov_icept_y_slope_x <- 0

cor_slope_x_slope_y <- 1
cov_slope_x_slope_y <- cor_slope_x_slope_y * sqrt(var_slope_x*var_slope_y)

### Define the models
lgcm_x <- lgcm(
    timepoints = 1:10,
    intercept.variance = var_icept_x,
    slope.variance = var_slope_x,
    residual.variance = error_var_x,
    intercept.slope.covariance = cov_icept_x_slope_x,
    intercept.mean = mean_icept_x,
    slope.mean = mean_slope_x
)

lgcm_y <- lgcm(
    timepoints = 1:10,
    intercept.variance = var_icept_y,
    slope.variance = var_slope_y,
    residual.variance = error_var_y,
    intercept.slope.covariance = cov_icept_y_slope_y,
    intercept.mean = mean_icept_y,
    slope.mean = mean_slope_y
)

lgcm_xy <- bivariateLgcm(lgcm_x, lgcm_y,
                         icept.x.icept.y.covariance = cov_icept_x_icept_y,
                         slope.x.slope.y.covariance = cov_slope_x_slope_y)

### Lemma 1 Appendix A
# mean_x <- mean_icept_x + 1/2*mean_slope_x
# mean_y <- mean_icept_y + 1/2*mean_slope_y

var_x <- var_icept_x + 1/3*var_slope_x + cov_icept_x_slope_x + 
    1/12*mean_slope_x**2 + error_var_x
var_y <- var_icept_y + 1/3*var_slope_y + cov_icept_y_slope_y + 
    1/12*mean_slope_y**2 + error_var_y

cov_x_y <- cov_icept_x_icept_y + 
    1/2*(cov_icept_x_slope_y + cov_icept_y_slope_x) +
    1/3*cov_slope_x_slope_y +
    1/12*mean_slope_x*mean_slope_y

### Compute SOS
sos <- 1 - 12*var_y*(mean_slope_y*var_x - mean_slope_x*cov_x_y)**2 /
    ((var_y*var_x - cov_x_y**2) * (12*var_x - mean_slope_x**2)*mean_slope_y**2)

### Make an OpenMx model
mod1 <- toOpenMx(lgcm_xy)

### Simulate some data
dat <- simulateData(mod1, 20)



###
cor_vec <- seq(0,1,0.01)
sos_vec <- rep(NA, 101)
for (i in 1:101) {
    cor_slope_x_slope_y <- cor_vec[i]
    cov_slope_x_slope_y <- cor_slope_x_slope_y * sqrt(var_slope_x*var_slope_y)
   
    var_x <- var_icept_x + 1/3*var_slope_x + cov_icept_x_slope_x + 
        1/12*mean_slope_x**2 + error_var_x
    var_y <- var_icept_y + 1/3*var_slope_y + cov_icept_y_slope_y + 
        1/12*mean_slope_y**2 + error_var_y
    
    cov_x_y <- cov_icept_x_icept_y + 
        1/2*(cov_icept_x_slope_y + cov_icept_y_slope_x) +
        1/3*cov_slope_x_slope_y +
        1/12*mean_slope_x*mean_slope_y
    
    sos <- 1 - 12*var_y*(mean_slope_y*var_x - mean_slope_x*cov_x_y)**2 /
        ((var_y*var_x - cov_x_y**2) * (12*var_x - mean_slope_x**2)*mean_slope_y**2)
    
    sos_vec[i] <- sos
}

plot(cor_vec, sos_vec, type = "l")



### Plot the simulated data for two variables similar to Figure 1 in the paper
library(tidyverse)

dat_x <- as.data.frame(dat[, 1:10])
dat_y <- as.data.frame(dat[, 11:20])

names(dat_x) <- paste0("x", 1:10)
names(dat_y) <- paste0("y", 1:10)

dat_x_long <- gather(dat_x)
dat_y_long <- gather(dat_y)

dat_x_long$id <- rep(1:20, 10)
dat_y_long$id <- rep(1:20, 10)

names(dat_x_long)[1] <- "timepoint"
names(dat_y_long)[1] <- "timepoint"

dat_x_long$timepoint <- factor(rep(1:10, each = 20))
dat_y_long$timepoint <- factor(rep(1:10, each = 20))

plt_x <- ggplot(dat_x_long, aes(x = timepoint, y = value)) +
    geom_line(aes(group = id))
plt_y <- ggplot(dat_y_long, aes(x = timepoint, y = value)) +
    geom_line(aes(group = id))


set.seed(42)

dat_x_one_per_person <- data.frame(matrix(NA, nrow = 20, ncol = 3))
for (i in 1:20) {
    one_timepoint <- dat_x_long %>% filter(id == i) %>% sample_n(1)
    dat_x_one_per_person[i, ] <- one_timepoint
}

dat_y_one_per_person <- data.frame(matrix(NA, nrow = 20, ncol = 3))
for (i in 1:20) {
    one_timepoint <- dat_y_long %>% filter(id == i) %>% sample_n(1)
    dat_y_one_per_person[i, ] <- one_timepoint
}

names(dat_x_one_per_person) <- c("timepoint", "value", "id")
names(dat_y_one_per_person) <- c("timepoint", "value", "id")

plt_x_scatter <- ggplot(dat_x_one_per_person, aes(x = timepoint, y = value)) +
                            geom_point()

plt_y_scatter <- ggplot(dat_y_one_per_person, aes(x = timepoint, y = value)) +
                            geom_point()
