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