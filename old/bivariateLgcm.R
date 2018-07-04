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