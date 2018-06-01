dat_transform <- function(data, sample_n, timepoints_n) {
    dat_x <- as.data.frame(data[, 1:(sample_n/2)])
    dat_y <- as.data.frame(data[, (sample_n/2+1):sample_n])
    
    names(dat_x) <- paste0("x", 1:timepoints_n)
    names(dat_y) <- paste0("y", 1:timepoints_n)
    
    dat_x_long <- gather(dat_x)
    dat_y_long <- gather(dat_y)
    
    dat_x_long$id <- rep(1:sample_n, timepoints_n)
    dat_y_long$id <- rep(1:sample_n, timepoints_n)
    
    names(dat_x_long)[1] <- "timepoint"
    names(dat_y_long)[1] <- "timepoint"
    
    dat_x_long$timepoint <- factor(rep(1:timepoints_n, each = sample_n))
    dat_y_long$timepoint <- factor(rep(1:timepoints_n, each = sample_n))
    
    # dat_long <- rbind(dat_x_long, dat_y_long)
    # dat_long$var <- factor(c(rep("x", nrow(dat_x_long)), 
    #                          rep("y", nrow(dat_y_long))))
    # 
    # dat_long
    
    list(dat_x_long, dat_y_long)
    
}

dat_sample <- function(data, n_sample, seed = sample(0:100, 1)) {
    set.seed(seed)
    
    dat_long <- rbind(data[[1]], data[[2]])
    
    dat_x_one_per_person <- data.frame(matrix(NA, nrow = n_sample, ncol = 3))
    for (i in 1:n_sample) {
        one_timepoint <- dat_long[1:nrow(dat_long)/2, ] %>% filter(id == i) %>% sample_n(1)
        dat_x_one_per_person[i, ] <- one_timepoint
    }
    
    dat_y_one_per_person <- data.frame(matrix(NA, nrow = n_sample, ncol = 3))
    for (i in 1:n_sample) {
        one_timepoint <- dat_long[nrow(dat_long)/2+1:nrow(dat_long), ] %>% filter(id == i) %>% sample_n(1)
        dat_y_one_per_person[i, ] <- one_timepoint
    }
    
    names(dat_x_one_per_person) <- c("timepoint", "value", "id")
    names(dat_y_one_per_person) <- c("timepoint", "value", "id")
    
    # dat_one_per_person <- rbind(dat_x_one_per_person, dat_y_one_per_person)
    # dat_one_per_person$var <- factor(c(
    #     rep("x", nrow(dat_one_per_person)/2), 
    #     rep("y", nrow(dat_one_per_person)/2)))
    # 
    # dat_one_per_person
    
    list(dat_x_one_per_person, dat_y_one_per_person)
}