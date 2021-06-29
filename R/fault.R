#' Calibrate the algorithm
#' @export
#' @param database csv file with the profiles
#' @param dif_per_1 lower percentile to remove in the difference stage
#' @param dif_per_2 upper percentile to remove in the difference stage
#' @param dif_per_3 percentile used in fine tuning (should be higher than per_2)

calib <- function(database, dif_per_1, dif_per_2, dif_per_3) {

  df <- database[,-c(1,2)] * -0.01

  db <- numeric()
  dif <- matrix(ncol=(ncol(df)-4), nrow=nrow(df))
  q1 <- numeric()
  q3 <- numeric()
  iqr <- numeric()
  x_fit <- seq(1:2044)
  x2 <- x_fit^2

  for (j in 1:nrow(df)) {

    profile <- as.numeric(df[j,])
    profile <- profile - stats::median(profile)

    for (i in 3:(length(profile)-2)) {
      db[i-2]=profile[i]-profile[i-1]
    }
    dif[j,] <- db
    print(100*j/nrow(df))
  }

  for (j in 1:nrow(database)) {

    profile <- as.numeric(df[j,])
    profile <- profile - stats::median(profile)
    q1[j] <- stats::quantile(x = profile, probs = 0.25)
    q3[j] <- stats::quantile(x = profile, probs = 0.75)
    print(100*j/nrow(database))
  }

  iq <- q3 - q1
  t1 <- mean(q1) - 3*mean(iq)
  t3 <- mean(q3) + 3*mean(iq)

  a <- c(stats::quantile(dif, probs = c(dif_per_1, dif_per_2, dif_per_3)), t1, t3)

  return (a)
}

#' Calibrate the algorithm
#' @export
#' @param database csv file with the profiles
#' @param dt_1 lower differential threshold
#' @param dt_2 upper differential threshold
#' @param dt_3 fine tuned differential threshold
#' @param t1 upper boxplot threshold
#' @param t3 lower boxplot threshold
denoise <- function(database, dt_1, dt_2, dt_3, t1, t3) {

  df <- database[,-c(1,2)] * -0.01

  int=nrow(df)
  pb <- progress::progress_bar$new(total = int)
  na_counter <- numeric()
  dl <- numeric()
  db <- numeric()
  f_database <- matrix(ncol=ncol(df), nrow=int)

  for (j in 1:int) {

    profile <- as.numeric(df[j,])
    profile <- profile - median(profile)

    #Difference with lag
    for (i in 3:(length(profile)-2)) {
      db[i-2]=profile[i]-profile[i-1]
    }

    #Difference with lead
    for (i in 3:(length(profile)-2)) {
      dl[i-2]=profile[i]-profile[i+1]
    }

    #BoxPlot Outlier Removal
    for (i in 1:length(profile)) {
      if (profile[i] > t3 | profile[i] < t1 ) {
        profile[i]=NA
      }
    }

    #Flatline Detection
    for (i in 3:(length(profile)-2)) {
      if (db[i-2]==0) {
        profile[i]=NA
      }
    }

    #Spike Detection
    for (i in 3:(length(profile)-2)) {
      if ((db[i-2]>dt_2 & dl[i-2]>dt_2) | (db[i-2]< dt_1 & dl[i-2]< dt_1) ) {
        profile[i]=NA
      }
    }

    #Spike + Flatline Detection
    for (i in 3:(length(profile)-2)) {
      if ((db[i-2]>150 & dl[i-2]==0) | (db[i-2]< dt_1 & dl[i-2]==0) ) {
        profile[i]=NA
      }
    }

    #Fine fine tune spike Detection
    for (i in 3:(length(profile)-2)) {
      if ( i == 3 | i ==  2046) {
        next()
      } else {
        if ( (db[i-1] >= dt_3 & dl[i-2] <= -dt_3) | (db[i-2] <= -dt_3 & dl[i-3] >= dt_3) |
             (db[i-2] >= dt_3 & dl[i-3] <= -dt_3) ) {
          profile[i]=NA
        } }
    }

    #Spike in between flatline Detection
    for (i in 3:(length(profile)-2)) {
      if (is.na(profile[i])) {
        next()
      } else {
        if ( is.na(profile[i+1]) & is.na(profile[i-1]) ) {
          profile[i]=NA
        }
      }
    }

    # End Point Pre-imputation
    #This is because if the endpoints are flatlines they will not be removed based on the for loop iteration limits
    for (i in 3:(length(profile)-2)) {
      if (is.na(profile[i] & is.na(profile[i+1]) & i==3) | (is.na(profile[i]) & is.na(profile[i-1]) & i==2046)) {
        profile[i]= median(x = profile, na.rm = T)
      }
    }

    counter <- sum(is.na(profile))*100/length(profile)

    ##Imputation##
    #This part of the code will impute the profile by linear interpolation
    k <- imputeTS::na_interpolation(profile)

    ##Detrending##
    #This part of the code detrends the profiles
    fit <- lm(k ~ x_fit + x2)
    y_hat <- unname(fit$coefficients[[1]]) + unname(fit$coefficients[[2]])*x_fit + unname(fit$coefficients[[3]])*x2

    k_d <- k - y_hat

    na_counter[j] <- counter
    f_database[j,] <- k_d

    pb$tick()
    Sys.sleep(1 / int)
  }

  f_database <- as.data.frame(f_database)
  cf_database <- f_database[,-c(2046,2045,1,2),]
  cf_database <- cbind(na_counter,cf_database)

  return(cf_database)
}


#' @import "stats"
#' @import "tidyverse"
#' @import "imputeTS"
#' @import "pracma"
#' @import "progress"
