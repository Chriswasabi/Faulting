#' Calibrate the algorithm
#' @export
#' @param database csv file with the profiles
#' @param dif_per_1 lower percentile to remove in the difference stage
#' @param dif_per_2 upper percentile to remove in the difference stage
#' @param dif_per_3 percentile used in fine tuning (should be higher than per_2)

calib <- function(database, dif_per_1, dif_per_2, dif_per_3) {

  df <- database[,-c(1,2)] * -0.01

  tot = nrow(df)
  pb <- progress::progress_bar$new( format = "  Calibrating difference thresholds [:bar] :percent eta: :eta",
                                    total = tot, clear = FALSE, width= 140)
  db <- numeric()
  dif <- matrix(ncol=(ncol(df)-4), nrow=nrow(df))
  q1 <- numeric()
  q3 <- numeric()
  iqr <- numeric()

  for (j in 1:tot) {

    profile <- as.numeric(df[j,])
    profile <- profile - stats::median(profile)

    for (i in 3:(length(profile)-2)) {
      db[i-2]=profile[i]-profile[i-1]
    }
    dif[j,] <- db
    pb$tick()
    Sys.sleep(1/tot)
  }

  pb <- progress::progress_bar$new( format = "  Calibrating boxplot thresholds [:bar] :percent eta: :eta",
                                    total = tot, clear = FALSE, width= 140)

  for (j in 1:tot) {

    profile <- as.numeric(df[j,])
    profile <- profile - stats::median(profile)
    q1[j] <- stats::quantile(x = profile, probs = 0.05)
    q3[j] <- stats::quantile(x = profile, probs = 0.95)

    pb$tick()
    Sys.sleep(1/tot)
  }

  iq <- q3 - q1
  t1 <- mean(q1) - 3*stats::quantile(x = iq, probs = 0.95)
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
  pb <- progress::progress_bar$new( format = "  Denoising Profiles [:bar] :percent eta: :eta",
                                    total = int, clear = FALSE, width= 140)
  na_counter <- numeric()
  dl <- numeric()
  db <- numeric()
  f_database <- matrix(ncol=ncol(df), nrow=int)
  x_fit <- seq(1:2048)
  x2 <- x_fit^2

  for (j in 1:int) {

    profile <- as.numeric(df[j,])
    fit1 <- lm(profile ~ x_fit)
    y_hat1 <- unname(fit1$coefficients[[1]]) + unname(fit1$coefficients[[2]])*x_fit
    profile <- profile - y_hat1

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
      if ((db[i-2]>dt_2 & dl[i-2]==0) | (db[i-2]< dt_1 & dl[i-2]==0) ) {
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

#' Detect joint
#' @export
#' @param database csv file with the profiles
joint_detect <- function(database) {

  df=database
  tot = nrow(df)
  sample_size <- ncol(df)-1
  jc <- numeric()

  for (int in 1:tot) {

    profile <- as.numeric(df[int,][-1])
    profile[profile>-1] = -1


    d <- numeric()
    w <- numeric()
    t <- numeric(); kk=1; k=1
    tt <- numeric(); j=1


    for (i in 1:(sample_size-1)) {
      if((profile[i]<(-1) & profile[i+1]==(-1)) | (profile[i]==(-1) & profile[i+1]<(-1))) {
        d[i] <- 1

      } else {
        d[i] <- 0
      }

    }
    d[sample_size] <- 0

    for (i in 1:sample_size) {
      if(i==1) {

        w[i] <- 0

      } else { if ( d[i]==0  ) {

        w[i] <- w[i-1] + 1

      } else {

        w[i] <- 0

      }
      }
    }

    for (i in 1:sample_size) {
      if(i==1) {
        next()

      } else { if ( w[i]==0  ) {

        t[k] <- kk
        k=k+1
        kk=1

      } else { if (i==2044) {

        t[k] <- kk + 1

      } else {

        kk = kk + 1

      }
      }
      }
    }

    for (i in 1:sample_size) {
      if(i==1) {

        tt[i] = t[j]

      } else { if ( w[i]==0  ) {

        j=j+1
        tt[i] <- t[j]

      } else {

        tt[i] <- t[j]
      }
      }
    }

    z <- as.data.frame(profile)
    z$width = tt

    z1 <- z %>% mutate(joint_test = ifelse(test = profile < -3 & width > 30, 1, 0))

    z2 <- z1 %>% distinct(width, joint_test) %>% filter(joint_test==1)
    jw = max(z2$width)

    z3 <- z1 %>% mutate(joint = ifelse(test = profile < -3 & width == jw, 1, 0))

    joint_center = ifelse(max(z3$joint)==1,
                          (min(which(z3$joint %in% 1)) + max(which(z3$joint %in% 1)))/2,0)
    joint  =ifelse(joint_center==0, 0, which(z1$joint %in% 1))


    jc[int] = joint_center

    print(100*int/tot)
  }

  f_database <- cbind(jc, df)
  return(f_database)
}


#' @import "stats"
#' @import "tidyverse"
#' @import "imputeTS"
#' @import "progress"
