#' @import "stats"
#' @import "tidyverse"


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

  a <- c(stats::quantile(dif, probs = c(dif_per_1, dif_per_2, dif_per_3)), mean(q1), mean(q3), mean(iq))

  return (a)
}

