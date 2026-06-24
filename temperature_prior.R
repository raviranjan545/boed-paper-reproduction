# If you are running this independently, 
# Load MASS and Matrix packages first
prior <- function(B){
  # MVN parameters for (log a, log b, tmax, tmin)
  mean_a <- 0.002; sd_a <- 0.0002
  mean_b <- 0.025; sd_b <- 0.005
  mean_tmax <- 35;  sd_tmax <- 2
  mean_tmin <- 2;   sd_tmin <- 3
  
  # Convert to lognormal params for a, b
  sdlog_a   <- sqrt(log(1 + (sd_a/mean_a)^2))
  meanlog_a <- log(mean_a) - 0.5 * sdlog_a^2
  sdlog_b   <- sqrt(log(1 + (sd_b/mean_b)^2))
  meanlog_b <- log(mean_b) - 0.5 * sdlog_b^2
  
  # Correlations
  rho_ab     <- -0.6
  rho_a_tmax <- -0.5
  rho_a_tmin <-  0.3
  
  # Build correlation matrix
  C <- diag(4)
  dimnames(C) <- list(c("loga","logb","tmax","tmin"), c("loga","logb","tmax","tmin"))
  C["loga","logb"] <- C["logb","loga"] <- rho_ab
  C["loga","tmax"] <- C["tmax","loga"] <- rho_a_tmax
  C["loga","tmin"] <- C["tmin","loga"] <- rho_a_tmin
  C <- as.matrix(nearPD(C, corr = TRUE)$mat)
  
  # MVN on (log a, log b, tmax, tmin)
  mu_vec <- c(meanlog_a, meanlog_b, mean_tmax, mean_tmin)
  sd_vec <- c(sdlog_a,   sdlog_b,   sd_tmax,  sd_tmin)
  Sigma  <- diag(sd_vec) %*% C %*% diag(sd_vec)
  
  # Draw from multivariate normal
  draws <- mvrnorm(B, mu = mu_vec, Sigma = Sigma)
  
  # Transform a and b back to original scale
  a    <- exp(draws[,1])
  b    <- exp(draws[,2])
  tmax <- draws[,3]
  tmin <- draws[,4]
  
  # Draw sigma and square it
  sig <- rlnorm(n = B, meanlog = -2.3, sdlog = 0.1)
  sig2 <- sig^2
  
  # Combine and return
  out <- cbind(a, b, tmax, tmin, sig2)
  colnames(out) <- c("a", "b", "tmax", "tmin", "sig2")
  return(out)
}