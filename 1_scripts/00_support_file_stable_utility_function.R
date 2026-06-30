utilitynlmTemp <- function (formula, prior, desvars, criterion = c("D", "A", "E", "SIG", "NSEL"), 
                        method = c("quadrature", "MC"), nrq) 
{
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  
  if(missing(nrq)) {
    nrq <- switch(method,
                  quadrature = c(2, 8),
                  NULL)
  }  
  nr <- nrq[[1]]
  nq <- nrq[[2]]
  
  
  allvars <- all.vars(formula)
  
  
  paravars <- setdiff(allvars, desvars)
  p <- length(paravars)
  k <- length(desvars)
  
  DD <- deriv(expr = formula, namevec = paravars)
  
  aDD <- as.character(DD)
  grad <- function() {
  }
  gradtext <- "grad<-function(d,paras){ \n"
  gradtext <- paste(gradtext, desvars[1], "<-d[,1]", " \n ", 
                    sep = "")
  if (k > 1) {
    for (j in 2:k) {
      gradtext <- paste(gradtext, desvars[j], "<-d[,", 
                        j, "]", " \n ", sep = "")
    }
  }
  gradtext <- paste(gradtext, paravars[1], "<-paras[,1]", " \n ", 
                    sep = "")
  if (p > 1) {
    for (j in 2:p) {
      gradtext <- paste(gradtext, paravars[j], "<-paras[,", 
                        j, "]", " \n ", sep = "")
    }
  }
  gradtext <- paste(gradtext, substr(x = aDD, start = 2, stop = nchar(aDD)), 
                    sep = "")
  eval(parse(text = gradtext))
  if (criterion == "SIG") {
    inte <- function(d, B) {
      n1 <- dim(d)[1]
      B2 <- B * 2
      sam <- prior(B2)
      d3 <- matrix(0, ncol = k, nrow = B2 * n1)
      for (i in 1:k) {
        d3[, i] <- rep(d[, i], B2)
      }
      sam3 <- matrix(0, ncol = p, nrow = B2 * n1)
      for (i in 1:p) {
        sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                   colnames(sam) == paravars[i]], each = n1)
      }
      mu1 <- matrix(grad(d = d3, paras = sam3)[1:(B2 * 
                                                    n1)], ncol = n1, byrow = TRUE)
      
      mu2 <- matrix(mu1[-(1:B), ],nrow=B)
      mu1 <- matrix(mu1[1:B, ],nrow=B)
      
      y <- matrix(rnorm(n = n1 * B, mean = as.vector(mu1), 
                        sd = rep(sqrt(sam[1:B, colnames(sam) == "sig2"]), 
                                 n1)), ncol = n1)
      
      as.vector(SIGnlmcpp_stable(y, mu1, mu2, 
                                 sam[1:B, colnames(sam) == "sig2"], 
                                 sam[-(1:B), colnames(sam) == "sig2"]))
    }
  }
  
  if (criterion == "NSEL") {
    inte <- function(d, B) {
      n1 <- dim(d)[1]
      B2 <- B * 2
      sam <- prior(B2)
      d3 <- matrix(0, ncol = k, nrow = B2 * n1)
      for (i in 1:k) {
        d3[, i] <- rep(d[, i], B2)
      }
      sam3 <- matrix(0, ncol = p, nrow = B2 * n1)
      for (i in 1:p) {
        sam3[, paravars == paravars[i]] <- rep(sam[, 
                                                   colnames(sam) == paravars[i]], each = n1)
      }
      mu1 <- matrix(grad(d = d3, paras = sam3)[1:(B2 * 
                                                    n1)], ncol = n1, byrow = TRUE)
      mu2 <- matrix(mu1[-(1:B), ],nrow=B)
      mu1 <- matrix(mu1[1:B, ],nrow=B)
      y <- matrix(rnorm(n = n1 * B, mean = as.vector(mu1), 
                        sd = rep(sqrt(sam[1:B, colnames(sam) == "sig2"]), 
                                 n1)), ncol = n1)
      
      
      
      as.vector(NSELnlmcpp_stable(y, mu2, 
                                       sam[-(1:B), colnames(sam) == "sig2"], 
                                       sam[1:B, colnames(sam) != "sig2"], 
                                       sam[-(1:B), colnames(sam) != "sig2"]))
    }
  }
  
  output <- list(utility = inte)
  output
}
