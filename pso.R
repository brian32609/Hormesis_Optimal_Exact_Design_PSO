library(globpso)
source("HEOD_PSO.R")

#Initialize PSO settings
psoinfo <- psoinfo_setting()


# Hunt-Bowman d-optimal
hb1_parms <- hb_parms(c1 = 170, tau = 0.04, b0 = 1.46, b1 = 40) # Set Hunt-Bowman model parameters.
pso_test <- hb_pso(nPoints = 4, parms = hb1_parms, psoinfo = psoinfo) # Set PSO settings
pso_test$par 
pso_test$val 

hb_doptimal(d = c(0, 0, 0.02, 0.02, 0.02, 0.04, 0.04, 0.04, 0.0991, 0.0991, 0.0991, 
                  0, 0.02, 0.04, 0.0991, 0, 0.02, 0.04, 0.0991), 
            loc = hb1_parms)

hb_res <- hb_pso_rep(mRep = 10, nPoints = 5, parms = hb1_parms, psoinfo = psoinfo)

hb_res$efficiency
hb_res$best_point

hb_results <- list()
psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)

for(i in 16:20){
  hb_results[[i]] <- hb_pso_rep(nRep = 10, nPoints = i, 
                                          parms = hb1_parms, psoinfo = psoinfo)
  print(i)
  print(hb_results[[i]]$result$best_eff)
  print(hb_results[[i]]$result$best_points)
}

saveRDS(hb_results, "doptimal.rds")

dop <- sapply(4:20, function(x) hb_results[[x]]$result$best_eff^(1/x))
plot(x = c(4:20), y = dop, type = "l", xlab = "Number of Design Points", ylab = "Efficiency")
points(x = c(4:20), y = dop, pch = 19)


# exp-log h-optimal
el_par <- exp_log_params(0.15, 89, 3.2, 41)
1 / exp_log_hoptimal(d = c(0.0124, 0.0620, 0.1243), el_par)



exp_log_1$val
exp_log_1$par
plot(x = c(1:1001), y = exp_log_1$history, type = "l", 
     xlab = "iterations", ylab = "criterion value", 
     main = "Optimal value from each iteration.")

exp_log_6 <- exp_log_pso_rep(nRep = 3, nPoints = 6, parms = el_par, psoinfo = psoinfo)
exp_log_6$efficiency
exp_log_6$design_points
exp_log_6$result$best_eff
exp_log_6$result$best_points


exp_log_results <- list()
psoinfo <- psoinfo_setting(nSwarms = 512, Iters = 30000)

for(i in 15:17){
  exp_log_results[[i]] <- exp_log_pso_rep(nRep = 10, nPoints = i, 
                                          parms = el_par, psoinfo = psoinfo)
  print(i)
}


psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 2000)

for(i in 4:10){
  exp_log_results[[i]] <- exp_log_pso_rep(nRep = 1, nPoints = i, 
                                          parms = el_par, psoinfo = psoinfo)
  
  print(i)
  print(exp_log_results[[i]]$result$best_eff)
  print(exp_log_results[[i]]$result$best_points)
}

psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 30000)

for(i in 15:17){
  exp_log_results[[i]] <- exp_log_pso_rep(nRep = 10, nPoints = i, 
                                          parms = el_par, psoinfo = psoinfo)
  print(i)
}

exp_log_results[[15]]$efficiency
exp_log_results[[15]]$design_points
exp_log_results[[17]]$result$best_eff
exp_log_results[[17]]$result$best_points

saveRDS(exp_log_results, "hoptimal.rds")


ndesigns <- matrix(ncol = 4, nrow = 17)
rownames(ndesigns) <- c(4:20)
colnames(ndesigns) <- c(0.000, 0.0108, 0.0526, 0.1187)
for(j in 4:20){
  ndesigns[j-3, 1] = .371 * j
  ndesigns[j-3, 2] = .501 * j
  ndesigns[j-3, 3] = .087 * j
  ndesigns[j-3, 4] = .041 * j
}


h_round <- list(points = matrix(ncol = 4, nrow = 17), eff = c(0,0,0))
rownames(h_round$points) <- c(4:20)
colnames(h_round$points) <- c(0.000, 0.0108, 0.0526, 0.1187)
h_round$eff[4] <- 1 / exp_log_hoptimal(d = c(rep(0, 1), rep(0.0108, 2), rep(0.0526, 1), rep(0.1187, 0)), 
                 loc = el_par)
h_round$points[4-3,] <- c(1, 2, 1, 0)
saveRDS(h_round, "hround.rds")

ndesigns[5-3,]
1 / exp_log_hoptimal(d = c(rep(0, 1), rep(0.0108, 2), rep(0.0526, 1), rep(0.1187, 0)), 
                     loc = el_par)


h_best <- list(points = matrix(ncol = 4, nrow = 17), eff = c(0,0,0))
rownames(h_best$points) <- c(4:20)
colnames(h_best$points) <- c(0.000, 0.0108, 0.0526, 0.1187)
for(i in 4:20){
  res = 0
  points = c(0,0,0,0)
  for(a in 0:i){
    for(b in 0:(i-a)){
      for(c in 0:(i-a-b)){
        eff <- 1 / exp_log_hoptimal(d = c(rep(0, a), rep(0.0108, b), rep(0.0526, c), rep(0.1187, i-a-b-c)), 
                                    loc = el_par)
        if(eff > res){
          res = eff
          points <- c(a, b, c, i-a-b-c)
        }
      }
    }
  }
  h_best$points[i-3,] <- points
  h_best$eff[i] <- res
}
h_best$points
h_best$eff
saveRDS(h_best, "hbest.rds")

exp_log_loc_mat <- function(d, c0, c1, b0, b1){
  d <- c(0, d)
  n <- length(d)
  w <- c(0.371, 0.501, 0.087, 0.041)
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1))
  
  #print(mat_list)
  M <- Reduce("+", mat_list)
  M
}

exp_log_loc_hoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 999999999
  
  n <- length(d) + 1
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  M <- exp_log_loc_mat(d, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  } 
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res
}

efb <- exp_log_loc_hoptimal(d = c(0.0108, 0.0526, 0.1187), el_par)
efb

psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 2000)
exp_log_tau <- exp_log_tau_pso(4, el_par, psoinfo)
exp_log_tau$points
exp_log_tau$weight

exp_log_tau_res <- list()

for(j in 4:10){
  exp_log_tau_res[[j]] <- exp_log_loc_tau_pso_rep(3, j, el_par, psoinfo)
  print(j)
  print(exp_log_tau_res[[j]]$design_points)
  print(exp_log_tau_res[[j]]$weight)
}

psoinfo <- psoinfo_setting(nSwarms = 512, Iters = 10000)
exp_log_tau_res1 <- list()
for(j in 7:10){
  exp_log_tau_res1[[j]] <- exp_log_tau_pso_rep(10, j, el_par, psoinfo)
  print(j)
  print(exp_log_tau_res1[[j]]$design_points)
  print(exp_log_tau_res1[[j]]$val)
}

tauop <- c()
tauop[17] <- exp_log_tauoptimal(d = c(rep(0,10), rep(0.04197718,10)), loc = c(el_par, 0.04197718))
tauoptimal <- exp_log_tauoptimal(d = c(rep(0,2), rep(0.04197718,2)), loc = c(el_par, 0.04197718)) / tauop
saveRDS(tauoptimal, "tauop.rds")

tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
               c0 = 0.15, c1 = 89, b0 = 3.2, b1 = 41)$root
tau

tauoptimal <- sapply(1:17, function(x) tauoptimal[x]^(1/(x+3)))
plot(x = c(4:20), y = tauoptimal, xlab = "Number of Points", ylab = "Efficiency", type = "l")
points(x = c(4:20), y = tauoptimal, pch = 19, cex = .75)

h_op <- list()
for (j in 16:20){
  if(j <= 10){
    psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 2000)
    h_op[[j]] <- exp_log_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else if (j <= 15){
    psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)
    h_op[[j]] <- exp_log_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else {
    psoinfo <- psoinfo_setting(nSwarms = 1024, Iters = 30000)
    h_op[[j]] <- exp_log_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  }
  print(psoinfo$nSwarm)
  print(j)
}

par(mfrow = c(2,2))
plot(x = c(0:2000), y = h_op[[5]]$history^(1/5), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "5-Point Exact Design")
plot(x = c(0:2000), y = h_op[[10]]$history^(1/10), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "10-Point Exact Design")
plot(x = c(0:10000), y = h_op[[15]]$history^(1/15), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "15-Point Exact Design")
plot(x = c(0:30000), y = h_op[[20]]$history^(1/20), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "20-Point Exact Design")




h_cpu <- lapply(4:20, function(x) h_op[[x]]$cputime)
plot(x = c(4:20), y = h_cpu, type = "l",
                xlab = "Number of Design Points", ylab = "CPU Runtime(sec.)")

d_op <- list()
for (j in 4:20){
  if(j <= 10){
    psoinfo <- psoinfo_setting(nSwarms = 64, Iters = 1000)
    d_op[[j]] <- hb_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else if (j <= 15){
    psoinfo <- psoinfo_setting(nSwarms =128, Iters = 5000)
    d_op[[j]] <- hb_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  } 
  else {
    psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 10000)
    d_op[[j]] <- hb_pso(nPoints = j, parms = el_par, psoinfo = psoinfo)
  }
  print(psoinfo$nSwarm)
  print(j)
}

par(mfrow = c(2,2))
plot(x = c(0:2000), y = h_op[[5]]$history^(1/5), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "5-Point Exact Design")
plot(x = c(0:2000), y = h_op[[10]]$history^(1/10), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "10-Point Exact Design")
plot(x = c(0:10000), y = h_op[[15]]$history^(1/15), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "15-Point Exact Design")
plot(x = c(0:30000), y = h_op[[20]]$history^(1/20), type = "l", 
     xlab = "Number of Iterations", ylab = "Efficiency", main = "20-Point Exact Design")




d_cpu <- lapply(4:20, function(x) d_op[[x]]$cputime)
plot(x = c(4:20), y = d_cpu, type = "l",
     xlab = "Number of Design Points", ylab = "CPU Runtime(sec.)")


logistic_doptimal <- function(d, loc){
  a <- loc[1]
  b <- loc[2]
  n <- length(d)
  
  theta <- a + b * d
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d * omega), sum(d * omega), sum(d^2 * omega)), nrow = 2)
  -det(inf_mat) / n^2
}

logistic_pso <- function(npoints, parms, psoinfo){
  lb <- rep(-10, npoints)
  ub <- rep(10, npoints)
  
  pso_res <- globpso(objFunc = logistic_doptimal, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

psoinfo <- psoinfo_setting(nSwarms = 128, Iters = 2000)
logistic_base <- logistic_pso(2, parms = c(2, 1), psoinfo = psoinfo)
sort(2.145 + 0.658 * logistic_best$par)
logistic_base$val
logistic_best_list <- list(points = list(), val = c())

for(k in 2:20){
  logistic_best <- logistic_pso(k, parms = c(2, 1), psoinfo = psoinfo)
  logistic_best_list$points[[k]] <- sort(2 + 1 * logistic_best$par)
  logistic_best_list$val <- c(logistic_best_list$val, (logistic_best$val/logistic_base$val)^(1/k))
}

logistic_best_list$val
plot(x = c(2:20), y = logistic_best_list$val, type = "l", 
     xlab = "Number of Design Points", ylab = "Efficiency", 
     main = "Efficiency to number of design points for logistic model.\nalpha = 2, beta = 1")

qlogsitic_mat <- function(){
  
}

qlogistic_doptimal1 <- function(d, loc){
  a <- loc[1]
  b <- loc[2]
  mu <- loc[3]
  n <- length(d)
  
  eta <- a + b * (x - mu)^2
  pi <- exp(eta) / (1 + exp(eta))
  v <- pi * (1 - pi)
  
  
  
}

qlogistic_doptimal2 <- function(d, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  n <- length(d)
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  
  inf_mat <- matrix(c(sum(omega), sum(d * omega), sum(d^2 * omega), 
                      sum(d * omega), sum(d^2 * omega), sum(d^3 * omega), 
                      sum(d^2 * omega), sum(d^3 * omega), sum(d^4 * omega)), 
                    nrow = 3)
  #print(inf_mat)
  
  -det(inf_mat) / n^3
}

qlogistic_doptimal3 <- function(d, loc){
  a <- loc[1]
  b1 <- loc[2]
  b2 <- loc[3]
  d <- d
  n <- length(d)
  
  theta <- a + b1 * d + b2 * d^2
  omega <- exp(theta) / (1 + exp(theta))^2
  w <- c(0.297, 0.203, 0.203, 0.297)
  
  inf_mat <- matrix(c(sum(w * omega), sum(w * d * omega), sum(w * d^2 * omega), 
                      sum(w * d * omega), sum(w * d^2 * omega), sum(w * d^3 * omega), 
                      sum(w * d^2 * omega), sum(w * d^3 * omega), sum(w * d^4 * omega)), 
                    nrow = 3)
  #print(inf_mat)
  
  -det(inf_mat)
}

qlogistic_pso1 <- function(npoints, parms, psoinfo){
  lb <- rep(-3, npoints-1)
  ub <- rep(3, npoints-1)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal1, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

qlogistic_pso2 <- function(npoints, parms, psoinfo){
  lb <- rep(-3, npoints)
  ub <- rep(3, npoints)
  
  pso_res <- globpso(objFunc = qlogistic_doptimal2, lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, loc = parms)
  
  pso_res
}

psoinfo <- psoinfo_setting(nSwarms = 256, Iters = 2000)

qlogistic_best1 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso1(j, parms = c(3, 0, -1), psoinfo = psoinfo)
  qlogistic_best1$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best1$val <- c(qlogistic_best1$val, qlogistic_best$val)
}  

qlogistic_best1$points
qlogistic_best1$val

qlogistic_best2 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso2(j, parms = c(3, 0 ,-1), psoinfo = psoinfo)
  qlogistic_best2$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best2$val <- c(qlogistic_best2$val, qlogistic_best$val)
}  

qlogistic_best2$points
qlogistic_best2$val
qlogistic_base <- qlogistic_doptimal3(d = c(-2.061, -1.324, 1.324, 2.061), loc = c(3, 0, -1))
eff2 <- qlogistic_best2$val[2:9] / qlogistic_base
eff2

plot(x = c(3:10), y = sapply(3:10, function(x) eff2[x-2]^(1/x)), type = "l", 
     xlab = "Number of Design Points", ylab = "Efficiency", 
     main = "Efficiency of D-optimal exact design for quadratic loistic model
     alpha = 3, beta = 0, mu = -1")


qlogistic_best01 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso1(j, parms = c(0, 0 ,-1), psoinfo = psoinfo)
  qlogistic_best01$points[[j]] <- sort(c(qlogistic_best$par,0)) |> round(4)
  qlogistic_best01$val <- c(qlogistic_best01$val, qlogistic_best$val)
}  

qlogistic_best01$points
qlogistic_best01$val



qlogistic_best02 <- list(points = list(), val = c())
for(j in 2:10){
  qlogistic_best <- qlogistic_pso2(j, parms = c(0, 0 ,-1), psoinfo = psoinfo)
  qlogistic_best02$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best02$val <- c(qlogistic_best02$val, qlogistic_best$val)
}  

qlogistic_best02$points
qlogistic_best02$val



qlogistic_best11 <- list(points = list(), val = c())
for(j in 3:10){
  qlogistic_best <- qlogistic_pso1(j, parms = c(-3, 0, -1), psoinfo = psoinfo)
  qlogistic_best11$points[[j]] <- sort(c(qlogistic_best$par,0)) |> round(4)
  qlogistic_best11$val <- c(qlogistic_best11$val, qlogistic_best$val)
}  
qlogistic_best11$points
qlogistic_best11$val



qlogistic_best12 <- list(points = list(), val = c())
for(j in 3:10){
  qlogistic_best <- qlogistic_pso2(j, parms = c(-3, 0, -1), psoinfo = psoinfo)
  qlogistic_best12$points[[j]] <- sort(qlogistic_best$par) |> round(4)
  qlogistic_best12$val <- c(qlogistic_best12$val, qlogistic_best$val)
}  
qlogistic_best <- qlogistic_pso2(3, parms = c(-3, 0, -1), psoinfo = psoinfo)
qlogistic_best$par
qlogistic_best12$points
qlogistic_best12$val
qlogistic_base1 <- qlogistic_best12$val[1]
eff1 <- qlogistic_best12$val / qlogistic_base1
sapply(3:10, function(x) eff1[x-2]^(1/x))
plot(x = c(3:10), y = sapply(3:10, function(x) eff1[x-2]^(1/x)), type = "l", 
     xlab = "Number of Design Points", ylab = "Efficiency", 
     main = "Efficiency of D-optimal exact design for quadratic loistic model
     alpha = -3, beta = 0, mu = -1")

(qlogistic_doptimal2(d = c(qlogistic_best$par, 1.237986), loc = c(-3, 0, -1)) / qlogistic_base1)^(1/3)

