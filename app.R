library(shiny)
library(globpso)

# Define UI for application that draws a histogram
ui <- fluidPage(

  withMathJax(), 
  tabsetPanel(
    selected = "Find Hormesis Exact Optimal Design", 
    type = "tabs", 
    id = "mainpanel", 
    
    tabPanel(
      "Find Hormesis Exact Optimal Design", 
      
      tabsetPanel(
        selected = "Hunt-Bowman D-optimal Exact Design", 
        type = "pills", 
        id = "HBD", 
        tabPanel("Hunt-Bowman D-optimal Exact Design", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("hb_swarm", 
                                  "Number of Swarms", 
                                  value = 64, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("hb_iter", 
                                  "Number of Iterations", 
                                  value = 1000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ), 
                     numericInput("hb_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 50, 
                                  step = 1
                     ), 
                     numericInput("hb_rep", 
                                  "Number of Replications", 
                                  value = 10, 
                                  min = 1, 
                                  max = 50, 
                                  step = 1
                     )
                   ), 
                   mainPanel(
                     tags$h3("Hunt-Bowman Model D-optimal Exact Design"), 
                     tags$p(""),
                     fluidRow(
                       column(8, 
                              "Hunt-Bowman Model Parameters", 
                              numericInput("hb_c1", 
                                           "\\(c_1\\)", 
                                           value = 170, 
                                           min = 1, 
                                           max = Inf), 
                              numericInput("hb_tau", 
                                           "\\(\\tau\\)", 
                                           value = 0.04, 
                                           min = 0, 
                                           max = 0.15), 
                              numericInput("hb_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 1.46, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("hb_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 40, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     actionButton("hb_pso", "Find D-optimal Exact Design"), 
                     verbatimTextOutput("hb_out")
                   )
                 )
        ), 
        
        tabPanel("exp-log \\(\\tau\\)-Optiaml Exact Design", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("eltau_swarm", 
                                  "Number of Swarms", 
                                  value = 128, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("eltau_iter", 
                                  "Number of Iterations", 
                                  value = 2000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ), 
                     numericInput("eltau_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 50, 
                                  step = 1
                     ), 
                     numericInput("eltau_rep", 
                                  "Number of Replications", 
                                  value = 10, 
                                  min = 1, 
                                  max = 50, 
                                  step = 1
                     )
                   ), 
                   mainPanel(
                     tags$h3("exp-log Model \\(\\tau\\)-optimal Exact Design"), 
                     fluidRow(
                       column(8, 
                              "exp-log Model Parameters", 
                              numericInput("eltau_c0", 
                                           "\\(c_0\\)", 
                                           value = 0.15, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("eltau_c1", 
                                           "\\(c_1\\)", 
                                           value = 89, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("eltau_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 3.2, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("eltau_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 41, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     actionButton("eltau_pso", "Find \\(\\tau\\)-optimal Exact Design"), 
                     verbatimTextOutput("eltau_out")
                   )
                 )
        ), 
        tabPanel("exp-log h-Optiaml Exact Design", 
                 sidebarLayout(
                   sidebarPanel(
                     tags$h3("PSO parameters"), 
                     numericInput("elh_swarm", 
                                  "Number of Swarms", 
                                  value = 128, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1), 
                     numericInput("elh_iter", 
                                  "Number of Iterations", 
                                  value = 2000, 
                                  min = 1, 
                                  max = Inf, 
                                  step = 1
                     ), 
                     numericInput("elh_dim", 
                                  "Number of Design Points", 
                                  value = 4, 
                                  min = 1, 
                                  max = 50, 
                                  step = 1
                     ), 
                     numericInput("elh_rep", 
                                  "Number of Replications", 
                                  value = 10, 
                                  min = 1, 
                                  max = 50, 
                                  step = 1
                     )
                   ), 
                   mainPanel(
                     tags$h3("exp-log Model h-optimal Exact Design"), 
                     fluidRow(
                       column(8, 
                              "exp-log Model Parameters", 
                              numericInput("elh_c0", 
                                           "\\(c_0\\)", 
                                           value = 0.15, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("elh_c1", 
                                           "\\(c_1\\)", 
                                           value = 89, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("elh_b0", 
                                           "\\(\\beta_0\\)", 
                                           value = 3.2, 
                                           min = 0, 
                                           max = Inf), 
                              numericInput("elh_b1", 
                                           "\\(\\beta_1\\)", 
                                           value = 41, 
                                           min = 0, 
                                           max = Inf), 
                              
                       )
                     ), 
                     actionButton("elh_pso", "Find h-optimal Exact Design"), 
                     verbatimTextOutput("elh_out")
                   )
                 )
        )
      )
      
    )
  )
  
    
  )


# Define server logic required to draw a histogram
server <- function(input, output) {
  values <- reactiveValues()
  values$hbd <- list(eff = numeric(), points = c(), process = 0)
  values$eltau <- list(eff = numeric(), points = c())
  values$elh <- list(eff = numeric(), points = c())
  
  observeEvent(input$hb_pso, {
    values$hbd$process = 1
    c1 = input$hb_c1
    tau = input$hb_tau
    b0 = input$hb_b0
    b1 = input$hb_b1
    nswarm = input$hb_swarm
    iter = input$hb_iter
    ndp = input$hb_dim
    nrep = input$hb_rep
    
    hb_par <- hb_parms(c1 = c1, tau = tau, b0 = b0, b1 = b1)
    psoinfo_hb <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    hb_res <- hb_pso_rep(nRep = nrep, nPoints = ndp, parms = hb_par, psoinfo = psoinfo_hb)
    
    values$hbd$eff = hb_res$result$best_eff
    values$hbd$points = hb_res$result$best_points
    values$hbd$process = 2
  }
               )
  
  output$hb_out <- renderPrint({
    cat("Efficiency:", values$hbd$eff, 
        "\nDesign Points:", values$hbd$points, "\n")
  })
  
  
  observeEvent(input$eltau_pso, {
    c0 = input$eltau_c0
    c1 = input$eltau_c1
    b0 = input$eltau_b0
    b1 = input$eltau_b1
    nswarm = input$eltau_swarm
    iter = input$eltau_iter
    ndp = input$eltau_dim
    nrep = input$eltau_rep
    
    el_par <- exp_log_params(c0 = c0, c1 = c1, b0 = b0, b1 = b1)
    psoinfo_el <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    eltau_res <- exp_log_tau_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    
    values$eltau$eff = eltau_res$result$best_eff
    values$eltau$points = eltau_res$result$best_points
  }
  )
  
  output$eltau_out <- renderPrint({
    cat("Efficiency:", values$eltau$eff, 
        "\nDesign Points:", values$eltau$points, "\n")
  })
  
  observeEvent(input$elh_pso, {
    c0 = input$elh_c0
    c1 = input$elh_c1
    b0 = input$elh_b0
    b1 = input$elh_b1
    nswarm = input$elh_swarm
    iter = input$elh_iter
    ndp = input$elh_dim
    nrep = input$elh_rep
    
    el_par <- exp_log_params(c0 = c0, c1 = c1, b0 = b0, b1 = b1)
    psoinfo_el <- psoinfo_setting(nSwarms = nswarm, Iters = iter)
    elh_res <- exp_log_h_pso_rep(nRep = nrep, nPoints = ndp, parms = el_par, psoinfo = psoinfo_el)
    
    values$elh$eff = elh_res$result$best_eff
    values$elh$points = elh_res$result$best_points
  }
  )
  
  output$elh_out <- renderPrint({
    cat("Efficiency:", values$elh$eff, 
        "\nDesign Points:", values$elh$points, "\n")
  })
}



## Initalize PSO settings

psoinfo_setting <- function(nSwarms = 64, Iters = 1000){
  getPSOInfo(nSwarm = nSwarms, maxIter = Iters, w0 = 0.9, w1 = 0.4, w_var = 1)
}

## Hunt-Bowman model

# Hunt-Bowman parameters
hb_parms <- function(c1, tau, b0, b1){
  c(c1, tau, b0, b1)
}

# Hunt-Bowman information matrix of 1 point
hb_mat <- function(d, c1, tau, b0, b1){
  
  if (d <= tau){
    f1 <- d^2 - tau * d
    f2 <- -c1 * d
    f3 <- -exp(b0) / (1 + exp(b0))^2
    f4 <- 0
    f <- matrix(c(f1, f2, f3, f4))
    mat <- f %*% t(f)
  } 
  else if (d > tau){
    f1 <- 0
    f2 <- -b1 * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f3 <- -exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f4 <- (d - tau) * exp(b0 - b1 * (d - tau)) / (1 + exp(b0 - b1 * (d - tau)))^2
    f <- matrix(c(f1, f2, f3, f4))
    mat <- f %*% t(f)
  }
  
  mat
}

# Hunt-Bowman d-optimality criterion value efficiency
hb_doptimal <- function(d, loc){
  
  # Hunt-Bowman parameters
  c1 <- loc[1]
  tau <- loc[2]
  b0 <- loc[3]
  b1 <- loc[4]
  eff_base <- 9.938331e-14 #Criterion value of 4-point uniform design
  
  # Number of experiment points
  d <- c(0, d)
  w <- length(d)
  
  # Evaluate d-optimality criterion value
  mat_list <- lapply(d, function(x) 1/w * hb_mat(x, c1, tau, b0, b1))
  inf_mat <- Reduce("+", mat_list)
  - (det(inf_mat) / eff_base)
}

# Find the exact design points that minimizes the d-optimality efficiency
# of hunt-bowman model by PSO
hb_pso <- function(nPoints, parms, psoinfo){
  
  # set the lower and upper bounds for PSO
  lb <- rep(0, nPoints-1)
  ub <- rep(0.15, nPoints-1)
  
  pso_res <- globpso(objFunc = hb_doptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- (-1 * pso_res$val)^(1/nPoints) |> round(4)
  pso_res$par <- c(0, pso_res$par) |> sort() |> round(4)
  pso_res$history <- pso_res$history * -1
  pso_res
}

# Replicate m PSO results of Hunt-Bowman model d-optimality
hb_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- hb_pso(nPoints, parms, psoinfo))
  
  hb_list <- list(efficiency = c(), design_points = c(), history = c(), 
                  result = list(best_eff = 0, best_points = c()))
  
  # Efficiency of each replication.
  hb_list$efficiency <- sapply(1:nRep, function(x) pso_results[[x]]$val) 
  # Design points found in each replication.
  hb_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par) 
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  hb_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history) 
  
  best_idx <- which.max(hb_list$efficiency)
  # Best efficiency among all replications.
  hb_list$result$best_eff <- hb_list$efficiency[best_idx]
  # Design points correspond to the best efficiency.
  hb_list$result$best_points <- hb_list$design_points[,best_idx]
  
  hb_list
}

## exp-log model

# exp-log model parameters
exp_log_params <- function(c0, c1, b0, b1){
  c(c0, c1, b0, b1)
}

# exp-log model information matrix
exp_log_f <- function(d, c0, c1, b0, b1){
  f1 <- exp(-c1 * d)
  f2 <- -c0 * d * exp(-c1 * d)
  f3 <- -exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f4 <- d * exp(b0 - b1 * d) / (1 + exp(b0 - b1 * d))^2
  f <- matrix(c(f1, f2, f3, f4))
  f %*% t(f)
}

exp_log_mat <- function(d, c0, c1, b0, b1){
  #d <- c(0, d)
  n <- length(d)
  mat_list <- lapply(d, function(x) 1/n * exp_log_f(x, c0, c1, b0, b1))
  
  #print(mat_list)
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  M
}

# exp-log model h-optimality criterion value
exp_log_hoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  pen = 999999999 # Penalty value
  eff_base = 243738.4 # Criterion value of 4-point locally optimal design
  
  n <- length(d)
  h1 <- -c1
  h2 <- -c0
  h3 <- b1 * exp(b0) * (1 - exp(b0)) / (exp(b0) + 1)^3
  h4 <- exp(b0) / (exp(b0) + 1)^2 
  h <- matrix(c(h1, h2, h3, h4))
  
  M <- exp_log_mat(d, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    res = t(h) %*% M_inv %*% h
  }
  
  res / eff_base
}

# Find the exact design points that minimizes the h-optimality efficiency
# of exp-log model by PSO
exp_log_pso <- function(nPoints, parms, psoinfo){
  
  lb <- rep(0, nPoints)
  ub <- rep(0.15, nPoints)
  
  pso_res <- globpso(objFunc = exp_log_hoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = parms)
  
  pso_res$val <- (1 / pso_res$val)^(1/nPoints)  |> round(4)
  pso_res$par <- c(pso_res$par) |> sort() |> round(4)
  pso_res$history <- 1 / pso_res$history
  pso_res
}

# Replicate m PSO results of exp-log model h-optimality
exp_log_h_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(efficiency = c(), design_points = c(), history = c(), 
                       result = list(best_eff = 0, best_points = c()))
  
  # Efficiency of each replication.
  exp_log_list$efficiency <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  exp_log_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.max(exp_log_list$efficiency)
  # Best efficiency among all replications.
  exp_log_list$result$best_eff <- exp_log_list$efficiency[best_idx]
  # Design points correspond to the best efficiency.
  exp_log_list$result$best_points <- exp_log_list$design_points[,best_idx]
  
  exp_log_list
}


## tau-optimal 

tau_func <- function(d, c0, c1, b0, b1){
  c0 * exp(-c1*d) + 1 / (1 + exp(b0 - b1*d)) - (c0 + 1 / (1 + exp(b0)))
}

exp_log_b <- function(d, c0, c1, b0, b1){
  h <- ((b1 * exp(b1*d + b0)) / (exp(b1*d) + exp(b0))^2) - c0 * c1 * exp(-c1*d)
  h1 <- -1*(exp(-c1*d) - 1) / h
  h2 <- -1*(-c0 * d * exp(-c1*d)) / h
  h3 <- -1*(exp(b0) * (1 / (exp(b0) + 1)^2 - exp(b1*d) / (exp(b1*d) + exp(b0))^2 )) / h
  h4 <- -1*(d * exp(b1*d + b0) / (exp(b1 * d) + exp(b0))^2) / h
  
  b <- matrix(c(h1, h2, h3, h4))
}

exp_log_loc_mat <- function(d, w, c0, c1, b0, b1){
  n <- length(d)
  mat_list <- lapply(1:n, function(x) w[x] * exp_log_f(d[x], c0, c1, b0, b1))
  
  #print(mat_list)
  M <- Reduce("+", mat_list)
  diag(M) <- diag(M) + 1e-10
  M
}

exp_log_loc_tauoptimal <- function(x, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  tau = loc[5]
  pen = 999999999
  
  n <- round(length(x) / 2)
  d <- x[1:n]
  w <- x[-(1:n)]
  w[n] = 1 - sum(w)
  
  #tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, c0 = c0, c1 = c1, b0 = b0, b1 = b1)$root
  #print(tau)
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  M <- exp_log_loc_mat(d, w, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  } 
  else if (w[n] <= 0){
    res = pen
  }
  else {
    M_inv <- solve(M)
    #print(M_inv)
    #print(b)
    res = t(b) %*% M_inv %*% b
  }
  
  res
}

exp_log_loc_tau_pso <- function(nPoints, parms, psoinfo){
  
  lb <- c(rep(0, nPoints), rep(0, nPoints-1))
  ub <- c(rep(0.15, nPoints), rep(1, nPoints-1))
  tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                 c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
  
  pso_res <- globpso(objFunc = exp_log_loc_tauoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = c(parms, tau))
  
  pso_res$val <- pso_res$val
  pso_res$points <- pso_res$par[1:nPoints] |> round(4)
  pso_res$weight <- c(pso_res$par[-c(1:nPoints)], 1 - sum(pso_res$par[-c(1:nPoints)])) |> round(4)
  pso_res
}

exp_log_loc_tau_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_loc_tau_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(design_points = c(), weight = c(), val = c())
  
  # Efficiency of each replication.
  exp_log_list$weight <- sapply(1:nRep, function(x) pso_results[[x]]$weight)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$points)
  exp_log_list$val <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  
  exp_log_list
}

exp_log_tauoptimal <- function(d, loc){
  c0 = loc[1]
  c1 = loc[2]
  b0 = loc[3]
  b1 = loc[4]
  tau = loc[5]
  pen = 999999999
  
  n <- length(d)
  w <- rep(1/n, n)
  
  b <- exp_log_b(tau, c0, c1, b0, b1)
  
  M <- exp_log_loc_mat(d, w, c0, c1, b0, b1)
  if (rcond(M) < 2.220446e-16){
    res = pen
  }
  else {
    M_inv <- solve(M)
    #print(M_inv)
    #print(b)
    res = t(b) %*% M_inv %*% b
  }
  
  #print(tau)
  res / 0.1157493
}

exp_log_tau_pso <- function(nPoints, parms, psoinfo){
  
  lb <- c(rep(0, nPoints))
  ub <- c(rep(0.15, nPoints))
  tau <- uniroot(tau_func, c(0.00001, 0.15), tol = 1e-10, 
                 c0 = parms[1], c1 = parms[2], b0 = parms[3], b1 = parms[4])$root
  
  pso_res <- globpso(objFunc = exp_log_tauoptimal, 
                     lower = lb, upper = ub, 
                     PSO_INFO = psoinfo, verbose = F, 
                     loc = c(parms, tau))
  
  pso_res$val <- (1 / pso_res$val)^(1/nPoints) |> round(4)
  pso_res$par <- c(pso_res$par) |> sort() |> round(4)
  pso_res$history <- 1 / pso_res$history
  pso_res
}

exp_log_tau_pso_rep <- function(nRep, nPoints = 4, parms, psoinfo){
  
  # Create nRep replicates of nPoints exact design result found by PSO.
  pso_results <- list()
  pso_results <- lapply(1:nRep, function(x) pso_results[[x]] <- exp_log_tau_pso(nPoints, parms, psoinfo))
  
  exp_log_list <- list(design_points = c(), weight = c(), val = c())
  
  # Efficiency of each replication.
  exp_log_list$efficiency <- sapply(1:nRep, function(x) pso_results[[x]]$val)
  # Design points found in each replication.
  exp_log_list$design_points <- sapply(1:nRep, function(x) pso_results[[x]]$par)
  # Best result of each iteration in each replication. (For finding the appropriate nRep.)
  exp_log_list$history <- sapply(1:nRep, function(x) pso_results[[x]]$history)
  
  best_idx <- which.max(exp_log_list$efficiency)
  # Best efficiency among all replications.
  exp_log_list$result$best_eff <- exp_log_list$efficiency[best_idx]
  # Design points correspond to the best efficiency.
  exp_log_list$result$best_points <- exp_log_list$design_points[,best_idx]
  
  exp_log_list
}

# Run the application 
shinyApp(ui = ui, server = server)
