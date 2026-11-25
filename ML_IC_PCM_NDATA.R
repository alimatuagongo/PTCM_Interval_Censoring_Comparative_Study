ptm <- proc.time()
library(e1071)
library(xtable)
library(survival)
library(dplyr)
library(writexl)
library(readxl)
library(PRROC)
library(MASS)
library(caTools)
library(splines)
library(rpart)
library(randomForest)
library(xgboost)
library(neuralnet)
library(pROC)
library(thregI)


# Creating interval censoring
LR_int=function(y1,len1,l1){ # function to find the left and right end points that contain the event time, provided the event takes place
  if (len1 < 0) stop("len1 must be positive.")
  if (y1 < 0) stop("y1 must be positive.")
  if(y1>0 & y1<=l1){
    a=c(.Machine$double.eps,l1)
  }else{
    k=as.integer((y1-l1)/len1)+1 # i changed from as.integer to round
    a=c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}


#LOGIT
em.LOGIT.ic <- function(TimeL.train, TimeR.train, Status.train, X.train, Z.train,TimeL.test, TimeR.test,X.test, Z.test,b, alpha, beta,emmax = 1000,eps = 1e-3) {
  
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  
  while (convergence > eps && iter < emmax) {
    
    ## ===================================================
    ## Current incidence estimate
    ## Z includes the intercept
    ## ===================================================
    
    eta_inc.train <- drop(Z.train %*% b)
    uncureprob.train <- plogis(eta_inc.train)
    uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
    
    ## ===================================================
    ## Current Weibull latency survival
    ## X includes the intercept in column 1
    ## beta corresponds to X without the intercept
    ## ===================================================
    
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    ## ===================================================
    ## E-step
    ## ===================================================
    
    M.train <- Status.train - survival.train * log1p(-uncureprob.train)
    
    ## ===================================================
    ## M-step: incidence/LOGIT component
    ## ===================================================
    Q1 <- function(par) {
      eta1 <- drop(Z.train %*% par)
      uncure1 <- plogis(eta1)
      uncure1 <- pmin(pmax(uncure1, 1e-6), 1 - 1e-6)
      theta1 <- -log1p(-uncure1)
      theta1 <- pmax(theta1, 1e-10)
      out <- sum(M.train * log(theta1)) + sum(log1p(-uncure1))
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_b <- optim(par = b, fn = Q1, method = "Nelder-Mead")$par
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    #convergence <- sum(c(update_b - b, update_latency - c(alpha, beta))^2)
    convergence <- sum((update_b - b)^2) + sum((update_latency - c(alpha, beta))^2) 
    
    ## Update parameters
    b <- update_b
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_inc.train <- drop(Z.train %*% b)
  uncureprob.train <- plogis(eta_inc.train)
  uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob.train)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob.train)) / uncureprob.train
  ## =====================================================
  ## Test predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_inc.test <- drop(Z.test %*% b)
  uncureprob.test <- plogis(eta_inc.test)
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    b = b,
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob.train, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter
  ))
  
}
smcure.LOGIT.ic <- function(formula, cureform,alpha.init, beta.init,train_data, test_data, na.action = na.omit,Var = F,emmax = 1000, eps = 1e-3,nboot = 1) {
  
  
  ## =====================================================
  ## Prepare data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nb <- ncol(Z.train)
  nbeta <- ncol(X.train) - 1
  ## Check beta.init length
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  ## Initial incidence coefficients
  b <- rep(0.5, nb)
  ## =====================================================
  ## Fit EM
  ## =====================================================
  emfit <- em.LOGIT.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    b = b,
    alpha = alpha.init,
    beta = beta.init,
    emmax = emmax,
    eps = eps
  )
  latency <- emfit$latencyfit
  b.est <- emfit$b
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap variance for latency parameters
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    ## ---------------------------------------------------
    ## Store the ORIGINAL row-level quantities
    ## We bootstrap row indices, then subset X and Z
    ## directly. This keeps the code universal.
    ## ---------------------------------------------------
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) {
      stop("Bootstrap requires both event and censored observations.")
    }
    boot_count <- 1
    while (boot_count <= nboot) {
      ## -------------------------------------------------
      ## Stratified bootstrap:
      ## resample D = 1 and D = 0 separately
      ## -------------------------------------------------
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      ## -------------------------------------------------
      ## Bootstrap survival data
      ## -------------------------------------------------
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      ## -------------------------------------------------
      ## Bootstrap generic X and Z matrices
      ## -------------------------------------------------
      bootX <- X.train[boot_id, , drop = FALSE]
      bootZ <- Z.train[boot_id, , drop = FALSE]
      ## -------------------------------------------------
      ## Bootstrap EM fit
      ## -------------------------------------------------
      bootfit <- tryCatch(
        em.LOGIT.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          b = b.est,
          alpha = alpha.est,
          beta = beta.est,
          emmax = emmax,
          eps = eps
        ),
        error = function(e) {
          NULL
        }
      )
      ## -------------------------------------------------
      ## Keep only converged bootstrap fits
      ## -------------------------------------------------
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    ## ===================================================
    ## Bootstrap variance and SD
    ## ===================================================
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.LOGIT.ic"
  fit$call <- call
  fit$b <- b.est
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  ## Store names automatically
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  # cat(" done.\n")
  return(fit)
  
  
}

#SPLINE
em.SPLINE.ic <- function(TimeL.train, TimeR.train, Status.train,X.train, Z.train,TimeL.test, TimeR.test,X.test, Z.test,b, alpha, beta,emmax = 1000, eps = 1e-3) {
  
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## Current incidence estimate
    ## Z is the SPLINE-expanded incidence design matrix
    ## and includes the intercept.
    ## ===================================================
    eta_inc.train <- drop(Z.train %*% b)
    uncureprob.train <- plogis(eta_inc.train)
    uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## Current Weibull latency survival
    ## X includes the intercept in column 1.
    ## beta corresponds to X without the intercept.
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    ## ===================================================
    ## E-step
    ## ===================================================
    M.train <- Status.train - survival.train * log1p(-uncureprob.train)
    ## ===================================================
    ## M-step: SPLINE incidence component
    ## ===================================================
    Q1 <- function(par) {
      eta1 <- drop(Z.train %*% par)
      uncure1 <- plogis(eta1)
      uncure1 <- pmin(pmax(uncure1, 1e-6), 1 - 1e-6)
      theta1 <- -log1p(-uncure1)
      theta1 <- pmax(theta1, 1e-10)
      out <- sum(M.train * log(theta1)) + sum(log1p(-uncure1))
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_b <- optim(par = b, fn = Q1, method = "Nelder-Mead")$par
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    #convergence <- sum(c(update_b - b, update_latency - c(alpha, beta))^2)
    convergence <- sum((update_b - b)^2) + sum((update_latency - c(alpha, beta))^2)  
    
    ## ===================================================
    ## Update parameters
    ## ===================================================
    b <- update_b
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_inc.train <- drop(Z.train %*% b)
  uncureprob.train <- plogis(eta_inc.train)
  uncureprob.train <- pmin(pmax(uncureprob.train, 1e-6), 1 - 1e-6)
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob.train)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob.train)) / uncureprob.train
  ## =====================================================
  ## Test predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_inc.test <- drop(Z.test %*% b)
  uncureprob.test <- plogis(eta_inc.test)
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    b = b,
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob.train, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter
  ))
}

smcure.SPLINE.ic <- function(formula, cureform, alpha.init, beta.init,train_data, test_data, df_spline = 2,na.action = na.omit,Var = TRUE,emmax = 1000,eps = 1e-3,nboot = 100) {
  call <- match.call()
  
  ## =====================================================
  ## Prepare data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Latency model: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Incidence variables
  ## =====================================================
  cvars <- all.vars(cureform)
  ## =====================================================
  ## Build SPLINE incidence matrices: Z
  ## Numeric variables with > 2 unique values:
  ## natural SPLINE
  ## Binary / categorical variables:
  ## ordinary model-matrix terms
  ## =====================================================
  Z.train.parts <- list()
  Z.test.parts <- list()
  spline_bases <- list()
  for (v in cvars) {
    x.train <- train_data[[v]]
    x.test <- test_data[[v]]
    ## -----------------------------------------------------
    ## Continuous variable -> natural spline
    ## -----------------------------------------------------
    if (is.numeric(x.train) && length(unique(x.train)) > 2) {
      basis.train <- splines::ns(x.train, df = df_spline)
      basis.test <- predict(basis.train, x.test)
      colnames(basis.train) <- paste0(v, "_s", seq_len(ncol(basis.train)))
      colnames(basis.test) <- colnames(basis.train)
      Z.train.parts[[v]] <- basis.train
      Z.test.parts[[v]] <- basis.test
      spline_bases[[v]] <- basis.train
    } else {
      ## ---------------------------------------------------
      ## Binary / categorical variable -> ordinary terms
      ## ---------------------------------------------------
      temp_formula <- as.formula(paste0("~ ", v))
      temp.train <- model.matrix(temp_formula, data = train_data)[, -1, drop = FALSE]
      temp.test <- model.matrix(temp_formula, data = test_data)[, -1, drop = FALSE]
      Z.train.parts[[v]] <- temp.train
      Z.test.parts[[v]] <- temp.test
    }
  }
  Z.train <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, Z.train.parts)))
  Z.test <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, Z.test.parts)))
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nb <- ncol(Z.train)
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  b <- rep(0.5, nb)
  ## =====================================================
  ## Fit EM
  ## =====================================================
  emfit <- em.SPLINE.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    b = b,
    alpha = alpha.init,
    beta = beta.init,
    emmax = emmax,
    eps = eps
  )
  latency <- emfit$latencyfit
  b.est <- emfit$b
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) {
      stop("Bootstrap requires both event and censored observations.")
    }
    boot_count <- 1
    while (boot_count <= nboot) {
      ## ===================================================
      ## Stratified bootstrap
      ## ===================================================
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      bootdata <- train_data[boot_id, , drop = FALSE]
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      ## ===================================================
      ## Bootstrap latency matrix
      ## ===================================================
      boot_mf <- model.frame(formula, data = bootdata)
      bootX <- model.matrix(attr(boot_mf, "terms"), boot_mf)
      ## ===================================================
      ## Bootstrap SPLINE incidence matrices
      ## ===================================================
      boot.Z.parts <- list()
      test.Z.parts <- list()
      for (v in cvars) {
        x.boot <- bootdata[[v]]
        x.test <- test_data[[v]]
        if (is.numeric(x.boot) && length(unique(x.boot)) > 2) {
          basis.boot <- splines::ns(x.boot, df = df_spline)
          basis.test <- predict(basis.boot, x.test)
          colnames(basis.boot) <- paste0(v, "_s", seq_len(ncol(basis.boot)))
          colnames(basis.test) <- colnames(basis.boot)
          boot.Z.parts[[v]] <- basis.boot
          test.Z.parts[[v]] <- basis.test
        } else {
          temp_formula <- as.formula(paste0("~ ", v))
          temp.boot <- model.matrix(temp_formula, data = bootdata)[, -1, drop = FALSE]
          temp.test <- model.matrix(temp_formula, data = test_data)[, -1, drop = FALSE]
          boot.Z.parts[[v]] <- temp.boot
          test.Z.parts[[v]] <- temp.test
        }
      }
      bootZ <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, boot.Z.parts)))
      bootZ.test <- as.matrix(cbind("(Intercept)" = 1, do.call(cbind, test.Z.parts)))
      ## ===================================================
      ## Bootstrap EM
      ## ===================================================
      bootfit <- tryCatch(
        em.SPLINE.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = bootZ.test,
          b = b.est,
          alpha = alpha.est,
          beta = beta.est,
          emmax = emmax,
          eps = eps
        ),
        error = function(e) NULL
      )
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.SPLINE.ic"
  fit$call <- call
  fit$b <- b.est
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$df_spline <- df_spline
  fit$cure_covariates <- cvars
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  # cat(" done.\n")
  return(fit)
}


#SUPPORT VECTORE MACHINE
em.SVM.ic <- function(TimeL.train, TimeR.train, Status.train,X.train, Z.train,TimeL.test, TimeR.test, X.test, Z.test,uncureprob, alpha, beta,best_gamma, best_cost,emmax = 1000,eps = 1e-3, M.imp = 5) {
  
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  n <- length(Status.train)
  m <- nrow(Z.test)
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  final_models <- vector("list", M.imp)
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    ## Keep incidence probabilities away from 0 and 1
    uncureprob <- pmin(pmax(uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## E-step
    ## ===================================================
    w <- Status.train + (1 - Status.train) * (1 - (1 - uncureprob)^survival.train)
    w <- pmin(pmax(w, 1e-6), 1 - 1e-6)
    M.train <- Status.train - survival.train * log1p(-uncureprob)
    ## ===================================================
    ## Posterior-based class weights
    ## ===================================================
    n_uncured_eff <- sum(w)
    n_cured_eff <- sum(1 - w)
    class_weights <- c("-1" = n / (2 * n_cured_eff), "1" = n / (2 * n_uncured_eff))
    ## ===================================================
    ## M-step: SVM incidence component
    ## ===================================================
    update_probs <- matrix(NA_real_, nrow = n, ncol = M.imp)
    current_models <- vector("list", M.imp)
    for (j in seq_len(M.imp)) {
      ## -------------------------------------------------
      ## Impute latent uncure indicator
      ## -------------------------------------------------
      imputed_label <- rbinom(n = n, size = 1, prob = w)
      yj <- factor(2 * imputed_label - 1, levels = c(-1, 1))
      if (length(unique(yj)) < 2) {
        update_probs[, j] <- uncureprob
        next
      }
      ## -------------------------------------------------
      ## Fit SVM using generic Z covariates
      ## -------------------------------------------------
      modj <- tryCatch(
        e1071::svm(
          x = Z.train[, -1, drop = FALSE],
          y = yj,
          type = "C-classification",
          kernel = "radial",
          gamma = best_gamma,
          cost = best_cost,
          class.weights = class_weights,
          probability = TRUE,
          scale = FALSE
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(modj)) {
        update_probs[, j] <- uncureprob
        next
      }
      ## -------------------------------------------------
      ## Training incidence probabilities
      ## -------------------------------------------------
      predj <- tryCatch(
        predict(modj, Z.train[, -1, drop = FALSE], probability = TRUE),
        error = function(e) {
          NULL
        }
      )
      if (is.null(predj)) {
        update_probs[, j] <- uncureprob
        next
      }
      probj <- attr(predj, "probabilities")
      if (!is.null(probj) && "1" %in% colnames(probj)) {
        update_probs[, j] <- probj[, "1"]
      } else {
        update_probs[, j] <- uncureprob
      }
      current_models[[j]] <- modj
    }
    ## ---------------------------------------------------
    ## Average incidence probability over imputations
    ## ---------------------------------------------------
    update_uncureprob <- rowMeans(update_probs, na.rm = TRUE)
    bad_prob <- !is.finite(update_uncureprob)
    if (any(bad_prob)) {
      update_uncureprob[bad_prob] <- uncureprob[bad_prob]
    }
    update_uncureprob <- pmin(pmax(update_uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    #convergence <- sum(c(mean(update_uncureprob) - mean(uncureprob), update_latency - c(alpha, beta))^2)
    convergence <- mean((update_uncureprob - uncureprob)^2) + sum((update_latency - c(alpha, beta))^2)
    
    ## ===================================================
    ## Update parameters
    ## ===================================================
    uncureprob <- update_uncureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    ## Keep models from the CURRENT EM iteration
    final_models <- current_models
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob)) / uncureprob
  ## =====================================================
  ## Test incidence predictions
  ## =====================================================
  test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
  for (j in seq_len(M.imp)) {
    if (is.null(final_models[[j]])) {
      next
    }
    pred.test <- tryCatch(
      predict(final_models[[j]], Z.test[, -1, drop = FALSE], probability = TRUE),
      error = function(e) {
        NULL
      }
    )
    if (is.null(pred.test)) {
      next
    }
    prob.test <- attr(pred.test, "probabilities")
    if (!is.null(prob.test) && "1" %in% colnames(prob.test)) {
      test_probs[, j] <- prob.test[, "1"]
    }
  }
  if (all(is.na(test_probs))) {
    stop("No valid SVM model available for test prediction.")
  }
  uncureprob.test <- rowMeans(test_probs, na.rm = TRUE)
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Test latency predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter,
    models = final_models,
    best_gamma = best_gamma,
    best_cost = best_cost
  ))
}
smcure.SVM.ic <- function(formula, cureform,
                          alpha.init, beta.init,
                          train_data, test_data,
                          na.action = na.omit,
                          Var = FALSE,
                          emmax = 1000,
                          eps = 1e-3,
                          nboot = 100,
                          gamma_grid = 2^(-4:-6),
                          cost_grid = 2^(4:6),
                          M.imp = 5) {
  call <- match.call()
  
  
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  ## =====================================================
  ## Initial SVM uncure probability
  ## =====================================================
  init_label <- factor(2 * Status.train - 1, levels = c(-1, 1))
  if (length(unique(init_label)) < 2) {
    stop("Initial SVM labels contain only one class. Cannot fit SVM.")
  }
  ## =====================================================
  ## Tune gamma and cost
  ## Z[, -1] is used because SVM does not need
  ## an explicit intercept column.
  ## =====================================================
  tune_dat <- data.frame(y = init_label, Z.train[, -1, drop = FALSE])
  obj <- e1071::tune(
    e1071::svm,
    y ~ .,
    data = tune_dat,
    ranges = list(gamma = gamma_grid, cost = cost_grid),
    tunecontrol = e1071::tune.control(sampling = "fix"),
    type = "C-classification",
    kernel = "radial",
    probability = TRUE,
    scale = FALSE
  )
  best_gamma.init <- obj$best.parameters$gamma
  best_cost.init <- obj$best.parameters$cost
  ## =====================================================
  ## Initial SVM fit
  ## =====================================================
  mod.init <- e1071::svm(
    x = Z.train[, -1, drop = FALSE],
    y = init_label,
    type = "C-classification",
    kernel = "radial",
    gamma = best_gamma.init,
    cost = best_cost.init,
    probability = TRUE,
    scale = FALSE
  )
  pred.init <- predict(mod.init, Z.train[, -1, drop = FALSE], probability = TRUE)
  proba.init <- attr(pred.init, "probabilities")
  if (is.null(proba.init) || !("1" %in% colnames(proba.init))) {
    stop("Initial SVM did not return probability for uncured class 1.")
  }
  uncureprob.init <- proba.init[, "1"]
  uncureprob.init <- pmin(pmax(uncureprob.init, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Run EM algorithm
  ## =====================================================
  emfit <- em.SVM.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    best_gamma = best_gamma.init,
    best_cost = best_cost.init,
    emmax = emmax,
    eps = eps,
    M.imp = M.imp
  )
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    ## ---------------------------------------------------
    ## Stratified bootstrap indices
    ## ---------------------------------------------------
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) {
      stop("Bootstrap requires both event and censored observations.")
    }
    boot_count <- 1
    while (boot_count <= nboot) {
      ## -------------------------------------------------
      ## Resample event and censored observations
      ## separately
      ## -------------------------------------------------
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      ## -------------------------------------------------
      ## Generic bootstrap quantities
      ## -------------------------------------------------
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      bootX <- X.train[boot_id, , drop = FALSE]
      bootZ <- Z.train[boot_id, , drop = FALSE]
      ## =================================================
      ## Initial SVM for bootstrap sample
      ## =================================================
      boot_label <- factor(2 * boot_Status - 1, levels = c(-1, 1))
      if (length(unique(boot_label)) < 2) {
        next
      }
      boot_mod.init <- tryCatch(
        e1071::svm(
          x = bootZ[, -1, drop = FALSE],
          y = boot_label,
          type = "C-classification",
          kernel = "radial",
          gamma = best_gamma.init,
          cost = best_cost.init,
          probability = TRUE,
          scale = TRUE
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_mod.init)) {
        next
      }
      boot_pred.init <- tryCatch(
        predict(boot_mod.init, bootZ[, -1, drop = FALSE], probability = TRUE),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_pred.init)) {
        next
      }
      boot_proba.init <- attr(boot_pred.init, "probabilities")
      if (is.null(boot_proba.init) || !("1" %in% colnames(boot_proba.init))) {
        next
      }
      boot_uncure.init <- boot_proba.init[, "1"]
      boot_uncure.init <- pmin(pmax(boot_uncure.init, 1e-6), 1 - 1e-6)
      ## =================================================
      ## Bootstrap EM
      ## =================================================
      bootfit <- tryCatch(
        em.SVM.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          best_gamma = best_gamma.init,
          best_cost = best_cost.init,
          emmax = emmax,
          eps = eps,
          M.imp = M.imp
        ),
        error = function(e) {
          NULL
        }
      )
      ## =================================================
      ## Keep converged bootstrap fits
      ## =================================================
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    ## ===================================================
    ## Bootstrap variance and SD
    ## ===================================================
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.SVM.ic"
  fit$call <- call
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$models <- emfit$models
  fit$best_gamma <- emfit$best_gamma
  fit$best_cost <- emfit$best_cost
  fit$initial_model <- mod.init
  fit$initial_gamma <- best_gamma.init
  fit$initial_cost <- best_cost.init
  ## =====================================================
  ## Store generic covariate names
  ## =====================================================
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  # cat(" done.\n")
  return(fit)
}



#DECISION TREE
em.DT.ic <- function(TimeL.train, TimeR.train, Status.train,
                     X.train, Z.train,
                     TimeL.test, TimeR.test,
                     X.test, Z.test,
                     uncureprob, alpha, beta,
                     best_cp, best_minsplit,
                     emmax = 1000, eps = 1e-3, M.imp = 5) {
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  n <- length(Status.train)
  m <- nrow(Z.test)
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  final_models <- vector("list", M.imp)
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## Current Weibull latency survival
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    uncureprob <- pmin(pmax(uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## E-step
    ## ===================================================
    w <- Status.train + (1 - Status.train) * (1 - (1 - uncureprob)^survival.train)
    w <- pmin(pmax(w, 1e-6), 1 - 1e-6)
    M.train <- Status.train - survival.train * log1p(-uncureprob)
    ## ===================================================
    ## M-step: DT incidence component
    ## ===================================================
    update_probs <- matrix(NA_real_, nrow = n, ncol = M.imp)
    current_models <- vector("list", M.imp)
    for (j in seq_len(M.imp)) {
      imputed_label <- rbinom(n = n, size = 1, prob = w)
      yj <- factor(imputed_label, levels = c(0, 1))
      if (length(unique(yj)) < 2) {
        update_probs[, j] <- uncureprob
        next
      }
      datj <- data.frame(y = yj, Z.train[, -1, drop = FALSE])
      modj <- tryCatch(
        rpart::rpart(
          y ~ .,
          data = datj,
          method = "class",
          control = rpart::rpart.control(cp = best_cp, minsplit = best_minsplit, xval = 0)
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(modj)) {
        update_probs[, j] <- uncureprob
        next
      }
      probj <- tryCatch(
        predict(modj, newdata = data.frame(Z.train[, -1, drop = FALSE]), type = "prob"),
        error = function(e) {
          NULL
        }
      )
      if (is.null(probj)) {
        update_probs[, j] <- uncureprob
        next
      }
      if ("1" %in% colnames(probj)) {
        update_probs[, j] <- probj[, "1"]
      } else {
        update_probs[, j] <- uncureprob
      }
      current_models[[j]] <- modj
    }
    update_uncureprob <- rowMeans(update_probs, na.rm = TRUE)
    bad_prob <- !is.finite(update_uncureprob)
    if (any(bad_prob)) {
      update_uncureprob[bad_prob] <- uncureprob[bad_prob]
    }
    update_uncureprob <- pmin(pmax(update_uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    #convergence <- sum(c(mean(update_uncureprob) - mean(uncureprob), update_latency - c(alpha, beta))^2)
    convergence <- mean((update_uncureprob - uncureprob)^2) + sum((update_latency - c(alpha, beta))^2)
    
    uncureprob <- update_uncureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    final_models <- current_models
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob)) / uncureprob
  ## =====================================================
  ## Test incidence predictions
  ## =====================================================
  test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
  for (j in seq_len(M.imp)) {
    if (is.null(final_models[[j]])) {
      next
    }
    prob.test <- tryCatch(
      predict(final_models[[j]], newdata = data.frame(Z.test[, -1, drop = FALSE]), type = "prob"),
      error = function(e) {
        NULL
      }
    )
    if (is.null(prob.test)) {
      next
    }
    if ("1" %in% colnames(prob.test)) {
      test_probs[, j] <- prob.test[, "1"]
    }
  }
  if (all(is.na(test_probs))) {
    stop("No valid decision-tree model available for test prediction.")
  }
  uncureprob.test <- rowMeans(test_probs, na.rm = TRUE)
  bad_test <- !is.finite(uncureprob.test)
  if (any(bad_test)) {
    uncureprob.test[bad_test] <- mean(uncureprob)
  }
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Test survival predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter,
    models = final_models,
    best_cp = best_cp,
    best_minsplit = best_minsplit
  ))
}

smcure.DT.ic <- function(formula, cureform,
                         alpha.init, beta.init,
                         train_data, test_data,
                         na.action = na.omit,
                         Var = FALSE,
                         emmax = 1000,
                         eps = 1e-3,
                         nboot = 100,
                         cp_grid = c( 0.001, 0.0025, 0.005, 0.01, 0.02),
                         minsplit_grid = c(2, 3, 5, 10, 15, 20, 30, 40),
                         M.imp = 5) {
  call <- match.call()
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  ## =====================================================
  ## Initial DT uncure probability
  ## =====================================================
  init_label <- factor(Status.train, levels = c(0, 1))
  if (length(unique(init_label)) < 2) {
    stop("Initial DT labels contain only one class. Cannot fit DT.")
  }
  tune_dat <- data.frame(y = init_label, Z.train[, -1, drop = FALSE])
  ## =====================================================
  ## Tune cp and minsplit
  ## =====================================================
  best_xerror <- Inf
  best_cp.init <- NULL
  best_minsplit.init <- NULL
  for (cp_val in cp_grid) {
    for (minsplit_val in minsplit_grid) {
      fit_tmp <- tryCatch(
        rpart::rpart(
          y ~ .,
          data = tune_dat,
          method = "class",
          control = rpart::rpart.control(cp = cp_val, minsplit = minsplit_val, xval = 10)
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(fit_tmp)) {
        next
      }
      tmp_xerror <- min(fit_tmp$cptable[, "xerror"], na.rm = TRUE)
      if (is.finite(tmp_xerror) && tmp_xerror < best_xerror) {
        best_xerror <- tmp_xerror
        best_cp.init <- cp_val
        best_minsplit.init <- minsplit_val
      }
    }
  }
  if (is.null(best_cp.init) || is.null(best_minsplit.init)) {
    stop("DT tuning failed to select cp and minsplit.")
  }
  ## =====================================================
  ## Initial DT fit
  ## =====================================================
  mod.init <- rpart::rpart(
    y ~ .,
    data = tune_dat,
    method = "class",
    control = rpart::rpart.control(cp = best_cp.init, minsplit = best_minsplit.init, xval = 0)
  )
  proba.init <- predict(mod.init, newdata = data.frame(Z.train[, -1, drop = FALSE]), type = "prob")
  if (!("1" %in% colnames(proba.init))) {
    stop("Initial DT did not return probability for uncured class 1.")
  }
  uncureprob.init <- proba.init[, "1"]
  uncureprob.init <- pmin(pmax(uncureprob.init, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Run EM algorithm
  ## =====================================================
  emfit <- em.DT.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    best_cp = best_cp.init,
    best_minsplit = best_minsplit.init,
    emmax = emmax,
    eps = eps,
    M.imp = M.imp
  )
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) {
      stop("Bootstrap requires both event and censored observations.")
    }
    boot_count <- 1
    while (boot_count <= nboot) {
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      bootX <- X.train[boot_id, , drop = FALSE]
      bootZ <- Z.train[boot_id, , drop = FALSE]
      ## =================================================
      ## Initial DT for bootstrap sample
      ## =================================================
      boot_label <- factor(boot_Status, levels = c(0, 1))
      if (length(unique(boot_label)) < 2) {
        next
      }
      boot_dat <- data.frame(y = boot_label, bootZ[, -1, drop = FALSE])
      boot_mod.init <- tryCatch(
        rpart::rpart(
          y ~ .,
          data = boot_dat,
          method = "class",
          control = rpart::rpart.control(cp = best_cp.init, minsplit = best_minsplit.init, xval = 0)
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_mod.init)) {
        next
      }
      boot_proba.init <- tryCatch(
        predict(boot_mod.init, newdata = data.frame(bootZ[, -1, drop = FALSE]), type = "prob"),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_proba.init)) {
        next
      }
      if (!("1" %in% colnames(boot_proba.init))) {
        next
      }
      boot_uncure.init <- boot_proba.init[, "1"]
      boot_uncure.init <- pmin(pmax(boot_uncure.init, 1e-6), 1 - 1e-6)
      ## =================================================
      ## Bootstrap EM
      ## =================================================
      bootfit <- tryCatch(
        em.DT.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          best_cp = best_cp.init,
          best_minsplit = best_minsplit.init,
          emmax = emmax,
          eps = eps,
          M.imp = M.imp
        ),
        error = function(e) {
          NULL
        }
      )
      ## =================================================
      ## Keep converged bootstrap fits
      ## =================================================
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    ## ===================================================
    ## Bootstrap variance and SD
    ## ===================================================
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.DT.ic"
  fit$call <- call
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$models <- emfit$models
  fit$best_cp <- emfit$best_cp
  fit$best_minsplit <- emfit$best_minsplit
  fit$initial_model <- mod.init
  fit$initial_cp <- best_cp.init
  fit$initial_minsplit <- best_minsplit.init
  ## =====================================================
  ## Store generic covariate names
  ## =====================================================
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}

#RANDOM FOREST
em.RF.ic <- function(TimeL.train, TimeR.train, Status.train,
                     X.train, Z.train,
                     TimeL.test, TimeR.test,
                     X.test, Z.test,
                     uncureprob, alpha, beta,
                     best_nodesize, best_ntree,
                     emmax = 1000, eps = 1e-3, M.imp = 5) {
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  n <- length(Status.train)
  m <- nrow(Z.test)
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  final_models <- vector("list", M.imp)
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## Current Weibull latency survival
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    uncureprob <- pmin(pmax(uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## E-step
    ## ===================================================
    w <- Status.train + (1 - Status.train) * (1 - (1 - uncureprob)^survival.train)
    w <- pmin(pmax(w, 1e-6), 1 - 1e-6)
    M.train <- Status.train - survival.train * log1p(-uncureprob)
    ## ===================================================
    ## M-step: RF incidence component
    ## ===================================================
    update_probs <- matrix(NA_real_, nrow = n, ncol = M.imp)
    current_models <- vector("list", M.imp)
    for (j in seq_len(M.imp)) {
      imputed_label <- rbinom(n = n, size = 1, prob = w)
      yj <- factor(imputed_label, levels = c(0, 1))
      if (length(unique(yj)) < 2) {
        update_probs[, j] <- uncureprob
        next
      }
      datj <- data.frame(y = yj, Z.train[, -1, drop = FALSE])
      modj <- tryCatch(
        randomForest::randomForest(
          y ~ .,
          data = datj,
          ntree = best_ntree,
          nodesize = best_nodesize
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(modj)) {
        update_probs[, j] <- uncureprob
        next
      }
      probj <- tryCatch(
        predict(modj, newdata = data.frame(Z.train[, -1, drop = FALSE]), type = "prob"),
        error = function(e) {
          NULL
        }
      )
      if (is.null(probj)) {
        update_probs[, j] <- uncureprob
        next
      }
      if ("1" %in% colnames(probj)) {
        update_probs[, j] <- probj[, "1"]
      } else {
        update_probs[, j] <- uncureprob
      }
      current_models[[j]] <- modj
    }
    update_uncureprob <- rowMeans(update_probs, na.rm = TRUE)
    bad_prob <- !is.finite(update_uncureprob)
    if (any(bad_prob)) {
      update_uncureprob[bad_prob] <- uncureprob[bad_prob]
    }
    update_uncureprob <- pmin(pmax(update_uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    #convergence <- sum(c(mean(update_uncureprob) - mean(uncureprob), update_latency - c(alpha, beta))^2)
    convergence <- mean((update_uncureprob - uncureprob)^2) + sum((update_latency - c(alpha, beta))^2)
    
    uncureprob <- update_uncureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    final_models <- current_models
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob)) / uncureprob
  ## =====================================================
  ## Test incidence predictions
  ## =====================================================
  test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
  for (j in seq_len(M.imp)) {
    if (is.null(final_models[[j]])) {
      next
    }
    prob.test <- tryCatch(
      predict(final_models[[j]], newdata = data.frame(Z.test[, -1, drop = FALSE]), type = "prob"),
      error = function(e) {
        NULL
      }
    )
    if (is.null(prob.test)) {
      next
    }
    if ("1" %in% colnames(prob.test)) {
      test_probs[, j] <- prob.test[, "1"]
    }
  }
  if (all(is.na(test_probs))) {
    stop("No valid random-forest model available for test prediction.")
  }
  uncureprob.test <- rowMeans(test_probs, na.rm = TRUE)
  bad_test <- !is.finite(uncureprob.test)
  if (any(bad_test)) {
    uncureprob.test[bad_test] <- mean(uncureprob)
  }
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Test survival predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter,
    models = final_models,
    best_nodesize = best_nodesize,
    best_ntree = best_ntree
  ))
}
smcure.RF.ic <- function(formula, cureform,
                         alpha.init, beta.init,
                         train_data, test_data,
                         na.action = na.omit,
                         Var = FALSE,
                         emmax = 1000,
                         eps = 1e-3,
                         nboot = 100,
                         nodesize_grid =  c(1, 2, 3, 5, 8, 10, 15, 20, 30),
                         ntree_grid = c(500, 800, 1000, 1500),
                         M.imp = 5) {
  call <- match.call()
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  ## =====================================================
  ## Initial RF uncure probability
  ## =====================================================
  init_label <- factor(Status.train, levels = c(0, 1))
  if (length(unique(init_label)) < 2) {
    stop("Initial RF labels contain only one class. Cannot fit RF.")
  }
  init_dat <- data.frame(y = init_label, Z.train[, -1, drop = FALSE])
  ## =====================================================
  ## Tune nodesize and ntree
  ## =====================================================
  best_oob.init <- -Inf
  tuning_results <- data.frame(nodesize = integer(), ntree = integer(), OOB_AUC = numeric())
  best_nodesize.init <- NULL
  best_ntree.init <- NULL
  for (nodesize_val in nodesize_grid) {
    for (ntree_val in ntree_grid) {
      fit_tmp <- tryCatch(
        randomForest::randomForest(
          y ~ .,
          data = init_dat,
          ntree = ntree_val,
          nodesize = nodesize_val
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(fit_tmp)) {
        next
      }
      # Each row's vote is calculated using trees where that row was OOB.
      oob_prob <- fit_tmp$votes[, "1"]
      valid <- is.finite(oob_prob)
      tmp_oob <- NA_real_
      if (sum(valid) >= 2 && length(unique(Status.train[valid])) == 2) {
        tmp_oob <- tryCatch(
          as.numeric(pROC::auc(
            response = Status.train[valid],
            predictor = oob_prob[valid],
            levels = c(0, 1),
            direction = "<",
            quiet = TRUE
          )),
          error = function(e) NA_real_
        )
      }
      tuning_results <- rbind(
        tuning_results,
        data.frame(nodesize = nodesize_val,
                   ntree = ntree_val,
                   OOB_AUC = tmp_oob)
      )
      if (is.finite(tmp_oob) && tmp_oob > best_oob.init) {
        best_oob.init <- tmp_oob
        best_nodesize.init <- nodesize_val
        best_ntree.init <- ntree_val
      }
    }
  }
  if (is.null(best_nodesize.init) || is.null(best_ntree.init)) {
    stop("RF OOB-AUC tuning failed to select nodesize and ntree.")
  }
  ## =====================================================
  ## Initial RF fit
  ## =====================================================
  mod.init <- randomForest::randomForest(
    y ~ .,
    data = init_dat,
    ntree = best_ntree.init,
    nodesize = best_nodesize.init
  )
  proba.init <- predict(mod.init, newdata = data.frame(Z.train[, -1, drop = FALSE]), type = "prob")
  if (!("1" %in% colnames(proba.init))) {
    stop("Initial RF did not return probability for uncured class 1.")
  }
  uncureprob.init <- proba.init[, "1"]
  uncureprob.init <- pmin(pmax(uncureprob.init, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Run EM algorithm
  ## =====================================================
  emfit <- em.RF.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    best_nodesize = best_nodesize.init,
    best_ntree = best_ntree.init,
    emmax = emmax,
    eps = eps,
    M.imp = M.imp
  )
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) {
      stop("Bootstrap requires both event and censored observations.")
    }
    boot_count <- 1
    while (boot_count <= nboot) {
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      bootX <- X.train[boot_id, , drop = FALSE]
      bootZ <- Z.train[boot_id, , drop = FALSE]
      ## =================================================
      ## Initial RF for bootstrap sample
      ## =================================================
      boot_label <- factor(boot_Status, levels = c(0, 1))
      if (length(unique(boot_label)) < 2) {
        next
      }
      boot_dat <- data.frame(y = boot_label, bootZ[, -1, drop = FALSE])
      boot_mod.init <- tryCatch(
        randomForest::randomForest(
          y ~ .,
          data = boot_dat,
          ntree = best_ntree.init,
          nodesize = best_nodesize.init
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_mod.init)) {
        next
      }
      boot_proba.init <- tryCatch(
        predict(boot_mod.init, newdata = data.frame(bootZ[, -1, drop = FALSE]), type = "prob"),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_proba.init)) {
        next
      }
      if (!("1" %in% colnames(boot_proba.init))) {
        next
      }
      boot_uncure.init <- boot_proba.init[, "1"]
      boot_uncure.init <- pmin(pmax(boot_uncure.init, 1e-6), 1 - 1e-6)
      ## =================================================
      ## Bootstrap EM
      ## =================================================
      bootfit <- tryCatch(
        em.RF.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          best_nodesize = best_nodesize.init,
          best_ntree = best_ntree.init,
          emmax = emmax,
          eps = eps,
          M.imp = M.imp
        ),
        error = function(e) {
          NULL
        }
      )
      ## =================================================
      ## Keep converged bootstrap fits
      ## =================================================
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    ## ===================================================
    ## Bootstrap variance and SD
    ## ===================================================
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.RF.ic"
  fit$call <- call
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$models <- emfit$models
  fit$best_nodesize <- emfit$best_nodesize
  fit$best_ntree <- emfit$best_ntree
  fit$initial_model <- mod.init
  fit$initial_nodesize <- best_nodesize.init
  fit$initial_ntree <- best_ntree.init
  fit$initial_oob <- best_oob.init  # Retained for compatibility; now stores OOB AUC
  fit$initial_auc <- best_oob.init
  fit$tuning_results <- tuning_results[order(-tuning_results$OOB_AUC), ]
  ## =====================================================
  ## Store generic covariate names
  ## =====================================================
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}



#XGBOOST
em.XGB.ic <- function(TimeL.train, TimeR.train, Status.train,
                      X.train, Z.train,
                      TimeL.test, TimeR.test,
                      X.test, Z.test,
                      uncureprob, alpha, beta,
                      best_eta, best_nrounds,
                      emmax = 1000, eps = 1e-3, M.imp = 5) {
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  n <- length(Status.train)
  m <- nrow(Z.test)
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  convergence <- Inf
  iter <- 1
  final_models <- vector("list", M.imp)
  Xgb.train <- as.matrix(Z.train[, -1, drop = FALSE])
  Xgb.test <- as.matrix(Z.test[, -1, drop = FALSE])
  dtest <- xgboost::xgb.DMatrix(data = Xgb.test)
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## Current Weibull latency survival
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    uncureprob <- pmin(pmax(uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## E-step
    ## ===================================================
    w <- Status.train + (1 - Status.train) * (1 - (1 - uncureprob)^survival.train)
    w <- pmin(pmax(w, 1e-6), 1 - 1e-6)
    M.train <- Status.train - survival.train * log1p(-uncureprob)
    ## ===================================================
    ## M-step: XGBoost incidence component
    ## ===================================================
    update_probs <- matrix(NA_real_, nrow = n, ncol = M.imp)
    current_models <- vector("list", M.imp)
    for (j in seq_len(M.imp)) {
      imputed_label <- rbinom(n = n, size = 1, prob = w)
      if (length(unique(imputed_label)) < 2) {
        update_probs[, j] <- uncureprob
        next
      }
      dtrainj <- xgboost::xgb.DMatrix(data = Xgb.train, label = imputed_label)
      modj <- tryCatch(
        xgboost::xgb.train(
          params = list(
            objective = "binary:logistic",
            eval_metric = "logloss",
            learning_rate = best_eta,
            max_depth = 4,
            min_child_weight = 0.8,
            gamma = 0.1,
            subsample = 0.7,
            colsample_bytree = 1,
            lambda = 5,
            alpha = 0
          ),
          data = dtrainj,
          nrounds = best_nrounds,
          verbose = 0
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(modj)) {
        update_probs[, j] <- uncureprob
        next
      }
      predj <- tryCatch(
        predict(modj, dtrainj),
        error = function(e) {
          rep(NA_real_, n)
        }
      )
      if (length(predj) != n || any(!is.finite(predj))) {
        update_probs[, j] <- uncureprob
        next
      }
      update_probs[, j] <- predj
      current_models[[j]] <- modj
    }
    update_uncureprob <- rowMeans(update_probs, na.rm = TRUE)
    bad_prob <- !is.finite(update_uncureprob)
    if (any(bad_prob)) {
      update_uncureprob[bad_prob] <- uncureprob[bad_prob]
    }
    update_uncureprob <- pmin(pmax(update_uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) {
        return(1e20)
      }
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) {
        return(1e20)
      }
      return(-out)
    }
    update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    ## ===================================================
    ## Convergence
    ## ===================================================
    #convergence <- sum(c(mean(update_uncureprob) - mean(uncureprob), update_latency - c(alpha, beta))^2)
    convergence <- mean((update_uncureprob - uncureprob)^2) + sum((update_latency - c(alpha, beta))^2)
    
    uncureprob <- update_uncureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    final_models <- current_models
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob)) / uncureprob
  ## =====================================================
  ## Test incidence predictions
  ## =====================================================
  test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
  for (j in seq_len(M.imp)) {
    if (is.null(final_models[[j]])) {
      next
    }
    pred.test <- tryCatch(
      predict(final_models[[j]], dtest),
      error = function(e) {
        rep(NA_real_, m)
      }
    )
    if (length(pred.test) == m && all(is.finite(pred.test))) {
      test_probs[, j] <- pred.test
    }
  }
  if (all(is.na(test_probs))) {
    stop("No valid XGBoost model available for test prediction.")
  }
  uncureprob.test <- rowMeans(test_probs, na.rm = TRUE)
  bad_test <- !is.finite(uncureprob.test)
  if (any(bad_test)) {
    uncureprob.test[bad_test] <- mean(uncureprob)
  }
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Test survival predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter,
    models = final_models,
    best_eta = best_eta,
    best_nrounds = best_nrounds
  ))
}
smcure.XGB.ic <- function(formula, cureform,
                          alpha.init, beta.init,
                          train_data, test_data,
                          na.action = na.omit,
                          Var = FALSE,
                          emmax = 1000,
                          eps = 1e-3,
                          nboot = 100,
                          eta_grid = c(0.01, 0.05, 0.10, 0.20),
                          nrounds_grid = c(50, 100, 200, 500),
                          M.imp = 5) {
  call <- match.call()
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) {
    stop("train_data must contain status column D.")
  }
  if (!("D" %in% names(test_data))) {
    stop("test_data must contain status column D.")
  }
  ## =====================================================
  ## Build latency model matrices: X
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) {
    stop("Training response must be a survival object.")
  }
  if (!inherits(Y.test, "Surv")) {
    stop("Test response must be a survival object.")
  }
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  ## =====================================================
  ## Build incidence model matrices: Z
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- train_data$D
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) {
    stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  }
  ## =====================================================
  ## Initial XGBoost uncure probability
  ## =====================================================
  init_label <- Status.train
  if (length(unique(init_label)) < 2) {
    stop("Initial XGBoost labels contain only one class. Cannot fit XGBoost.")
  }
  Xgb.train <- as.matrix(Z.train[, -1, drop = FALSE])
  dtrain.init <- xgboost::xgb.DMatrix(data = Xgb.train, label = init_label)
  ## =====================================================
  ## Tune learning rate and nrounds
  ## =====================================================
  best_logloss.init <- Inf
  best_eta.init <- NULL
  best_nrounds.init <- NULL
  for (eta_val in eta_grid) {
    for (nrounds_val in nrounds_grid) {
      cvfit <- tryCatch(
        xgboost::xgb.cv(
          params = list(
            objective = "binary:logistic",
            eval_metric = "logloss",
            learning_rate = best_eta,
            max_depth = 4,
            min_child_weight = 0.8,
            gamma = 0.1,
            subsample = 0.7,
            colsample_bytree = 1,
            lambda = 5,
            alpha = 0),
          data = dtrain.init,
          nrounds = nrounds_val,
          nfold = 10,
          verbose = 0
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(cvfit)) {
        next
      }
      tmp_logloss <- min(cvfit$evaluation_log$test_logloss_mean, na.rm = TRUE)
      if (is.finite(tmp_logloss) && tmp_logloss < best_logloss.init) {
        best_logloss.init <- tmp_logloss
        best_eta.init <- eta_val
        best_nrounds.init <- nrounds_val
      }
    }
  }
  if (is.null(best_eta.init) || is.null(best_nrounds.init)) {
    best_eta.init <- 0.05
    best_nrounds.init <- 100
  }
  ## =====================================================
  ## Initial XGBoost fit
  ## =====================================================
  mod.init <- xgboost::xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      learning_rate = best_eta.init,
      max_depth = 3
    ),
    data = dtrain.init,
    nrounds = best_nrounds.init,
    verbose = 0
  )
  uncureprob.init <- predict(mod.init, dtrain.init)
  uncureprob.init <- pmin(pmax(uncureprob.init, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Run EM algorithm
  ## =====================================================
  emfit <- em.XGB.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    best_eta = best_eta.init,
    best_nrounds = best_nrounds.init,
    emmax = emmax,
    eps = eps,
    M.imp = M.imp
  )
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    event_idx <- which(Status.train == 1)
    cens_idx <- which(Status.train == 0)
    n1 <- length(event_idx)
    n0 <- length(cens_idx)
    if (n1 == 0 || n0 == 0) {
      stop("Bootstrap requires both event and censored observations.")
    }
    boot_count <- 1
    while (boot_count <= nboot) {
      id1 <- sample(event_idx, size = n1, replace = TRUE)
      id0 <- sample(cens_idx, size = n0, replace = TRUE)
      boot_id <- c(id1, id0)
      boot_TimeL <- TimeL.train[boot_id]
      boot_TimeR <- TimeR.train[boot_id]
      boot_Status <- Status.train[boot_id]
      bootX <- X.train[boot_id, , drop = FALSE]
      bootZ <- Z.train[boot_id, , drop = FALSE]
      ## =================================================
      ## Initial XGBoost for bootstrap sample
      ## =================================================
      if (length(unique(boot_Status)) < 2) {
        next
      }
      bootXgb <- as.matrix(bootZ[, -1, drop = FALSE])
      boot_dtrain.init <- xgboost::xgb.DMatrix(data = bootXgb, label = boot_Status)
      boot_mod.init <- tryCatch(
        xgboost::xgb.train(
          params = list(
            objective = "binary:logistic",
            eval_metric = "logloss",
            learning_rate = best_eta.init,
            max_depth = 3
          ),
          data = boot_dtrain.init,
          nrounds = best_nrounds.init,
          verbose = 0
        ),
        error = function(e) {
          NULL
        }
      )
      if (is.null(boot_mod.init)) {
        next
      }
      boot_uncure.init <- tryCatch(
        predict(boot_mod.init, boot_dtrain.init),
        error = function(e) {
          rep(NA_real_, length(boot_Status))
        }
      )
      if (length(boot_uncure.init) != length(boot_Status) || any(!is.finite(boot_uncure.init))) {
        next
      }
      boot_uncure.init <- pmin(pmax(boot_uncure.init, 1e-6), 1 - 1e-6)
      ## =================================================
      ## Bootstrap EM
      ## =================================================
      bootfit <- tryCatch(
        em.XGB.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          best_eta = best_eta.init,
          best_nrounds = best_nrounds.init,
          emmax = emmax,
          eps = eps,
          M.imp = M.imp
        ),
        error = function(e) {
          NULL
        }
      )
      ## =================================================
      ## Keep converged bootstrap fits
      ## =================================================
      if (!is.null(bootfit) && is.finite(bootfit$tau) && bootfit$tau < eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    ## ===================================================
    ## Bootstrap variance and SD
    ## ===================================================
    latency_var <- apply(latency_boot, 2, var, na.rm = TRUE)
    latency_sd <- sqrt(latency_var)
  }
  ## =====================================================
  ## Return fitted object
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.XGB.ic"
  fit$call <- call
  fit$latency <- latency
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$models <- emfit$models
  fit$best_eta <- emfit$best_eta
  fit$best_nrounds <- emfit$best_nrounds
  fit$best_logloss <- best_logloss.init
  fit$initial_model <- mod.init
  ## =====================================================
  ## Store generic covariate names
  ## =====================================================
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  ## =====================================================
  ## Bootstrap inference
  ## =====================================================
  if (Var) {
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}


#NEURAL NETWORK
em.NN.ic <- function(TimeL.train, TimeR.train, Status.train,
                     X.train, Z.train,
                     TimeL.test, TimeR.test,
                     X.test, Z.test,
                     uncureprob, alpha, beta,
                     best_params,
                     emmax = 1000, eps = 1e-3, M.imp = 5,
                     nn_stepmax = 5e5, nn_threshold = 0.05) {
  ## =====================================================
  ## Basic quantities
  ## =====================================================
  n <- length(Status.train)
  m <- nrow(Z.test)
  TimeY.train <- ifelse(Status.train == 0, TimeL.train, TimeR.train)
  NN.train <- as.data.frame(Z.train[, -1, drop = FALSE])
  NN.test <- as.data.frame(Z.test[, -1, drop = FALSE])
  colnames(NN.train) <- paste0("Input_", seq_len(ncol(NN.train)))
  colnames(NN.test) <- colnames(NN.train)
  nn_formula <- as.formula(paste("y ~", paste(colnames(NN.train), collapse = " + ")))
  convergence <- Inf
  iter <- 1
  final_models <- vector("list", M.imp)
  final_test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
  ## =====================================================
  ## EM algorithm
  ## =====================================================
  while (convergence > eps && iter < emmax) {
    ## ===================================================
    ## Current Weibull latency survival
    ## ===================================================
    eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
    H.train <- (TimeY.train^alpha) * exp(eta_lat.train)
    survival.train <- exp(-H.train)
    uncureprob <- pmin(pmax(uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## E-step
    ## ===================================================
    w <- Status.train + (1 - Status.train) * (1 - (1 - uncureprob)^survival.train)
    w <- pmin(pmax(w, 1e-6), 1 - 1e-6)
    M.train <- Status.train - survival.train * log1p(-uncureprob)
    ## ===================================================
    ## M-step: NN incidence component
    ## ===================================================
    train_probs <- matrix(NA_real_, nrow = n, ncol = M.imp)
    current_models <- vector("list", M.imp)
    current_test_probs <- matrix(NA_real_, nrow = m, ncol = M.imp)
    for (j in seq_len(M.imp)) {
      imputed_label <- rbinom(n = n, size = 1, prob = w)
      if (length(unique(imputed_label)) < 2) {
        train_probs[, j] <- uncureprob
        next
      }
      datj <- data.frame(y = imputed_label, NN.train, check.names = FALSE)
      modj <- tryCatch(
        neuralnet::neuralnet(
          formula = nn_formula,
          data = datj,
          hidden = best_params,
          linear.output = FALSE,
          err.fct = "ce",
          lifesign = "none",
          stepmax = nn_stepmax,
          threshold = nn_threshold
        ),
        error = function(e) NULL
      )
      if (is.null(modj)) {
        train_probs[, j] <- uncureprob
        next
      }
      ## =================================================
      ## Check training prediction
      ## =================================================
      pred_train_j <- tryCatch(
        as.numeric(neuralnet::compute(modj, covariate = NN.train)$net.result),
        error = function(e) rep(NA_real_, n)
      )
      if (length(pred_train_j) != n || any(!is.finite(pred_train_j))) {
        train_probs[, j] <- uncureprob
        next
      }
      ## =================================================
      ## Check test prediction before accepting model
      ## =================================================
      pred_test_j <- tryCatch(
        as.numeric(neuralnet::compute(modj, covariate = NN.test)$net.result),
        error = function(e) rep(NA_real_, m)
      )
      if (length(pred_test_j) != m || any(!is.finite(pred_test_j))) {
        train_probs[, j] <- uncureprob
        next
      }
      ## =================================================
      ## Store valid model and predictions
      ## =================================================
      pred_train_j <- pmin(pmax(pred_train_j, 1e-6), 1 - 1e-6)
      pred_test_j <- pmin(pmax(pred_test_j, 1e-6), 1 - 1e-6)
      train_probs[, j] <- pred_train_j
      current_test_probs[, j] <- pred_test_j
      current_models[[j]] <- modj
    }
    ## ===================================================
    ## Preserve previous valid models when current fit fails
    ## ===================================================
    for (j in seq_len(M.imp)) {
      if (!is.null(current_models[[j]])) {
        final_models[[j]] <- current_models[[j]]
        final_test_probs[, j] <- current_test_probs[, j]
      }
    }
    valid_models <- sum(!vapply(final_models, is.null, logical(1)))
    if (valid_models == 0) {
      stop("No valid neural-network model has been obtained during the EM algorithm.")
    }
    ## ===================================================
    ## Update incidence probabilities
    ## ===================================================
    update_uncureprob <- rowMeans(train_probs, na.rm = TRUE)
    bad_train <- !is.finite(update_uncureprob)
    if (any(bad_train)) update_uncureprob[bad_train] <- uncureprob[bad_train]
    update_uncureprob <- pmin(pmax(update_uncureprob, 1e-6), 1 - 1e-6)
    ## ===================================================
    ## M-step: Weibull latency component
    ## ===================================================
    Q2 <- function(par) {
      alpha_new <- par[1]
      beta_new <- par[-1]
      if (!is.finite(alpha_new) || alpha_new <= 0) return(1e20)
      eta2 <- drop(X.train[, -1, drop = FALSE] %*% beta_new)
      HY <- (TimeY.train^alpha_new) * exp(eta2)
      S.Y <- exp(-HY)
      S.Y <- pmin(pmax(S.Y, 1e-12), 1)
      out1 <- sum((M.train - Status.train) * log(S.Y))
      event_idx <- which(Status.train == 1)
      if (length(event_idx) > 0) {
        HL <- (TimeL.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        HR <- (TimeR.train[event_idx]^alpha_new) * exp(eta2[event_idx])
        S.L <- exp(-HL)
        S.R <- exp(-HR)
        diff_SR <- pmax(S.L - S.R, 1e-12)
        out2 <- sum(log(diff_SR))
      } else {
        out2 <- 0
      }
      out <- out1 + out2
      if (!is.finite(out)) return(1e20)
      return(-out)
    }
    #update_latency <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")$par
    #if (any(!is.finite(update_latency))) stop("The Weibull latency update produced non-finite parameters.")
    latency_opt <- optim(par = c(alpha, beta), fn = Q2, method = "Nelder-Mead")
    update_latency <- latency_opt$par
    
    ## ===================================================
    ## Convergence
    ## ===================================================
    convergence <- sum(c(mean(update_uncureprob) - mean(uncureprob), update_latency - c(alpha, beta))^2)
    #convergence <- mean((update_uncureprob - uncureprob)^2) + sum((update_latency - c(alpha, beta))^2)
    
    #fresh_models <- sum(!vapply(current_models, is.null, logical(1)))
   # incidence_rms <- sqrt(mean((update_uncureprob - uncureprob)^2))
   # latency_change <- sqrt(sum((update_latency - c(alpha, beta))^2))
    #message(sprintf("EM %d | fresh NN fits: %d/%d | incidence RMS: %.5f | latency change: %.5f | tau: %.6g | optim code: %d", iter, fresh_models, M.imp, incidence_rms, latency_change, convergence, latency_opt$convergence))
    
    
    uncureprob <- update_uncureprob
    alpha <- update_latency[1]
    beta <- update_latency[-1]
    iter <- iter + 1
  }
  ## =====================================================
  ## Training predictions
  ## =====================================================
  TimeK.train <- ifelse(is.infinite(TimeR.train), TimeL.train, TimeR.train)
  Time.mid.train <- (TimeL.train + TimeK.train) / 2
  eta_lat.train <- drop(X.train[, -1, drop = FALSE] %*% beta)
  H.mid.train <- (Time.mid.train^alpha) * exp(eta_lat.train)
  S.train <- exp(-H.mid.train)
  Sp.train <- (1 - uncureprob)^(1 - S.train)
  S1.train <- (Sp.train - (1 - uncureprob)) / uncureprob
  ## =====================================================
  ## Test incidence predictions
  ## =====================================================
  valid_cols <- which(colSums(is.finite(final_test_probs)) == m)
  if (length(valid_cols) == 0) {
    stop("No valid neural-network model was able to predict the test data.")
  }
  uncureprob.test <- rowMeans(final_test_probs[, valid_cols, drop = FALSE])
  uncureprob.test <- pmin(pmax(uncureprob.test, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Test survival predictions
  ## =====================================================
  TimeK.test <- ifelse(is.infinite(TimeR.test), TimeL.test, TimeR.test)
  Time.mid.test <- (TimeL.test + TimeK.test) / 2
  eta_lat.test <- drop(X.test[, -1, drop = FALSE] %*% beta)
  H.mid.test <- (Time.mid.test^alpha) * exp(eta_lat.test)
  S.test <- exp(-H.mid.test)
  Sp.test <- (1 - uncureprob.test)^(1 - S.test)
  S1.test <- (Sp.test - (1 - uncureprob.test)) / uncureprob.test
  ## =====================================================
  ## Return
  ## =====================================================
  return(list(
    latencyfit = c(alpha, beta),
    train = list(Time.mid = Time.mid.train, Uncureprob = uncureprob, Sp = Sp.train, S1 = S1.train, S = S.train),
    test = list(Time.mid = Time.mid.test, Uncureprob = uncureprob.test, Sp = Sp.test, S1 = S1.test, S = S.test),
    tau = convergence,
    iter = iter,
    converged = is.finite(convergence) && convergence <= eps,
    models = final_models,
    best_params = best_params
  ))
}
smcure.NN.ic <- function(formula, cureform,
                         alpha.init, beta.init,
                         train_data, test_data,
                         na.action = na.omit,
                         Var = FALSE,
                         emmax = 1000,
                         eps = 1e-3,
                         nboot = 100,
                         hidden_grid = list(c(2), c(3), c(4), c(5), c(6), c(8), c(10), c(12), c(16), c(3, 2), c(4, 2), c(5, 3), c(6, 3), c(8, 4), c(10, 3), c(10, 5), c(12, 4), c(12, 6), c(16, 8)),
                         M.imp = 5,
                         nn_stepmax = 5e5,
                         nn_threshold = 0.05) {
  call <- match.call()
  ## =====================================================
  ## Prepare training and test data
  ## =====================================================
  train_data <- na.action(train_data)
  test_data <- na.action(test_data)
  if (!("D" %in% names(train_data))) stop("train_data must contain status column D.")
  if (!("D" %in% names(test_data))) stop("test_data must contain status column D.")
  ## =====================================================
  ## Build latency model matrices
  ## =====================================================
  mf.train <- model.frame(formula, data = train_data)
  mf.test <- model.frame(formula, data = test_data)
  Y.train <- model.extract(mf.train, "response")
  Y.test <- model.extract(mf.test, "response")
  if (!inherits(Y.train, "Surv")) stop("Training response must be a survival object.")
  if (!inherits(Y.test, "Surv")) stop("Test response must be a survival object.")
  X.train <- model.matrix(attr(mf.train, "terms"), mf.train)
  X.test <- model.matrix(attr(mf.test, "terms"), mf.test)
  if (!identical(colnames(X.train), colnames(X.test))) stop("Training and test latency design matrices have different columns.")
  ## =====================================================
  ## Build incidence matrices
  ## =====================================================
  Z.train <- model.matrix(cureform, data = train_data)
  Z.test <- model.matrix(cureform, data = test_data)
  if (!identical(colnames(Z.train), colnames(Z.test))) stop("Training and test incidence design matrices have different columns.")
  ## =====================================================
  ## Survival quantities
  ## =====================================================
  TimeL.train <- Y.train[, 1]
  TimeR.train <- Y.train[, 2]
  Status.train <- as.numeric(train_data$D)
  TimeL.test <- Y.test[, 1]
  TimeR.test <- Y.test[, 2]
  if (any(!Status.train %in% c(0, 1))) stop("Training status D must contain only 0 and 1.")
  ## =====================================================
  ## Parameter dimensions
  ## =====================================================
  nbeta <- ncol(X.train) - 1
  if (length(beta.init) != nbeta) stop("Length of beta.init must equal the number of latency covariates: ", nbeta)
  if (!is.finite(alpha.init) || alpha.init <= 0) stop("alpha.init must be positive and finite.")
  ## =====================================================
  ## Construct NN data
  ## =====================================================
  NN.train <- as.data.frame(Z.train[, -1, drop = FALSE])
  NN.test <- as.data.frame(Z.test[, -1, drop = FALSE])
  nn_names <- paste0("Input_", seq_len(ncol(NN.train)))
  colnames(NN.train) <- nn_names
  colnames(NN.test) <- nn_names
  nn_formula <- as.formula(paste("y ~", paste(nn_names, collapse = " + ")))
  init_label <- Status.train
  if (length(unique(init_label)) < 2) stop("Initial NN labels contain only one class. The NN incidence model cannot be fitted.")
  init_data <- data.frame(y = init_label, NN.train, check.names = FALSE)
  ## =====================================================
  ## Tune architecture once using training log-loss
  ## =====================================================
  best_loss.init <- Inf
  best_params.init <- NULL
  best_model.init <- NULL
  for (hidden_val in hidden_grid) {
    fit_tmp <- tryCatch(neuralnet::neuralnet(formula = nn_formula, data = init_data, hidden = hidden_val, linear.output = FALSE, err.fct = "ce", lifesign = "none", stepmax = nn_stepmax, threshold = nn_threshold), error = function(e) NULL)
    if (is.null(fit_tmp)) next
    pred_tmp <- tryCatch(as.numeric(neuralnet::compute(fit_tmp, covariate = NN.train)$net.result), error = function(e) rep(NA_real_, nrow(NN.train)))
    if (length(pred_tmp) != nrow(NN.train)) next
    if (any(!is.finite(pred_tmp))) next
    pred_tmp <- pmin(pmax(pred_tmp, 1e-6), 1 - 1e-6)
    loss_tmp <- -mean(init_label * log(pred_tmp) + (1 - init_label) * log1p(-pred_tmp))
    if (is.finite(loss_tmp) && loss_tmp < best_loss.init) {
      best_loss.init <- loss_tmp
      best_params.init <- hidden_val
      best_model.init <- fit_tmp
    }
  }
  if (is.null(best_model.init) || is.null(best_params.init)) stop("All initial NN architecture fits failed. Increase nn_stepmax or use a simpler hidden_grid.")
  ## =====================================================
  ## Initial uncured probabilities
  ## =====================================================
  best_params <- best_params.init
  mod.init <- best_model.init
  uncureprob.init <- tryCatch(as.numeric(neuralnet::compute(mod.init, covariate = NN.train)$net.result), error = function(e) rep(NA_real_, nrow(NN.train)))
  if (length(uncureprob.init) != nrow(NN.train) || any(!is.finite(uncureprob.init))) stop("The selected initial NN produced invalid uncured probabilities.")
  uncureprob.init <- pmin(pmax(uncureprob.init, 1e-6), 1 - 1e-6)
  ## =====================================================
  ## Run EM
  ## =====================================================
  emfit <- em.NN.ic(
    TimeL.train = TimeL.train,
    TimeR.train = TimeR.train,
    Status.train = Status.train,
    X.train = X.train,
    Z.train = Z.train,
    TimeL.test = TimeL.test,
    TimeR.test = TimeR.test,
    X.test = X.test,
    Z.test = Z.test,
    uncureprob = uncureprob.init,
    alpha = alpha.init,
    beta = beta.init,
    best_params = best_params,
    emmax = emmax,
    eps = eps,
    M.imp = M.imp,
    nn_stepmax = nn_stepmax,
    nn_threshold = nn_threshold
  )
  latency <- emfit$latencyfit
  alpha.est <- latency[1]
  beta.est <- latency[-1]
  conv <- emfit$tau
  ## =====================================================
  ## Bootstrap latency variance
  ## =====================================================
  if (Var) {
    latency_boot <- matrix(NA_real_, nrow = nboot, ncol = nbeta + 1)
    event_index <- which(Status.train == 1)
    cens_index <- which(Status.train == 0)
    if (length(event_index) == 0 || length(cens_index) == 0) stop("Bootstrap requires both event and censored observations in the training data.")
    boot_count <- 1
    attempt <- 0
    max_attempts <- 10 * nboot
    while (boot_count <= nboot && attempt < max_attempts) {
      attempt <- attempt + 1
      boot_event <- sample(event_index, size = length(event_index), replace = TRUE)
      boot_cens <- sample(cens_index, size = length(cens_index), replace = TRUE)
      boot_index <- c(boot_event, boot_cens)
      boot_TimeL <- TimeL.train[boot_index]
      boot_TimeR <- TimeR.train[boot_index]
      boot_Status <- Status.train[boot_index]
      bootX <- X.train[boot_index, , drop = FALSE]
      bootZ <- Z.train[boot_index, , drop = FALSE]
      if (length(unique(boot_Status)) < 2) next
      bootNN <- as.data.frame(bootZ[, -1, drop = FALSE])
      colnames(bootNN) <- nn_names
      boot_init_data <- data.frame(y = boot_Status, bootNN, check.names = FALSE)
      boot_initial_model <- tryCatch(neuralnet::neuralnet(formula = nn_formula, data = boot_init_data, hidden = best_params, linear.output = FALSE, err.fct = "ce", lifesign = "none", stepmax = nn_stepmax, threshold = nn_threshold), error = function(e) NULL)
      if (is.null(boot_initial_model)) next
      boot_uncure.init <- tryCatch(as.numeric(neuralnet::compute(boot_initial_model, covariate = bootNN)$net.result), error = function(e) rep(NA_real_, length(boot_index)))
      if (length(boot_uncure.init) != length(boot_index) || any(!is.finite(boot_uncure.init))) next
      boot_uncure.init <- pmin(pmax(boot_uncure.init, 1e-6), 1 - 1e-6)
      bootfit <- tryCatch(
        em.NN.ic(
          TimeL.train = boot_TimeL,
          TimeR.train = boot_TimeR,
          Status.train = boot_Status,
          X.train = bootX,
          Z.train = bootZ,
          TimeL.test = TimeL.test,
          TimeR.test = TimeR.test,
          X.test = X.test,
          Z.test = Z.test,
          uncureprob = boot_uncure.init,
          alpha = alpha.est,
          beta = beta.est,
          best_params = best_params,
          emmax = emmax,
          eps = eps,
          M.imp = M.imp,
          nn_stepmax = nn_stepmax,
          nn_threshold = nn_threshold
        ),
        error = function(e) NULL
      )
      if (is.null(bootfit)) next
      if (!is.finite(bootfit$tau) || any(!is.finite(bootfit$latencyfit))) next
      if (bootfit$tau <= eps) {
        latency_boot[boot_count, ] <- bootfit$latencyfit
        boot_count <- boot_count + 1
      }
    }
    valid_boot <- complete.cases(latency_boot)
    nboot_valid <- sum(valid_boot)
    if (nboot_valid < 2) {
      warning("Fewer than two valid NN bootstrap fits were obtained. Bootstrap standard deviations cannot be computed.")
      latency_var <- rep(NA_real_, nbeta + 1)
      latency_sd <- rep(NA_real_, nbeta + 1)
    } else {
      latency_var <- apply(latency_boot[valid_boot, , drop = FALSE], 2, var)
      latency_sd <- sqrt(latency_var)
    }
  } else {
    latency_boot <- NULL
    latency_var <- NULL
    latency_sd <- NULL
    nboot_valid <- NULL
  }
  ## =====================================================
  ## Return
  ## =====================================================
  fit <- list()
  class(fit) <- "smcure.NN.ic"
  fit$call <- call
  fit$latency <- latency
  fit$alpha <- alpha.est
  fit$beta <- beta.est
  fit$tau <- conv
  fit$iter <- emfit$iter
  fit$converged <- emfit$converged
  fit$train <- emfit$train
  fit$test <- emfit$test
  fit$models <- emfit$models
  fit$best_params <- best_params
  fit$best_hidden <- paste(best_params, collapse = "-")
  fit$best_loss <- best_loss.init
  fit$initial_model <- mod.init
  fit$initial_uncureprob <- uncureprob.init
  fit$incidence_names <- colnames(Z.train)
  fit$latency_names <- c("alpha", colnames(X.train)[-1])
  if (Var) {
    fit$latency_boot <- latency_boot
    fit$nboot_valid <- nboot_valid
    fit$latency_var <- latency_var
    fit$latency_sd <- latency_sd
    fit$latency_zvalue <- fit$latency / fit$latency_sd
    fit$latency_pvalue <- 2 * (1 - pnorm(abs(fit$latency_zvalue)))
  }
  return(fit)
}




#Data
data("hdsd", package = "thregI")
H <- hdsd[-c(1, 2), ]
H$x1<-H$Age
H$x2<- H$Sex
H$x3<- H$TR360
H$x4<-H$Noadyn
H$L<- H$left
H$R<-H$right
H$D <- as.numeric(!is.infinite(H$R))
H$x1 = (H$x1- mean(H$x1))/(sd(H$x1))
#H$x3 = (H$x3- mean(H$x3))/(sd(H$x3))


#my_data <- H[sample(nrow(H), size = 500, replace = TRUE), ]
#H <- my_data

## =====================================================
## 10-fold Cross Validation
## =====================================================


set.seed(2026)

K <- 10
n <- nrow(H)
fold_id <- rep(NA_integer_, n)

for (dval in unique(H$D)) {
  idx <- which(H$D == dval)
  idx <- sample(idx)
  fold_id[idx] <- rep(1:K, length.out = length(idx))
}

surv_formula <- Surv(L, R, type = "interval2") ~ x1 + x2 + x3 + x4
cure_formula <- ~ x1 + x2 + x3 + x4

alpha.init <- 1.5
beta.init  <- rep(0, 4)

#methods <- c ( "LOGIT", "SPLINE", "DT")
              
methods <- c("LOGIT", "SPLINE", "DT", "SVM",  "NN",  "XGB", "RF")#, "RF", "XGB", "NN"

## ---- 1. Out-of-fold predictions, once per method ----

method_seeds <-c(LOGIT = 2022, SPLINE = 2022, SVM = 2022, DT = 2022, RF = 2022, XGB = 2022, NN = 2022)

cv_results <- list()
cv_fits <- list()

for (method_name in methods) {
  set.seed(method_seeds[method_name])
  cv_fits[[method_name]] <- vector("list", K)
  
  pred_uncure <- rep(NA_real_, n)
  pred_Sp <- rep(NA_real_, n)
  pred_S1 <- rep(NA_real_, n)
  pred_S <- rep(NA_real_, n)
  pred_time <- rep(NA_real_, n)
  
  func <- get(paste0("smcure.", method_name, ".ic"))
  
  for (k in 1:K) {
    
    train_data <- H[fold_id != k, , drop = FALSE]
    test_data  <- H[fold_id == k, , drop = FALSE]
    test_index <- which(fold_id == k)
    
    fit <- func(
      formula     = surv_formula,
      cureform    = cure_formula,
      alpha.init  = alpha.init,
      beta.init   = beta.init,
      train_data  = train_data,
      test_data   = test_data,
      Var         = F,
      emmax       = 1000,
      eps         = 1e-3
    )
    cv_fits[[method_name]][[k]] <- fit
    
    pred_uncure[test_index] <- fit$test$Uncureprob
    pred_Sp[test_index]     <- fit$test$Sp
    pred_S1[test_index]     <- fit$test$S1
    pred_S[test_index]      <- fit$test$S
    pred_time[test_index]   <- fit$test$Time.mid
  }
  
  cv_results[[method_name]] <- list(
    Uncureprob = pred_uncure,
    Sp         = pred_Sp,
    S1         = pred_S1,
    S          = pred_S,
    Time.mid   = pred_time
  )
}
## ---- 2. Multiple imputation of the latent cure status, once ----

set.seed(2032)
B <- 500
imputed_status <- list()

for (method_name in methods) {
  
  pred_pi <- cv_results[[method_name]]$Uncureprob
  pred_S  <- cv_results[[method_name]]$S
  
  w_post <- H$D + (1 - H$D) * (1 - (1 - pred_pi)^pred_S)
  w_post <- pmin(pmax(w_post, 1e-6), 1 - 1e-6)
  
  V_imp <- matrix(NA_integer_, nrow = n, ncol = B)
  cens_idx <- which(H$D == 0)
  
  for (b in 1:B) {
    V_imp[, b] <- H$D
    V_imp[cens_idx, b] <- rbinom(length(cens_idx), size = 1, prob = w_post[cens_idx])
  }
  
  imputed_status[[method_name]] <- list(w_post = w_post, V_imp = V_imp)
}

## ---- 3. ROC / AUC, once, printed once ----

roc_results <- list()
auc_results <- data.frame(Method = methods, Mean_AUC = NA_real_, SD_AUC = NA_real_)

for (method_name in methods) {
  
  pred_pi <- cv_results[[method_name]]$Uncureprob
  V_imp   <- imputed_status[[method_name]]$V_imp
  
  auc_vec  <- rep(NA_real_, B)
  roc_list <- vector("list", B)
  
  for (b in 1:B) {
    roc_b <- roc(response = V_imp[, b], predictor = pred_pi, quiet = TRUE, direction = "<")
    auc_vec[b] <- as.numeric(auc(roc_b))
    roc_list[[b]] <- roc_b
  }
  
  roc_results[[method_name]] <- list(roc = roc_list, auc = auc_vec)
  
  auc_results[auc_results$Method == method_name, "Mean_AUC"] <- mean(auc_vec, na.rm = TRUE)
  auc_results[auc_results$Method == method_name, "SD_AUC"]   <- sd(auc_vec, na.rm = TRUE)
}

print(auc_results)



## ---- 4. Plot mean ROC curves with AUC ----

method_colors <- c(LOGIT = "black", SPLINE = "orange", SVM = "blue", DT = "green4", RF = "purple", XGB = "brown", NN = "red")

fpr_grid <- seq(0, 1, length.out = 1001)

plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "False Positive Rate", ylab = "True Positive Rate", main = "", xaxs = "i", yaxs = "i", asp = 1)

abline(a = 0, b = 1, lty = 2, col = "gray")

for (method_name in methods) {
  
  roc_list <- roc_results[[method_name]]$roc
  tpr_matrix <- matrix(NA_real_, nrow = length(fpr_grid), ncol = length(roc_list))
  
  for (b in seq_along(roc_list)) {
    
    roc_b <- roc_list[[b]]
    fpr <- 1 - roc_b$specificities
    tpr <- roc_b$sensitivities
    
    tpr_matrix[, b] <- approx(x = fpr, y = tpr, xout = fpr_grid, ties = max, rule = 2)$y
  }
  
  mean_tpr <- rowMeans(tpr_matrix, na.rm = TRUE)
  
  lines(c(0, fpr_grid, 1), c(0, mean_tpr, 1), col = method_colors[method_name], lwd = 2)
}

legend_labels <- sprintf("%s (AUC = %.4f)", auc_results$Method, auc_results$Mean_AUC)

legend("bottomright", legend = legend_labels, col = method_colors[auc_results$Method], lwd = 2, bty = "n", cex = 0.85)
## =====================================================
## MIDPOINT KAPLAN-MEIER WITH 95% CONFIDENCE INTERVALS
## =====================================================



## 1. Prepare midpoint event times
km_data <- H

## Identify observed events and right-censored observations
km_data$event <- as.integer(is.finite(km_data$R))

## Midpoint for interval-censored events
## Last observed time for right-censored observations
km_data$time <- ifelse(km_data$event == 1, (km_data$L + km_data$R) / 2, km_data$L)

## Check the constructed data
print(table(km_data$event))
print(summary(km_data$time))

## 2. Fit Kaplan-Meier estimator
fit_km <- survival::survfit(Surv(time, event) ~ 1, data = km_data, conf.int = 0.95, conf.type = "log-log")

## 3. Extract survival estimates and confidence intervals
plot_data <- data.frame(time = fit_km$time, survival = fit_km$surv, lower = fit_km$lower, upper = fit_km$upper)

## Add starting point at time zero
plot_data <- rbind(data.frame(time = 0, survival = 1, lower = 1, upper = 1), plot_data)



## 4. Plot Kaplan-Meier curve with 95% confidence intervals
km_plot <- ggplot(plot_data, aes(x = time, y = survival)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#E99689", alpha = 0.22, na.rm = TRUE) +
  geom_step(color = "#E57368", linewidth = 0.85, direction = "hv") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Time", y = "Survival probability", title = "") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(size = 12), axis.text = element_text(size = 11), axis.line = element_line(linewidth = 0.5), axis.ticks = element_line(linewidth = 0.5))

print(km_plot)



## =====================================================
## 4. Full-data Weibull latency estimation
##    Bootstrap SE, Z-values and P-values
## =====================================================

set.seed(2040)

nboot <- 100

latency_fits <- list()
latency_results <- list()

for (method_name in methods) {
  
  cat("\n====================================\n")
  cat("Latency analysis:", method_name, "\n")
  cat("Bootstrap replications:", nboot, "\n")
  cat("====================================\n")
  
  set.seed(method_seeds[method_name])
  
  func <- get(paste0("smcure.", method_name, ".ic"))
  
  fit_lat <- tryCatch(
    func(formula = surv_formula, cureform = cure_formula, alpha.init = alpha.init, beta.init = beta.init, train_data = H, test_data = H, Var = F, nboot = nboot, emmax = 1000, eps = 1e-3),
    error = function(e) {
      warning(paste("Latency fit failed for", method_name, ":", conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(fit_lat)) next
  
  latency_fits[[method_name]] <- fit_lat
  
  estimates <- as.numeric(fit_lat$latency)
  se <- as.numeric(fit_lat$latency_sd)
  
  if (is.null(fit_lat$latency_names)) {
    parameter_names <- c("alpha", paste0("beta", seq_along(estimates[-1])))
  } else {
    parameter_names <- fit_lat$latency_names
  }
  
  zvalue <- estimates / se
  pvalue <- 2 * pnorm(-abs(zvalue))
  
  latency_results[[method_name]] <- data.frame(
    Method = method_name,
    Parameter = parameter_names,
    Estimate = estimates,
    SE = se,
    Z_value = zvalue,
    P_value = pvalue,
    row.names = NULL
  )
  
  print(latency_results[[method_name]])
}

## =====================================================
## 5. Combined latency comparison table
## =====================================================

latency_table <- do.call(rbind, latency_results)

if (!is.null(latency_table)) {
  rownames(latency_table) <- NULL
  print(latency_table)
}

## =====================================================
## 6. Save results
## =====================================================

if (!is.null(latency_table)) {
  write.csv(latency_table, "Latency_Comparison.csv", row.names = FALSE)
}

saveRDS(latency_fits, "Latency_Fitted_Models.rds")
saveRDS(latency_results, "Latency_Results.rds")

## =====================================================
## 7. Combined AUC and latency results
## =====================================================

print(auc_results)
print(latency_table)





time_taken <- proc.time() - ptm
time_taken