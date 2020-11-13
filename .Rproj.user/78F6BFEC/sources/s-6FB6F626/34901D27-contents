#' @import dplyr 
#' @import rlist
#' @import progress 
#' @import purrr
#' @importFrom truncnorm rtruncnorm
NULL


#' Clustering Variying Coefficients with a Pitman-Yor Process
#' 
#' @param ID Subjects' ID vector
#' @param W A part of design matrix for varying coefficients
#' @param X A part of design matrix for overall fixed effects 
#' @param Z A part of design matrix for random effects
#' @param t A time vector
#' @param Y A binary response vector
#' @param num_knots The number of knots for splines
#' @param K.max The upper bound of the number of clusters for the stick-breaking process
#' @param num_iters The number of iterations for the partially collapsed Gibbs sampler
#' @param sampling_unit Thinning unit
#' 
#' @export
cvarpyp <- function(ID, W, X=NULL, Z, t, Y, 
                    num_knots=10, K.max, num_iters=10000, sampling_unit=5,
                    nu=1, lambda=0.5, Psi=NULL, kappa=10, 
                    sd_nu = 1, sd_lambda = 1,
                    pi.a=1, pi.b=1,
                    stick1 = c(1,10), stick2 = c(1,10)){

  temp <- preprocess_data(ID, W, X, Z, t, Y)
  ID = temp$ID
  W = temp$W
  X = temp$X
  Z = temp$Z
  t = temp$t
  Y = temp$Y
  rm(temp)
  
  num_obs_per_sbj <- table(ID)
  NSbj <- max(ID)
  one_obsr <- which(num_obs_per_sbj==1)
  
  # translate to matrix values
  p <- ncol(W);  r <- ncol(Z) ; q <- ncol(X)
  
  Z_sbj <- map(split(Z %>% data.matrix, f=ID), ~ matrix(.x, ncol=r))
  X_sbj <- map(split(X %>% data.matrix, f=ID), ~ matrix(.x, ncol=q))
  XSum <- map(1:NSbj, ~ t(X_sbj[[.x]])%*%X_sbj[[.x]]) %>% reduce(`+`)
  
  
  lower.tn = log(Y)
  upper.tn = -log(1-Y)
  
  # beta-binomial
  a_pi <-pi.a ; b_pi <- pi.b  # noninformative
  nu_param = 1

  # rinvGamma
  n_g_k <- stick1[1] ; n_h_k <- stick1[2] # noninformative
  i_g_k <- stick2[2] ; i_h_k <- stick2[2]
  
  knots <- quantile(t, ppoints(num_knots, a=0))
  

  # Check Psi 
  if (is.null(Psi)) {
    Psi <- diag(r)
  }else {
    if (is.matrix(Psi)) {
      if (((ncol(Psi)==nrow(Psi))==r)) {
        Psi <- Psi 
      } else {
        print(paste0("Psi should be a positive definite matrix with","(",r,",",r,")"))
      }
    }else {
      print(paste0("Psi should be a positive definite matrix with","(",r,",",r,")"))
    }
  }
  Psi_0 <- diag(r)
  beta <- rep(0,q)

  
  bi <- map(1:NSbj, ~ MASS::mvrnorm(1, mu=rep(0,r), Sigma=Psi))
  Li <- map_dbl(Y, function(sign) if(sign>0){rtruncnorm(1, a=0, mean=2, sd=0.5)}else{rtruncnorm(1,b=0,mean=-2,sd=0.5)}) %>% as.matrix()
  L_sbj <- split(Li, f=ID)
  
  pu <- r ; D <- diag(r)
  
  g_vec <- rep(n_g_k, K.max); h_vec <- rep(n_h_k, K.max)
  tau <- MCMCpack::rinvgamma(n=K.max, shape=n_g_k/2, scale=n_h_k/2)
  gamma <- map(1:K.max, function(x) {
    temp <- matrix(rbinom(n = p*(num_knots + 2), size=1, prob=1/(num_knots+2)),nrow=p)
    temp[,1]<-1 
    temp
  })
  
  v <- c(map_dbl(1:(K.max-1), ~ rbeta(1, 1, 1)),1)
  Diri_p <- get_DP_pi(v = v, K = K.max)
  Clusters <- sample(1:K.max, size=NSbj, replace=T)
  
  active_ <- 1:K.max %in% Clusters  
  active <- which(active_==T) 
  inactive <- which(active_==F)
  
  theta_k <- map(1:K.max, function(k) MASS::mvrnorm(1, mu=rep(0,sum(gamma[[k]])), Sigma = diag(sum(gamma[[k]]))))
  
  W_star <- get_basis_data(tData = W %>% data.matrix, t = t, knots = knots, sd = sd(t))
  C_star <- do.call(cbind, W_star)
  rm(W_star)
  
  C_star_sbj <- split(C_star %>% data.matrix, f=ID) %>% map( ~ matrix(.x, ncol=dim(C_star)[2]) )
  
  temp <- byproducts(num_obs_per_sbj, K.max, NSbj, C_star_sbj, gamma, X_sbj, beta, Z_sbj, L_sbj, Psi, Psi_0, ID-1, Clusters, active, flatten_gamma = flatten_gamma_with_fixC)
  InvAdj <- temp[[1]]
  InvAdj_0 <- temp[[2]]
  XXi_k <- temp[[3]]
  Xi_k <- temp[[4]]
  
  temp <- MakeRk(K.max, tau, active, C_star_sbj, InvAdj, InvAdj_0, gamma, NSbj, flatten_gamma =flatten_gamma_with_fixC)
  Rk <- temp[[1]]

  updated_theta <-list()
  updated_beta <- list()
  updated_gamma <- list()
  updated_cluster <- list()
  updated_tau <-list()
  updated_Psi <- list()
  update_indi <- 1
  
  # nu_store <- c()
  # nu_accept_store <- c()
  # lambda_store <- c()
  # lambda_accept_store <- c()
  
  pb <- progress_bar$new(
    format = "  running pcg [:bar] :percent eta: :eta",
    clear = FALSE, total = num_iters)
  
  
  for (iter in 1:num_iters) {

    pb$tick()
    
    # ================== activation_clusters,boolean ===================== #
    active_ <- 1:K.max %in% Clusters  
    active <- which(active_==T) 
    inactive <- which(active_==F)
    
    gamma <- sample_gamma(ID, K.max, p, active, num_knots, a_pi, b_pi, NSbj, gamma, tau, C_star_sbj, InvAdj, L_sbj, X_sbj, beta, Clusters)
    
    # Draw stick breaking beta
    num_of_ones_K_sbjs <- get_sbj_clusters(K.max, Clusters)
    v <-  c(map_dbl(1:(K.max-1), ~ rbeta(1, 1-lambda+num_of_ones_K_sbjs[.x], nu + .x*lambda+sum(num_of_ones_K_sbjs[(.x+1):K.max]))) ,1)
    Diri_p <- get_DP_pi(v, K.max)
    
    # Draw nu
    nu <- sample_nu(nu, lambda, v, sd_nu, K.max, nu_param)
    # nu <- 1.0
    # nu_accept_store <- append(nu_accept_store, nu)

    # Draw lambda
    lambda <- sample_lambda(lambda, nu, v, K.max, sd_lambda)
    # lambda <- 0.1
    #lambda_accept_store <- append(lambda_accept_store,proposed_lambda)
    temp <- byproducts(num_obs_per_sbj, K.max, NSbj, C_star_sbj, gamma, X_sbj, beta, Z_sbj, L_sbj, Psi, Psi_0, ID-1, Clusters, active, flatten_gamma = flatten_gamma_with_fixC)
    InvAdj <- temp[[1]]
    InvAdj_0 <- temp[[2]]
    XXi_k <- temp[[3]]
    Xi_k <- temp[[4]]
    
    temp <- MakeRk(K.max, tau, active, C_star_sbj, InvAdj, InvAdj_0, gamma, NSbj, flatten_gamma =flatten_gamma_with_fixC)
    Rk <- temp[[1]]
    
    C_star_sbj_k <- map(1:NSbj, ~ C_star_sbj[[.x]][,which(flatten_gamma_with_fixC(gamma[[Clusters[.x]]])==TRUE)])
    
    # Draw theta
    theta_k <- sample_theta(K.max, active, XXi_k, Xi_k, Rk, tau)
    theta_sbj_k <- map(1:NSbj, ~ theta_k[[Clusters[[.x]]]])
    
    # Draw tau
    num_of_ones_vars <- map_dbl(1:K.max, ~ sum(gamma[[.x]]))
    g_vec <- map_dbl(1:K.max, ~ ifelse((.x%in%active), n_g_k, i_g_k))
    h_vec <- map_dbl(1:K.max, ~ ifelse((.x%in%active), n_h_k, i_h_k))
    tau <- map_dbl(1:K.max, ~ MCMCpack::rinvgamma(1, shape=(num_of_ones_vars[.x]+g_vec[.x])/2, scale=h_vec[.x]/2 + (theta_k[[.x]] %*% Rk[[.x]] %*% theta_k[[.x]])/2))
    
    # Draw bi 
    temp <- sample_randomeffects(1,Z_sbj,InvAdj,L_sbj,C_star_sbj,gamma,X_sbj, beta, Psi,theta_k,Clusters,r, flatten_gamma = flatten_gamma_with_fixC)
    bi <- temp[[1]]
    bi_mats <- temp[[2]]

    # Draw Psi
    proposed_Psi <- MCMCpack::riwish(pu+NSbj, D+bi_mats)
    probs_Psi <- MakeRkPsi(num_obs_per_sbj, K.max, Rk, tau, active, C_star_sbj, Z_sbj, proposed_Psi,gamma,theta_k,NSbj,flatten_gamma = flatten_gamma_with_fixC)
    if(runif(1) < sum(unlist(probs_Psi[[1]]))) Psi <- proposed_Psi
    
    # Draw Li
    L_sbj <- sample_L(ID, C_star_sbj_k, theta_sbj_k, Z_sbj, X_sbj, bi, beta, lower.tn, upper.tn)
    
    # Draw cluster
    out_mat <- ObtainOutput(NSbj, K.max, Diri_p , L_sbj,C_star_sbj, X_sbj, beta, Z_sbj,bi,gamma,theta_k,flatten_gamma = flatten_gamma_with_fixC)
    Clusters <- sample_cluster(out_mat, K.max)
    
    # Draw beta 
    if (!is.null(X)) {
      beta <- sample_beta(kappa, q, XSum, NSbj, X_sbj, L_sbj, C_star_sbj_k, theta_sbj_k, Z_sbj, bi)
    } else {
      beta <- 0
    }
    
    # ======================= store variables ======================#
    if(iter %% sampling_unit == 0){
      updated_gamma<-list.append(updated_gamma, gamma)
      updated_cluster<-list.append(updated_cluster, Clusters)
      updated_theta<-list.append(updated_theta, theta_k)
      updated_tau<-list.append(updated_tau, tau)
      updated_Psi <- list.append(updated_Psi, Psi)
      if (!is.null(X)) updated_beta <- list.append(updated_beta, beta)
    }
  }
  list(gamma.sample = updated_gamma, 
       cluster.sample = updated_cluster,
       theta.sample = updated_theta,
       tau.sample = updated_tau, 
       Psi.sample = updated_Psi,
       beta.sample = updated_beta,
       knots = knots,
       time.range = c(min(t),max(t)),
       p = p,
       q = ifelse(!is.null(X), dim(X)[2],0),
       r = r
       ) 
}
