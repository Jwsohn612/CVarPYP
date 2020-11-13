#' @import dplyr 
#' @import rlist 
#' @import progress 
#' @import purrr
NULL

#' Draw the posterior results of random effects
#' 
#' @export
summarise_Psi <- function(Psi.sample, row.position, col.position, burn=NULL){
  len <- length(Psi.sample)
  if (is.null(burn)) burn = round(len/2)
  post.sample <- map_dbl(Psi.sample, ~ .x[row.position, col.position])
  post.sample <- post.sample[(burn+1):len]

  par(mfrow=c(1,2))
  temp.hist=hist(post.sample, breaks=10, probability=T,main="",xlab="",ylab="",col=gray(0:9/9)[8],border=gray(0:9/9)[8]);box()
  lines(temp.hist$breaks,c(temp.hist$density,0),type="s")
  abline(v=quantile(post.sample,c(0.025,0.975)),lty=2)
  abline(v=median(post.sample))
  mtext("Posterior density",side=2,line=2.3,cex=0.9)
  
  acf(post.sample,main="")
} 

#' Draw the posterior results of fixed effects
#' 
#' @export
summarise_beta <- function(beta.sample, position, burn=NULL) {
  len <- length(beta.sample)
  if (is.null(burn)) burn = round(len/2)
  post.sample <- map_dbl(beta.sample, ~ .x[position])
  post.sample <- post.sample[(burn+1):len]
  
  par(mfrow=c(1,2))
  temp.hist=hist(post.sample, breaks=10, probability=T,main="",xlab="",ylab="",col=gray(0:9/9)[8],border=gray(0:9/9)[8]);box()
  lines(temp.hist$breaks,c(temp.hist$density,0),type="s")
  abline(v=quantile(post.sample,c(0.025,0.975)),lty=2)
  abline(v=median(post.sample))
  mtext("Posterior density",side=2,line=2.3,cex=0.9)
  
  acf(post.sample,main="")
}

#' Draw the posterior results of clusters
#' 
#' @export
summarise_cluster <- function(cluster.sample, burn=NULL) {
  len <- length(cluster.sample)
  if (is.null(burn)) burn = round(len/2)
  cluster.matrix <- do.call(rbind, cluster.sample)
  apply(cluster.matrix[(burn+1):len,], 2, function(x) {
    mode.value <- DescTools::Mode(x)
    mode.value[1]
  }) %>% unlist()
}

get_theta_gamma_matrix <- function(p, unit_gamma, unit_theta){
  positions <- c(1,apply(unit_gamma, 1, sum))
  indexing <- map(1:p, function(i) {
    start <- positions[1:i] %>% sum
    end <- start + positions[i+1] - 1
    c(start,end)
  })
  
  do.call(rbind, map(1:p, function(i) {
    index <- indexing[[i]]
    temp_gamma <- unit_gamma[i,]
    temp_gamma[temp_gamma==1] <- unit_theta[index[1]:index[2]]
    temp_gamma
  }))
}

#' Get the posterior results of varying coefficients
#' 
#' @export
get_post_summary_vc <- function(time_range, knots, gamma.sample, theta.sample, grids=50, burn=NULL) {
  
  len <- length(gamma.sample)
  if (is.null(burn)) burn = round(len/2)
  
  time.domain <- seq(from=time_range[1], to=time_range[2], l=grids)
  
  sd.t <- sd(t)
  K <- length(gamma.sample[[1]])
  p <- dim(gamma.sample[[1]][[1]])[[1]]
  

  basis_functions <- map(time.domain, ~ matrix(rep(B(.x, knots, sd.t),p),nrow=p, byrow=T )) 
  
  # Time => MCMC => Cluster
  Bgamma <- map(1:length(time.domain), function(time) {
    map(gamma.sample[(burn+1):len], function(gamma)  map(gamma, function(each_cls) {
      bgam.mat <- each_cls * basis_functions[[time]] # p times knots matrices
      
    }))
  })
  
  Tgamma <-  map((burn+1):len, function(MC) {
    map(1:K, function(k) get_theta_gamma_matrix(p, gamma.sample[[MC]][[k]], theta.sample[[MC]][[k]]) %>% t)  
  })
  
  vc_summary <-  map(1:K, function(k) {
    
         map(1:length(time.domain), function(time) {
           
            temp_mat<- do.call(rbind, 
                               map(1:length((burn+1):len), function(MC) {
                                diag(Bgamma[[time]][[MC]][[k]] %*% Tgamma[[MC]][[k]])
                                 }))
            
            apply(temp_mat,2,function(x) quantile(x, c(0.025,0.5,0.975)))
            
                          })
        })
  
  list(vc_summary = vc_summary,
       time.domain = time.domain)
}

#' Draw the posterior results of varying coefficients
#' 
#' @export
summarise_varying_coefficient <- function(vc_object, cluster_number, variable_number){
  vc_summary <- vc_object$vc_summary
  time.domain <- vc_object$time.domain
  
  vc_mat <- do.call(rbind,map(vc_summary[[cluster_number]], ~.x[,variable_number]))
  L <- vc_mat[,1]
  M <- vc_mat[,2]
  U <- vc_mat[,3]
  
  poly_range <- c(time.domain, rev(time.domain))
  poly_coef_UL <- c(L, rev(U))
  
  plot(NULL,type="l",ylim=c(min(vc_mat),max(vc_mat)),xlim=c(min(time.domain),max(time.domain)), xlab="",ylab="")
  polygon(poly_range, poly_coef_UL,col=gray(0:9/9)[8],border=F)
  lines(time.domain, M,lty=1)
  
  # mtext(paste0('cluster',which_clusters,';varying',j),side=2,line=2.3,cex=0.9)
  # mtext("t",side=1,line=2.3,cex=0.9)
}

# Time => Cluster => VC
# vc_object <- get_post_summary_vc(from=0, to=1, knots, gamma.sample, theta.sample, grids=50,  burn=NULL)
# summarise.varying.coefficient(vc_object, cluster_number = 1, variable_number = 3)

#' Draw the posterior results of latent location parameters
#' 
#' @export
summarise_latent_location <- function(time_range, knots, gamma.sample, cluster_number, variable_number, burn=NULL) {
  
  len <- length(gamma.sample)
  if (is.null(burn)) burn = round(len/2)
    
  location_mat <- do.call(rbind, map(gamma.sample[(burn+1):len], ~ .x[[cluster_number]][variable_number,]))
  proba <- apply(location_mat,2,mean)
  
  df <- data.frame(knots=c(time_range[1],knots), freq=proba[2:(length(knots)+2)])

  plot(df$knots[2:(length(knots)+1)], df$freq[2:(length(knots)+1)], type='h', ylim=c(0,1), xlim=c(time_range[1],time_range[2]), lwd=3,xlab="",ylab="")
  points(0,df$freq[1],pch= 1)
  segments(0,0,df$knots[1],df$freq[1],lty=2)
  # mtext("klp",side=2,line=2.5,cex=0.9)
  # mtext("t",side=1,line=2.5,cex=0.9)
}

# summarise.latent.location(time_range=c(0,1), gamma.sample, cluster_number=1, variable_number=2, burn=NULL)

