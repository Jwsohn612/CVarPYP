remove.packages('cvarpyp')

devtools::install_github('jwsohn612/cvarpyp')

library(cvarpyp)

data <- get_simulation_data(N=600, K=3, M=10)

# Input 
ID <- data$dataset$ID
W <- data$dataset[c('var1','var2','var3','var4')]
X <- data$dataset[c('fix1','fix2')]
Z <- data$dataset[c('re1','re2')]
t <- data$dataset$t
Y <- data$dataset$Y

model.output <- cvarpyp(ID=ID,
                        W=W,
                        X=X,
                        Z=Z,
                        t=t,
                        Y=Y, 
                        K.max=6,
                        num_knots = 10,
                        num_iters = 1000)

summarise_Psi(Psi.sample=model.output$Psi.sample, 
              row.position=1, col.position=2)

summarise_beta(beta.sample=model.output$beta.sample, 
               position=2, burn=NULL)

summarise_cluster(cluster.sample = model.output$cluster.sample)

vc_object <- get_post_summary_vc(time_range = model.output$time.range, 
                                 knots = model.output$knots, 
                                 gamma.sample = model.output$gamma.sample, 
                                 theta.sample = model.output$theta.sample, grids=50,  burn=NULL)

summarise_varying_coefficient(vc_object, cluster_number = 5, variable_number = 1)

summarise_latent_location(time_range = model.output$time.range, 
                          knots = model.output$knots,
                          gamma.sample = model.output$gamma.sample, 
                          cluster_number = 1, variable_number = 2, burn=NULL)

