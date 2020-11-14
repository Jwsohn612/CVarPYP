library(cvarpyp)
library(dplyr)

data <- get_simulation_data(N=300, M=10)

# Input 
ID <- data %>% pull(ID)
W <- data %>% select(c('var1','var2','var3','var4'))
X <- data %>% select(c('fix1','fix2'))
Z <- data %>% select(c('re1','re2'))
t <- data %>% pull(t)
Y <- data %>% pull(Y)

output <- cvarpyp(ID=ID,
                  W=W,
                  X=X,
                  Z=Z,
                  t=t,
                  Y=Y, 
                  K_max=5,
                  num_knots = 10,
                  num_iters = 100)

plot(output$random_effect, row_position=1, col_position=1, burn=NULL)

plot(output$fixed_effect, position=2, burn=NULL)

plot(output$varying_coefficient, cluster_number = 5, variable_number = 1)

plot(output$latent_location, cluster_number = 5, variable_number = 1, time_range=output$time_range, knot_position=output$knot_position, burn=NULL)

# plot(output$cluster)

