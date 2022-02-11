###############################################################################
# Sensitivity analysis, and plotting covarage when changing the prior on b0. 
###############################################################################

set.seed(seed)


prior.means=b0+rep(sort(rep(c(-1.8,-1.2,-0.6,0,0.6,1.2,1.8),6)),3) # Vary the mean of the prior for beta_0.
prior.sd=rep(rep(c(1.2,1,0.8,0.6,0.4,0.2),7),3) # Vary the standard deviation of the prior for beta_0.
v.index=sort(rep(c(1,3,4),42)) # Vary the strength of the under-reporting covariate.


start_time_mcmc <- Sys.time()

sim.code.mcmc=nimbleCode({
  for(i in 1:N){
    lambda[i] <- exp(a0+a1*x[i]+phi[i]) 
    pi[i] <- ilogit(b0+b1*v[i]) 
    z[i] ~ dpois(lambda[i]*pi[i])
  }
  phi[1:N] ~ dcar_normal(adj=adj[1:l_adj],num=n_adj[1:N],tau=tau,zero_mean = 1)
  a0 ~ dnorm(0,sd=10)
  a1 ~ dnorm(0,sd=10)
  b0 ~ dnorm(prior_mean,sd=prior_sd)
  b1 ~ dnorm(0,sd=10) 
  tau ~ dgamma(1, 0.00005)
  #epsilon ~ T(dnorm(0,sd=1),0,) #REMOVE THIS, NOT USED IN MODEL
}) 



# Setup NIMBLE.

sim.constants=list(N=sim$N,x=sim$x,adj=sim$adj,l_adj=length(sim$adj),n_adj=sim$n_adj)
sim.data=list(z=sim$z,prior_sd=prior.sd[1],prior_mean=prior.means[1],v=sim$v[,v.index[1]])
sim.inits=list(a0=0,a1=0,b0=0,b1=0, phi=rep(0,sim$N), tau = 1) 

sim.model <- nimbleModel(sim.code.mcmc, sim.constants,sim.data,sim.inits)
sim.compiled.model<-compileNimble(sim.model,resetFunctions = TRUE)

sim.mcmc.conf <- configureMCMC(sim.model,
                               monitors=c('a0','a1','b0','b1','tau', 'lambda', 'pi'),
                               useConjugacy = FALSE)
sim.mcmc.conf$removeSamplers(c('a0','a1','b0','b1', 'tau')) 
sim.mcmc.conf$addSampler(target=c('a0','a1','b0','b1','tau'),type='AF_slice',control=list(adaptInterval=1000))

sim.mcmc<-buildMCMC(sim.mcmc.conf)
sim.compiled.mcmc<-compileNimble(sim.mcmc, project = sim.model,resetFunctions = TRUE)



print("Running MCMC and calculating the coverage when changing the prior distribution on beta for the three examples used in the specialisation project)")

sim.mcmc.list=list()
for(i in c(1:(length(v.index)))){
  print(i)
  sim.compiled.model$prior_mean=prior.means[i]
  sim.compiled.model$prior_sd=prior.sd[i]
  sim.compiled.model$v=sim$v[,v.index[i]]
  sim.samples<-runMCMC(sim.compiled.mcmc,inits=sim.inits,
                       nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                       summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)
  sim.mcmc.list[[i]]=as.mcmc(as.mcmc.list(sim.samples))
  
}

end_time_mcmc <- Sys.time()
print(c("Time taken to run",length(v.index)," iterations of the MCMC algorithm:", difftime(end_time_mcmc,  start_time_mcmc, units = "secs")))
print("Finished:")
print(Sys.time())


#####################################################################################
# Generating y samples from MCMC, computing coverage, and computing lambda error and bias.
#####################################################################################

# Produce y samples.
sim.y=lapply(sim.mcmc.list,function(x){
  lambda=x[,5:(4+sim$N)]
  pi=x[,(5+sim$N):(4+2*sim$N)]
  y=matrix(nrow=dim(lambda)[1],ncol=sim$N)
  for(j in 1:sim$N){
    y[,j]=rpois(dim(lambda)[1],lambda[,j])
  }
  return(y)
})


# Compute coverage of 95% prediction intervals.
coverage=unlist(lapply(sim.y,function(x){
  covered=numeric(sim$N)
  for(j in 1:sim$N){
    covered[j]=sim$y[j]>=quantile(x[,j],0.025)&sim$y[j]<=quantile(x[,j],0.975)
  }
  return(mean(covered))
}))

# Compute mean error of log(lambda).
lambda.bias=unlist(lapply(sim.mcmc.list,function(x){
  log.lambda.means=apply(log(x[,5:(4+sim$N)]),2,mean)
  return(mean(log.lambda.means-log(sim$lambda)))
}))

# Compute root mean squared error of log(lambda).
lambda.error=unlist(lapply(sim.mcmc.list,function(x){
  log.lambda.means=apply(log(x[,5:(4+sim$N)]),2,mean)
  return(sqrt(mean((log.lambda.means-log(sim$lambda))^2)))
}))

mcmc.data=data.frame(x=prior.means,y=prior.sd,z=v.index,cov=cor(sim$v)[v.index],
                     c=coverage,l=lambda.error,b=lambda.bias)


##################################################
# Setting up and running inlabru
##################################################

run.inlabru.coverage <- function(means, sd, v.index, dump = FALSE){
  print("If dump == TRUE, then all output from running inlabru here is dumped to file.")
  print("Then, it will not show up when running R")
  print(c("Here, dump ==", dump))
  
  #Dump into file
  if(dump){
    sink(file = "inlabru_output.txt")
    sink(stdout(),type="message")
    sink(stdout(),type= "output")
  }
  coverage.tot <- rep(0, length(v.index))
  inla.lambda.bias.tot <- rep(0, length(v.index))
  inla.lambda.error.tot <- rep(0, length(v.index))
  start_time_inla <- Sys.time()
  for (i in 1:(length(v.index))){ 
    print("------------------------------------------------------")
    print(c("Run no: ",i, "Mean: ", means[i], "SD: ",sd[i]))
    
    df = data.frame(yy = sim$z, true = sim$y, loc = 1:100,
                    x = sim$x, v = sim$v[,v.index[i]]) #This is the correct place to change v, depending on correlation with w
    
    Q = sim$A
    dd  = diag(apply(Q,1,sum))
    Q = dd-Q
    func = function(v,beta0,beta1)
    {
      aa = beta0 + beta1 * v
      return(-log(1+exp(-aa)))
    }
    
    cmp = ~ 
      covariate(x) +
      Intercept(1)+
      beta0(main = 1, model = "linear",
            mean.linear = means[i],
            prec.linear = 1/(sd[i])^2) +
      beta1(main = 1, model = "linear",
            mean.linear = 0,
            prec.linear = 0.1^2)+ 
      spatial(loc, model = "besag", graph = Q,constr = T,
              hyper = list(prec = list(param = c(1, 0.00005)))) 
    
    formula = yy   ~ 
      Intercept +
      covariate  + 
      spatial+ 
      func(v,beta0, beta1)
    
    lik = like("poisson",
               formula = formula,
               data= df)
    
    bru_options_set(bru_verbose = TRUE, bru_max_iter = 30)
    
    
    
    fit <- bru(components = cmp,  
               lik,
               options = list(verbose = F,
                              num.threads = "1:1"))
    
    
    #Generate lambda
    sim$post_lambda <- generate(fit, data = df, formula =~ exp(Intercept + covariate + spatial), n.samples = 1000)
    
    
    #Generate pi 
    sim$post_pi <- generate(fit, data = df, formula =~ expit(beta0 + beta1 * v), n.samples = 1000)
    
    #Generate sorted pi. 
    sim$post_pi_sorted <- generate(fit, data = df, formula =~ expit(beta0 + beta1 * sort(v)), n.samples = 1000)
    
    
    generate.y.samples = function(lambda, pi = FALSE){
      if (pi != FALSE){
        y <- sim$z + rpois(length(lambda), lambda*(1-pi)) 
      } else{
        y <- rpois(length(lambda), lambda)
      }
      dim(y) <- dim(lambda)
      return (y) 
    }
    
    sim$post.y <- generate.y.samples(sim$post_lambda)
    
    
    
    #Create an estimated 95% PI for y.
    #Compute upper and lower limit of the 95% prediction interval for y.
    post.pred.int <- rowQuantiles(sim$post.y, rows = NULL, cols = NULL, probs = c(0.025, 0.975))
    #post_pred_int
    # From inlabru, Compute coverage using a 95% prediction interval.
    
    covered = numeric(sim$N)
    covered=sim$y>=post.pred.int[,1]&sim$y<=post.pred.int[,2]
    print(covered)
    print(c("Coverage:",mean(covered)))
    coverage.tot[i] <- mean(covered)
    
    
    
    #Computing the mean error of log(lambda) and the root mean squared error of log(lambda)
    
    inla.lambda.bias.tot[i] <- mean(rowMeans(log(sim$post_lambda)) - log(sim$lambda))
    inla.lambda.error.tot[i] <- sqrt( mean((rowMeans(log(sim$post_lambda)) - log(sim$lambda))^2))
    
    
  } 
  
  end_time_inla <- Sys.time()
  #Stopping the dumping of all output to file.
  if(dump){
    print(c("Time taken for inlabru to run the entire sensitivity analysis for all three examples:", difftime(end_time_inla,  start_time_inla, units = "secs"), "seconds"))
    sink(type = "output")
    sink(type = "message")
  }
  print(c("Time taken for inlabru to run the entire sensitivity analysis for all three examples:", difftime(end_time_inla,  start_time_inla, units = "secs"), "seconds"))
  print("Finished:")
  print(Sys.time())
  df <- list(coverage.tot, inla.lambda.error.tot, inla.lambda.bias.tot)
  return(df)
}

print("Running the sensitivity analysis using INLA (inlabru) for the three examples used in the specialisation project)")

inlabru.coverage.sim <-run.inlabru.coverage(prior.means,prior.sd, v.index, dump = TRUE)


inlabru.data <- data.frame(x=prior.means,y=prior.sd,z=v.index,cov=cor(sim$v)[v.index],
                           c=inlabru.coverage.sim[[1]],
                           l=inlabru.coverage.sim[[2]],
                           b=inlabru.coverage.sim[[3]])



#####################################################################################
# Generating  plots for v_{s,1}, with correlation 1 to w_s. From MCMC and inlabru
#####################################################################################

# Coverage plot when correlation is 1.
ggplot(filter(mcmc.data,z==1))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='MCMC, correlation 1'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_coverage_mcmc_corr1.pdf',device='pdf',width=4.5,height=3)

# Coverage plot when correlation is 1.
ggplot(filter(inlabru.data,z==1))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='inlabru, correlation 1'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_coverage_inla_corr1.pdf',device='pdf',width=4.5,height=3)

# Lambda error plot when correlation is 1. 
ggplot(filter(mcmc.data,z==1))+
  geom_tile(aes(x=x,y=y,fill=round(l,2)))+
  geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='MCMC, lambda error, correlation 1'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_lambda_error_mcmc_corr1.pdf',device='pdf',width=4.5,height=3)

ggplot(filter(inlabru.data,z==1))+
  geom_tile(aes(x=x,y=y,fill=round(l,2)))+
  geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='inlabru, lambda error, correlation 1'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_lambda_error_inla_corr1.pdf',device='pdf',width=4.5,height=3)

#####################################################################################
# Generating  plots for v_{s,3}, with correlation 0.6 to w_s. 
#####################################################################################

# Coverage plot when correlation is 0.6.
ggplot(filter(mcmc.data,z==3))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='MCMC, correlation 0.6'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_coverage_mcmc_corr06.pdf',device='pdf',width=4.5,height=3)

ggplot(filter(inlabru.data,z==3))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='inlabru, correlation 0.6'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_coverage_inla_corr06.pdf',device='pdf',width=4.5,height=3)

# Lambda error plot when correlation is 0.6. 
ggplot(filter(mcmc.data,z==3))+
  geom_tile(aes(x=x,y=y,fill=round(l,2)))+
  geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='MCMC, lambda error, correlation 0.6'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_lambda_error_mcmc_corr06.pdf',device='pdf',width=4.5,height=3)

ggplot(filter(inlabru.data,z==3))+
  geom_tile(aes(x=x,y=y,fill=round(l,2)))+
  geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='inlabru, lambda error, correlation 0.6'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_lambda_error_inla_corr06.pdf',device='pdf',width=4.5,height=3)


#####################################################################################
# Generating  plots for v_{s,4}, with correlation 0.4 to w_s. 
#####################################################################################

# Coverage plot when correlation is 0.4.
ggplot(filter(mcmc.data,z==4))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='MCMC, correlation 0.4'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_coverage_mcmc_corr04.pdf',device='pdf',width=4.5,height=3)

ggplot(filter(inlabru.data,z==4))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='inlabru, correlation 0.4'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_coverage_inla_corr04.pdf',device='pdf',width=4.5,height=3)

# Lambda error plot when correlation is 0.4. 
ggplot(filter(mcmc.data,z==4))+
  geom_tile(aes(x=x,y=y,fill=round(l,2)))+
  geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='MCMC, lambda error, correlation 0.4'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_lambda_error_mcmc_corr04.pdf',device='pdf',width=4.5,height=3)

ggplot(filter(inlabru.data,z==4))+
  geom_tile(aes(x=x,y=y,fill=round(l,2)))+
  geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='inlabru, lambda error, correlation 0.4'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour="none",fill="none")+
  scale_fill_distiller(palette = col.pal) +scale_colour_gradientn(colours = c(col.pal.green, col.pal.seagreen, col.pal.blue, col.pal.violet))
ggsave('sim_lambda_error_inla_corr04.pdf',device='pdf',width=4.5,height=3)
