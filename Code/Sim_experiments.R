##################################################
# Simulation experiments, using inlabru and MCMC, 
# and plotting the results to compare. 
##################################################


set.seed(seed)
prior.means=b0+rep(0.6,3) # The mean of the prior for beta_0.
prior.sd=rep(0.6,3) # The standard deviation of the prior for beta_0.
v.index= c(1,3,4) # The strength of the under-reporting covariate.

##################################################
# Setting up and running MCMC
##################################################

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



print("Running MCMC for the three examples used in the specialisation project)")

sim.mcmc.list=list()
for(i in c(1:3)){
  sim.compiled.model$prior_mean=prior.means[i]
  sim.compiled.model$prior_sd=prior.sd[i]
  sim.compiled.model$v=sim$v[,v.index[i]]
  sim.samples<-runMCMC(sim.compiled.mcmc,inits=sim.inits,
                         nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                         summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)
  sim.mcmc.list[[i]]=as.mcmc(as.mcmc.list(sim.samples))
}

end_time_mcmc <- Sys.time()
print(c("Time taken to run three iterations of the MCMC algorithm:", difftime(end_time_mcmc,  start_time_mcmc, units = "secs")))

##################################################
# Setting up and running inlabru
##################################################
set.seed(seed)

run_inlabru <- function(means, sd, v.index){
  sim.inla.list <- list()
  for (i in 1:(length(v.index))){ 
    start_time_inla <- Sys.time()
    
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
    
    sim$post_y <- generate.y.samples(sim$post_lambda)
    
    
    sim.inla.list[[i]] <- list(fit, sim$post_lambda, sim$post_pi, sim$post_y, sim$post_pi_sorted)
    end_time_inla <- Sys.time()
    print(c("Time taken to run one iteration of the inla algorithm:", difftime(end_time_inla,  start_time_inla, units = "secs"), "seconds"))
    
  } 
  return(sim.inla.list)
}

print("Running INLA (inlabru) for the three examples used in the specialisation project)")

inlabru.results.sim <-run_inlabru(prior.means,prior.sd, v.index) 


 #Under here: generating plots for all three examples. 

##################################################################################
# Generating results and plots for v_{s,1}, correlation 1
##################################################################################

index=1 
lambda.1=sim.mcmc.list[[index]][,5:(4+sim$N)]
pi.1=sim.mcmc.list[[index]][,(5+sim$N):(4+2*sim$N)]
y.1=matrix(nrow=dim(lambda.1)[1],ncol=sim$N)
for(j in 1:sim$N){
  y.1[,j]=rpois(dim(lambda.1)[1],lambda.1[,j])
}
pi.sorted.1=expit(sim.mcmc.list[[index]][,3]+sim.mcmc.list[[index]][,4]%*%t(sort(sim$v[,1])))
phi.1=log(lambda.1)-sim.mcmc.list[[index]][,1]-sim.mcmc.list[[index]][,2]%*%t(sim$x)

inla.fit.1 <- inlabru.results.sim[[index]][[1]]
inla.lambda.1 <- inlabru.results.sim[[index]][[2]]
inla.pi.1 <- inlabru.results.sim[[index]][[3]]
inla.y.1 <-  inlabru.results.sim[[index]][[4]]
inla.pi.sorted.1 <- inlabru.results.sim[[index]][[5]]




a0.1=ggplot(data=data.frame(a0_mcmc=as.numeric(sim.mcmc.list[[index]][,1]),
                              x=seq(2.5,5,length=1000),prior=dnorm(seq(2.5,5,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  geom_density(mapping=aes(x=a0_mcmc, colour = "MCMC"), adjust = 3, fill=NA,alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.1$marginals.fixed), aes(x = Intercept.x, y = Intercept.y, colour = "inlabru"), alpha = 0.4, fill = NA, size = 0.6) +
  geom_vline(aes(xintercept = a0,colour="True Value"), show.legend = F) +
  labs(
    y='Density',
    x=expression(alpha[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+scale_x_continuous(limits = c(2.5,5))+scale_y_continuous(limits = c(0,5.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))
a0.1
ggsave('mcmc_inla_a0_corr1.pdf',device='pdf',width=4.5,height=3)

a1.1=ggplot(data=data.frame(a1_mcmc=as.numeric(sim.mcmc.list[[index]][,2]),
                              x=seq(0.5,2,length=1000),prior=dnorm(seq(0.25,1.75,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=a1),adjust=3,fill=col.pal.red,alpha=0.4, colour = "black")+
  geom_density(mapping=aes(x=a1_mcmc, colour = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.1$marginals.fixed), aes(x = covariate.x, y = covariate.y, colour = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(alpha[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = a1,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(0.5,2))+scale_y_continuous(limits = c(0,6.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(shape = NA,colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
a1.1
ggsave('mcmc_inla_a1_corr1.pdf',device='pdf',width=4.5,height=3)

b0.1=ggplot(data=data.frame(b0_mcmc=as.numeric(sim.mcmc.list[[index]][,3]),
                              x=seq(-1.7,3,length=1000),
                              prior=dnorm(seq(-1.7,3,length=1000),0.6,0.6)))+
  geom_line(aes(x=x,y=prior, colour = "Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=b0),adjust=3,fill=col.pal.blue,alpha=0.4)+
  geom_density(mapping=aes(x=b0_mcmc, color = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.1$marginals.fixed), aes(x = beta0.x, y = beta0.y, color = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(beta[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = b0,colour='True Value'), show.legend = F
  )+scale_x_continuous(limits = c(-2.5,3))+scale_y_continuous(limits = c(0,1.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
b0.1
ggsave('mcmc_inla_b0_corr1.pdf',device='pdf',width=4.5,height=3)

b1.1=ggplot(data=data.frame(b1_mcmc=as.numeric(sim.mcmc.list[[index]][,4]),
                              x=seq(-0.5,4,length=1000),
                              prior=dnorm(seq(-0.5,4,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=b1),adjust=3,fill=col.pal.blue,alpha=0.4)+
  geom_density(mapping=aes(x=b1_mcmc, color = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.1$marginals.fixed), aes(x = beta1.x, y = beta1.y, color = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(beta[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = b1,colour='True Value'), show.legend = F
  )+scale_x_continuous(limits = c(-0.5,4))+scale_y_continuous(limits = c(0,1.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
b1.1
ggsave('mcmc_inla_b1_corr1.pdf',device='pdf',width=4.5,height=3)


mcmc.pi.plot.1=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,1],index.return=TRUE)$ix],x=sort(sim$v[,1]),
                              l=apply(pi.sorted.1,2,quantile,probs=0.025),
                              m=apply(pi.sorted.1,2,quantile,probs=0.5),
                              u=apply(pi.sorted.1,2,quantile,probs=0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=col.pal.green)+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.4,fill=col.pal.green)+
  geom_line(aes(x=x,y=m),colour=col.pal.green,size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,1']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_y_continuous(limits = c(0,1))
mcmc.pi.plot.1
ggsave('mcmc_pi_corr1.pdf',device='pdf',width=4.5,height=3)

inla.pi.plot.1=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,1],index.return=TRUE)$ix],x=sort(sim$v[,1]),
                                 l=rowQuantiles(inla.pi.sorted.1, rows = NULL, cols = NULL, probs = 0.025),
                                 m=rowQuantiles(inla.pi.sorted.1, rows = NULL, cols = NULL, probs = 0.5),
                                 u=rowQuantiles(inla.pi.sorted.1, rows = NULL, cols = NULL, probs = 0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=col.pal.green)+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.4, fill = col.pal.green)+
  geom_line(aes(x=x,y=m),colour=col.pal.green,size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,1']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_y_continuous(limits = c(0,1))
inla.pi.plot.1
ggsave('inla_pi_corr1.pdf',device='pdf',width=4.5,height=3)


mcmc.phi.plot.1=ggplot(data=data.frame(x=sim$phi,y=apply(phi.1,2,mean)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=col.pal.green)+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.5,1))+scale_y_continuous(limits = c(-1.5,1.5))
mcmc.phi.plot.1
ggsave('mcmc_phi_corr1.pdf',device='pdf',width=4.5,height=3)

inla.phi.plot.1 =ggplot(data=data.frame(x=sim$phi,y=inla.fit.1$summary.random$spatial$mean))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=col.pal.green)+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.5,1))+scale_y_continuous(limits = c(-1.5,1.5))
inla.phi.plot.1
ggsave('inla_phi_corr1.pdf',device='pdf',width=4.5,height=3)


mcmc.y.plot.1=ggplot(data.frame(x=sim$y,
                           l=apply(y.1,2,quantile,0.025),
                           u=apply(y.1,2,quantile,0.975),
                           m=apply(y.1,2,mean)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=col.pal.blue)+
  geom_point(mapping=aes(x=x,y=u),colour=col.pal.green)+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+ scale_y_continuous(limits = c(0,400)) +scale_x_continuous(limits = c(0,300))
mcmc.y.plot.1
ggsave('mcmc_pred_y_corr1.pdf',device='pdf',width=4.5,height=3)


inla.y.plot.1 = ggplot(data.frame(x=sim$y,
                                     l=rowQuantiles(inla.y.1, rows = NULL, cols = NULL, probs = 0.025),
                                     m=rowQuantiles(inla.y.1, rows = NULL, cols = NULL, probs = 0.5),
                                     u=rowQuantiles(inla.y.1, rows = NULL, cols = NULL, probs = 0.975)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=col.pal.blue)+
  geom_point(mapping=aes(x=x,y=u),colour=col.pal.green)+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+ scale_y_continuous(limits = c(0,400))  +scale_x_continuous(limits = c(0,300))
inla.y.plot.1
ggsave('inla_pred_y_corr1.pdf',device='pdf',width=4.5,height=3)



##################################################################################
# Generating results and plots for v_{s,3}, correlation 0.6
##################################################################################

index=2 
lambda.06=sim.mcmc.list[[index]][,5:(4+sim$N)]
pi.06=sim.mcmc.list[[index]][,(5+sim$N):(4+2*sim$N)]
y.06=matrix(nrow=dim(lambda.06)[1],ncol=sim$N)
for(j in 1:sim$N){
  y.06[,j]=rpois(dim(lambda.06)[1],lambda.06[,j])
}
pi.sorted.06=expit(sim.mcmc.list[[index]][,3]+sim.mcmc.list[[index]][,4]%*%t(sort(sim$v[,3])))
phi.06=log(lambda.06)-sim.mcmc.list[[index]][,1]-sim.mcmc.list[[index]][,2]%*%t(sim$x)

inla.fit.06 <- inlabru.results.sim[[index]][[1]]
inla.lambda.06 <- inlabru.results.sim[[index]][[2]]
inla.pi.06 <- inlabru.results.sim[[index]][[3]]
inla.y.06 <-  inlabru.results.sim[[index]][[4]]
inla.pi.sorted.06 <- inlabru.results.sim[[index]][[5]]


a0.06=ggplot(data=data.frame(a0_mcmc=as.numeric(sim.mcmc.list[[index]][,1]),
                              x=seq(2.5,5,length=1000),prior=dnorm(seq(2.5,5,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  geom_density(mapping=aes(x=a0_mcmc, colour = "MCMC"), adjust = 3, fill=NA,alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.06$marginals.fixed), aes(x = Intercept.x, y = Intercept.y, colour = "inlabru"), alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(alpha[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = a0,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(2.5,5))+scale_y_continuous(limits = c(0,3)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
a0.06
ggsave('mcmc_inla_a0_corr06.pdf',device='pdf',width=4.5,height=3)

a1.06=ggplot(data=data.frame(a1_mcmc=as.numeric(sim.mcmc.list[[index]][,2]),
                              x=seq(0.5,2,length=1000),prior=dnorm(seq(0.25,1.75,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=a1),adjust=3,fill=col.pal.red,alpha=0.4, colour = "black")+
  geom_density(mapping=aes(x=a1_mcmc, colour = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.06$marginals.fixed), aes(x = covariate.x, y = covariate.y, colour = "inlabru"), 
                 alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(alpha[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = a1,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(0.5,2))+scale_y_continuous(limits = c(0,6.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior","True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
a1.06
ggsave('mcmc_inla_a1_corr06.pdf',device='pdf',width=4.5,height=3)

b0.06=ggplot(data=data.frame(b0_mcmc=as.numeric(sim.mcmc.list[[index]][,3]),
                              x=seq(-1.7,3,length=1000),
                              prior=dnorm(seq(-1.7,3,length=1000),0.6,0.6)))+
  geom_line(aes(x=x,y=prior, colour = "Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=b0),adjust=3,fill=col.pal.blue,alpha=0.4)+
  geom_density(mapping=aes(x=b0_mcmc, color = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.06$marginals.fixed), aes(x = beta0.x, y = beta0.y, color = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(beta[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = b0,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(-2.5,3))+scale_y_continuous(limits = c(0,1)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2,"True Value" = 1))
b0.06
ggsave('mcmc_inla_b0_corr06.pdf',device='pdf',width=4.5,height=3)

b1.06=ggplot(data=data.frame(b1_mcmc=as.numeric(sim.mcmc.list[[index]][,4]),
                              x=seq(-0.5,4,length=1000),
                              prior=dnorm(seq(-0.5,4,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=b1),adjust=3,fill=col.pal.blue,alpha=0.4)+
  geom_density(mapping=aes(x=b1_mcmc, color = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.06$marginals.fixed), aes(x = beta1.x, y = beta1.y, color = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(beta[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = b1,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(-0.5,4))+scale_y_continuous(limits = c(0,1.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior","True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
b1.06
ggsave('mcmc_inla_b1_corr06.pdf',device='pdf',width=4.5,height=3)


mcmc.pi.plot.06=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,3],index.return=TRUE)$ix],x=sort(sim$v[,3]),
                              l=apply(pi.sorted.06,2,quantile,probs=0.025),
                              m=apply(pi.sorted.06,2,quantile,probs=0.5),
                              u=apply(pi.sorted.06,2,quantile,probs=0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=col.pal.green)+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.4,fill=col.pal.green)+
  geom_line(aes(x=x,y=m),colour=col.pal.green,size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,3']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.35,1.35))+scale_y_continuous(limits = c(0,1))
mcmc.pi.plot.06
ggsave('mcmc_pi_corr06.pdf',device='pdf',width=4.5,height=3)

inla.pi.plot.06=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,3],index.return=TRUE)$ix],x=sort(sim$v[,3]),
                                      l=rowQuantiles(inla.pi.sorted.06, rows = NULL, cols = NULL, probs = 0.025),
                                      m=rowQuantiles(inla.pi.sorted.06, rows = NULL, cols = NULL, probs = 0.5),
                                      u=rowQuantiles(inla.pi.sorted.06, rows = NULL, cols = NULL, probs = 0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=col.pal.green)+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.4, fill = col.pal.green)+
  geom_line(aes(x=x,y=m),colour=col.pal.green,size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,3']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.35,1.35))+scale_y_continuous(limits = c(0,1))
inla.pi.plot.06
ggsave('inla_pi_corr06.pdf',device='pdf',width=4.5,height=3)

mcmc.phi.plot.06=ggplot(data=data.frame(x=sim$phi,y=apply(phi.06,2,mean)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=col.pal.green)+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.5,1))+scale_y_continuous(limits = c(-1.5,1.5))
mcmc.phi.plot.06
ggsave('mcmc_phi_corr06.pdf',device='pdf',width=4.5,height=3)

inla.phi.plot.06 =ggplot(data=data.frame(x=sim$phi,y=inla.fit.06$summary.random$spatial$mean))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=col.pal.green)+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.5,1))+scale_y_continuous(limits = c(-1.5,1.5))
inla.phi.plot.06
ggsave('inla_phi_corr06.pdf',device='pdf',width=4.5,height=3)

mcmc.y.plot.06=ggplot(data.frame(x=sim$y,l=apply(y.06,2,quantile,0.025),u=apply(y.06,2,quantile,0.975),
                         m=apply(y.06,2,mean)))+geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=col.pal.blue)+
  geom_point(mapping=aes(x=x,y=u),colour=col.pal.green)+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+
  scale_y_continuous(limits = c(0,400)) + scale_x_continuous(limits = c(0,300))
mcmc.y.plot.06
ggsave('mcmc_pred_y_corr06.pdf',device='pdf',width=4.5,height=3)

inla.y.plot.06 = ggplot(data.frame(x=sim$y,
                                  l=rowQuantiles(inla.y.06, rows = NULL, cols = NULL, probs = 0.025),
                                  m=rowQuantiles(inla.y.06, rows = NULL, cols = NULL, probs = 0.5),
                                  u=rowQuantiles(inla.y.06, rows = NULL, cols = NULL, probs = 0.975)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=col.pal.blue)+
  geom_point(mapping=aes(x=x,y=u),colour=col.pal.green)+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+ scale_y_continuous(limits = c(0,400))  +scale_x_continuous(limits = c(0,300))
inla.y.plot.06
ggsave('inla_pred_y_corr06.pdf',device='pdf',width=4.5,height=3)


##################################################################################
# Generating results and plots for v_{s,4}, correlation 0.4
##################################################################################

index=3
lambda.04=sim.mcmc.list[[index]][,5:(4+sim$N)]
pi.04=sim.mcmc.list[[index]][,(5+sim$N):(4+2*sim$N)]
y.04=matrix(nrow=dim(lambda.04)[1],ncol=sim$N)
for(j in 1:sim$N){
  y.04[,j]=rpois(dim(lambda.04)[1],lambda.04[,j])
}
pi.sorted.04=expit(sim.mcmc.list[[index]][,3]+sim.mcmc.list[[index]][,4]%*%t(sort(sim$v[,1])))
phi.04=log(lambda.04)-sim.mcmc.list[[index]][,1]-sim.mcmc.list[[index]][,2]%*%t(sim$x)

inla.fit.04 <- inlabru.results.sim[[index]][[1]]
inla.lambda.04 <- inlabru.results.sim[[index]][[2]]
inla.pi.04 <- inlabru.results.sim[[index]][[3]]
inla.y.04 <-  inlabru.results.sim[[index]][[4]]
inla.pi.sorted.04 <- inlabru.results.sim[[index]][[5]]


a0.04=ggplot(data=data.frame(a0_mcmc=as.numeric(sim.mcmc.list[[index]][,1]),
                              x=seq(2.5,5,length=1000),prior=dnorm(seq(2.5,5,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  geom_density(mapping=aes(x=a0_mcmc, colour = "MCMC"), adjust = 3, fill=NA,alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.04$marginals.fixed), aes(x = Intercept.x, y = Intercept.y, colour = "inlabru"), alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(alpha[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = a0,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(2.5,5))+scale_y_continuous(limits = c(0,3)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
a0.04
ggsave('mcmc_inla_a0_corr04.pdf',device='pdf',width=4.5,height=3)

a1.04=ggplot(data=data.frame(a1_mcmc=as.numeric(sim.mcmc.list[[index]][,2]),
                              x=seq(0,2,length=1000),prior=dnorm(seq(0.25,1.75,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=a1),adjust=3,fill=col.pal.red,alpha=0.4, colour = "black")+
  geom_density(mapping=aes(x=a1_mcmc, colour = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.04$marginals.fixed), aes(x = covariate.x, y = covariate.y, colour = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(alpha[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = a1,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(0,2))+scale_y_continuous(limits = c(0,4)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
a1.04
ggsave('mcmc_inla_a1_corr04.pdf',device='pdf',width=4.5,height=3)

b0.04=ggplot(data=data.frame(b0_mcmc=as.numeric(sim.mcmc.list[[index]][,3]),
                              x=seq(-1.7,3,length=1000),
                              prior=dnorm(seq(-1.7,3,length=1000),0.6,0.6)))+
  geom_line(aes(x=x,y=prior, colour = "Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=b0),adjust=3,fill=col.pal.blue,alpha=0.4)+
  geom_density(mapping=aes(x=b0_mcmc, color = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.04$marginals.fixed), aes(x = beta0.x, y = beta0.y, color = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(beta[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = b0,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(-2.5,3))+scale_y_continuous(limits = c(0,1.0)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
b0.04
ggsave('mcmc_inla_b0_corr04.pdf',device='pdf',width=4.5,height=3)

b1.04=ggplot(data=data.frame(b1_mcmc=as.numeric(sim.mcmc.list[[index]][,4]),
                              x=seq(-0.5,4,length=1000),
                              prior=dnorm(seq(-0.5,4,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior,colour="Prior", linetype = "Prior"))+
  #stat_density(mapping=aes(x=b1),adjust=3,fill=col.pal.blue,alpha=0.4)+
  geom_density(mapping=aes(x=b1_mcmc, color = "MCMC"), adjust = 3, fill=NA, alpha=0.4, size = 0.6)+
  geom_area(data = data.frame(inla.fit.04$marginals.fixed), aes(x = beta1.x, y = beta1.y, color = "inlabru"), 
            alpha = 0.4, fill = NA, size = 0.6) +
  labs(
    y='Density',
    x=expression(beta[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.position = "bottom"
  )+geom_vline(aes(xintercept = b1,colour="True Value"), show.legend = F
  )+scale_x_continuous(limits = c(-0.5,4))+scale_y_continuous(limits = c(0,1.5)
  )+scale_color_manual(" " ,limits=c("MCMC", "inlabru", "Prior", "True Value"), values = c(col.pal.red,col.pal.blue, "black", col.pal.green)
  )+guides(colour = guide_legend(override.aes = list( values =c(col.pal.red,col.pal.blue, "black", col.pal.green))))+
  scale_linetype_manual(" ",values=c("MCMC"=1,"inlabru"=1, "Prior" = 2, "True Value" = 1))
b1.04
ggsave('mcmc_inla_b1_corr04.pdf',device='pdf',width=4.5,height=3)



mcmc.pi.plot.04=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,4],index.return=TRUE)$ix],x=sort(sim$v[,4]),
                              l=apply(pi.sorted.04,2,quantile,probs=0.025),
                              m=apply(pi.sorted.04,2,quantile,probs=0.5),
                              u=apply(pi.sorted.04,2,quantile,probs=0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=col.pal.green)+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.4,fill=col.pal.green)+
  geom_line(aes(x=x,y=m),colour=col.pal.green,size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,4']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_y_continuous(limits = c(0,1))+scale_x_continuous(limits = c(-1.35,1.35))
mcmc.pi.plot.04
ggsave('mcmc_pi_corr04.pdf',device='pdf',width=4.5,height=3)

inla.pi.plot.04=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,4],index.return=TRUE)$ix],x=sort(sim$v[,4]),
                                       l=rowQuantiles(inla.pi.sorted.04, rows = NULL, cols = NULL, probs = 0.025),
                                       m=rowQuantiles(inla.pi.sorted.04, rows = NULL, cols = NULL, probs = 0.5),
                                       u=rowQuantiles(inla.pi.sorted.04, rows = NULL, cols = NULL, probs = 0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=col.pal.green)+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.4, fill = col.pal.green)+
  geom_line(aes(x=x,y=m),colour=col.pal.green,size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,4']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_y_continuous(limits = c(0,1))
inla.pi.plot.04
ggsave('inla_pi_corr04.pdf',device='pdf',width=4.5,height=3)


mcmc.phi.plot.04=ggplot(data=data.frame(x=sim$phi,y=apply(phi.04,2,mean)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=col.pal.green)+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.5,1))+scale_y_continuous(limits = c(-1.5,1.5))
mcmc.phi.plot.04
ggsave('mcmc_phi_corr04.pdf',device='pdf',width=4.5,height=3)

inla.phi.plot.04 =ggplot(data=data.frame(x=sim$phi,y=inla.fit.04$summary.random$spatial$mean))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=col.pal.green)+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+scale_x_continuous(limits = c(-1.5,1))+scale_y_continuous(limits = c(-1.5,1.5))
inla.phi.plot.04
ggsave('inla_phi_corr04.pdf',device='pdf',width=4.5,height=3)


mcmc.y.plot.04=ggplot(data.frame(x=sim$y,l=apply(y.04,2,quantile,0.025),u=apply(y.04,2,quantile,0.975),
                         m=apply(y.04,2,mean)))+geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=col.pal.blue)+
  geom_point(mapping=aes(x=x,y=u),colour=col.pal.green)+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+ scale_y_continuous(limits = c(0,400)) +scale_x_continuous(limits = c(0,300))
mcmc.y.plot.04
ggsave('mcmc_pred_y_corr04.pdf',device='pdf',width=4.5,height=3)


inla.y.plot.04 = ggplot(data.frame(x=sim$y,
                                   l=rowQuantiles(inla.y.04, rows = NULL, cols = NULL, probs = 0.025),
                                   m=rowQuantiles(inla.y.04, rows = NULL, cols = NULL, probs = 0.5),
                                   u=rowQuantiles(inla.y.04, rows = NULL, cols = NULL, probs = 0.975)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=col.pal.blue)+
  geom_point(mapping=aes(x=x,y=u),colour=col.pal.green)+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+ scale_y_continuous(limits = c(0,400))  +scale_x_continuous(limits = c(0,300))
inla.y.plot.04
ggsave('inla_pred_y_corr04.pdf',device='pdf',width=4.5,height=3)

