#Plotting simulated values, 
#and the relationship between the observed counts and
#the process covariates. 


set.seed(seed)

sim=list()
sim$N=100 # Number of observations.

# Set-up adjacency structure for ICAR model.
sim$A=adjacency.matrix(10)
sim$n_adj=rowSums(sim$A)
sim$D=diag(sim$n_adj)
sim$adj=as.carAdjacency(sim$A)$adj

tau=4 # Spatial effect precision parameter.

# Note: the function to simulate ICAR fields calls eigen
# which we found produced numerically different results
# on different machines, so to reproduce the plots in the
# paper we include our simulated phi.
sim$phi <- ricar_simple(tau*(sim$D-sim$A)) 
sim$phi <- simulated_phi # Overwrite. From Stoner et al. which crated this simulated_phi. 

# True covariates.
sim$x=runif(sim$N,-1,1)
sim$w=runif(sim$N,-1,1) # Uniform(-1,1) has variance 1/3

# Under-reporting covariate noise.
sim$gamma=rnorm(sim$N,0,1)
sim$gamma=sim$gamma-mean(sim$gamma)

# Imperfect under-reporting covariates.
sim$v=matrix(nrow=sim$N,ncol=6)
sim$v[,1]=sim$w
sim$v[,2]=0.8*sim$w + sqrt((1-0.8^2)/3)*sim$gamma 
sim$v[,3]=0.6*sim$w + sqrt((1-0.6^2)/3)*sim$gamma
sim$v[,4]=0.4*sim$w + sqrt((1-0.4^2)/3)*sim$gamma
sim$v[,5]=0.2*sim$w + sqrt((1-0.2^2)/3)*sim$gamma
sim$v[,6]=sqrt(1/3)*sim$gamma

cor(sim$v) # Check correlation is as desired.

a0=4
a1=1
b0=0
b1=2

sim$lambda=exp(a0+a1*sim$x+sim$phi) # Poisson means.

sim$y=rpois(sim$N,sim$lambda) # True counts.

sim$pi=expit(b0+b1*sim$v[,1]) # Reporting probabilities. # do not change v here, simulate pi using the true under-reporting covariate w. 

sim$z=rbinom(sim$N,sim$y,sim$pi) # Observed counts.


################################################################################
# Plotting the simulated data. Not needed for running the rest of the scripts. 
################################################################################


s1=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=x,y=y))+geom_point(colour=col.pal.green)+
  labs(
    y=expression('True Count ('*y[s]*')'),
    x=expression('Process Covariate ('*x[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s1
#ggsave('sim_xy.pdf',device='pdf',width=4.5,height=2.5)

s2=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=c,y=y))+geom_point(colour=col.pal.green)+
  labs(
    y=expression('True Count ('*y[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,1']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s2
#ggsave('sim_v1y.pdf',device='pdf',width=4.5,height=2.5)

s3=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,3],p=sim$pi),
          mapping=aes(x=c,y=y))+geom_point(colour=col.pal.green)+
  labs(
    y=expression('True Count ('*y[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,3']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s3
#ggsave('sim_v3y.pdf',device='pdf',width=4.5,height=2.5)

s4=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=x,y=z))+geom_point(colour=col.pal.green)+
  labs(
    y=expression('Observed Count ('*z[s]*')'),
    x=expression('Process Covariate ('*x[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s4
#ggsave('sim_xz.pdf',device='pdf',width=4.5,height=2.5)

s5=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=c,y=z))+geom_point(colour=col.pal.green)+
  labs(
    y=expression('Observed Count ('*z[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,1']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s5
#ggsave('sim_v1z.pdf',device='pdf',width=4.5,height=2.5)

s6=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,3],p=sim$pi),
          mapping=aes(x=c,y=z))+geom_point(colour=col.pal.green)+
  labs(
    y=expression('Observed Count ('*z[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,3']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s6
#ggsave('sim_v3z.pdf',device='pdf',width=4.5,height=2.5)

