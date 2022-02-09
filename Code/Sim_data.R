#Plotting simulated values, 
#and the relationship between the observed counts and
#the process covariates. 


set.seed(seed)

# Storing the value for simulated phi here. This makes the simulation-code possible to run without loading in the Data-file, where Stoner
# Has stored the simulated phi, as it changes values depending on computer used to calculate ut, because of the eigen()-function. 
# Using the same simulated phi as in the paper by Stoner et al. 
simulated_phi <- c(0.2844821276, -0.0363941265, -0.0707305993, 0.0427993385, -0.0017689272, 0.0008966238, 0.5634107744,
          0.2764364953, -0.8024495523,  -0.8561177709, 0.5271112148, 0.3809533570, -0.1819777796, 0.2477208992, -0.1182326767,
          -0.1516608025, 0.2040801653, -0.2160254748, -0.8508768879, -0.8204159336, 0.5689268672, 0.3076277727, -0.3250235462,
          -0.5634853700, -0.4016911559,-0.7568827615, -0.3636256845, -0.8918397721, -0.7434928724, -0.3197348102,  0.6355634230,
          0.3739016964, 0.2692871375, -0.0825380532, -0.5627878853, -0.1719729603,-0.7043588282, -0.9596770810, -1.4120089445,
          -0.3691143419,  0.5539610822, 0.7728624530, 0.9414071883,  0.4371045211, -0.3517693934, -0.1157574608, -0.4380227330, 
          -0.8139820155, -0.8439393041, -0.3245777684,  0.4809039258, 0.5230876870,0.9001388411, 0.3650667645, 0.1361678012,
          0.0078659201, -0.4247900319, -0.5862446291, -0.1379619510, -0.6258588073, 0.5512059103, 0.5166236169,  0.4088752449,
          0.2283988913, 0.5150460738,  0.2408916793, 0.0893253367, -0.2754093517, 0.3040021363, -0.6754143825, 0.7185599346,
          0.0475291408,  0.1564015858,  0.2246614043, 0.5712582000, 0.1825522324, 0.2022921257, 0.5200029705, 0.3655308218,
          -0.1764287029,-0.0054663672,-0.0511434216, 0.0202323774,  0.7017659297, 0.8503653324, 0.3515035738, 0.1504637600,
          -0.2613447203,0.0146129158,-0.0431936051,-0.0684673858, 0.2725965583, 0.1471866894, 0.0136802858, 0.1200327661,
          0.2336880557,0.3809539810,0.0139818627,-0.1939428604, 0.2326140199)


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
#sim$phi <- ricar_simple(tau*(sim$D-sim$A)) 
sim$phi <- simulated_phi #From Stoner et al. which crated this simulated_phi. 

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
ggsave('sim_xy.pdf',device='pdf',width=4.5,height=2.5)

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
ggsave('sim_v1y.pdf',device='pdf',width=4.5,height=2.5)

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
ggsave('sim_v3y.pdf',device='pdf',width=4.5,height=2.5)

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
ggsave('sim_xz.pdf',device='pdf',width=4.5,height=2.5)

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
ggsave('sim_v1z.pdf',device='pdf',width=4.5,height=2.5)

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
ggsave('sim_v3z.pdf',device='pdf',width=4.5,height=2.5)

