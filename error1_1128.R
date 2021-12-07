set.seed(5)
options(digits=10)
# This is an error of taken wt freedom of 5 rather than 6
#
#
nu_right=6
nu=5

niter=50000

p=0.2
#para of p-prior
alpha_p <- 1
beta_p <- 1
#para of s-prior
lambda_s=c(1,2)
beta_s=c(1,2)
alpha_s=c(1,3)
#initial of mu and sigmasquare
mu=c(2,40)
sigmasquare=c(4,16)
#_chains are value of mcmc simulator
p_chain=rep(0.5,niter)
s_chain=matrix(0, nr= niter, nc = 2)

w_chain=rep(1,niter)
mu_chain=matrix(0, nr= niter, nc = 2)
mu_chain[1,]=mu
sigmasquare_chain=matrix(0, nr= niter, nc = 2)
sigmasquare_chain[1,]=sigmasquare


y_chain=rep(0,niter)
w0_chain=rep(0,niter)
#true p, mu1, mu2, sigmasquare1, sigmasuqare2
p0=0.3  
mu0=c(5,2)
sigmasquare0=c(9,1)
for(i in 1:niter){
  delta=rbinom(1,1,p0)+1
  w0_chain[i]=rchisq(1,nu)/6
  y_chain[i]<-rnorm(1,mean=mu0[delta],sqrt(sigmasquare0[delta]/w0_chain[i]))
}




#mcmc
for(iter in 2:niter){
  #sample s
  s=rbinom(1,1,p_chain[iter-1])
  s_chain[iter,1]=s
  s_chain[iter,2]=1-s
  #update conjugate p
  p_chain[iter]=rbeta(1,alpha_p+sum(s_chain[1:iter,1]),beta_p+iter-sum(s_chain[1:iter,1]))
  #sample w
  w_chain[iter]=rchisq(1,nu)/6
  # if(s==1){
  #   #sample y
  #   temp000<-rnorm(1,mean=mu_chain[iter,1],sqrt(sigmasquare_chain[iter,1]/w_chain[iter]))
  #   y_chain[iter]<-temp000
  # }else{
  #   y_chain[iter]<-rnorm(1,mean=mu_chain[iter,2],sqrt(sigmasquare_chain[iter,2]/w_chain[iter]))
  # }
  if(iter==1){
    n1y1=y_chain[1]*s_chain[1,1]
    s1_square=0
    n2y2=y_chain[2]*s_chain[1,2]
    s2_square=0
  }else{
    n1y1=sum(s_chain[1:iter,1]*y_chain[1:iter])
    s1_square=sum((y_chain[1:iter]-mean(y_chain[1:iter]))*s_chain[1:iter,1]*
                    (y_chain[1:iter]-mean(y_chain[1:iter]))*s_chain[1:iter,1])
    n2y2=sum(s_chain[1:iter,2]*y_chain[1:iter])
    s2_square=sum((y_chain[1:iter]-mean(y_chain[1:iter]))*s_chain[1:iter,2]*
                    (y_chain[1:iter]-mean(y_chain[1:iter]))*s_chain[1:iter,2])
  }
  if(s==1){
    temp1_0=lambda_s[1]*s1_square+beta_s[1]+s1_square-(lambda_s[1]*alpha_s[1]+n1y1)/(lambda_s[1]+sum(s_chain[1:iter,1]))
    temp1=rgamma(1,shape=((lambda_s[1]+sum(s_chain[1:iter,1])+3)/2),rate=(2/temp1_0))
    #update conjugate sigmasquare
    sigmasquare_chain[iter,1]=1/temp1
    sigmasquare_chain[iter,2]=sigmasquare_chain[iter-1,2]
    #update conjugate mu
    mu_chain[iter,1]=rnorm(1,(lambda_s[1]*alpha_s[1]+n1y1)/(lambda_s[1]+sum(s_chain[1:iter,1])),sqrt(sigmasquare_chain[iter,1]/(lambda_s[1]+sum(s_chain[1:iter,1]))))
    mu_chain[iter,2]=mu_chain[iter-1,2]
    
  }else{
    temp2_0=lambda_s[2]*s2_square+beta_s[2]+s2_square-(lambda_s[2]*alpha_s[2]+n2y2)/(lambda_s[2]+sum(s_chain[1:iter,2]))
    temp2=rgamma(1,shape=((lambda_s[2]+sum(s_chain[1:iter,2])+3)/2),rate=(2/temp2_0))
    #update conjugate sigmasquare
    sigmasquare_chain[iter,2]=1/temp2
    sigmasquare_chain[iter,1]=sigmasquare_chain[iter-1,1]
    #update conjugate mu
    mu_chain[iter,2]=rnorm(1,(lambda_s[2]*alpha_s[2]+n2y2)/(lambda_s[2]+sum(s_chain[1:iter,2])),sqrt(sigmasquare_chain[iter,2]/(lambda_s[2]+sum(s_chain[1:iter,2]))))
    mu_chain[iter,1]=mu_chain[iter-1,1]
    
  }
}


#_cons are value of marginal-conditional simulator
y_con=rep(20,niter)
s_con=matrix(0,nr=niter,nc=2)
p_con=rep(0.5,niter)
mu_con=matrix(NA,nr=niter,nc=2)
mu_con[1,]=mu
sigmasquare_con=matrix(NA,nr=niter,nc=2)
sigmasquare_con[1,]=sigmasquare



#marginal conditional
for(iter in 2:niter){
  p_con[iter]=rbeta(1,alpha_p,beta_p)
  s=rbinom(1,1,p_con[iter])
  s_con[iter,1]=s
  s_con[iter,2]=1-s
  if(s==1){
    mu_con[iter,1]=rnorm(1,lambda_s[1],sigmasquare[1]/lambda_s[1])
    sigmasquare_con[iter,1]=1/rgamma(1,0.5*(lambda_s[1]+3),2/beta_s[1])
    mu_con[iter,2]=mu_con[iter-1,2]
    sigmasquare_con[iter,2]=sigmasquare_con[iter-1,1]
    
    # temp_w=sigmasquare_con[iter,1]*rchisq(1,nu_right)/6
    # temp_e=rnorm(1,mu_con[iter,1],sqrt(sigmasquare_con[iter,1]/temp_w))
    # y_con[iter]=temp_e
  }else{
    mu_con[iter,2]=rnorm(1,lambda_s[2],sigmasquare[2]/lambda_s[2])
    sigmasquare_con[iter,2]=1/rgamma(1,0.5*(lambda_s[2]+3),2/beta_s[2])
    mu_con[iter,1]=mu_con[iter-1,1]
    sigmasquare_con[iter,1]=sigmasquare_con[iter-1,2]
    
    # temp_w=sigmasquare_con[iter,2]*rchisq(1,nu_right)/6
    # temp_e=rnorm(1,mu_con[iter,2],sqrt(sigmasquare_con[iter,2]/temp_w))
    # y_con[iter]=temp_e
  }
}


mcmc_mu1=mean(mu_chain[40000:niter,1])
mcmc_mu1_var=var(mu_chain[40000:niter,1])/10000
con_mu1=mean(mu_con[40000:niter,1])
con_mu1_var=var(mu_con[40000:niter,1])/10000
test_mu1=(mcmc_mu1-con_mu1)/sqrt(mcmc_mu1_var+con_mu1_var)

mcmc_mu2=mean(mu_chain[40000:niter,2])
mcmc_mu2_var=var(mu_chain[40000:niter,2])/10000
con_mu2=mean(mu_con[40000:niter,2])
con_mu2_var=var(mu_con[40000:niter,2])/10000
test_mu2=(mcmc_mu2-con_mu2)/sqrt(mcmc_mu2_var+con_mu2_var)

mcmc_sigmasquare1=mean(sigmasquare_chain[40000:niter,1])
mcmc_sigmasquare1_var=var(sigmasquare_chain[40000:niter,1])/10000
con_sigmasquare1=mean(sigmasquare_con[40000:niter,1])
con_sigmasquare1_var=var(sigmasquare_con[40000:niter,1])/10000
test_sigmasquare1=(mcmc_sigmasquare1-con_sigmasquare1)/sqrt(mcmc_sigmasquare1_var+con_sigmasquare1_var)


mcmc_sigmasquare2=mean(sigmasquare_chain[40000:niter,2])
mcmc_sigmasquare2_var=var(sigmasquare_chain[40000:niter,2])/10000
con_sigmasquare2=mean(sigmasquare_con[40000:niter,2])
con_sigmasquare2_var=var(sigmasquare_con[40000:niter,2])/10000
test_sigmasquare2=(mcmc_sigmasquare2-con_sigmasquare2)/sqrt(mcmc_sigmasquare2_var+con_sigmasquare2_var)

mcmc_p=mean(p_chain[40000:niter])
mcmc_p_var=var(p_chain[40000:niter])/10000
con_p=mean(p_con[40000:niter])
con_p_var=var(p_con[40000:niter])/1000
test_p=(mcmc_p-con_p)/sqrt(mcmc_p_var+con_p_var)

print(test_mu1)
print(test_mu2)
print(test_p)
print(test_sigmasquare1)
print(test_sigmasquare2)



