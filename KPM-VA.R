###EXAMPLE####
rm(list=ls())
library(fpc)
library(GenSA)
library(MASS)
library(glmnet)
library(fda)
library(splines)
library(mclust)
library(igraph)
library(igraph)
library(flexmix)
setwd("e:/heye/subgroup for panel data/code_two_dim/s1")

genedata<-function(sam)
{
  set.seed(sam)
  t0<-seq(0,1,length.out=m)
  beta_1 <- 3*(1+exp(-(t0-0.5)/0.1))^(-1)
  beta_2 <- 3*(2*t0-(6*t0^2)+(4*t0^3)+(1+exp(-((t0-0.7)/0.05)))^(-1))
  beta_3 <- 3*(4*t0-(8*t0^2)+4*t0^3+(1+exp(-(t0-0.6)/0.05))^(-1))
  beta_4 <- 3*(2*t0-4*t0^2+2*t0^3+(1+exp(-(t0-0.6)/0.1))^(-1))
  beta_5 <- 3*(t0-(3*t0^2)+(2*t0^3)+(1+exp(-(t0-0.7)/0.04))^(-1))
  beta_6 <- 3*(0.5*t0-(0.5*t0^2)+(1+exp(-(t0-0.4)/0.07))^(-1))
  noise<-matrix(rnorm(N,mu1,sigma1),n,m)
  X_data<-matrix(rnorm(2*N,mu2,sigma2),2*n,m)
  X1<-X_data[c(1:n),c(1:m)]
  X2<-X_data[c((n+1):(2*n)),c(1:m)]
  #X2<-matrix(rnorm(N,mu3,sigma3),n,m)
  X <- rbind(X1,X2)
  betax1 <-t(matrix(c(rep(beta_1,15),rep(beta_2,10),rep(beta_3,20)),m,n))####################n????????????
  betax2 <-t(matrix(c(rep(beta_4,15),rep(beta_5,10),rep(beta_6,20)),m,n))####################n???????????
  BETA_31 <- cbind(beta_1,beta_2,beta_3)
  BETA_32 <- cbind(beta_4,beta_5,beta_6)
  beta_3all <- t(rbind(BETA_31,BETA_32))
  betareal <- rbind(betax1,betax2)
  
  ui<-rnorm(n,0,1)
  u<-outer(ui,rep(1,m))
  Y<-u+X1*betax1+X2*betax2+noise
  return(list(Y=Y,X1=X1,X2=X2,betax1=betax1,betax2=betax2,t0=t0))
}


####the optim function######
f_optim <- function(beta_vec){ 
  khm <- 0
  beta_est<-matrix(beta_vec,n,2*d)
  beta_est1<-beta_est[,1:d]
  beta_est2<-beta_est[,-c(1:d)]
 rss=0
  for (j in 1:n)
  {
    rss<-rss+sum((ymean[j,]-x1mean[j,,]%*%beta_est1[j,]-x2mean[j,,]%*%beta_est2[j,])^2)/m
  }
    
  for (i in 1:K){
    khm <- khm + (apply((beta_est-rep(1,n)%o%gamma0[i,])^2,1,sum))^s######################################g
  }
  
  pen<-lam1*sum((khm/K)^(1/s))
  fun_value<-rss+pen 
  fun_value
} 

####the first gradient####
f_d<-function(beta_vec){
  beta_est<-matrix(beta_vec,n,2*d)
  beta_est1<-beta_est[,1:d]
  beta_est2<-beta_est[,-c(1:d)]
  
  dis<-matrix(0,nrow=1,ncol=K)
  f_d<-f_d1<-f_d2<-matrix(0,n,(2*d))
  WW <- matrix(0,nrow=(2*d),ncol=K)
  for(j in 1:n)
  { 
    xjcmean<-cbind(x1mean[j,,],x2mean[j,,])
    yjmean<-ymean[j,]
    f_d1[j,]<--2/m*(t(xjcmean)%*%(yjmean-xjcmean%*%beta_est[j,]))
    for(i in 1:K){
      dis[1,i]<-sum((beta_est[j,]-gam0[i,])^2)
    }
    W<-(1/K*(sum(dis^s)))^(1/s-1)
    for (l in 1:K)
    {
      WW[,l] <- W*dis[1,l]^(s-1)*2*(beta_est[j,]-gam0[l,])/K
    }
    f_d2[j,]<-apply(WW,1,sum)
    f_d[j,]<-f_d1[j,]+lam1*f_d2[j,]
  }
  return(as.vector(f_d))
    
  }
  


#update center####
gamhat<-function(betaest,gammaest){ 
  dis<-matrix(0,nrow=1,ncol=K)
  wei<-matrix(0,n,K)
  for (l in 1:n){
    for(i in 1:K){
      dis[1,i]<-sum((betaest[l,]-gammaest[i,])^2)
    }
    W<-((1/K)^(1/s))*((sum(dis^s)))^(1/s -1)
    wei[l,]<-t(W*(dis^(s-1)))
  }
  wei[is.nan (wei)] <- 0
  gamma_up<-t(wei)%*%betaest/(apply(wei,2,sum)%o%rep(1,2*d)) 
  return(gamma_up)
}

function_rss <- function(a2_1){
  k_v<-sort(unique(a2_1))
  rss<-0
  for (i in k_v){
    len_k<-sum(a2_1==i)
    if (len_k==1){
      yk<-ymean[which(a2_1==i),]
      xkt<-cbind(x1mean[which(a2_1==i),,],x2mean[which(a2_1==i),,])
      fit<-lm(yk~xkt-1)$residual
      rss<-rss+sum(fit^2)
    }else{
      yk<-as.vector(t(ymean[which(a2_1==i),]))
      xkt1<-x1mean[which(a2_1==i),,];xkt2<-x2mean[which(a2_1==i),,]
      xk1new<-do.call(rbind, lapply(1:dim(xkt1)[1], function(i) xkt1[i,,]))
      xk2new<-do.call(rbind, lapply(1:dim(xkt2)[1], function(i) xkt2[i,,]))
      res<-lm(yk~xk1new+xk2new-1)$residual
      rss<-rss+sum(res^2)
      }
  }
  return(rss=rss)
}



function_betahat <- function(a2_1){
  k_v<-sort(unique(a2_1))
  rss<-0
  for (i in k_v){
    len_k<-sum(a2_1==i)
    if (len_k==1){
      yk<-ymean[which(a2_1==i),]
      xkt<-cbind(x1mean[which(a2_1==i),,],x2mean[which(a2_1==i),,])
      beta_est<-lm(yk~xkt-1)$coef
      beta1k<-beta10[which(a2_1==i),]
      beta2k<-beta20[which(a2_1==i),]
      rss<-rss+sum((beta1k-B%*%beta_est[1:d])^2)+sum((beta2k-B%*%beta_est[-c(1:d)])^2)
    }else{
      yk<-as.vector(t(ymean[which(a2_1==i),]))
      beta1k<-beta10[which(a2_1==i),]
      beta2k<-beta20[which(a2_1==i),]
      xkt1<-x1mean[which(a2_1==i),,];xkt2<-x2mean[which(a2_1==i),,]
      xk1new<-do.call(rbind, lapply(1:dim(xkt1)[1], function(i) xkt1[i,,]))
      xk2new<-do.call(rbind, lapply(1:dim(xkt2)[1], function(i) xkt2[i,,]))
      beta_est<-lm(yk~xk1new+xk2new-1)$coef
      
      rss<-rss+sum((beta1k-t(matrix(rep(B%*%beta_est[1:d],len_k),nrow=m)))^2)
      +sum((beta2k-t(matrix(rep(B%*%beta_est[-c(1:d)],len_k),nrow=m)))^2)
      }
  }
  return(rss=rss)
}

n=45;m=50;N=n*m
mu1=mu2=0
sigma1=1####the variance of x
sigma2=1#####the variance of noise 
s=-0.5
epis<-10^(-8);emax=200;n0=100
a<-c(rep(1,15),rep(2,10),rep(3,20))

#d<-floor(N^(1/6))+4
k_sec<- c(1:5)####choose K
d_sec <- c(5:8)####choose d
lam_sec <- c(0.001,0.005,0.01,0.03,0.05,0.1,1)######choose lambda
d_len<-length(d_sec);k_len <- length(k_sec);lam_len <- length(lam_sec)

re<-matrix(0,n0,8)
colnames(re)=c("bic","ri","k","d","rss","rss0","ri_ini","rssmean")

for (sam in 8:n0)
{
  data<-genedata(sam)
  y<-data$Y
  x1<-data$X1
  x2<-data$X2
  t0<-data$t0
  beta10<-data$betax1
  beta20<-data$betax2
  BIC_RI_K<-matrix(0,d_len,6)
  colnames(BIC_RI_K)=c("bic","ri","k","rss","ri_ini","rssmean")
for (l in 1:d_len)
{
  d<-d_sec[l]
  Bsp = create.bspline.basis(rangeval=c(0,1),norder = 4,nbasis = d)
  B=eval.basis(t0, Bsp, returnMatrix=T)
  beta0<-matrix(0,n,2*d)
  
  ymean<-matrix(0,n,m)
  x1mean<-x2mean<-array(0,dim=c(n,m,d))
  for(j in 1:n)
  {
    yj<-y[j,]
    xj1<-diag(x1[j,])%*%B
    xj2<-diag(x2[j,])%*%B
    ymean[j,]<-yjt<-yj-mean(yj)
    x1mean[j,,]<-xjt1<-xj1-outer(rep(1,m),apply(xj1,2,mean))
    x2mean[j,,]<-xjt2<-xj2-outer(rep(1,m),apply(xj2,2,mean))
    beta0[j,]<-lm(yjt~xjt1+xjt2-1)$coef
  }
  RI<-BIC<-RSS0<-matrix(0,k_len,lam_len)
RI_ini<-matrix(0,k_len,2)
for(k_sam in 1:k_len)
{
  K<-k_sec[k_sam]
if(K==1){
       RI[1,1:lam_len]<-RI_ini[K,1]<-compare(rep(1,n),a,method="rand")
      yk<-as.vector(t(ymean))
      xk1new<-do.call(rbind, lapply(1:dim(x1mean)[1], function(i) x1mean[i,,]))
      xk2new<-do.call(rbind, lapply(1:dim(x2mean)[1], function(i) x2mean[i,,]))
      fit<-lm(yk~xk1new+xk2new-1)
      RSS<-sum(fit$residual^2)
      beta_est<-fit$coef
      BIC[K,1:lam_len]<-log(RSS/N)+(2*d)*((log(N))/N)*(K)
      RSS0[K,1:lam_len]<-RI_ini[K,2]<-sum((beta10-t(matrix(rep(B%*%beta_est[1:d],n),nrow=m)))^2)
      +sum((beta20-t(matrix(rep(B%*%beta_est[-c(1:d)],n),nrow=m)))^2)
}else{
  clu<-kmeans(beta0,K,nstart = 50)
  gamma0<-clu$center;a1<-clu$cluster
  RI_ini[K,1]<-compare(a1,a,method="rand")
  RI_ini[K,2]<-function_betahat(a1) 
  for (lam in 1:lam_len)
  {
    lam1<-lam_sec[lam]
  gam0 <- gamma0
  be0<-beta0
  e1<-1;q1=1
  while(e1>epis&q1<=emax){ 
    beta_vec0<-as.vector(be0)
    outcome <- optim(par=beta_vec0,fn=f_optim,gr=f_d,method = "BFGS",control = list(maxit=100))
    beta_vec1<-outcome$par####beta
    beta1<-matrix(beta_vec1,n,(2*d))
    q1<-q1+1
    gamma1<-gamhat(beta1,gam0)
    e1<-max((beta1-be0)^2)+max((gamma1-gam0)^2)
    gam0<-gamma1
    be0<-beta1
    a2 <- tryCatch(
      {a2<-kmeans(beta1,K,nstart = 50,centers=gam0)$cluster
      q1=q1+1
      a2},
      warning=function(w){
        print("warning")
        rep(1,n)},
      error=function(e){
        print("All sample in one cluster")
        rep(1,n)})
    if(length(unique(a2))==1){q1=emax+1}
    #print(e1)
  }
   RI[K,lam]<-compare(a2,a,method="rand")
    RSS<-function_rss(a2)
    BIC[K,lam]<-log(RSS/N)+(2*d*K)*((log(N))/N)
    RSS0[K,lam]<-function_betahat(a2)
  } ###else 
   }####lambda  
cat("this is the", K, "-th of the", l, "-th\n")
} ###d 
bic_pos<-which(BIC == min(BIC), arr.ind = TRUE)
BIC_RI_K[l,1]<-min(BIC)
BIC_RI_K[l,2]<-RI[bic_pos[1,1],bic_pos[1,2]]
BIC_RI_K[l,3]<-k_sec[bic_pos[1,1]]
BIC_RI_K[l,4]<-RSS0[bic_pos[1,1],bic_pos[1,2]]
BIC_RI_K[l,5]<-RI_ini[bic_pos[1,1],1]
BIC_RI_K[l,6]<-RI_ini[bic_pos[1,1],2]
  } ####K   
d_pos<-which.min(BIC_RI_K[,"bic"])
re[sam,"bic"]<-BIC_RI_K[d_pos,"bic"]
re[sam,"ri"]<-BIC_RI_K[d_pos,"ri"]
re[sam,"k"]<-BIC_RI_K[d_pos,"k"]
re[sam,"d"]<-d_sec[d_pos]
re[sam,"rss"]<-BIC_RI_K[d_pos,"rss"]
d<-d_sec[d_pos]
#re[sam,"rss0"]<-function_betahat(a)
re[sam,"ri_ini"]<-BIC_RI_K[d_pos,"ri_ini"]
re[sam,"rssmean"]<-BIC_RI_K[d_pos,"rssmean"]
cat("this is the", sam, "-th\n")
}

write.table(re,"KPM_n=45_t=50_s=-0.5.r")



