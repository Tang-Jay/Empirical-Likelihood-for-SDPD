#===========================================================
#                EL for Endogenous SDPD models             #
#===========================================================

# install.packages("spData",type='binary')
# install.packages("terra",type='binary')
# install.packages( "spdep", dependencies = TRUE)

rm(list = ls())
source('GlambdaChen.R')
library('sp')
library('sf')
library('spData')
library('spdep')

f<-function(Matrix,Vector){
  irow = nrow(Matrix)
  v = c(0)
  for(i in 2:irow){
    v[i] = Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
  }
  return(v)
}

Rou <- function(T){
  Rou = matrix(0,T+2,T+2)
  for(j in 2:(T+2)){
    Rou[j,1:(j-1)] = rou^seq((j-2),0,-1)
  }
  return(Rou)
}

#-------------------setting parameters---------------------#
nsim = 100              # nsim = 500 in article
size = c(7,10,13,16,20) 
m = 6                   # Pre-m period
beta = 1; p = 1
gama = 1; q = 1
sigma_zeta = 1; sigma_v = 1
Phi_zeta = (sigma_zeta/sigma_v)^2
dx <- deriv(y ~ (1-x^m)/(1-x), "x") # to generate derivative
a = 0.95
rou_lambda = matrix(NA, nrow=4, ncol=2)
colnames(rou_lambda) <- c("rou","lambda")
rou_lambda[,1] = c(-0.85,0.85,-0.15,0.15)
rou_lambda[,2] = c(-0.8 ,0.8 ,-0.2 ,0.2)

#-----------------Assignment calculation-------------------#
for(l in 1:4){
  # l means l-th parameters
  rou = rou_lambda[l,1]
  lambda = rou_lambda[l,2]
  x = rou^c(1,2)
  re = eval(dx)
  am = re[1]; am_pot = attr(re,"gradient")[1]
  bm = re[2]; bm_pot = attr(re,"gradient")[2]
  
  for(T in c(3,4)){
    # T means T periods 
    k = (T+1)*p+q+1
    cut = qchisq(a,p+q+k+4)
    
    cat('m',' ','rou',' ','lambda','   ','n*T','   ','EL','\n')
    for(j in 1:length(size)){
      m1 = size[j]
      n = m1*m1       # Sample size
      Wnb = cell2nb(m1,m1,type='queen')
      Ws = nb2listw(Wnb)
      Wn = listw2mat(Ws)
      In = diag(n)
      Bn = In - lambda*Wn
      Bni = solve(Bn)
      BtB = t(Bn)%*%Bn
      BtBni = solve(BtB)
      A = BtBni%*%(t(Wn)%*%Bn+t(Bn)%*%Wn)%*%BtBni
      z = rbinom(n, 1, 0.5)
      
      #------------Estimation equation preparation---------#
      n1 = n*(m+1)
      n2 = n*(T+1)
      n3 = n*(m+T+1)
      
      X = rnorm(n*T, 0, 2)
      Z = rep(z,T)
      x0 = rnorm(n, 0, 2)
      x = matrix(c(x0,X),n,(T+1)*p)
      X_wave = cbind(rep(1,n),x,z)
      zm = am*z
      X_s = matrix(NA,n2,(p+q+k))
      X_s[1:n,] = cbind(x0,zm,X_wave)
      X_s[(n+1):n2,] = cbind(X,Z,matrix(rep(0,n*T*k),n*T,k))
      
      F = kronecker(Rou(T-2), In)
      F_s = kronecker(Rou(T-1), In)
      
      vm = rnorm(n*m,0,sigma_v)
      zeta0 = kronecker(t(rou^seq(0,m-1,1)),Bni)%*%vm
      zeta = rnorm(n,0,sigma_zeta)
      u0 = zeta + zeta0
      eta_wave0 = X_wave%*%rep(1,k)+x0+zm
      y0 = eta_wave0+u0
      Y0 = kronecker(rou^seq(0,T-1,1), eta_wave0)
      u1_s = c(am_pot*z,F%*%X*beta+F%*%Z*gama+Y0)
      
      C1 = cbind(In,kronecker(t(c(rou^seq(m-1,0,-1))),Bni))
      C2 = kronecker(diag(T),Bni)
      C = matrix(rep(0,n2*n3),nrow=n2,ncol=n3)
      C[1:n,1:n1] = C1
      C[(n+1):n2,(n1+1):n3] = C2
      
      P = matrix(rep(0,n3^2),nrow=n3,ncol=n3)
      P[1:n,1:n] = kronecker(sqrt(Phi_zeta),In)
      P[(n+1):n3,(n+1):n3] = diag(n*(m+T))
      CP = C%*%P
      
      Ome_s = matrix(rep(0,n2*n2),nrow=n2,ncol=n2)
      Ome_s[1:n,1:n] = Phi_zeta*In+bm*BtBni
      Ome_s[(n+1):n2,(n+1):n2] = kronecker(diag(T),BtBni)
      Ome_s_ni = solve(Ome_s)
      
      Ome_s_rou = matrix(rep(0,n2*n2),nrow=n2,ncol=n2)
      Ome_s_rou[1:n,1:n] = bm_pot*BtBni
      Ome_s_lam = kronecker(diag(c(bm,rep(1,T))),A)
      Ome_s_phi = kronecker(diag(c(1,rep(0,T))),In)
      
      Prou_s = Ome_s_ni%*%Ome_s_rou%*%Ome_s_ni
      Plam_s = Ome_s_ni%*%Ome_s_lam%*%Ome_s_ni
      Pphi_s = Ome_s_ni%*%Ome_s_phi%*%Ome_s_ni
      
      H1 = t(CP)%*%Ome_s_ni%*%CP
      H2 = t(CP)%*%(F_s+t(F_s))%*%CP
      H3 = t(CP)%*%Prou_s%*%CP
      H4 = t(CP)%*%Plam_s%*%CP
      H5 = t(CP)%*%Pphi_s%*%CP
      
      h1 = diag(H1)
      h3 = diag(H3)
      h4 = diag(H4)
      h5 = diag(H5)
      
      a1 =   t(X_s )%*%Ome_s_ni%*%CP
      a2 = 2*t(u1_s)%*%Ome_s_ni%*%CP
      
      #----------------Starting simulation-----------------#
      f1 = 0
      for(i in 1:nsim){
        # i means i-th simulation
        
        # You should select an error distribution as following
        # e = rnorm(n3);sigma2=1
        # e = rt(n3,5);sigma2=5/3
        e = rchisq(n3,4)-4;sigma2=8
        
        # Score function
        g = matrix(NA,nrow=n3,ncol=p+q+k+4)
        g[,1:(p+q+k)] = t(a1)*e
        g[,p+q+k+1]   = h1*(e^2-sigma2) + 2*e*f(H1,e)
        g[,p+q+k+2]   = a2*e + 2*e*f(H2,e) + h3*(e^2-sigma2) + 2*e*f(H3,e)
        g[,p+q+k+3]   = h4*(e^2-sigma2) + 2*e*f(H4,e)
        g[,p+q+k+4]   = h5*(e^2-sigma2) + 2*e*f(H5,e)
        
        # Calculation of El value
        lam = lambdaChen(g)
        el = 2*sum( log(1+t(lam)%*%t(g)) )
        if(el<cut) f1=f1+1
      }
      
      # You will see the results on the screen.
      cat(m,' ',rou,' ',lambda,'  ',n,'*',T,' ',f1/nsim,'\n')
    }
    cat('\n')
  }
  Sys.sleep(10)
}

