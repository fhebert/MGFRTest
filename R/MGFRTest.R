ScoreTest = function(X,Y,U=NULL,Y0=NULL,N=1000){
  if(is.null(U)){
    U = matrix(rep(1,length(Y)),ncol=1)
  }else{
    U = cbind(1,U)
  }
  hatY0 = GLMProbs(U,Y)
  Z = as.vector(t(Y-hatY0)%*%X)
  s2Y = sum((Y-hatY0)^2)/length(Y)
  A = t(X)%*%X-t(X)%*%U%*%solve(t(U)%*%U)%*%t(U)%*%X
  Gamma = s2Y*A
  d = sqrt(diag(Gamma))
  Z = Z/d
  D = matrix(rep(d,ncol(X)),ncol=ncol(X))
  Sigma = Gamma/(D*t(D))
  if(is.null(Y0)){
    Y0 = sapply(1:N,function(i){sample(Y)})
  }
  hatY0 = PermGLMProbs(U,Y0)
  d0 = Y0-hatY0
  Z0 = t(d0)%*%X
  sY0 = sqrt((rep(1,length(Y))%*%(d0^2))[1,]/length(Y))
  dA = sqrt(diag(A))
  Z0 = Z0/(matrix(sY0,ncol=1)%*%matrix(dA,nrow=1))
  return(list(Z=Z,Sigma=Sigma,Z0=Z0))
}

OptT2 = function(vtt,lambda,mgamma) {
  m = length(lambda)
  mlambda = rep(1,nrow(mgamma))%*%t(lambda)
  tmp = sqrt(mgamma/mlambda)
  mtmp = rowSums(tmp)
  tmp = tmp/(mtmp%*%t(rep(1,m)))
  tmp = rowSums(tmp*mlambda*mgamma)
  rowSums(mgamma)%*%t(1/(2*vtt))+tmp%*%t(sign(vtt)*abs(m-(1/(2*vtt))*sum(1/lambda)))
}

MGFR = function(Z,Z0,Sigma=NULL,eigSigma=NULL,vtt=NULL){
  if(is.null(eigSigma)){
    eigSigma = eigen(Sigma)
    ind = which(eigSigma$values>10^(-8))
    lambda = eigSigma$values[ind]
    Omega = t(eigSigma$vectors[,ind,drop=FALSE])
  }else{
    lambda = eigSigma$values
    Omega = t(eigSigma$vectors)
  }
  gamma = ((Omega%*%Z)[,1])^2/lambda
  m = length(Z)
  Z0star = tcrossprod(Z0,Omega)
  mgamma = (Z0star^2)/(rep(1,nrow(Z0))%*%t(lambda))
  mlambda = rep(1,nrow(Z0))%*%t(lambda)
  mchi20 = mgamma
  mlambda = rep(1,nrow(mgamma))%*%t(lambda)
  if(is.null(vtt)){
    vtt = seq(-mean(1/lambda),mean(1/lambda),length=1000)
  }
  vpval = rep(0,length(vtt))
  mT20 = OptT2(vtt,lambda,mchi20)
  indNA = which(rowSums(is.na(mT20))>0)
  if(length(indNA)>0){
    mT20[indNA,] = 0
  }
  vT2 = OptT2(vtt,lambda,matrix(gamma,nrow=1))
  vpval = rep(0,length(vtt))
  mpval0 = matrix(0,nrow=nrow(mgamma),ncol=length(vtt))
  for (j in 1:length(vtt)) {
    vpval[j] = mean(abs(mT20[,j])>abs(vT2[j]))
    mpval0[,j] = 1-(order(order(abs(mT20[,j])))-1)/nrow(mT20)
  }
  minpval0 = apply(mpval0,1,min)
  minpval = min(vpval)
  p = mean(minpval0<=minpval)
  return(p)
}

MGFRTest = function(X,Y,U=NULL,Y0=NULL,N=1000){
  tmp = ScoreTest(X,Y,U,Y0,N)
  Z = tmp$Z
  Z0 = tmp$Z0
  Sigma = tmp$Sigma
  p = MGFR(Z,Z0,Sigma)
  return(p)
}
