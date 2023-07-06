library("ggplot2")
library("expm")
library("stats")
p=3;#the dimension
D=5;#the interpolation time points
d=100;#the number of observed time points
n=1*1000;#sample size
roundnumber=100;#number of Mento Carlo runs
t=c(1:d)/d;#the observed time points
subj=(c(1:D)-1)*(d/D)+1;#the dummy index of the interpolation points.
h=rep(0,(D-1));
for(i in 1:(D-1)){
  h[i]=t[subj[i+1]]-t[subj[i]]; 
}
timestart=proc.time();#time reckon


#preliminaries
#to generate data
phi<-function(t){
  return(2*pi*t);
}
theta<-function(t){
  return(0.5*sin(pi*t));
}
e1<-function(t){
  x1=cos(theta(t));
  x2=sin(theta(t))*sin(phi(t));
  x3=sin(theta(t))*cos(phi(t));
  return(c(x1,x2,x3));
}
e2<-function(t){
  x1=-sin(theta(t));
  x2=cos(theta(t))*sin(phi(t));
  x3=cos(theta(t))*cos(phi(t));
  return(c(x1,x2,x3));
}
e3<-function(t){
  x1=0;
  x2=cos(phi(t));
  x3=-sin(phi(t));
  return(c(x1,x2,x3));
}
ematrix<-function(t){  #the orthorgonal matrix for testing
  e=cbind(e1(t),e2(t),e3(t));
  return(e);
}


#decimal to binary
tentotwo=function(x){
  if(x==0){
    return(0);
  }
  else{
    y=c();
    while(x!=1){
      y=c(y,x%%2);
      x=floor(x/2);
    }
    return(c(1,rev(y)));
  }
}

#generate G
#G means the group and g means a matrix
G=array(0,c((2**p)*(factorial(p)),p,p));
G1=array(0,c(2**p,p,p));
G2=array(0,c(factorial(p),p,p));
#G1
for(i in 1:(2**p)){
  temp=tentotwo((2**p)-i);
  G1[i,,]=diag(c(array(0,p-length(temp)),temp)*2-1);
}
#G2
l=1;k=1;
temp=array(0,p);
while(k>0){
  flag=array(0,p+1);
  for(i in 1:(k-1)){
    flag[temp[i]]=1;
  }
  possi=p+1;
  for(i in 1:(p+1)){
    if((flag[i]==0)&(i>temp[k])){
      possi=i;
      break;
    }
  }
  if(possi<(p+1)){
    temp[k]=possi;
    if(k<p){
      k=k+1;
      next;
    }
    if(k==p){
      for(j in 1:p){
        G2[l,j,temp[j]]=1
      }
      l=l+1;
      temp[k]=0;
      k=k-1;
      next;
    }
  }else{
    temp[k]=0;
    k=k-1;
    next;
  }
}
#G
for(i in 1:(2**p)){
  for(j in 1:factorial(p)){
    G[(i-1)*factorial(p)+j,,]=G1[i,,]%*%G2[j,,];
  }
}


#find g such that d(A,Bg) reaches the minimum 
nearest=function(A,B,G,p){
  distance=100*p;
  bestg=0;
  I=diag(array(1,p));
  for(i in 1:((2**p)*(factorial(p)))){
    D=solve(A)%*%B%*%G[i,,];
    if((det(D)>0)&(norm((D-I),type='F')<0.9)){
      newdistance=norm(asantisym(logm(D)),type='F');
      if(distance>newdistance){
        bestg=i;
        distance=newdistance;
      }
    }
  }
  return(list(bestg=bestg,distance=distance));
}

#make A anti-symmetric
asantisym=function(A){
  return((A-t(A))/2);
}
#make A orthogonal
asorth=function(A){
  return(expm(asantisym(logm(A))));
}
#parallel transport, uA is the tangent vector on T_{A}M, output is on T_{I}M, output%*%B is on T_{B}M
para=function(A,B,u){
  sqrtbma=expm(0.5*asantisym(logm(B%*%solve(A))));
  sqrtamb=expm(0.5*asantisym(logm(A%*%solve(B))));
  return(sqrtbma%*%u%*%sqrtamb);
}
#u is the tangent vector on T_{G[j,,]}M, the output is on T_{G[1,,]}M
piecepara=function(G,j,u){
  startj=j;
  tanu=u;
  while(startj>1){
    tanu=para(G[startj,,],G[startj-1,,],tanu);
    startj=startj-1;
  }
  return(asantisym(tanu));
}
#u is the tangent vector on T_{G[1,,]}M, the output is on T_{G[j,,]}M
antipiecepara=function(G,j,u){
  startj=1;
  tanu=u;
  while(startj<j){
    tanu=para(G[startj,,],G[startj+1,,],tanu);
    startj=startj+1;
  }
  return(asantisym(tanu));
}
#u is the tangent vector on T_{G[j1,,]}M, the output is on T_{G[j2,,]}M
pieceparafield=function(G,j1,j2,u){
  startj=j1;
  tanu=u;
  if(j2==j1){return(asantisym(tanu));}
  if(j2<j1){
    while(startj>j2){
      tanu=para(G[startj,,],G[startj-1,,],tanu);
      startj=startj-1;
    }
    return(asantisym(tanu));
  }
  if(j2>j1){
    while(startj<j2){
      tanu=para(G[startj,,],G[startj+1,,],tanu);
      startj=startj+1;
    }
    return(asantisym(tanu));
  }
}
# Riemannian curvature on T_{I}M
Rie=function(A,B,C){
  return(((A%*%B-B%*%A)%*%C-C%*%(A%*%B-B%*%A))/4);
}
#integral z along gamma, output=initial+\int_{0}^{t}z(s)ds, all stored in T_{I}M
Integral=function(gamma,z,initial){
  output=array(0,c(d,p,p));
  output[1,,]=initial;
  for(j in 1:(d-1)){
    output[j+1,,]=pieceparafield(gamma,j,j+1,(output[j,,]+z[j,,]*(t[j+1]-t[j])));
  }
  return(output);
}



#######gradient-descent method, a directly smoothing spline method
msealphalist=c();
recorldlambda=c();
round=0;
#Mento Carlo starts
while(round<roundnumber){
cat("round=",round,"\n");
#a round of simulation
#generate data
x=array(0,c(n,d,p)); #the original data
for(i in 1:n){
  for(j in 1:floor(d/2)){
    em=ematrix(t[j]);
    for(k in 1:p){
      x[i,j,]=x[i,j,]+rnorm(1,0,100*k)*em[,k];
    }
  }
  for(j in (floor(d/2)+1):d){
    em=ematrix(t[j]);
    for(k in 1:p){
      x[i,j,]=x[i,j,]+rnorm(1,0,100*(p+1-k))*em[,k];
    }
  }
}
#obtain \hat{P}_{j}
P=array(0,c(d,p,p));
for(j in 1:d){
  P[j,,]=eigen(cov(x[,j,]))$vectors; #t(P)*var*P=diag
}
remove(x);
#modify P so that PG is the new representative
PG=array(0,c(d,p,p));
I=diag(array(1,p));
PG[1,,]=P[1,,]%*%G[nearest(ematrix(t[1]),P[1,,],G,p)$bestg,,];
nicePGflag=1;
for(j in 2:d){
  bestg=nearest(PG[j-1,,],P[j,,],G,p)$bestg;
  if(bestg==0){
    cat("wrong PG",j,"\n");
    nicePGflag=0;
    break;
  }
  PG[j,,]=P[j,,]%*%G[bestg,,];
}
if(nicePGflag==0){next;}

#K-fold method to chose gralambda, where the cost function is E=E1+gralambda*E2
gralambdanumber=8;
gralambdalist=exp(c(1:gralambdanumber))*(10**(-6));
Egralambda=array(0,gralambdanumber);
#tau is the length of each step
candidatenumber=5;
candidatetau=exp(c(1:candidatenumber))/(10**3);
#K-fold method
foldnum=5;
validationindex=array(0,c(foldnum,d/foldnum));
estimaindex=array(0,c(foldnum,(d-d/foldnum)));
#the index of the validation group is random
permu=sort(runif(d,0,1),index.return=TRUE)$ix;
for(fold in 1:foldnum){
  validationindex[fold,]=sort(permu[(fold-1)*d/foldnum+c(1:(d/foldnum))]);
}
for(fold in 1:foldnum){
  estimaindex[fold,]=c(1:d)[-validationindex[fold,]];
}
#search gralambda with minimum Egralambda (prediction error)
for(gralambdait in 1:gralambdanumber){
gralambda=gralambdalist[gralambdait];
for(fold in 1:foldnum){
  #cat("fold=",fold,"\n");
  #follow PG[,,], gammagraph records all results in iteration
  itermax=10;
  gammagraph=array(0,c(itermax,d,p,p));
  for(j in estimaindex[fold,]){#set initial curve
    gammagraph[1,j,,]=PG[j,,];
  };
  for(j in validationindex[fold,]){
    newj=estimaindex[fold,which.min(abs(estimaindex[fold,]-j))];
    gammagraph[1,j,,]=gammagraph[1,newj,,];
  }
  #iteration of deepest descent
  for(it in 1:(itermax-1)){
    #the tangent vector u is stored in T_{I}M, representing uA on T_{A}M
    #cat("iteration of deepest descent, it=",it,"\n");
    #G
    nu=array(0,c(d,p,p));
    for(i in estimaindex[fold,]){
      nu[i,,]=-asantisym(logm(PG[i,,]%*%solve(gammagraph[it,i,,])));
    }
    tildenu=array(0,c(d,d,p,p))
    for(i in estimaindex[fold,]){
      tildenu[i,i,,]=nu[i,,];
      if(i==1){
        for(j in (i+1):d){
          tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j-1,j,tildenu[i,j-1,,]);
        }
      }
      if(i==d){
        for(j in (i-1):1){
          tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildenu[i,j+1,,]);
        }
      }
      if((1<i)&&(i<d)){
        for(j in (i+1):d){
          tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j-1,j,tildenu[i,j-1,,]);
        }
        for(j in (i-1):1){
          tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildenu[i,j+1,,]);
        }
      }
    }
    git=array(0,c(d,d,p,p));
    for(i in estimaindex[fold,]){
      for(j in 1:d){
        if(i>=j){
          git[i,j,,]=(1+t[i]*t[j]+t[i]*t[j]*t[j]/2-t[j]*t[j]*t[j]/6)*tildenu[i,j,,];
          git[i,j,,]=asantisym(git[i,j,,]);
        }
        if(i<j){
          git[i,j,,]=(1+t[i]*t[j]+t[i]*t[i]*t[j]/2-t[i]*t[i]*t[i]/6)*tildenu[i,j,,];
          git[i,j,,]=asantisym(git[i,j,,]);
        }
      }
    }
    Git=array(0,c(d,p,p));
    for(i in estimaindex[fold,]){
      Git[,,]=Git[,,]+git[i,,,];
    }
    #Git=Git/length(estimaindex[fold,]);
    #H2
    z=array(0,c(d,p,p)); # first order derivative of gammagraph[it,,,]
    dz=array(0,c(d,p,p)); # second order derivative of gammagraph[it,,,]
    for(j in 1:(d-1)){
      z[j,,]=asantisym(logm(gammagraph[it,j+1,,]%*%solve(gammagraph[it,j,,])))/(t[j+1]-t[j]);
    }
    z[d,,]=pieceparafield(gammagraph[it,,,],d-1,d,z[d-1,,]);
    for(j in 1:(d-1)){
      dz[j,,]=(pieceparafield(gammagraph[it,,,],j+1,j,z[j+1,,])-z[j,,])/(t[j+1]-t[j]);
    }
    dz[d,,]=pieceparafield(gammagraph[it,,,],d-1,d,dz[d-1,,]);
    R=array(0,c(d,p,p));
    for(j in 1:d){
      R[j,,]=Rie(dz[j,,],z[j,,],z[j,,]);
    }
    zeromatrix=array(0,c(p,p));
    hatH23=Integral(gammagraph[it,,,],R,zeromatrix);
    hatH22=Integral(gammagraph[it,,,],hatH23,zeromatrix);
    hatH21=Integral(gammagraph[it,,,],hatH22,zeromatrix);
    hatH2=Integral(gammagraph[it,,,],hatH21,zeromatrix);
    tildeQ=array(0,c(d,p,p));
    tildeS=array(0,c(d,p,p));
    tildeQ[d,,]=hatH22[d,,];
    tildeS[d,,]=hatH23[d,,];
    for(j in (d-1):1){
      tildeQ[j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildeQ[j+1,,]);
      tildeS[j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildeS[j+1,,]);
    }
    H2=array(0,c(d,p,p));
    for(j in 1:d){
      H2[j,,]=hatH2[j,,]-t[j]*t[j]*t[j]*tildeS[j,,]/6-t[j]*t[j]*(tildeQ[j,,]-tildeS[j,,])/2-t[j]*(tildeQ[j,,]-tildeS[j,,])+tildeS[j,,];
      H2[j,,]=asantisym(H2[j,,]);
    }
    #H3
    H31=Integral(gammagraph[it,,,],dz,zeromatrix);
    H3=Integral(gammagraph[it,,,],H31,zeromatrix);
    #total gradient
    totalgra=Git+gralambda*(H2+H3);
    #update
    E1=array(0,candidatenumber);
    E2=array(0,candidatenumber);
    newgamma=array(0,c(candidatenumber,d,p,p))
    for(candi in 1:candidatenumber){
      tau=candidatetau[candi];
      for(j in 1:d){
        newgamma[candi,j,,]=expm(-tau*totalgra[j,,])%*%gammagraph[it,j,,];
      }
      #calculate E
      for(j in estimaindex[fold,]){
        E1[candi]=E1[candi]+(norm(asantisym(logm(PG[j,,]%*%solve(newgamma[candi,j,,]))),type="F"))**2;
      }
      #E1[candi]=E1[candi]/length(estimaindex[fold,]);
      znewgamma=array(0,c(d,p,p)); # first order derivative of newgamma[candi,,,]
      dznewgamma=array(0,c(d,p,p)); # second order derivative of newgamma[candi,,,]
      for(j in 1:(d-1)){
        znewgamma[j,,]=asantisym(logm(newgamma[candi,j+1,,]%*%solve(newgamma[candi,j,,])))/(t[j+1]-t[j]);
      }
      znewgamma[d,,]=pieceparafield(newgamma[candi,,,],d-1,d,znewgamma[d-1,,]);
      for(j in 1:(d-1)){
        dznewgamma[j,,]=(pieceparafield(newgamma[candi,,,],j+1,j,znewgamma[j+1,,])-znewgamma[j,,])/(t[j+1]-t[j]);
      }
      dznewgamma[d,,]=pieceparafield(newgamma[candi,,,],d-1,d,dznewgamma[d-1,,]);
      for(j in 1:d){
        E2[candi]=E2[candi]+norm(dznewgamma[j,,],type="F")**2;
      }
      E2[candi]=gralambda*E2[candi]/d;
    }
    E=E1+E2;
    #cat("candi=",which.min(E),"\n");
    gammagraph[it+1,,,]=newgamma[which.min(E),,,];
    #relative error
    numerator=0;
    denominator=0;
    for(j in 1:d){
      numerator=numerator+norm(gammagraph[it+1,j,,]-gammagraph[it,j,,],type="F")**2;
      denominator=denominator+norm(gammagraph[it,j,,],type="F")**2;
    }
    if((numerator/denominator)<0.001){break;}
  }#ends of the iteration procedure of the deepest-descent algorithm

  #goodness of fit
  preerrorgammagraph=0;
  for(j in validationindex[fold,]){
    preerrorgammagraph=preerrorgammagraph+norm(PG[j,,]-gammagraph[it+1,j,,],type="F")**2;
  }
  preerrorgammagraph=preerrorgammagraph/length(validationindex[fold,]);
  Egralambda[gralambdait]=Egralambda[gralambdait]+preerrorgammagraph;
  }#end of K-fold cross validation
}#end of selecting gralambda

#set gralambda
gralambda=gralambdalist[which.min(Egralambda)];
cat("gralambda=",gralambda,"\n");
recorldlambda=c(recorldlambda,gralambda);
#do smoothing spline for the whole data
rm(gammagraph,nu,tildenu,git,Git,z,dz,R,hatH23,hatH22,hatH21,hatH2,tildeQ,tildeS,H2,H31,H3,newgamma)
#follow PG[,,], gammagraph records all results in iteration
itermax=10;
gammagraph=array(0,c(itermax,d,p,p));
for(j in 1:d){#set initial curve
  gammagraph[1,j,,]=PG[j,,];
};
#iteration of deepest descent
for(it in 1:(itermax-1)){
  #the tangent vector u is stored in T_{I}M, representing uA on T_{A}M
  #cat("iteration of deepest descent, it=",it,"\n");
  #G
  nu=array(0,c(d,p,p));
  for(i in 1:d){
    nu[i,,]=-asantisym(logm(PG[i,,]%*%solve(gammagraph[it,i,,])));
  }
  tildenu=array(0,c(d,d,p,p))
  for(i in 1:d){
    tildenu[i,i,,]=nu[i,,];
    if(i==1){
      for(j in (i+1):d){
        tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j-1,j,tildenu[i,j-1,,]);
      }
    }
    if(i==d){
      for(j in (i-1):1){
        tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildenu[i,j+1,,]);
      }
    }
    if((1<i)&&(i<d)){
      for(j in (i+1):d){
        tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j-1,j,tildenu[i,j-1,,]);
      }
      for(j in (i-1):1){
        tildenu[i,j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildenu[i,j+1,,]);
      }
    }
  }
  git=array(0,c(d,d,p,p));
  for(i in 1:d){
    for(j in 1:d){
      if(i>=j){
        git[i,j,,]=(1+t[i]*t[j]+t[i]*t[j]*t[j]/2-t[j]*t[j]*t[j]/6)*tildenu[i,j,,];
        git[i,j,,]=asantisym(git[i,j,,]);
      }
      if(i<j){
        git[i,j,,]=(1+t[i]*t[j]+t[i]*t[i]*t[j]/2-t[i]*t[i]*t[i]/6)*tildenu[i,j,,];
        git[i,j,,]=asantisym(git[i,j,,]);
      }
    }
  }
  Git=array(0,c(d,p,p));
  for(i in 1:d){
    Git[,,]=Git[,,]+git[i,,,];
  }
  #Git=Git/d;
  #H2
  z=array(0,c(d,p,p)); # first order derivative of gammagraph[it,,,]
  dz=array(0,c(d,p,p)); # second order derivative of gammagraph[it,,,]
  for(j in 1:(d-1)){
    z[j,,]=asantisym(logm(gammagraph[it,j+1,,]%*%solve(gammagraph[it,j,,])))/(t[j+1]-t[j]);
  }
  z[d,,]=pieceparafield(gammagraph[it,,,],d-1,d,z[d-1,,]);
  for(j in 1:(d-1)){
    dz[j,,]=(pieceparafield(gammagraph[it,,,],j+1,j,z[j+1,,])-z[j,,])/(t[j+1]-t[j]);
  }
  dz[d,,]=pieceparafield(gammagraph[it,,,],d-1,d,dz[d-1,,]);
  R=array(0,c(d,p,p));
  for(j in 1:d){
    R[j,,]=Rie(dz[j,,],z[j,,],z[j,,]);
  }
  zeromatrix=array(0,c(p,p));
  hatH23=Integral(gammagraph[it,,,],R,zeromatrix);
  hatH22=Integral(gammagraph[it,,,],hatH23,zeromatrix);
  hatH21=Integral(gammagraph[it,,,],hatH22,zeromatrix);
  hatH2=Integral(gammagraph[it,,,],hatH21,zeromatrix);
  tildeQ=array(0,c(d,p,p));
  tildeS=array(0,c(d,p,p));
  tildeQ[d,,]=hatH22[d,,];
  tildeS[d,,]=hatH23[d,,];
  for(j in (d-1):1){
    tildeQ[j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildeQ[j+1,,]);
    tildeS[j,,]=pieceparafield(gammagraph[it,,,],j+1,j,tildeS[j+1,,]);
  }
  H2=array(0,c(d,p,p));
  for(j in 1:d){
    H2[j,,]=hatH2[j,,]-t[j]*t[j]*t[j]*tildeS[j,,]/6-t[j]*t[j]*(tildeQ[j,,]-tildeS[j,,])/2-t[j]*(tildeQ[j,,]-tildeS[j,,])+tildeS[j,,];
    H2[j,,]=asantisym(H2[j,,]);
  }
  #H3
  H31=Integral(gammagraph[it,,,],dz,zeromatrix);
  H3=Integral(gammagraph[it,,,],H31,zeromatrix);
  #total gradient
  totalgra=Git+gralambda*(H2+H3);
  #update
  E1=array(0,candidatenumber);
  E2=array(0,candidatenumber);
  newgamma=array(0,c(candidatenumber,d,p,p));
  for(candi in 1:candidatenumber){
    tau=candidatetau[candi];
    for(j in 1:d){
      newgamma[candi,j,,]=expm(-tau*totalgra[j,,])%*%gammagraph[it,j,,];
    }
    for(j in 1:d){
      E1[candi]=E1[candi]+(norm(asantisym(logm(PG[j,,]%*%solve(newgamma[candi,j,,]))),type="F"))**2;
    }
    #E1[candi]=E1[candi]/d;
    znewgamma=array(0,c(d,p,p)); # first order derivative of newgamma[candi,,,]
    dznewgamma=array(0,c(d,p,p)); # second order derivative of newgamma[candi,,,]
    for(j in 1:(d-1)){
      znewgamma[j,,]=asantisym(logm(newgamma[candi,j+1,,]%*%solve(newgamma[candi,j,,])))/(t[j+1]-t[j]);
    }
    znewgamma[d,,]=pieceparafield(newgamma[candi,,,],d-1,d,znewgamma[d-1,,]);
    for(j in 1:(d-1)){
      dznewgamma[j,,]=(pieceparafield(newgamma[candi,,,],j+1,j,znewgamma[j+1,,])-znewgamma[j,,])/(t[j+1]-t[j]);
    }
    dznewgamma[d,,]=pieceparafield(newgamma[candi,,,],d-1,d,dznewgamma[d-1,,]);
    for(j in 1:d){
      E2[candi]=E2[candi]+norm(dznewgamma[j,,],type="F")**2;
    }
    E2[candi]=gralambda*E2[candi]/d;
  }
  E=E1+E2;
  #cat("candi=",which.min(E),"\n");
  gammagraph[it+1,,,]=newgamma[which.min(E),,,];
  #relative error
  numerator=0;
  denominator=0;
  for(j in 1:d){
    numerator=numerator+norm(gammagraph[it+1,j,,]-gammagraph[it,j,,],type="F")**2;
    denominator=denominator+norm(gammagraph[it,j,,],type="F")**2;
  }
  if((numerator/denominator)<0.001){break;}
}#ends of the iteration procedure of the deepest-descent algorithm


#goodness of fit
errorgammagraph=0;
for(j in 1:d){
  errorgammagraph=errorgammagraph+norm(ematrix(t[j])-gammagraph[it+1,j,,],type="F")**2;
}
errorgammagraph=errorgammagraph/d;
#record result
round=round+1;
msealphalist=c(msealphalist,errorgammagraph);
} #end of Mento Carlo run
timeend=proc.time();


output=c("gradient-descent algorithm",round,d,n,
         mean(msealphalist),sd(msealphalist),
         mean(recorldlambda),sd(recorldlambda),
         (timeend-timestart)[3]/roundnumber);
output=as.data.frame(output);
print(output);
write.table(output,file=paste("gradient-descent_d",d,"_n",n,".txt",sep=""))
#write.table(output,file=paste("/home/shaolx/EDA/gradient-descent_d",d,"_n",n,".txt",sep=""))



 





