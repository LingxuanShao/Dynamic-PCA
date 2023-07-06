library("ggplot2")
library("expm")
library("stats")
n=1000;#sample size
m=100;#sampling rate
d=2;#dimension
swit=0;#1 for switch and 0 for non-switch
#######generate G
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
#G means the group and g means a matrix
G=array(0,c((2**d)*(factorial(d)),d,d));
G1=array(0,c(2**d,d,d));
G2=array(0,c(factorial(d),d,d));
#G1
for(i in 1:(2**d)){
  temp=tentotwo((2**d)-i);
  G1[i,,]=diag(c(array(0,d-length(temp)),temp)*2-1);
}
#G2
l=1;k=1;
temp=array(0,d);
while(k>0){
  flag=array(0,d+1);
  for(i in 1:(k-1)){
    flag[temp[i]]=1;
  }
  possi=d+1;
  for(i in 1:(d+1)){
    if((flag[i]==0)&(i>temp[k])){
      possi=i;
      break;
    }
  }
  if(possi<(d+1)){
    temp[k]=possi;
    if(k<d){
      k=k+1;
      next;
    }
    if(k==d){
      for(j in 1:d){
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
for(i in 1:(2**d)){
  for(j in 1:factorial(d)){
    G[(i-1)*factorial(d)+j,,]=G1[i,,]%*%G2[j,,];
  }
}

#find g such that d(A,Bg) reaches the minimum 
nearest=function(A,B,G,d){
  distance=100*d;
  bestg=0;
  I=diag(array(1,d));
  for(i in 1:((2**d)*(factorial(d)))){
    D=solve(A)%*%B%*%G[i,,];
    if((det(D)>0)&(norm((D-I),type='F')<0.9)){
      newdistance=norm(logm(D),type='F');
      if(distance>newdistance){
        bestg=i;
        distance=newdistance;
      }
    }
  }
  return(list(bestg=bestg,distance=distance));
}

#eigenvalue (lambda) for switch case
lambswitch=function(j,k){
  t=j/m;
  if(k==1){
    if((0<=t)&&(t<0.25)){
      return(1);
    }
    if((0.25<=t)&&(t<0.5)){
      return(2);
    }
    if((0.5<=t)&&(t<0.75)){
      return(1);
    }
    if((0.75<=t)&&(t<=1)){
      return(2);
    }
  }
  if(k==2){
    if((0<=t)&&(t<0.25)){
      return(2);
    }
    if((0.25<=t)&&(t<0.5)){
      return(1);
    }
    if((0.5<=t)&&(t<0.75)){
      return(2);
    }
    if((0.75<=t)&&(t<=1)){
      return(1);
    }
  }
  return(1);
}
#eigenvalue (lambda) for non-switch case
lambnonswitch=function(j,k){
  t=j/m;
  if(k==1){
    return(2);
  }
  if(k==2){
    return(1);
  }
}
#eigenvalue (lambda)
lamb=function(j,k,swit){
  #swit: 1 for switch and 0 for non-switch
  if(swit==1){
    return(lambswitch(j,k));
  }
  if(swit==0){
    return(lambnonswitch(j,k));
  }
}
    

#generate data
Y=array(0,c(n,m,d));#data
psi=array(0,c(n,m,d));#score
e=array(0,c(m,d,d));#eigen-frame
Var=array(0,c(m,d,d));#variance matrix
for(i in 1:n){
  for(j in 1:m){
    for(k in 1:d){
      psi[i,j,k]=rnorm(1,0,sqrt(lamb(j,k,swit)));
    }
  }
}
for(j in 1:m){
  t=j/m;
  #P(t_{j})=(e[j,,])^{t}
  e[j,1,1]=cos(t);
  e[j,1,2]=sin(t);
  e[j,2,1]=-sin(t);
  e[j,2,2]=cos(t);
}
for(i in 1:n){
  for(j in 1:m){
    for(k in 1:d){
      Y[i,j,]=Y[i,j,]+psi[i,j,k]*e[j,k,];
    }
  }
}
for(j in 1:m){
  Var[j,,]=var(Y[,j,]);
}

####################################method 1: one-step unrolling method
P=array(0,c(m,d,d));
hatD=array(0,c(m,d,d));
for(j in 1:m){
  P[j,,]=eigen(Var[j,,])$vectors;
  hatD[j,,]=diag(eigen(Var[j,,])$values);
}
#modify P such that PG is the new representative
PG=array(0,c(m,d,d));
DG=array(0,c(m,d,d));
I=diag(array(1,d));
PG[1,,]=P[1,,]%*%G[nearest(I,P[1,,],G,d)$bestg,,];
for(j in 2:m){
  bestg=nearest(PG[j-1,,],P[j,,],G,d)$bestg;
  if(bestg==0){
    print(j);break;
  }
  PG[j,,]=P[j,,]%*%G[bestg,,];
}
#DG is the estimate for eigenvalues
for(j in 1:m){
  DG[j,,]=t(PG[j,,])%*%Var[j,,]%*%PG[j,,];
}
#since S1 is flat, no parallel transport is needed 
#Euclidean smoothing spline can be directly applied
#calculate lambda by K-fold cross validation with K=5
foldnum=5;
validationindex=array(0,c(foldnum,m/foldnum));
estimaindex=array(0,c(foldnum,(m-m/foldnum)));
#the index of the validation group is random
permu=sort(runif(m,0,1),index.return=TRUE)$ix;
for(fold in 1:foldnum){
  validationindex[fold,]=sort(permu[(fold-1)*m/foldnum+c(1:(m/foldnum))]);
}
for(fold in 1:foldnum){
  estimaindex[fold,]=c(1:m)[-validationindex[fold,]];
}
lambdaselection=function(t,Y){
  lambdanumber=20;
  lambdalist=exp(c(1:lambdanumber))*(10**(-8));
  MSElambda=array(0,lambdanumber);
  for(lambdait in 1:lambdanumber){
    lambda=lambdalist[lambdait];
    for(fold in 1:foldnum){
      estimate=predict(smooth.spline(t[estimaindex[fold,]],Y[estimaindex[fold,]],lambda=lambda),t[validationindex[fold,]])
      MSElambda[lambdait]=MSElambda[lambdait]+sum((estimate$y-Y[validationindex[fold,]])**2);
    }
  }
  return(lambdalist[which.min(MSElambda)]);
}

roughhattheta1=array(0,m);
hattheta1=array(0,m);
for(j in 1:m){
  roughhattheta1[j]=acos(PG[j,1,1]);
}
#add K-fold
lambda=lambdaselection((c(1:m)/m),roughhattheta1);
cat(lambda,"\n")
hattheta1=smooth.spline((c(1:m)/m),roughhattheta1,lambda=lambda)$y;



###############################method 2: smoothing time-varying covariance method
#under HS norm, it is equivalent to smooth each coordinate separately
smoothVar=array(0,c(m,d,d));
smoothVar[,1,1]=smooth.spline((c(1:m)/m),Var[,1,1],cv=1)$y;
smoothVar[,1,2]=smooth.spline((c(1:m)/m),Var[,1,2],cv=1)$y;
smoothVar[,2,1]=smooth.spline((c(1:m)/m),Var[,2,1],cv=1)$y;
smoothVar[,2,2]=smooth.spline((c(1:m)/m),Var[,2,2],cv=1)$y;
hattheta2=array(0,m);
hatD2_1=array(0,m);
hatD2_2=array(0,m);
for(j in 1:m){
  temp=eigen(smoothVar[j,,]);
  hattheta2[j]=acos(abs((temp$vectors)[1,1]));
  hatD2_1[j]=temp$values[1];
  hatD2_2[j]=temp$values[2];
}



#plot hattheta1 and roughhattheta2
anchort=c(1:m)/m;
plotdata=data.frame(anchort,hattheta1,hattheta2);
ggplot(plotdata)+geom_path(aes(x=anchort,y=anchort),color="red",size=1)+
  geom_path(aes(x=anchort,y=hattheta1),color="blue",linetype="dashed",size=1)+
  geom_point(aes(x=anchort,y=hattheta2),color="black",size=1)+
  xlab(label="")+ylab(label="")+ylim(0,1.8);

#plot DG and hatD2
if(swit==1){
  anchort=c(1:m)/m;
  true=array(0,m);
  for(j in 1:m){
    true[j]=lamb(j,1,swit);
  }
  e1=DG[,1,1];
  e2=hatD2_1;
  plotdata=data.frame(anchort,e1,e2);
  ggplot(plotdata)+geom_path(aes(x=anchort,y=true),color="red",size=1)+
    geom_path(aes(x=anchort,y=e1),color="blue",linetype="dashed",size=1)+
    geom_point(aes(x=anchort,y=e2),color="black",size=1)+
    xlab(label="")+ylab(label="")+ylim(0.9,2.5);
}
if(swit==0){
  anchort=c(1:m)/m;
  true=array(0,m);
  for(j in 1:m){
    true[j]=lamb(j,1,swit);
  }
  e1=DG[,1,1];
  e2=hatD2_1;
  plotdata=data.frame(anchort,e1,e2);
  ggplot(plotdata)+geom_path(aes(x=anchort,y=true),color="red",size=1)+
    geom_path(aes(x=anchort,y=e1),color="blue",linetype="dashed",size=1)+
    geom_point(aes(x=anchort,y=e2),color="black",size=1)+
    xlab(label="")+ylab(label="")+ylim(0.9,2.5);
}

