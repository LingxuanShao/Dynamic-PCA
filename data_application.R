library("ggplot2")
library("expm")
library("stats")
library("ggpubr")
#set dimension p=3;
p=3;
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
    D=t(A)%*%B%*%G[i,,];
    if(p==2){
      if((det(D)>0)&&(D[1,1]>0)){
        newdistance=acos(D[1,1]);
        if(distance>newdistance){
           bestg=i;
           distance=newdistance;
        }
      }
    }
    if(p>2){
     if((det(D)>0)&(norm((D-I),type='F')<0.9)){
        newdistance=norm(logm(D),type='F');
        if(distance>newdistance){
          bestg=i;
          distance=newdistance;
        }
      }
    }
  }
  if(bestg==0){
  for(i in 1:((2**p)*(factorial(p)))){
    D=t(A)%*%B%*%G[i,,];
    if(det(D)>0){
      newdistance=norm(A-B%*%G[i,,],type='F');
      if(distance>newdistance){
        bestg=i;
        distance=newdistance;
      }
    }
  }}
  return(list(bestg=bestg,distance=distance));
}

#make A anti-symmetric
asantisym=function(A){
  return((A-t(A))/2);
}
#make A orth
asorth=function(A){
  return(expm(asantisym(logm(A))));
}
#parallel transport, uA is the tangent vector on T_{A}M
para=function(A,B,u){
  sqrtbma=expm(0.5*logm(B%*%solve(A)));
  sqrtamb=expm(0.5*logm(A%*%solve(B)));
  return(sqrtbma%*%u%*%sqrtamb);
}
#parallel transport u from G[j] to G[1]
piecepara=function(G,j,u){
  startj=j;
  tanu=u;
  while(startj>1){
    tanu=para(G[startj,,],G[startj-1,,],tanu);
    startj=startj-1;
  }
  return(tanu);
}
#parallel transport u from G[1] to G[j]
antipiecepara=function(G,j,u){
  startj=1;
  tanu=u;
  while(startj<j){
    tanu=para(G[startj,,],G[startj+1,,],tanu);
    startj=startj+1;
  }
  return(tanu);
}


##############################################
#import data
#Indicator is the flag showing the 208 countries and locations
#MIL is the military  expenditure
#GDI is gross capital formation
#IMP is the import
#EXP is the export
Indicator <- read.csv("Indicator.csv", header = 0)
API_MS.MIL.XPND.GD.ZS_DS2_en_csv_v2_2055780 <- read.csv("API_MS.MIL.XPND.GD.ZS_DS2_en_csv_v2_2055780.csv", header = 1)
API_NE.GDI.TOTL.ZS_DS2_en_csv_v2_2252056 <- read.csv("API_NE.GDI.TOTL.ZS_DS2_en_csv_v2_2252056.csv", header = 1)
API_NE.EXP.GNFS.ZS_DS2_en_csv_v2_2058633 <- read.csv("API_NE.EXP.GNFS.ZS_DS2_en_csv_v2_2058633.csv", header = 1)
API_NE.IMP.GNFS.ZS_DS2_en_csv_v2_2059780 <- read.csv("API_NE.IMP.GNFS.ZS_DS2_en_csv_v2_2059780.csv", header = 1)
Indicator=as.data.frame(Indicator);
MIL=as.data.frame(API_MS.MIL.XPND.GD.ZS_DS2_en_csv_v2_2055780)[which(Indicator[,3]==1),];
GDI=as.data.frame(API_NE.GDI.TOTL.ZS_DS2_en_csv_v2_2252056)[which(Indicator[,3]==1),];
EXP=as.data.frame(API_NE.EXP.GNFS.ZS_DS2_en_csv_v2_2058633)[which(Indicator[,3]==1),];
IMP=as.data.frame(API_NE.IMP.GNFS.ZS_DS2_en_csv_v2_2059780)[which(Indicator[,3]==1),];
n=dim(MIL)[1];# sample size
d=55;# sampling rate
D=5; #the interpolation time points
t=c(1:d)/d;# the observed time points
anchort=1959+c(1:d);# the year for each observed time points
subj=(c(1:D)-1)*(d/D)+1;#the dummy index of the interpolation points
h=rep(0,(D-1));
for(i in 1:(D-1)){
  h[i]=t[subj[i+1]]-t[subj[i]]; 
}
x=array(0,c(n,d,p)); # the data, IMP-EXP is the net import
for(i in 1:n){
  for(j in 1:d){
    x[i,j,1]=IMP[i,j]-EXP[i,j];
    x[i,j,2]=GDI[i,j];
    x[i,j,3]=MIL[i,j];
  }
}

#obtain \hat{P}_{j}
P=array(0,c(d,p,p));
cova=array(0,c(d,p,p));
initialindex=1;
for(j in initialindex:d){
  cova[j,,]=var(x[,j,],na.rm=TRUE);
  P[j,,]=eigen(cova[j,,])$vectors; #t(P)*var*P=diagonal
  print(j);
}
#modify P such that PG is the new representatives
PG=array(0,c(d,p,p));
DG=array(0,c(d,p,p));
I=diag(array(1,p));
bestgflag=array(0,d);
#PG[1,,]=P[1,,];
PG[1,,]=P[1,,]%*%G[nearest(I,P[1,,],G,p)$bestg,,];
for(j in 2:d){
  temp=nearest(PG[j-1,,],P[j,,],G,p);
  bestg=temp$bestg;
  distance=temp$distance;
  if(bestg==0){
    print(j);
    print("PG is not continuous!");
    break;
  }
  bestgflag[j]=bestg;
  PG[j,,]=P[j,,]%*%G[bestg,,];
}
#DG is the estimate for eigenvalues on the observed time points
for(j in 1:d){
  DG[j,,]=t(PG[j,,])%*%cova[j,,]%*%PG[j,,];
}

#plot DG
plotDG=function(DG){
  e1=DG[,1,1];
  e2=DG[,2,2];
  e3=DG[,3,3];
  plotdata=data.frame(anchort,e1,e2,e3);
  ggplot(plotdata)+
    geom_path(aes(x=anchort,y=e1),linetype="dashed",color="blue",size=0.5)+
    geom_path(aes(x=anchort,y=e2),linetype="solid",color="red",size=0.5)+
    geom_path(aes(x=anchort,y=e3),linetype="dotted",color="black",size=0.5)+
    geom_point(aes(x=anchort,y=e1),color="blue",size=1.5)+
    geom_point(aes(x=anchort,y=e2),color="red",size=1.5)+
    geom_point(aes(x=anchort,y=e3),color="black",size=1.5)+
    xlab(label="")+ylab(label="");
}
plotDG(DG);

#estimate the dynamic eigen-frame
#set the regular base curve r(t)
a=array(0,c(D-1,p,p));
for(i in 1:(D-1)){
  a[i,,]=logm(PG[subj[i+1],,]%*%solve(PG[subj[i],,]))/h[i];
  a[i,,]=asantisym(a[i,,]);
}
z=array(0,c(d,p,p));
for(j in 1:(subj[1]/4+3*subj[2]/4)){
  z[j,,]=a[1,,];
}
for(j in (subj[D-1]*3/4+subj[D]/4):d){
  z[j,,]=a[D-1,,];
}
for(i in 2:(D-2)){
  for(j in (subj[i]*3/4+subj[i+1]/4):(subj[i]/4+subj[i+1]*3/4)){
    z[j,,]=a[i,,]
  }
}
for(i in 2:(D-1)){
  leftj=(subj[i]*3/4+subj[i-1]/4);
  rightj=(subj[i]*3/4+subj[i+1]/4);
  for(j in leftj:rightj){
    z[j,,]=(a[i-1,,]*(rightj-j)+a[i,,]*(j-leftj))/(rightj-leftj);
  }
}
r=array(0,c(d,p,p));
r[1,,]=PG[1,,];
for(j in 1:(d-1)){
  r[j+1,,]=expm((t[j+1]-t[j])*z[j,,])%*%r[j,,];
}
#obtain \tilde{r}
tilder=array(0,c(d,p,p));
tilder[1,,]=array(0,c(p,p));
for(j in 1:(d-1)){
  tilder[j+1,,]=tilder[j,,]+asantisym(piecepara(r,j,z[j,,]))*(t[j+1]-t[j]);
}
#obtain \tilde{[Q]}
tildeq=array(0,c(d,p,p));
for(j in 1:d){
  tildeq[j,,]=tilder[j,,]+asantisym(piecepara(r,j,logm(PG[j,,]%*%solve(r[j,,]))));
}

#calculate lambda by K-fold cross validation with K=5
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

#obtain \tilde{\alpha}
tildealpha=array(0,c(d,p,p));
for(k in 1:p){
  for(l in k:p){
    if(k==l){next;}
    lambda=lambdaselection(t,tildeq[,k,l]);
    cat(k,l,lambda,"\n")
    tildealpha[,k,l]=smooth.spline(t,tildeq[,k,l],lambda=lambda)$y;
    tildealpha[,l,k]=-tildealpha[,k,l];
  }
}
#selected lambda
#k=1, l=2, lambda=0.0005987414 
#k=1, l=3, lambda=0.004424134 
#k=2, l=3, lambda=2.980958e-05 

#obtain \alpha
alpha=array(0,c(d,p,p));
for(j in 1:d){
  alpha[j,,]=expm(antipiecepara(r,j,(tildealpha[j,,]-tilder[j,,])))%*%r[j,,];
}

#plot each column of the estimate alpha 
tar=3;
coe1_alpha=alpha[,1,tar];
coe2_alpha=alpha[,2,tar];
coe3_alpha=alpha[,3,tar];
plotdata=data.frame(anchort,coe1_alpha,coe2_alpha,coe3_alpha);
ggplot(plotdata)+
  geom_path(aes(x=anchort,y=coe1_alpha),linetype="dashed",color="blue",size=1)+
  geom_path(aes(x=anchort,y=coe2_alpha),linetype="solid",color="red",size=1)+
  geom_path(aes(x=anchort,y=coe3_alpha),linetype="dotted",color="black",size=1)+
  xlab(label="")+ylab(label="");




