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
noisetype=2;#0 for no noise; 1 for normal, 2 for uniform, 3 for laplace.
h=rep(0,(D-1));
for(i in 1:(D-1)){
  h[i]=t[subj[i+1]]-t[subj[i]]; 
}
timestart=proc.time();#time reckon


#generate data
#e1 e2 e3 form the eigen-frame
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
ematrix<-function(t){  #the orthogonal matrix for testing
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
      newdistance=norm(logm(D),type='F');
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
#parallel u from G[j] to G[1]
piecepara=function(G,j,u){
  startj=j;
  tanu=u;
  while(startj>1){
    tanu=para(G[startj,,],G[startj-1,,],tanu);
    startj=startj-1;
  }
  return(tanu);
}
#parallel u from G[1] to G[j]
antipiecepara=function(G,j,u){
  startj=1;
  tanu=u;
  while(startj<j){
    tanu=para(G[startj,,],G[startj+1,,],tanu);
    startj=startj+1;
  }
  return(tanu);
}

#Plot vectors: the Mercator projection
mercator<-function(a,tra){
  longi=longitude(c(a[1],a[2]));
  if(longi>=tra){
    mercatorx=longi-tra;
  }else{
    mercatorx=longi+2*pi-tra;
  }
  mercatory=log(tan(pi/4+asin(a[3])/2));
  return(list(mercatorx=mercatorx,mercatory=mercatory));
}
#longitude ranges in (0,2pi)
longitude<-function(a){
  x=a[1]/(sqrt(sum(a*a)));
  y=a[2]/(sqrt(sum(a*a)));
  if(y==0){
    if(x>0){
      return(0);
    }
    if(x<0){
      return(pi);
    }
  }
  if(y>0){
    return(-asin(x)+0.5*pi);
  }
  if(y<0){
    return(asin(x)+1.5*pi);
  }
}


#calculate lambda
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
  lambdanumber=8;
  lambdalist=exp(c(1:lambdanumber))*(10**(-5));
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

#calculate mse by Mento Carlo method
msePGlist=c();
msealphalist=c();
recordlambda12=c();
recordlambda13=c();
recordlambda23=c();
#Mento Carlo run starts
round=0;
while(round<roundnumber){
cat("round=",round,"\n");
# a round of simulation
# generate data
x=array(0,c(n,d,p)); #the original data
for(i in 1:n){
  for(j in 1:d){
   if(j<=floor(d/2)){
      em=ematrix(t[j]);
      for(k in 1:p){
        x[i,j,]=x[i,j,]+rnorm(1,0,100*k)*em[,k];#the scores. In R: rnorm(n,mean,sd)
     }
   }
   if(j>=(floor(d/2)+1)){
     em=ematrix(t[j]);
     for(k in 1:p){
       x[i,j,]=x[i,j,]+rnorm(1,0,100*(p+1-k))*em[,k];#the scores 
     }
   }
  #the noise
    if(noisetype==1){x[i,j,]=x[i,j,]+rnorm(3,0,sqrt(28000/3));}
    if(noisetype==2){x[i,j,]=x[i,j,]+runif(3,-sqrt(28000),sqrt(28000));}
    if(noisetype==3){
      U=runif(3,-0.5,0.5);
      x[i,j,]=x[i,j,]+sqrt(28000/6)*sign(U)*log(1-2*abs(U));}
  }
}
#obtain \hat{P}_{j}
P=array(0,c(d,p,p));
for(j in 1:d){
  P[j,,]=eigen(cov(x[,j,]))$vectors; #t(P)*var*P=diag
}
#modify P such that PG is the new representative
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
if(nicePGflag==0){next;}#next run

# ######DG is the estimate for eigenvalues on the observed time points
# DG=array(0,c(d,p,p));
# for(j in 1:d){
#   DG[j,,]=t(PG[j,,])%*%cov(x[,j,])%*%PG[j,,];
# }
# #plot DG
# plotDG=function(DG){
#   anchort=t;
#   e1=DG[,1,1];
#   underlyinge1=array(0,d);
#   contaminatede1=array(0,d);
#   for(j in 1:d){
#     if(j<=floor(d/2)){
#       underlyinge1[j]=(100*1)**2;
#       contaminatede1[j]=(100*1)**2+28000/3;
#     }
#     if(j>=(floor(d/2)+1)){
#       underlyinge1[j]=(100*(p+1-1))**2;
#       contaminatede1[j]=(100*(p+1-1))**2+28000/3;
#     }
#   }
#   plotdata=data.frame(anchort,e1,underlyinge1,contaminatede1);
#   ggplot(plotdata)+
#     geom_path(aes(x=anchort,y=e1),linetype="dashed",color="blue",size=1)+
#     geom_path(aes(x=anchort,y=underlyinge1),linetype="dotted",color="black",size=1)+
#     geom_path(aes(x=anchort,y=contaminatede1),color="red",size=1)+
#     xlab(label="")+ylab(label="")+
#     scale_y_continuous(limits=c(10**4,11*(10**4)));
# }
# plotDG(DG);
# round=round+1;}

######################################
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
#obtain \tilde{\alpha}
tildealpha=array(0,c(d,p,p));
for(k in 1:p){
  for(l in k:p){
    if(k!=l){
      lambda=lambdaselection(t,tildeq[,k,l]);
      if((k==1)&&(l==2)){recordlambda12=c(recordlambda12,lambda);}
      if((k==1)&&(l==3)){recordlambda13=c(recordlambda13,lambda);}
      if((k==2)&&(l==3)){recordlambda23=c(recordlambda23,lambda);}
      tildealpha[,k,l]=smooth.spline(t,tildeq[,k,l],lambda=lambda)$y;
      tildealpha[,l,k]=-tildealpha[,k,l];
    }
  }
}
#obtain \alpha
alpha=array(0,c(d,p,p));
for(j in 1:d){
  alpha[j,,]=expm(antipiecepara(r,j,(tildealpha[j,,]-tilder[j,,])))%*%r[j,,];
}
#goodness of fit
msealpha=0;
for(j in 1:d){
  msealpha=msealpha+norm(ematrix(t[j])-alpha[j,,],type="F")**2;
}
msealpha=msealpha/d;
#record result
round=round+1;
msealphalist=c(msealphalist,msealpha);
} #Mento Carlo methods
timeend=proc.time();

output=c("One-step unrolling",round,d,n,
         mean(msealphalist),sd(msealphalist),
         mean(recordlambda12),sd(recordlambda12),
         mean(recordlambda13),sd(recordlambda13),
         mean(recordlambda23),sd(recordlambda23),
         (timeend-timestart)[3]/roundnumber);
output=as.data.frame(output);
print(output);
#write.table(output,file=paste("/home/shaolx/EDA/onestepunrolling_noise",noisetype,"_d",d,"_n",n,".txt",sep=""))
write.table(output,file=paste("onestepunrolling_noise",noisetype,"_d",d,"_n",n,".txt",sep=""))






