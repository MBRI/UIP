clear;
strFolder = '';
% strFolder = 'D:\My_Documents\Density_forecasts\Code\TVAR\files\';


addpath(strcat(strFolder,'functions'));
% Number of the first/last data file
minfile= 354;  % [PG 7.1.13: set to 0 to include 0:T00 sample]
maxfile=354; 
scalex=0.1;
normalise=1; %to normalise the shock=scalex in bothe regimes. Otherwise it
%is equal to scalex times SD of shock in each regime
% Gibbs sampling reps
REPS = 20500;
BURN = 20000;
Update = 10;
MaxTrys=1000;
% Update is the step (no of reps) at which the screen is updated.

% Forecast horizon, initial sample and training sample:
HORZ=12;
HORZIR=40;

% The 1st subsample is t=(1, ..., T00), and the first forecast is for
% t=T00+1.

% VAR models:
VarBench.L=13;                       %lag length of the VAR
VarBench.lamdaP=0.1;                %This controls the tightness of the priors on the first lag
VarBench.tauP=10*VarBench.lamdaP;   % this controls the tightness of the priors on sum of coefficients
VarBench.epsilonP=1;                % this controls tightness of the prior on the constant
VarBench.TarVAR=4;                % this specifies the column in Y that is the threshold variable
VarBench.TarD=12;                % this specifies the max delay
VarBench.TarVariance=10;                % this specifies the variance of the prior for threshold value
VarBench.TarScale=0.01;                % this specifies the variance of random walk proposal



dirname=strcat(pwd,filesep);  %directory for data
fdirname=strcat(strFolder,filesep,'forecasts',filesep,'TAR',filesep); %directory to save forecast densities

% Redefine paras (just to use shorter names):
L       = VarBench.L;
lamdaP  = VarBench.lamdaP;
tauP    = VarBench.tauP;
epsilonP= VarBench.epsilonP;
tarvar=VarBench.TarVAR;
tard=1:VarBench.TarD;
tarvariance=VarBench.TarVariance;
tarscale=VarBench.TarScale;
% Display progress:
disp(sprintf(' '));
disp(sprintf('============================================================'));
disp(sprintf('BENCHMARK MODEL: %s replications with %s burns.', num2str(REPS), num2str(BURN)));


for file=minfile:maxfile %recursive estimation
    
% Get input data    
fname=strcat(dirname,'data',num2str(file),'.mat');
fname1=strcat(fdirname,'forecast',num2str(file),'.mat'); %name of file to save forecast density
load(fname);
data=dataout;

Y=data;
N=cols(Y);
ncrit=(N*L+1);

%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];

%compute threshold variable
Ystar=lag0(Y(:,tarvar),tard(1));


Y=Y(max([L,tard(1)])+1:end,:);
X=X(max([L,tard(1)])+1:end,:);
Ystar=Ystar(max([L,tard(1)])+1:end,:);
tarmean=mean(Ystar);  %mean of the prior on the threshold is the mean value of the threshold variable

% Additional priors for VAR coefficients
muP=mean(Y)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if abs(btemp(1))>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
    e0=[e0 etemp];
end

%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);


T=rows(Y);


fsave=zeros(REPS-BURN,HORZ,N); 
jgibbs=1;

% Display empty line (to separate samples in the screen shot)
disp(sprintf(' ')) 


  Y0=[Y;yd];
  X0=[X;xd];
  
  sigma1=eye(N); %starting value for sigma
  sigma2=eye(N);
  beta0=vec(X0\Y0);
  beta01=beta0;
  beta02=beta0;
  tar=tarmean; %initial value of the threshold 
  tarold=tar;
  naccept=0;
  b0filter=[Y(1,:) X(1,1:cols(X)-(N+1))];
  p0filter=eye(cols(b0filter)).*0.00000001;
  
  naccept=0;
  igibbs=1;
  jgibbs=0;
while jgibbs<REPS-BURN
	
    
    
    %step 1: Seperate into two regimes
    e1=Ystar<=tar;
    e2=Ystar>tar;
    
    Y1=Y(e1,:);
    X1=X(e1,:);
    
    Y2=Y(e2,:);
    X2=X(e2,:);
    
    %step 2 Sample Coefficients and variance regime 1
    
    Y0=[Y1;yd];
    X0=[X1;xd];
  %conditional mean of the VAR coefficients
  mstar1=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx1=xx\eye(cols(xx));
   [ beta1,PROBLEM1] = getcoef( mstar1,sigma1,ixx1,MaxTrys,N,L );
     if PROBLEM1
         beta1=beta01;
     else
         beta01=beta1;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta1,N*L+1,N);
    scale=e'*e;
    sigma1=iwpq(rows(Y0),inv(scale));  
    
    
    %step 3 Sample Coefficients and variance in regime 2
    
     Y0=[Y2;yd];
    X0=[X2;xd];
  %conditional mean of the VAR coefficients
  mstar2=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx2=xx\eye(cols(xx));
   [ beta2,PROBLEM2] = getcoef( mstar2,sigma2,ixx2,MaxTrys,N,L );
     if PROBLEM2
         beta2=beta02;
     else
         beta02=beta2;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma2=iwpq(rows(Y0),inv(scale)); 
    
    
    
    %step 4 Sample Threshold via a Random Walk Metropolis Step
    
    tarnew=tarold+randn(1,1)*sqrt(tarscale);
     %compute conditional posterior at the old and new draw
%     postnew=kfilterTar(Y,beta1,beta2,sigma1,sigma2,L,b0filter,p0filter,tarnew,tarmean,tarvariance,Ystar,ncrit);
%     postold=kfilterTar(Y,beta1,beta2,sigma1,sigma2,L,b0filter,p0filter,tarold,tarmean,tarvariance,Ystar,ncrit);

      postnew=getvarpost(Y,X,beta1,beta2,sigma1,sigma2,L,tarnew,tarmean,tarvariance,Ystar,ncrit);
      postold=getvarpost(Y,X,beta1,beta2,sigma1,sigma2,L,tarold,tarmean,tarvariance,Ystar,ncrit);

     accept=exp(postnew-postold);
     u=rand(1,1);
     if u<accept
         tarold=tarnew;
         naccept=naccept+1;
     end
     tar=tarold;
     arate=naccept/igibbs;
     
     %step 5 sample delay parameter
     prob=[];
     for jj=1:length(tard)
         [yy,xx,yys]=prepare( data,L,tard(jj),tarvar );
         
         postx=getvarpost(yy,xx,beta1,beta2,sigma1,sigma2,L,tar,tarmean,tarvariance,yys,ncrit);
         prob=[prob postx];
     end
     prob=safeexp(prob);
     prob=prob/sum(prob);
     
     %draw delay
     index=discretesample(prob,1);
     %update data
     [Y,X,Ystar]=prepare( data,L,tard(index),tarvar );
     
     
     
     
    % Display progress:
    if mod(igibbs,Update)==0 
        disp(sprintf(' Sample %s    Replication %s of %s acceptance %s. delay is %s.', ... 
            num2str(file), num2str(igibbs), num2str(REPS),num2str(arate),num2str(tard(index)) ));
    end 
     
     
     
     if igibbs>100 && igibbs<5000
         if arate<0.35
             tarscale=tarscale*0.99;
         elseif arate>0.55
             tarscale=tarscale*1.01;
         end
     end
    
    

     if igibbs>BURN && ~PROBLEM1 &&~ PROBLEM2
         jgibbs=jgibbs+1;
         % Simulate and store forward paths:

  [r1yhat4]=get_irfTVARS(N, HORZIR, sum(e1), L, Y(e1,:), beta1, sigma1,beta2,sigma2,tar,tarvar,tard(index),50,scalex,normalise);
  [r2yhat4]=get_irfTVARS(N, HORZIR, sum(e2), L, Y(e2,:), beta1, sigma1,beta2,sigma2,tar,tarvar,tard(index),50,scalex,normalise);
fsave1(jgibbs,:,:)=r1yhat4;
fsave2(jgibbs,:,:)=r2yhat4;

 [r1yhat4m]=get_irfTVARS(N, HORZIR, sum(e1), L, Y(e1,:), beta1, sigma1,beta2,sigma2,tar,tarvar,tard(index),50,scalex*2,normalise);
  [r2yhat4m]=get_irfTVARS(N, HORZIR, sum(e2), L, Y(e2,:), beta1, sigma1,beta2,sigma2,tar,tarvar,tard(index),50,scalex*2,normalise);
fsave1m(jgibbs,:,:)=r1yhat4m;
fsave2m(jgibbs,:,:)=r2yhat4m;



   bsave1(jgibbs,:)=beta1';
   bsave2(jgibbs,:)=beta2';
   sigmaS1(jgibbs,:,:)=sigma1;
   sigmaS2(jgibbs,:,:)=sigma2;   
   tsave(jgibbs,:)=[tar tard(index)];
   regime(jgibbs,:)=e1;
     end
     igibbs=igibbs+1;
end

 save(fname1,'fsave');

end

% Display progress:
disp(sprintf(' '));
disp(sprintf('BENCHMARK MODEL - done'));

save results.mat;


figure(1)
tmp1=prctile(fsave1(:,:,1),[50 16 84],1);
tmp2=prctile(fsave2(:,:,1),[50 16 84],1);
h=0:39;
subplot(4,2,1);
plotx2(h,tmp1');
ylabel('IP');
title('Regime 1');

subplot(4,2,2);
plotx2(h,tmp2');
ylabel('IP');
title('Regime 2');
%
tmp1=prctile(fsave1(:,:,2),[50 16 84],1);
tmp2=prctile(fsave2(:,:,2),[50 16 84],1);
h=0:39;
subplot(4,2,3);
plotx2(h,tmp1');
ylabel('R');
title('Regime 1');

subplot(4,2,4);
plotx2(h,tmp2');
ylabel('R');
title('Regime 2');

%
tmp1=prctile(fsave1(:,:,3),[50 16 84],1);
tmp2=prctile(fsave2(:,:,3),[50 16 84],1);
h=0:39;
subplot(4,2,5);
plotx2(h,tmp1');
ylabel('\pi');
title('Regime 1');

subplot(4,2,6);
plotx2(h,tmp2');
ylabel('\pi');
title('Regime 2');

%
tmp1=prctile(fsave1(:,:,4),[50 16 84],1);
tmp2=prctile(fsave2(:,:,4),[50 16 84],1);
h=0:39;
subplot(4,2,7);
plotx2(h,tmp1');
ylabel('FCI');
title('Regime 1');

subplot(4,2,8);
plotx2(h,tmp2');
ylabel('FCI');
title('Regime 2');




figure(2)
tmp1=prctile(fsave1m(:,:,1),[50 16 84],1);
tmp2=prctile(fsave2m(:,:,1),[50 16 84],1);
h=0:39;
subplot(4,2,1);
plotx2(h,tmp1');
ylabel('IP');
title('Regime 1');

subplot(4,2,2);
plotx2(h,tmp2');
ylabel('IP');
title('Regime 2');
%
tmp1=prctile(fsave1m(:,:,2),[50 16 84],1);
tmp2=prctile(fsave2m(:,:,2),[50 16 84],1);
h=0:39;
subplot(4,2,3);
plotx2(h,tmp1');
ylabel('R');
title('Regime 1');

subplot(4,2,4);
plotx2(h,tmp2');
ylabel('R');
title('Regime 2');

%
tmp1=prctile(fsave1m(:,:,3),[50 16 84],1);
tmp2=prctile(fsave2m(:,:,3),[50 16 84],1);
h=0:39;
subplot(4,2,5);
plotx2(h,tmp1');
ylabel('\pi');
title('Regime 1');

subplot(4,2,6);
plotx2(h,tmp2');
ylabel('\pi');
title('Regime 2');

%
tmp1=prctile(fsave1m(:,:,4),[50 16 84],1);
tmp2=prctile(fsave2m(:,:,4),[50 16 84],1);
h=0:39;
subplot(4,2,7);
plotx2(h,tmp1');
ylabel('FCI');
title('Regime 1');

subplot(4,2,8);
plotx2(h,tmp2');
ylabel('FCI');
title('Regime 2');




figure(5)
T=1974+(3/12):(1/12):2012+(8/12);
z1=mean(regime==0,1);
names={'1-S','IP','R','\pi','financial indicator'};
subplot(1,1,1)
[AX,H1,H2] = plotregx( T,z1,Y,names );
title('TVAR regimes using Financial Conditions index');