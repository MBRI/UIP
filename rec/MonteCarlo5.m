% Versoin 1
%--------------------------------------------------------------------------
%
% Author:  Pedram Davoudi
%  p.davoudi@mbri.ac.ir
%
%--------------------------------------------------------------------------
function [Res]=MonteCarlo5(Dt,vU,vR,SameDatabase,alpha,sims,NPer,HasCons,constrA,InsPos,TrgPos)
%NPer=300;% Horizon of impulse response
% Dt is the IRIS database coressponde to VAR estimated model, vU and vR
% vU is the IRIS unrestricted VAR Results
% vR is the IRIS Rstricted VAR Results
% alpha is significance level
% sims is the number of simulation
% NPer is the horizon of IRFs
% HasCons the model has Constant term or not
% constrA is the Constraint imposed on restricted model
% Hsit is the estimataion strathist:endhist
% InsPos is the position of instrument in the endogenous variable.
% ChnPos is the position of Channel in the endogenous variable.
% TrgPos is the position of Target in the endogenous variable(s).
%% Get the estimated coefficients in the transition matrix and bootstrap the irfs
% convergency bound
CB=0.999;
burning=150;

AU = get(vU,'A*'); % Transition
KU = get(vU,'K'); % Const
OmgU = get(vU,'Omega');
lagsU=get(vU,'order');

AR = get(vR,'A*'); %?A%?
KR = get(vR,'K'); %?K?
lagsR=get(vR,'order');
if ~SameDatabase
    OmgR = get(vR,'Omega');
else
    OmgR =OmgU;
end
% 1. Estimating the Parameters of the Model
Evar=get(vU,'yList');
K=length(Evar);

X=Dt.(Evar{1}).data;
% for j=2:K
%     X=[X Dt.(Evar{j}).data]; %#ok<AGROW>
% end

T=length(X);

if T<800
    T=800;
end
%K=K%-1; %Constant Term\
%T=T-lagsU;% Remove lagsU
%per=T-lagsU;%+1;

%% 2. Generate Bootstrap Residuals ****  re-centred  ***** not autocorrelated
% **** contemporaneously correlated
mu = zeros(1,K);
sigma = OmgU;
rng default  % For reproducibility
U=nan(T,K,sims);
for i=1:sims
    r = mvnrnd(mu,sigma,T);
    % **  re-centred  **
    U(:,:,i)=r;%-repmat(mean(r,1),T,1);
end
if ~SameDatabase
    mu = zeros(1,K);
    sigma = OmgR;
    UR=nan(T,K,sims);
    for i=1:sims
        r = mvnrnd(mu,sigma,T);
        % **  re-centred  **
        UR(:,:,i)=r-repmat(mean(r,1),T,1);
    end
else
    UR=U;
end
%clear mu sigma r i
%% 3. Creating bootstrap datasets *** recursively *****

% 3-1. Created monteCarlo matrix for each simulation for unrestricted matrix

% unconditional mean comp
% Created tranistion coeficient
A=reshape(AU,K,[]);
for i=1:lagsU-1
    A=[A;zeros(K,(i-1)*K),eye(K,K),zeros(K,(lagsU-i)*K)];
end
KU=[KU; zeros((lagsU-1)*K, 1)];
% Starting Value
X = (eye(lagsU*K)-A)^-1*KU;

% burning X

%fill Simlation with nan
SimU=nan(T,K,sims);


for s=1:sims
    for i=1:burning
        X=KU+A*X+[ mvnrnd(mu,OmgU).';zeros((lagsU-1)*K, 1)];
    end
    
    for t=1: T
        X=KU+A*X+[U(t,:,s).';zeros((lagsU-1)*K, 1)];
        SimU(t,:,s)=X(1:K);
    end
end
% 3-2. Created monteCarlo matrix for each simulation for Restricted matrix

% unconditional mean comp
KR=[KR; zeros((lagsR-1)*K, 1)];
A=reshape(AR,K,[]);
for i=1:lagsR-1
    A=[A;zeros(K,(i-1)*K),eye(K,K),zeros(K,(lagsR-i)*K)];
end
X = (eye(lagsR*K)-A)^-1*KR;


%fill Simlation with nan
SimR=nan(T,K,sims);


for s=1:sims
    % burning X
    for i=1:burning
        X=KR+A*X+[ mvnrnd(mu,OmgR).';zeros((lagsR-1)*K, 1)];
    end
    
    
    
    for t=1: T
        X=KR+A*X+[UR(t,:,s).';zeros((lagsR-1)*K, 1)];
        SimR(t,:,s)=X(1:K);
    end
end


clear U UR Q1 Q2
%% 4. Estimating the IRF for Each Bootstrap Dataset
% 4-1. Unrestricted
AsU=nan(K,lagsU*K+1,sims); % one for Constant
OmgU=nan(K,K,sims);
IRFsU=nan(NPer,length(TrgPos),sims);
DtT=Dt;
for s=1:sims
    % Set Dataset by Simulation data
    for j=1:K
        DtT.(Evar{j}).data=SimU(:,j,s);
    end
       startHist = get(DtT.(Evar{1}),'start');
    endHist = get(DtT.(Evar{1}),'end');
    
    [v,~] = estimate(vU,DtT,startHist:endHist, 'order=',lagsU,'const=',HasCons,'covParameters=',true);
    if max(abs(eig(v)))<CB % Remove Divergent IRFs
        [Resp,~] = srf(SVAR(v),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
        for tt=1:length(TrgPos)
            IRFsU(:,tt,s)=Resp.(Evar{TrgPos(tt)}).data; %
        end
        
        
        % Store Coef
        AsU(:,:,s)=[get(v,'K'),reshape(get(v,'A*'),K,[])];
        OmgU(:,:,s)=get(v,'Omega');
    end
    %    As(s,:,:)=AA.';
end

% 4-2. Restricted
AsR=nan(K,lagsR*K+1,sims); % one for Constant
OmgR=nan(K,K,sims);
IRFsR=nan(NPer,length(TrgPos),sims);
DtT=Dt;
for s=1:sims
    % Set Dataset by Simulation data
    for j=1:K
        DtT.(Evar{j}).data=SimR(:,j,s);
    end
    startHist = get(DtT.(Evar{1}),'start');
    endHist = get(DtT.(Evar{1}),'end');
    [v,~] = estimate(vR,DtT,startHist:endHist, 'order=',lagsR,'const=',HasCons,'covParameters=',true,'A=',constrA);
    if max(abs(eig(v)))<CB % Remove Divergent IRFs
        [Resp,~] = srf(SVAR(v),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
        for tt=1:length(TrgPos)
            IRFsR(:,tt,s)=Resp.(Evar{TrgPos(tt)}).data; %
        end
        
        
        % Store Coef
        AsR(:,:,s)=[get(v,'K'),reshape(get(v,'A*'),K,[])];
        OmgR(:,:,s)=get(v,'Omega');
    end
end

% 4-3. Refine IRFs
% find the oulier irf in the last perid
[ind1,~]=find(isnan(squeeze(IRFsU(1,1,:))));
[ind2,~]=find(isnan(squeeze(IRFsR(1,1,:))));
ind=unique([ind1;ind2]);
IRFsU(:,:,ind)=[];
IRFsR(:,:,ind)=[];
AsR(:,:,ind)=[];
OmgR(:,:,ind)=[];
AsU(:,:,ind)=[];
OmgU(:,:,ind)=[];
% a=quantile(IRFsU(:,end),0.999999,1);
% IRFsU(IRFsU(:,end)>=a,:)=[];
% a=quantile(IRFsU(:,end),1-0.999999,1);
% IRFsU(IRFsU(:,end)<=a,:)=[];
%
% a=quantile(IRFsR(:,end),0.999999,1);
% IRFsR(IRFsR(:,end)>=a,:)=[];
% a=quantile(IRFsR(:,end),1-0.999999,1);
% IRFsR(IRFsR(:,end)<=a,:)=[];

%% 5. Confidence interval & Statistical Tests

% 5-1. IRFs statistics
for tt=1:length(TrgPos)
    IU=squeeze(IRFsU(:,tt,:));
    IRFU.Simul=IU;
    IRFU.UBound=quantile(IU,1-alpha/2,2);
    IRFU.LBound=quantile(IU,alpha/2,2);
    IRFU.Median=quantile(IU,1/2,2);
    IRFU.Mean=mean(IU,2);
    IRFU.std=std(IU,0,2);
    
    IR=squeeze(IRFsR(:,tt,:));
    IRFR.Simul=IR;
    IRFR.UBound=quantile(IR,1-alpha/2,2);
    IRFR.LBound=quantile(IR,alpha/2,2);
    IRFR.Median=quantile(IR,1/2,2);
    IRFR.Mean=mean(IR,2);
    IRFR.std=std(IR,0,2);
    
    
    DI=(IU-IR).^2;% abs Be or not to be
    Dif.Simul=DI;
    Dif.UBound=quantile(DI,1-alpha/2,2);%Uper bound
    Dif.LBound=quantile(DI,alpha/2,2);%lower bound
    Dif.Median=quantile(DI,1/2,2);%median
    Dif.Mean=mean(DI,2);
    Dif.std=std(DI,0,2);
    
    % Chi2
    
    [RespU,~] = srf(SVAR(vU),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
    [RespR,~] = srf(SVAR(vR),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
    
    RespU=RespU.(Evar{TrgPos(tt)}).data;
    RespR=RespR.(Evar{TrgPos(tt)}).data;
    
    %T=size(IRFsU,1);
    %RespU=RespU-mean(IRFsU,2)./std(IRFsU,2);
    %RespR=RespR-mean(IRFsR,2)./std(IRFsR,2);
    CI=(IU-IR);% abs Be or not to be
    Resp=RespU-RespR;
    %CI=(CI-repmat(mean(CI,1),T,1))./repmat(std(CI,1),T,1);
    Resp=Resp-mean(CI,2)./std(CI,0,2); % if the std ==0 then Resp=nan
    Test.CHI2.Simul=Resp;
    % CHI2.stat=sum(Resp.^2,2);
    Test.CHI2.stat=sum(Resp(~isnan(Resp)).^2);
    Test.CHI2.prob=chi2cdf(Test.CHI2.stat,NPer-sum(isnan(Resp))); % the Nan is not removed
    
    % Ttest
    % pairwise t test
    [h,p,~,stats]=ttest2(IU,IR,'Dim',2,'Alpha',alpha);
    h(isnan(h))=0;
    Test.t.Hypo=h; % one: the null hyphothesis is rejected, so the means are not equal
    Test.t.prob=p;
    Test.t.stats=stats.tstat;
    Test.t.df=stats.df;
    Test.t.sd=stats.sd;
    % Store in the output
    Res.(['Trg' num2str(TrgPos(tt))]).IRFsU=IRFU;
    Res.(['Trg' num2str(TrgPos(tt))]).IRFsR=IRFR;
    Res.(['Trg' num2str(TrgPos(tt))]).Diff=Dif;
    Res.(['Trg' num2str(TrgPos(tt))]).Test=Test;
    clear  IRFU IRFR Dif Test
end
% 5-2. Coef Statistics
for i=1:K
    for j=1:HasCons+lagsU*K
        AU=reshape(AsU(i,j,:),[],1);
        Res.Coef.U.(['C' num2str(i) '_' num2str(j)])=AU;
        %  end
        
        %    for j=1:HasCons+lagsR*K
        AR=reshape(AsR(i,j,:),[],1);
        Res.Coef.R.(['C' num2str(i) '_' num2str(j)])=AR;
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',alpha);
        
        Res.Coef.tTest.Hypo(i,j)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Coef.tTest.prob(i,j)=p.';
        Res.Coef.tTest.stats(i,j)=stats.tstat.';
        Res.Coef.tTest.df(i,j)=stats.df.';
        Res.Coef.tTest.sd(i,j)=stats.sd.';
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.1);
        Res.Coef.tTest10.Hypo(i,j)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Coef.tTest10.prob(i,j)=p.';
        Res.Coef.tTest10.stats(i,j)=stats.tstat.';
        Res.Coef.tTest10.df(i,j)=stats.df.';
        Res.Coef.tTest10.sd(i,j)=stats.sd.';
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.05);
        Res.Coef.tTest05.Hypo(i,j)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Coef.tTest05.prob(i,j)=p.';
        Res.Coef.tTest05.stats(i,j)=stats.tstat.';
        Res.Coef.tTest05.df(i,j)=stats.df.';
        Res.Coef.tTest05.sd(i,j)=stats.sd.';
        
    end
    
    for j=1:K
        AU=reshape(OmgU(i,j,:),[],1);
        Res.Omg.U.(['C' num2str(i) '_' num2str(j)])=AU;
        AR=reshape(OmgR(i,j,:),[],1);
        Res.Omg.R.(['C' num2str(i) '_' num2str(j)])=AR;
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',alpha);
        Res.Omg.tTest.Hypo(i,j)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Omg.tTest.prob(i,j)=p.';
        Res.Omg.tTest.stats(i,j)=stats.tstat.';
        Res.Omg.tTest.df(i,j)=stats.df.';
        Res.Omg.tTest.sd(i,j)=stats.sd.';
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.1);
        Res.Omg.tTest10.Hypo(i,j)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Omg.tTest10.prob(i,j)=p.';
        Res.Omg.tTest10.stats(i,j)=stats.tstat.';
        Res.Omg.tTest10.df(i,j)=stats.df.';
        Res.Omg.tTest10.sd(i,j)=stats.sd.';
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.05);
        Res.Omg.tTest05.Hypo(i,j)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Omg.tTest05.prob(i,j)=p.';
        Res.Omg.tTest05.stats(i,j)=stats.tstat.';
        Res.Omg.tTest05.df(i,j)=stats.df.';
        Res.Omg.tTest05.sd(i,j)=stats.sd.';
    end
end

end
