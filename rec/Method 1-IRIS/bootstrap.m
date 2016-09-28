function Res=bootstrap(vU,vUd,vR,vRd,HasCons,constrA,startHist,endHist,alpha,Sims,NPer,Evar,InsPos,TrgPos)

rng(0);
if endHist-startHist<500
endHist=endHist+1000-length(startHist:endHist); 
end% this may Correct the Confidence Intervals
%% ******************** Unrestricted Model *******************************
pU = get(vU,'order'); %?order?
bdU = resample(vU,vUd,startHist+pU:endHist,Sims,'wild=',true); %?resample?

disp('Bootstrapped Unrestricted Model');


% Estimate VAR Models from Bootstrapped Data 


yList = get(vU,'yList');
K=length(yList);
vUs = VAR(yList); %?emptyVAR?

vUs = estimate(vUs,bdU,startHist:endHist,'order=',pU,'const=',HasCons); %?estimate?

inxU = isstationary(vUs); %?nonstIndex?

disp('Remove explosive parameterisations');
vUs = vUs(inxU); %?removeNonst?

% Store Original and Bootstrapped Transition Matrices
AU = get(vU,'A*');
AUS = get(vUs,'A*');
OmgU=get(vUs,'omg');
%% ******************** Restricted Model *******************************
pR = get(vR,'order'); %?order?
bdR = resample(vR,vRd,startHist+pR:endHist,Sims,'wild=',true); %?resample?

disp('Bootstrapped Restricted Model');


%% Estimate 500 VAR Models from Bootstrapped Data 

yList = get(vR,'yList');
vRs = VAR(yList); %?emptyVAR?

vRs = estimate(vRs,bdR,startHist:endHist,'order=',pR,'const=',HasCons,'A=',constrA); %?estimate?

inxR = isstationary(vRs); %?nonstIndex?

disp('Remove explosive parameterisations');
vRs = vRs(inxR);

%% Compare Original and Bootstrapped Transition Matrices
% 
AR = get(vR,'A*'); %?getA?
ARS = get(vRs,'A*');
OmgR=get(vRs,'omg');
%% Store Bootstrapped IRFs
[RespU,~] = srf(SVAR(vUs),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
[RespR,~] = srf(SVAR(vRs),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
%% 5. Confidence interval & Statistical Tests
clear tt
% 5-1. IRFs statistics
for tt=1:length(TrgPos)
    IU=squeeze(RespU.(Evar{TrgPos(tt)}).data);
    IRFU.Simul=IU;
    IRFU.UBound=quantile(IU,1-alpha/2,2);
    IRFU.LBound=quantile(IU,alpha/2,2);
    IRFU.Median=quantile(IU,1/2,2);
    IRFU.Mean=mean(IU,2);
    IRFU.std=std(IU,0,2);
    
    IR=squeeze(RespR.(Evar{TrgPos(tt)}).data);
    IRFR.Simul=IR;
    IRFR.UBound=quantile(IR,1-alpha/2,2);
    IRFR.LBound=quantile(IR,alpha/2,2);
    IRFR.Median=quantile(IR,1/2,2);
    IRFR.Mean=mean(IR,2);
    IRFR.std=std(IR,0,2);
    
    
%{ 
   DI=(IU-IR).^2;% abs Be or not to be
    Dif.Simul=DI;
    Dif.UBound=quantile(DI,1-alpha/2,2);%Uper bound
    Dif.LBound=quantile(DI,alpha/2,2);%lower bound
    Dif.Median=quantile(DI,1/2,2);%median
    Dif.Mean=mean(DI,2);
    Dif.std=std(DI,0,2);
   %} 
    % Chi2
    
    [RespUU,~] = srf(SVAR(vU),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
    [RespRR,~] = srf(SVAR(vR),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
    
    RespUU=RespUU.(Evar{TrgPos(tt)}).data;
    RespRR=RespRR.(Evar{TrgPos(tt)}).data;
    
    %T=size(IRFsU,1);
    %RespU=RespU-mean(IRFsU,2)./std(IRFsU,2);
    %RespR=RespR-mean(IRFsR,2)./std(IRFsR,2);
    %CI=(IU-IR);% abs Be or not to be
    Resp=((RespUU-IRFU.Mean)./IRFU.std)-((RespRR-IRFR.Mean)./IRFR.std);
    %CI=(CI-repmat(mean(CI,1),T,1))./repmat(std(CI,1),T,1);
    %Resp=Resp-mean(CI,2)./std(CI,0,2); % if the std ==0 then Resp=nan
    Test.CHI2.Simul=Resp;
    % CHI2.stat=sum(Resp.^2,2);
    Test.CHI2.stat=sum(Resp(~isnan(Resp)).^2);
    Test.CHI2.prob=chi2cdf(Test.CHI2.stat,NPer-sum(isnan(Resp))); % the Nan is not removed
    %}
    % Ttest
    % pairwise t test
    [h,p,~,stats]=ttest2(IU,IR,'Dim',2,'Alpha',alpha,'Vartype','unequal');
    h(isnan(h))=0;
    Test.t.Hypo=h; % one: the null hyphothesis is rejected, so the means are not equal
    Test.t.prob=p;
    Test.t.stats=stats.tstat;
    Test.t.df=stats.df;
    Test.t.sd=stats.sd;
    % Store in the output
    Res.(Evar{TrgPos(tt)}).IRFsU=IRFU;
    Res.(Evar{TrgPos(tt)}).IRFsR=IRFR;
%    Res.(['Trg' num2str(TrgPos(tt))]).Diff=Dif;
    Res.(Evar{TrgPos(tt)}).Test=Test;
    clear  IRFU IRFR Dif Test
end
Res.Coef.U=AUS;
Res.Coef.R=ARS;
Res.Omg.U=OmgU;
Res.Omg.R=OmgR;
% 5-2. Coef Statistics
for i=1:K % for equations
    for j=1:K % for endo var
        for l=1:pU % for lags
        AU=squeeze(AUS(i,j,l,:));
       % Res.Coef.U.(['C_eq' num2str(i) '_Var' num2str(j) '_lag' num2str(l)])=AU;
        %  end
        
        %    for j=1:HasCons+lagsR*K
        AR=squeeze(ARS(i,j,l,:));
        %Res.Coef.R.(['C_eq' num2str(i) '_Var' num2str(j) '_lag' num2str(l)])=AR;
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',alpha,'Vartype','unequal');
        
        Res.Coef.tTest.Hypo(i,j,l)=h; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Coef.tTest.prob(i,j,l)=p;
        Res.Coef.tTest.stats(i,j,l)=stats.tstat;
        Res.Coef.tTest.df(i,j,l)=stats.df;
       Res.Coef.tTest.sd(i,j,l,1:2)=stats.sd;
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.1,'Vartype','unequal');
        Res.Coef.tTest10.Hypo(i,j,l)=h; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Coef.tTest10.prob(i,j,l)=p;
        Res.Coef.tTest10.stats(i,j,l)=stats.tstat;
        Res.Coef.tTest10.df(i,j,l)=stats.df;
        %Res.Coef.tTest10.sd(i,j,l)=stats.sd;
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.05,'Vartype','unequal');
        Res.Coef.tTest05.Hypo(i,j,l)=h.'; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Coef.tTest05.prob(i,j,l)=p.';
        Res.Coef.tTest05.stats(i,j,l)=stats.tstat.';
        Res.Coef.tTest05.df(i,j,l)=stats.df.';
        %Res.Coef.tTest05.sd(i,j,l)=stats.sd.';
        end

        
        
        
        % Var-Cov equality test
        AU=squeeze(OmgU(i,j,:));
        %Res.Omg.U.(['C' num2str(i) '_' num2str(j)])=AU;
        AR=squeeze(OmgR(i,j,:));
        %Res.Omg.R.(['C' num2str(i) '_' num2str(j)])=AR;
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',alpha,'Vartype','unequal');
        Res.Omg.tTest.Hypo(i,j)=h; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Omg.tTest.prob(i,j)=p;
        Res.Omg.tTest.stats(i,j)=stats.tstat;
        Res.Omg.tTest.df(i,j)=stats.df;
        Res.Omg.tTest.sd(i,j,1:2)=stats.sd;
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.1);
        Res.Omg.tTest10.Hypo(i,j)=h; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Omg.tTest10.prob(i,j)=p;
        Res.Omg.tTest10.stats(i,j)=stats.tstat;
        Res.Omg.tTest10.df(i,j)=stats.df;
        Res.Omg.tTest10.sd(i,j,1:2)=stats.sd;
        %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        [h,p,~,stats]=ttest2(AU,AR,'Alpha',0.05);
        Res.Omg.tTest05.Hypo(i,j)=h; % one: the null hyphothesis is rejected, so the means are not equal
        Res.Omg.tTest05.prob(i,j)=p;
        Res.Omg.tTest05.stats(i,j)=stats.tstat;
        Res.Omg.tTest05.df(i,j)=stats.df;
        Res.Omg.tTest05.sd(i,j,1:2)=stats.sd;
    end
end
%}
disp('Test''s Accomplished');
end

