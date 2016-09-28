%% ------------------------ Versoin 1 -------------------------------------
% Author:  Pedram Davoudi
% email: p.davoudi@mbri.ac.ir
% This script run the montecarlo simulation for irf asymptotic distribution
% All the estimations take place with the IRIS (iristoolbox.codeplex.com)
%% --------------------  To use the code  ---------------------------------
% you must enter your data in dt.csv available in the current folder.
% be aware of IRIs format for CSV files. you could find it in IRIS manual.
% uncomment line's 14 and 15 ******************************************
% this two line must be evaluated just once and comment them again for the
% next run's till matlab terminated.
% change the IRIS Pack Address which is copied in your Hard drive.
%% ******************** IRIS Pack Address & Startup ************************
%addpath('F:\Work\MBRI\CEPREMAP2015Jan - Original\IRIS_Tbx_20141219') %
%irisstartup  %IRIS Package path
%% ******************* Change Model Specification **************************
% Specify your model in the section 1
% Run the Code
% if there were any bug, don't hesitate to call me.
%
%% 1. Model Specification
tic
clear
clc
close all
% Change the default line styles
set(groot,'DefaultAxesColorOrder',[0.1 0.1 0.1], 'DefaultAxesLineStyleOrder','- | -- | : | -. | . | o | * | + | x | s | d | ^ | v | > | < | p | h'); % 'DefaultAxesColorOrder',[0 0 0]
set(groot, 'defaultTextInterpreter','latex');
% Use Same Database
SameDatabase=true;
% number of lags
pLags=4;
% Has Constant Term
HasCons=true; %true or fslse
% Degree
alpha=0.01;
% Number of Simulation and Horizo of impulse response
Sims=1000;
NPer=100;
% variable in estimation (all endogenuis)
% my var: lmb	lngdp	lp	fcn	lfc	lg	lm	lnda	lnfa	lryo	gdp	pi	lel	lgdp	lol	lpi	lva
Ins={'lnda'}; % don't confuse me, Just one var
Chn={'lol'}; % don't confuse me, Just one var
Oth={}; % Other Variables, {} means nothing
TrG={'lngdp' , 'lp'}; %

if isempty(Oth)>0
    Evar=[TrG Chn Ins];
else
    Evar=[TrG Oth Chn Ins];
end
% Restrictions
[~,InsPos]=find(strcmp(Evar,Ins)); % Find the Instrument Position
[~,ChnPos]=find(strcmp(Evar,Chn));% Find the Chanel Position

[~,TrgPos]=find(strcmp(Evar,TrG{1}));% Find the Terget Position
for i=2:length(TrG)
    [~,a]=find(strcmp(Evar,TrG{i}));% Find the Terget Position
    TrgPos=[TrgPos;a];
end
%ChnPos=strfind(Evar,Chn{1});
% Restrict channel cof to zero
constrA = nan(length(Evar),length(Evar),pLags);
constrA(ChnPos,InsPos,:) = 0;

%Load DataBase
Dt=dbload('dt.csv','dateFormat','YYYY-MM-01','freq',4);%

startHist = get(Dt.(Evar{1}),'start');
endHist = get(Dt.(Evar{1}),'end');


% 2. Model Estimation
%Estimated Unrestricted VAR'
v = VAR(Evar);
[vU,~] = estimate(v,Dt,startHist:endHist, ...
    'order=',pLags,'const=',HasCons, ...
    'covParameters=',true);

% % Get the estimated coefficients in the transition matrix (which is a lag
% % polynomial) <?A?> and the constant vector <?K?>.
%
% AU = get(vU,'A*'); %?A%?
% KU = get(vU,'K'); %?K?
% OmgU = get(vU,'Omega');

%Estimated Restricted VAR
%
[vR,~] = estimate(v,Dt,startHist:endHist, ...
    'order=',pLags,'const=',HasCons, ...
    'covParameters=',true,'A=',constrA);
% AR = get(vR,'A*'); %?A%?
% KR = get(vR,'K'); %?K?
% OmgR = get(vR,'Omega');
%}
%S=SVAR(vU);
%[Resp,~] = srf(S,1:20,'select=','res_X');
%Resp.(Evar{3});%{:,1}
%% Bootstraping
Res=MonteCarlo5(Dt,vU,vR,SameDatabase,alpha,Sims,NPer,HasCons,constrA,InsPos,TrgPos);
% Find
%% IRFs plot

close all
% Plot IRFs
[RespU,~] = srf(SVAR(vU),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator
[RespR,~] = srf(SVAR(vR),1:NPer,'select=',['res_' Evar{InsPos}]); % 1 must remove with the instruent indicator

for i=1:length(TrgPos)
    figure();
    
    subplot(2,2,1)
    hold on
    line([0 NPer],[0 0],'LineWidth',1,'Color',[.2 0 .1]);
    h1=plot(RespU.(Evar{TrgPos(i)}).data);
    h2=plot(Res.(['Trg' num2str(i)]).IRFsU.UBound);
    h3=plot(Res.(['Trg' num2str(i)]).IRFsU.LBound);
    h4=plot(Res.(['Trg' num2str(i)]).IRFsU.Mean);
    h5=plot(Res.(['Trg' num2str(i)]).IRFsU.Median);
    legend([h1 h2 h3 h4 h5],{'IRF', 'UpperBound', 'LowerBound', 'Mean', 'Median'});
    title(['Response of ' TrG{i} ' to impulse of ' Evar{InsPos} ' in Unrestricted Model']);
    %  plot(IRFsU.Mode)
    hold off
    
    
    subplot(2,2,2)
    hold on
    line([0 NPer],[0 0],'LineWidth',1,'Color',[.2 0 .1]);
    h1=plot(RespR.(Evar{TrgPos(i)}).data);
    h2=plot(Res.(['Trg' num2str(i)]).IRFsR.UBound);
    h3=plot(Res.(['Trg' num2str(i)]).IRFsR.LBound);
    h4=plot(Res.(['Trg' num2str(i)]).IRFsR.Mean);
    h5=plot(Res.(['Trg' num2str(i)]).IRFsR.Median);
    % plot(IRFsR.Mode)
    legend([h1 h2 h3 h4 h5],{'IRFs', 'UpperBound', 'LowerBound', 'Mean', 'Median'});
    title(['Response of ' TrG{i} ' to impulse of ' Evar{InsPos} ' in Restricted Model']);
    hold off
    
    %     subplot(2,2,3)
    %     hold on
    %     line([0 NPer],[0 0],'LineWidth',1,'Color',[.2 0 .1]);
    %      h1=plot((RespU.(Evar{TrgPos(i)}).data-RespR.(Evar{TrgPos(i)}).data).^2);
    %      h2=plot(Res.(['Trg' num2str(i)]).Diff.UBound);
    %      h3=plot(Res.(['Trg' num2str(i)]).Diff.LBound);
    %      h4=plot(Res.(['Trg' num2str(i)]).Diff.Mean);
    %      h5=plot(Res.(['Trg' num2str(i)]).Diff.Median);
    %     % plot(CI.Mode)
    %     legend([h1 h2 h3 h4 h5],{'IRFs', 'UpperBound', 'LowerBound', 'Mean', 'Median'});
    %     title(['The square of the difference in IRFs of ' TrG{i}]);
    %     %axis([mX MX 0 NPer]);
    %     hold off
    %
    subplot(2,2,4)
    hold on
    
    
    % *********-------------- plot pacthes
    X1=Res.(['Trg' num2str(i)]).IRFsU.Mean;
    X2=Res.(['Trg' num2str(i)]).IRFsR.Mean;
    X=[X1;X2;abs(X2-X1)];
    mX=min(X);
    MX=max(X);
    ss=Res.(['Trg' num2str(i)]).Test.t.Hypo;
    ss=ss-lagmatrix(ss,1);
    
    Pk=find(ss==1);
    tr=find(ss==-1);
    a=size(Pk,1)-size(tr,1);
    if a==1
        tr(end+1,1)=size(ss,1)+1;
    elseif a==-1
        Pk(end+1,1)=2;
        Pk=sort(Pk);
    end
    z=[Pk,tr-1];
    mm=repmat([mX, MX , MX,mX],size(z,1),1);
    zz=[z(:,1),z(:,1),z(:,2),z(:,2)];
    grbkgrnd = [.8 .8 .8];
    h = patch(zz.',mm.',grbkgrnd);
    set(h,'linestyle','none')
    
    %'-------------- IRFs
    line([0 NPer],[0 0],'LineWidth',1,'Color',[.8 0 .8]);
    h1=plot(X1);%Res.(['Trg' num2str(i)]).IRFsU.Mean);%RespU.(Evar{TrgPos(i)}).data);
    h2=plot(X2);%Res.(['Trg' num2str(i)]).IRFsR.Mean);%RespR.(Evar{TrgPos(i)}).data);
    h3=plot(abs(X1-X2));%Res.(['Trg' num2str(i)]).IRFsU.Mean-Res.(['Trg' num2str(i)]).IRFsR.Mean));
    
    legend([h1 h2 h3],{'Restricted\_Mean', 'Unrestricted\_Mean','Mean Dif.'});
    title(['Comparison of IRFs of ' TrG{i}]);
    %axis([mX MX 0 NPer]);
    hold off
    
end
toc
clear a constrA endHist alpha HasCons h h1 h2 h3 h4 h5 i Ins Chn ChnPos Dt grbkgrnd mm mX MX Oth Pk pLags RespR RespU Sims ss startHist tr X z zz v SameDatabase
set(groot,'DefaultAxesColorOrder','default', 'DefaultAxesLineStyleOrder','default');