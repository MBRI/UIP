function [BX,B2,Prc]=inf_hor(B,StartingPoint,Err,Tt,Var_names)
% bootstrap infinite horizon UIP
% B is the estimated coef. matrix without dummy and constant
% Starting point to run the bootstrap
% Err original error term
% Tt the length of the time series we use in estimation
% Var_names endogenus variables name
FireCount=500+Tt; % the size of samples
BurnCount=100; % the number of generated sapmles
EE_Position=strcmp(Var_names,'EE1');
% Step1. Get the estimation data and errors

K=size(B,1);
X=nan(FireCount,K,BurnCount);
BX=nan(K,K,BurnCount);
B2=nan(K,BurnCount);
for b=1:BurnCount;
% Step2. Generate Psedu Sample: We draw the error terms with equal
% probability from the vector of error terms estimated VAR

% Gener Randome seed to select Error 
ee =fix(random('uni',1,Tt,FireCount,1));
% Gener Startinf point 
X(1,:,b)=StartingPoint+Err(ee(1),:);
for f=2:FireCount
     X(f,:,b)=X(f-1,:,b)*B+Err(ee(f),:);
end
% Step3. Estimate model with the last t observation
[BX(:,:,b),B2(:,b)]=est(X(end-Tt-1:end,:,b),EE_Position);
% Step4. Store the results

end
% Step5. Descript the results
BX(:,:,isnan(B2(1,:)))=[];
B2(:,isnan(B2(1,:)))=[];

Prc=prctile(B2,[5,50,95],2);
for k=1:K
    figure;
    ksdensity(B2(k,:))
end
end
function [BX,B2,B1,Er]=est(X,EE_Position)
%% infinie  horizon
Y=lagmatrix(X,-1); % remove constant and dummy
inan=any(isnan(X),2) | any(isnan(Y),2);
X(inan,:)=[];
Y(inan,:)=[];
% remove invariant variables except constant
% indd=var(X)==0;
% 
% X(:,indd)=[];
% Y(:,indd)=[];
%Ind_var2=Ind_var;Ind_var2(:,indd)=[];dum_var2=dum_var();dum_var2(:,indd)=[];

% BX(~indd,~indd)=(X.'*X)\X.'*Y;
BX=(X.'*X)\X.'*Y;
B1=BX(EE_Position,:);
EigenValues=abs(eig(BX));
if any(EigenValues>=1)
   % warning(':: infinie  horizon Explosive ::');
    K=size(BX,1);
    BX=nan(K);
    B2=nan(K,1);
    B11=B2;
    Er=nan(size(X));
else
B2=(B1/(eye(size(BX,2))-BX)).';
Er=Y-X*BX;
end
end