function [BX,B2,CI,X]=BootStrap_inf_hor(B,StartingPoint,Err,Tt,Var_names,ahat)
% bootstrap infinite horizon UIP
% B is the estimated coef. matrix without dummy and constant
% Starting point to run the bootstrap
% Err original error term
% Tt the length of the time series we use in estimation
% Var_names endogenus variables name
% ahat is jacknife acceletion index contian a1 and a2
FireCount=500+Tt; % the size of samples
BurnCount=3000; % the number of generated sapmles
alpha=0.05;% significant level
EE_Position=strcmp(Var_names,'EE1');
% Step1. Get the estimation data and errors
Borig=B;
K=size(B,1);
B(EE_Position,:)=zeros(1,K);% impose Restriction
warning('off')

X=nan(FireCount,K,BurnCount);
BX=nan(K,K,BurnCount);
B2=nan(K,BurnCount);% infinite horizon bootstraped data
B1=nan(K,BurnCount);% EE1 eq. coef.
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
[BX(:,:,b),B2(:,b),B1(:,b)]=est(X(end-Tt-1:end,:,b),EE_Position);%(,:,b)
% Step4. Store the results
% if abs(B2(1,b))>4
%    disp('me'); 
% end
end
warning('on')
% Step5. Descript the results
BX(:,:,isnan(B2(1,:)))=[];
B1(:,isnan(B2(1,:)))=[];
B2(:,isnan(B2(1,:)))=[];

Prc1=prctile(B1,100*[alpha,0.5,1-alpha],2);
Prc2=prctile(B2,100*[alpha,0.5,1-alpha],2);
% Prc.Prc1=Prc1;
% Prc.Prc2=Prc2;
% BCa
B1O=B(EE_Position,:);% original B1
B2O=(B1O/(eye(size(B,2))-B)).';% original B2
B1O=B1O.';
Tk=size(B2,2);
if Tk<100
   warning(['inefficient sampling\\\ ' num2str(Tk) ' from ' num2str(BurnCount) ])
end
z0HatB2=norminv(sum(B2<repmat(B2O,1,Tk),2)/Tk);

a1B2=normcdf(z0HatB2+(z0HatB2+norminv(alpha))./(1-ahat.a2.*(z0HatB2+norminv(alpha))));
a2B2=normcdf(z0HatB2+(z0HatB2+norminv(0.5))./(1-ahat.a2.*(z0HatB2+norminv(0.5))));
a3B2=normcdf(z0HatB2+(z0HatB2+norminv(1-alpha))./(1-ahat.a2.*(z0HatB2+norminv(1-alpha))));
%BCa2=prctile(B2,[100*a1,100*a2,100*a3],2);
Tk=size(B1,2);
z0HatB1=norminv(sum(B1<repmat(B1O,1,Tk),2)/Tk);
a1B1=normcdf(z0HatB1+(z0HatB1+norminv(alpha))./(1-ahat.a1.*(z0HatB1+norminv(alpha))));
a2B1=normcdf(z0HatB1+(z0HatB1+norminv(0.5))./(1-ahat.a1.*(z0HatB1+norminv(0.5))));
a3B1=normcdf(z0HatB1+(z0HatB1+norminv(1-alpha))./(1-ahat.a1.*(z0HatB1+norminv(1-alpha))));
%----
BCa1=nan(K,3);
BCa2=BCa1;
for k=1:K
%     figure;
%     ksdensity(B2(k,:))
BCa1(k,1:3)=prctile(B1(k,:),[100*a1B1(k),100*a2B1(k),100*a3B1(k)],2);
BCa2(k,1:3)=prctile(B2(k,:),[100*a1B2(k),100*a2B2(k),100*a3B2(k)],2);
% [ci,bootstat]=bootci(100,{@mean,B2(:,k)},'type','bca');
% k
% ci
end
% Prc.BCa1=BCa1;
% Prc.BCa2=BCa2;
B1O=Borig(EE_Position,:);% original B1
B2O=(B1O/(eye(size(Borig,2))-Borig)).';% original B2
B1O=B1O.';
% BCa Stored
B2CI=mat2dataset([BCa2, B2O],'varnames',{'lower' 'Mod' 'Upper' 'Observed'},'obsnames',Var_names);
B1CI=mat2dataset([BCa1, B1O],'varnames',{'lower' 'Mod' 'Upper' 'Observed'},'obsnames',Var_names);
B2CI.Sig=B2CI.lower<=B2O & B2CI.Upper>=B2O;
B1CI.Sig=B1CI.lower<=B1O & B1CI.Upper>=B1O;
CI.BCa.B1=B1CI;
CI.BCa.B2=B2CI;
% percentile Stored
B2CI=mat2dataset([Prc2, B2O],'varnames',{'lower' 'Mod' 'Upper' 'Observed'},'obsnames',Var_names);
B1CI=mat2dataset([Prc1, B1O],'varnames',{'lower' 'Mod' 'Upper' 'Observed'},'obsnames',Var_names);
B2CI.Sig=B2CI.lower<=B2O & B2CI.Upper>=B2O;
B1CI.Sig=B1CI.lower<=B1O & B1CI.Upper>=B1O;

CI.Prc.B1=B1CI;
CI.Prc.B2=B2CI;
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
% if any(any(isnan(BX)))
%     EigenValues=1;
% else
% EigenValues=abs(eig(BX));
% end
% if any(EigenValues>=0.96)
%    % warning(':: infinie  horizon Explosive ::');
%     K=size(BX,1);
%     BX=nan(K);
%     B2=nan(K,1);
%     B1=B2;
%     Er=nan(size(X));
% else
B2=(B1/(eye(size(BX,2))-BX)).';
Er=Y-X*BX;
% end
end