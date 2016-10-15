function [ahat]=Jackknife_inf_hor(B,StartingPoint,Err,Tt,endo_names,Exo_names)
% bootstrap infinite horizon UIP
% B is the estimated coef. matrix without dummy and constant
% Starting point to run the bootstrap
% Err original error term
% Tt the length of the time series we use in estimation
% Var_names endogenus variables name
FireCount=size(Err,1)+1; % the size of samples
BurnCount=size(Err,1)+1; % the number of generated sapmles
alpha=0.05;% significant level
EE_Position=strcmp(endo_names,'EE1');
xx_Position=cmpr(endo_names,Exo_names);
% Step1. Get the estimation data and errors
XoX=StartingPoint(:,xx_Position);
B0X=B(xx_Position,:);
  StartingPoint(:,xx_Position)=[];% Remove exogenius variables
B(EE_Position,:)=0;%zeros(1,K);% impose Restriction
B(xx_Position,:)=[];

K0=size(B,2);
K=size(B,1);

X=nan(FireCount,K,BurnCount);
BX=nan(K,K0,BurnCount);
B2=nan(K0,BurnCount);% infinite horizon bootstraped data
B1=nan(K0,BurnCount);% EE1 eq. coef.
for b=1:BurnCount
    % Step2. Generate Psedu Sample: We exclude one error term in each turn

    % Gener Startinf point
    X(1,:,b)=StartingPoint(1,:);
  
    for f=2:FireCount
        if f~=b
            X(f,:,b)=X(f-1,:,b)*B+XoX(f-1,:)*B0X+Err(f-1,:);
        else
            X(f,:,b)=X(f-1,:,b)*B+XoX(f-1,:)*B0X;
        end
        
    end
    % Step3. Estimate model with the last t observation
    [BX(:,:,b),B2(:,b),B1(:,b)]=est(X(:,:,b),XoX,EE_Position);%(,:,b)
    % Step4. Store the results
    % if abs(B2(1,b))>4
    %    disp('me');
    % end
end
% Step5. Descript the results


BX(:,:,isnan(B2(1,:)))=[];
B1(:,isnan(B2(1,:)))=[];
B2(:,isnan(B2(1,:)))=[];
K=size(B1,2);
% ahat2
ahat1=sum((repmat(B1(:,1),1,K-1)-B1(:,2:end)).^3,2)./sum((repmat(B1(:,1),1,K-1)-B1(:,2:end)).^2,2).^(3/2)/6;
ahat2=sum((repmat(B2(:,1),1,K-1)-B2(:,2:end)).^3,2)./sum((repmat(B2(:,1),1,K-1)-B2(:,2:end)).^2,2).^(3/2)/6;
ahat.a1=ahat1;
ahat.a2=ahat2;
%{
Prc1=prctile(B1,100*[alpha,0.5,1-alpha],2);
Prc2=prctile(B2,100*[alpha,0.5,1-alpha],2);
Prc.Prc1=Prc1;
Prc.Prc2=Prc2;
% BCa
B1O=B(EE_Position,:);% original B1
B2O=(B1O/(eye(size(B,2))-B)).';% original B2
B1O=B1O.';
Tk=size(B2,2);
if Tk<100
    warning(['inefficient sampling\\\ ' num2str(Tk) ' from ' num2str(BurnCount) ])
end
z0HatB2=norminv(sum(B2<repmat(B2O,1,Tk),2)/Tk);
ahat=0;
a1B2=normcdf(z0HatB2+(z0HatB2+norminv(alpha))./(1-ahat.*(z0HatB2+norminv(alpha))));
a2B2=normcdf(z0HatB2+(z0HatB2+norminv(0.5))./(1-ahat.*(z0HatB2+norminv(0.5))));
a3B2=normcdf(z0HatB2+(z0HatB2+norminv(1-alpha))./(1-ahat.*(z0HatB2+norminv(1-alpha))));
%BCa2=prctile(B2,[100*a1,100*a2,100*a3],2);
Tk=size(B1,2);
z0HatB1=norminv(sum(B1<repmat(B1O,1,Tk),2)/Tk);
a1B1=normcdf(z0HatB1+(z0HatB1+norminv(alpha))./(1-ahat.*(z0HatB1+norminv(alpha))));
a2B1=normcdf(z0HatB1+(z0HatB1+norminv(0.5))./(1-ahat.*(z0HatB1+norminv(0.5))));
a3B1=normcdf(z0HatB1+(z0HatB1+norminv(1-alpha))./(1-ahat.*(z0HatB1+norminv(1-alpha))));
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
Prc.BCa1=BCa1;
Prc.BCa2=BCa2;
%}
end
function [BX,B2,B1,Er]=est(X,Exo,EE_Position)
% X endo genius time series
% Exo genuis Time series
% EE_
%% infinie  horizon
Y=lagmatrix(X,-1); % remove constant and dummy

X=[X,Exo];
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
BX0=(X.'*X)\X.'*Y;
B1=BX0(EE_Position,:);
BX=BX0(1:end-size(Exo,2),:);
EigenValues=abs(eig(BX));
if any(EigenValues>=0.96)
    % warning(':: infinie  horizon Explosive ::');
    K=size(BX,1);
    BX=nan(K);
    B2=nan(K,1);
    B1=B2;
    Er=nan(size(X));
else
    B2=(B1/(eye(size(BX,2))-BX)).';
    Er=Y-X*BX0;
end
end