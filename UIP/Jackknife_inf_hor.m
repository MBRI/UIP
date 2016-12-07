function [ahat]=Jackknife_inf_hor(B,StartingPoint,Err,Var_names,Exo_names)
% bootstrap infinite horizon UIP
% B is the estimated coef. matrix without dummy and constant
% Starting point to run the bootstrap
% Err original error term
% Tt the length of the time series we use in estimation
% Var_names endogenus variables name
FireCount=size(Err,1)+1; % the size of samples
BurnCount=size(Err,1)+1; % the number of generated sapmles
% alpha=0.05;% significant level
xx_Position=cellfun(@(x) ~isempty(x),regexpi(Var_names,strjoin(Exo_names,'|')));
endoVar=Var_names(~xx_Position);
EE_Position=cellfun(@(x) ~isempty(x),regexpi(endoVar,'EE[\w*]'));
dp_Position=cellfun(@(x) ~isempty(x),regexpi(endoVar,'dp[\w*]'));

% Step1. Get the estimation data and errors
XoX=StartingPoint(:,xx_Position);
if size(StartingPoint,1)<FireCount
    XoX(end+1,:)=XoX(end,:);
end
B0X=B(xx_Position,:);
% exo removed
B(xx_Position,:)=[];

  StartingPoint(:,xx_Position)=[];% Remove exogenius variables
  K=size(B,1);
% impose Restriction
% EE restriction
B(:,EE_Position)=zeros(1,K);

% 
% dp restriction
delta1=zeros(K,1);
delta1(strcmp(endoVar,'rer'),1)=1; %real exchange rate multipler
delta3=zeros(K,1);
delta3(cellfun(@(x) ~isempty(x),regexpi(endoVar,'dr[\w*]')),1)=1; % nOTE

B(:,dp_Position)=B(:,strcmp(endoVar,'rer'))+delta3-delta1;



% error term restriction
EC=eye(size(Err,2));
EC(EE_Position,EE_Position)=0;% zero sel error
EC(cellfun(@(x) ~isempty(x),regexpi(endoVar,'dp[\w*]')),EE_Position)=-1;% zero sel error
EC(cellfun(@(x) ~isempty(x),regexpi(endoVar,'rer')),EE_Position)=1;% zero sel error
% Constant Term restriction
CC=eye(size(Err,2));

CC(dp_Position,dp_Position)=0;% zero sel error
CC(cellfun(@(x) ~isempty(x),regexpi(endoVar,'EE[\w*]')),dp_Position)=-1;% zero sel error
CC(cellfun(@(x) ~isempty(x),regexpi(endoVar,'rer')),dp_Position)=1;% zero sel error

K0=size(B,2);
K=size(B,1);

X=nan(FireCount,K,BurnCount);
BX=nan(K,K0,BurnCount);
B2=nan(K0,BurnCount);% infinite horizon bootstraped data
B1=nan(K0,BurnCount);% EE1 eq. coef.
for b=1:BurnCount
    % Step2. Generate Psedu Sample: We exclude one error term in each turn

    % Gener Startinf point
%     if b>1
    X(1,:,b)=StartingPoint(1,:);
%     else
%         X(1,:,b)=StartingPoint(1,:)-Err(1,:);
%     end
    for f=2:FireCount
        if f~=b
            X(f,:,b)=X(f-1,:,b)*B+XoX(f-1,:)*B0X*CC+Err(f-1,:)*EC;
        else
            X(f,:,b)=X(f-1,:,b)*B+XoX(f-1,:)*B0X*CC;
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

% if K<100
%    warning(['inefficient sampling\\\ ' num2str(K) ' from ' num2str(BurnCount) ])
% end
% ahat2
ahat1=sum((repmat(B1(:,1),1,K-1)-B1(:,2:end)).^3,2)./sum((repmat(B1(:,1),1,K-1)-B1(:,2:end)).^2,2).^(3/2)/6;
ahat2=sum((repmat(B2(:,1),1,K-1)-B2(:,2:end)).^3,2)./sum((repmat(B2(:,1),1,K-1)-B2(:,2:end)).^2,2).^(3/2)/6;
ahat.a1=ahat1;
ahat.a2=ahat2;

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
BX=BX0(1:end-size(Exo,2),:);
B1=BX(:,EE_Position).';

% if any(any(isnan(BX)))
%     EigenValues=1;
% else
% EigenValues=abs(eig(BX));
% end
% if any(EigenValues>=0.96)
%     % warning(':: infinie  horizon Explosive ::');
%     K=size(BX,1);
%     BX=nan(K);
%     B2=nan(K,1);
%     B1=B2;
%     Er=nan(size(X));
% else
    B2=(B1/(eye(size(BX,2))-BX)).';
    Er=Y-X*BX0;
% end
end