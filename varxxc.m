%load Data
%Data=dataset('xls','DataBaseu.xlsx','sheet','Sheet3');
% estimate ols: EE? ~ rer EE1(-1) dr0 dpi NFA_MB Y_Y_  Openness tot dum
%Ind_var={'Con' 'rer' 'EE0' 'dr0' 'dpi0' 'NFA_MB' 'Y_YU'  'Openness' 'tot' 'dum' 'Rs'}; % independent Variables
%Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'dum' 'Con' 'NFA_MB' 'Openness' 'tot' 'Rs' 'Rh'}; % independent Variables 
Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'Con' 'NFA_MB'  'Rs' 'tot' 'Y_YU' 'Rh'}; % it is good
%Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'Con'}; % 

dum_var=cmpr(Ind_var,{'dum' 'Con' 'dum_time'});

B=nan(length(Ind_var),40);
X=double(Data(:,Ind_var));
% dumi and constant problem
X(:,~dum_var)=SA(X(:,~dum_var),4);
X(1:4,:)=[]; % remove from the first Data.Date>1992 ? <2006

for i=1:40
    X2=X;
    Dep_var={['EE' num2str(i)] }; % dependent Variables
    %X=[ones(size(X,1),1), X];
    Y=SA(double(Data(:,Dep_var)),4);
  
    Y=lagmatrix(Y,-i);
    %'******************
    Y(1:4,:)=[]; % remove from the first
    %'******************
    inan=any(isnan(X),2) | any(isnan(Y),2);
    X2(inan,:)=[];
    Y(inan,:)=[];
    % remove invariant variables except constant
    indd=var(X2)==0;
    indd(strcmp(Ind_var,'Con'))=0;
    X2(:,indd)=[];
    if isempty(X2)
        continue;
    end
    B(~indd,i)=(X2.'*X2)\X2.'*Y;
end
%% infinie  horizon
Y=lagmatrix(X(:,~dum_var),-1); % remove constant and dummy
X2=X;
inan=any(isnan(X2),2) | any(isnan(Y),2);
X2(inan,:)=[];
Y(inan,:)=[];
% remove invariant variables except constant
indd=var(X2)==0;
indd(strcmp(Ind_var,'Con'))=0;
X2(:,indd)=[];
Ind_var2=Ind_var;Ind_var2(:,indd)=[];dum_var2=dum_var();dum_var2(:,indd)=[];

BX=(X2.'*X2)\X2.'*Y;
BX=BX(~dum_var2,:).';
B1=BX(strcmp(Ind_var2,'EE1'),:)%B(~dum_var,1).';
EigenValues=abs(eig(BX));
if any(EigenValues>=1)
    warning(':: infinie  horizon Explosive ::');
end
B2=(B1/(eye(size(BX,2))-BX)).'
endo_var=Ind_var2(~dum_var2).'

[bsBX,bsB2,Prc]=inf_hor(BX,mean(Y),Y-X2*((X2.'*X2)\(X2.'*Y)),size(Y,1),endo_var);
Prc
 return
%% Set model
% Spec = vgxset();


% 'Series'

% Model Orders
%
% Name	Value
% n	:  A positive integer specifying the number of time series. The default is 1.
% nAR : A nonnegative integer specifying the number of AR lags. The default is 0.
% nX : A nonnegative integer specifying the number regression parameters. The default is 0.

% Model Parameters
% Constant
Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'NFA_MB'  'Rs' 'tot' 'Y_YU' 'Rh'}; % 
Spec = vgxset('Series',Ind_var,'nAR',1,'nX',0,'Constant',true);%,'ARsolve', repmat({ logical(eye(2)) }, 2, 1)));
% Estimate VARX model parameters
%Y=rand(100,2);

Y=SA(double(Data(:,Ind_var)),4);
%X0=SA(double(Data(:,'tot')),4);
[EstSpecX, EstStdErrorsX] = vgxvarx(Spec,Y,X0);

%Display Results
vgxdisp(EstSpecX);