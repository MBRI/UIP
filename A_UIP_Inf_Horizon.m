%load Data
%Data=dataset('xls','DataBaseu.xlsx','sheet','Sheet3');
% estimate ols: EE? ~ rer EE1(-1) dr0 dpi NFA_MB Y_Y_  Openness tot dum
%Ind_var={'Con' 'rer' 'EE0' 'dr0' 'dpi0' 'NFA_MB' 'Y_YU'  'Openness' 'tot' 'dum' 'Rs'}; % independent Variables
%Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'dum' 'Con' 'NFA_MB' 'Openness' 'tot' 'Rs' 'Rh'}; % independent Variables 
Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'Con' 'NFA_MB'  'Rs' 'tot' 'Y_YU' 'Rh'}; % it is good
% Ind_var={'rer' 'EE1'  'Con' 'dr0' 'dpi0'}; % 

dum_var=cmpr(Ind_var,{'dum' 'Con' 'dum_time'});

B=nan(length(Ind_var),40);
X=double(Data(:,Ind_var));
% dumi and constant problem
X(:,~dum_var)=SA(X(:,~dum_var),4);
% X(1:4,:)=[]; % remove from the first Data.Date>1992 ? <2006
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

BX0=(X2.'*X2)\X2.'*Y;
BX=BX0(~dum_var2,:).';
B1=BX(strcmp(Ind_var2,'EE1'),:);%B(~dum_var,1).';
EigenValues=abs(eig(BX));
if any(EigenValues>=1)
    warning(':: infinie  horizon Explosive ::');
end
B2=(B1/(eye(size(BX,2))-BX)).';
endo_var=Ind_var2(~dum_var2).';
B1=B1.';
X3=X(~any(isnan(X),2),:);
ahat=Jackknife_inf_hor(BX0,X3,Y-X2*BX0,size(Y,1),Ind_var2,{'dum' 'Con' 'dum_time'});
[bsBX,bsB2,Prc]=BootStrap_inf_hor(BX,zeros(1,size(Y,2)),Y-X2*((X2.'*X2)\(X2.'*Y)),size(Y,1),endo_var,ahat);%mean(Y)
B2CI=mat2dataset([Prc.BCa2, B2],'varnames',{'lower' 'Mod' 'Upper' 'Observed'},'obsnames',endo_var);
B1CI=mat2dataset([Prc.BCa1, B1],'varnames',{'lower' 'Mod' 'Upper' 'Observed'},'obsnames',endo_var);
B2CI.Sig=B2CI.lower<=B2 & B2CI.Upper>=B2;
B1CI.Sig=B1CI.lower<=B1 & B1CI.Upper>=B1;
%'B2'
B2CI
%'B1'
B1CI