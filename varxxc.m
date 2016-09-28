%load Data
Data=dataset('xls','DataBaseu.xlsx','sheet','Sheet2');
% estimate ols: EE? ~ rer EE1(-1) dr0 dpi NFA_MB Y_Y_  Openness tot dum
%Ind_var={'Con' 'rer' 'EE0' 'dr0' 'dpi0' 'NFA_MB' 'Y_Y_'  'Openness' 'tot' 'dum'}; % independent Variables
Ind_var={ 'rer' 'EE0' 'dr0' 'dpi0' 'dum' 'Con'}; % independent Variables

B=nan(length(Ind_var),40);
for i=1:40
Dep_var=['EE' num2str(i)]; % dependent Variables
X=double(Data(:,Ind_var));
%X=[ones(size(X,1),1), X];
Y=double(Data(:,Dep_var));
inan=any(isnan(X),2) | any(isnan(Y),2);
X(inan,:)=[];
Y(inan,:)=[];
B(:,i)=(X.'*X)\X.'*Y;
end
%% Set model 
Spec = vgxset();


% 'Series'

% Model Orders
% 
% Name	Value
% n	:  A positive integer specifying the number of time series. The default is 1.
% nAR : A nonnegative integer specifying the number of AR lags. The default is 0.
% nX : A nonnegative integer specifying the number regression parameters. The default is 0. 

% Model Parameters
% Constant
Spec = vgxset('Series',{'Y','X'},'n',2,'nAR',1,'nX',0,'Constant',true);
%% Estimate VARX model parameters
Y=rand(100,2);
[EstSpecX, EstStdErrorsX] = vgxvarx(Spec,Y);

%Display Results
vgxdisp(EstSpec);