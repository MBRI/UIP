%load Data
%Data=dataset('xls','DataBaseu.xlsx','sheet','Sheet3');
% estimate ols: EE? ~ rer EE1(-1) dr0 dpi NFA_MB Y_YU  Openness tot dum
%Ind_var={'Con' 'rer' 'drer' 'drt0' 'EE0' 'dr0' 'dpi0' 'NFA_MB' 'Y_YU'  'Openness' 'tot' 'dum' 'Rs'}; % independent Variables
% Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'dum' 'Con'  'Openness' 'tot' 'Rs' 'Rh'}; % independent Variables
% Ind_var={'rer' 'EE1' 'dr0' 'dpi0' 'Con' 'dum2'  }; % it is good
load('data5.mat')
for j=0:1
    Ind_var={ 'rer' ['EE' num2str(j)]  'Con' ['dr' num2str(j)] ['dp' num2str(j)]  'dum1' 'dum2' 'dum3' 'dum4' 'NFA_MB' ['rs' num2str(j)] ['rh' num2str(j)]}; %'NFA_MB'  'tot' 'Rs' 'Rh'
%     Ind_var={ 'rer' ['EE' num2str(j)]  'Con' ['dr' num2str(j)] ['dp' num2str(j)]  'dum1' 'dum2' 'dum3' 'dum4'}; %'NFA_MB'  'tot' 'Rs' 'Rh'
    % find dummies position
    % pay attention all dummies start with dum
    % the constant term must be con
    dum_var=cellfun(@(x) ~isempty(x),regexpi(Ind_var,'dum[\w*]|con'));
    regimVar=Data.dum1;
    Regims=unique(regimVar);
    
    %B=nan(length(Ind_var),40);
    X=double(Data(:,Ind_var));
    % dumi and constant problem
    %X(:,~dum_var)=SA(X(:,~dum_var),4);
    % X(1:4,:)=[]; % remove from the first Data.Date>1992 ? <2006
    %% infinie  horizon
    for i=0:length(Regims)
        if i==0
            X2=X;%(regimVar==Regims(i),:);
        else
            X2=X(regimVar==Regims(i),:);
        end
        Y=lagmatrix(X(:,~dum_var),-1); % remove constant and dummy
        if i~=0
            Y=Y(regimVar==Regims(i),:);
        end
        inan=any(isnan(X2),2) | any(isnan(Y),2);
        X2(inan,:)=[];
        Y(inan,:)=[];
        % remove invariant variables except constant
        indd=var(X2)==0;
        indd(strcmp(Ind_var,'Con'))=0;
        X2(:,indd)=[];
        Ind_var2=Ind_var;Ind_var2(:,indd)=[];dum_var2=dum_var;dum_var2(:,indd)=[];
        
        BX0=(X2.'*X2)\X2.'*Y;
        BX=BX0(~dum_var2,:).';
        % B1=BX(strcmp(Ind_var2,['EE' numstr(j)]),:);%B(~dum_var,1).';
        EigenValues=abs(eig(BX));
        if any(EigenValues>=0.99)
            warning([':: infinie  horizon, Found Explosion in Regime ' num2str(i) '  ::']);
        end
        % B2=(B1/(eye(size(BX,2))-BX)).';
        endo_var=Ind_var2(~dum_var2).';
        % B1=B1.';
        
        if i==0
            X3=X(~any(isnan(X),2),~indd);
        else
            X3=X(~any(isnan(X),2) & regimVar==Regims(i),~indd);
        end
        ahat=Jackknife_inf_hor(BX0,X3,Y-X2*BX0,Ind_var2,Ind_var2(dum_var2));
        [bsBX,bsB2,CI]=BootStrap_inf_hor(BX.',zeros(1,size(Y,2)),Y-X2*((X2.'*X2)\(X2.'*Y)),size(Y,1),endo_var,ahat);%mean(Y)
        
        Res.(['H' num2str(j) '_R' num2str(i)])=struct('CI',CI,'BX',bsBX,'B2',bsB2);
        if j>0
            Res.(['H' num2str(j) '_R' num2str(i)]).CI.BCa=rmfield(Res.(['H' num2str(j) '_R' num2str(i)]).CI.BCa,'B2');
            Res.(['H' num2str(j) '_R' num2str(i)]).CI.Prc=rmfield(Res.(['H' num2str(j) '_R' num2str(i)]).CI.Prc,'B2');
            
        end
        disp(['Regime ' num2str(i) ' is done']);
        
    end
end

clearvars -except Regims Res
exprt(Res,Regims)
