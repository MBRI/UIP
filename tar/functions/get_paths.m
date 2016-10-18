function yhat = get_paths(N, HORZ, T, L, Y, beta2, hlast,iamat,gmat)



%compute forecast
hhat=zeros(HORZ+L,N);
hhat(L,:)=log(hlast(end,:));
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y(T-L+1:T,:);
Q=diag(sqrt(gmat));
for fi=L+1:HORZ+L
    hhat(fi,:)=hhat(fi-1,:)+randn(1,N)*Q;
    epsilon=zeros(1,N);
        for jj=1:N
            epsilon(1,jj)=randn(1,1).*sqrt(exp(hhat(fi,jj)));
        end
    xhat=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
    end
    xhat=[xhat 1];
    
   
    yhat(fi,:) = xhat*reshape(beta2,N*L+1,N) + epsilon*iamat';

end
