%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histograms of PIT:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of buckets for histograms:
kk = round(T2/25);

TVTP.pit0 = TVTP.pit;
TVTP.pit0(:,4) = NaN; 

figure('Name', 'PITs')
for ii = 1:N
    subplot(3, N, ii), hist(VAR.pit(:, ii), kk), axis square, title(vnamesvar(ii)), ylabel('VAR')
    subplot(3, N, N+ii), hist(TAR.pit(:, ii), kk), axis square, title(vnamesvar(ii)), ylabel('TAR')
    subplot(3, N, 2*N+ii), hist(TVTP.pit0(:, ii), kk), axis square, title(vnamesvar(ii)), ylabel('MSVAR')
end

% figure('Name', 'PITs, inverse N')
% for ii = 1:N
%     subplot(2, N, ii), hist(VAR.pitin(:, ii), kk), axis square, title('VAR')
%     subplot(2, N, N+ii), hist(TAR.pitin(:, ii), kk), axis square, title('TAR')
% end

clear kk