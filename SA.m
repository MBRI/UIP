function [dt]=SA(x,freq)
% Seasonaly adjustement time series data
% x data
% date: time indicator
% freq 12 monthly, 4 frequently, 1 yearly
%% Step 1.

[T,K ]= size(x);
%% Step 2. Apply a 13-term moving average.
% Smooth the data using a 13-term moving average.
% To prevent observation loss, repeat the first and last smoothed values six(for monthly) times.
% Subtract the smoothed series from the original series to detrend the data.
% Add the moving average trend estimate to the observed time series plot.
if ~exist('freq','var')
    freq=4;
end
dt=nan(T,K);
%sW13 = [1/24;repmat(1/12,11,1);1/24];
sW13 = [1/(2*freq);repmat(1/freq,freq-1,1);1/(2*freq)];
for k=1:K
    y=x(:,k);
    yS = conv(y,sW13,'same');
    % w(1) = u(1)*v(1)
    % w(2) = u(1)*v(2)+u(2)*v(1)
    % w(3) = u(1)*v(3)+u(2)*v(2)+u(3)*v(1)
    % ...
    % w(n) = u(1)*v(n)+u(2)*v(n-1)+ ... +u(n)*v(1)
    % ...
    % w(2*n-1) = u(n)*v(n)
    yS(1:freq/2) = yS(freq/2+1); yS(T-(freq/2-1):T) = yS(T-freq/2);
    
    xt = y-yS;
    
    %% Step 3. Create seasonal indices.
    sidx = cell(freq,1);
    for i = 1:freq
        sidx{i,1} = i:freq:T;
    end
    
    %% Step 4. Apply a stable seasonal filter.
    
    %Apply a stable seasonal filter to the detrended series, xt. Using the indices constructed in Step 3, average the detrended data corresponding to each period. That is, average all of the January values (at indices 1, 13, 25,...,61), and then average all of the February values (at indices 2, 14, 26,...,62), and so on for the remaining months. Put the smoothed values back into a single vector.
    
    % Center the seasonal estimate to fluctuate around zero.
    
    sst = cellfun(@(x) mean(xt(x),'omitnan'),sidx);
    
    % Put smoothed values back into a vector of length N
    nc = floor(T/freq); % no. complete years
    rm = mod(T,freq); % no. extra perids
    sst = [repmat(sst,nc,1);sst(1:rm)];
    
    % Center the seasonal estimate (additive)
    sBar = mean(sst); % for centering
    sst = sst-sBar;
    
    %% Step 5. Deseasonalize the series.
    
    % Subtract the estimated seasonal component from the original data.
    
    dt(:,k) = y - sst;
%     plot(y,'r')
%     hold on
%     plot(yS,'b')
%     plot( dt(:,k),'g')
%     hold off
end