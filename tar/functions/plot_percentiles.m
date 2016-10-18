%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Percentiles against data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pick horizon (in terms of months):
hh = 12;

% Choose a pair of percentiles (ordered as in 'percentiles'):
temppct = [1];

% Tickmarks for dates:
% datesind = 1:12:nobsGR;
datesind = 2:48:342;

% NOTHING TO BE CHANGED FROM HERE ON.

% Get position of the horizon in the results data
hhp = find(horizons==hh);
if isempty(hhp)
    error('There is no such horizon')
end

% Define data and un-do IP and CPI transformation (see arrangedata):
tempData1   = datavar;
tempData1   = tempData1 ./ repmat([12 1 12 1], T,1);

% Get cumulative data
tempData2 = tempData1;
for tt = T0+hh:size(tempData1, 1) 
    tempData2(tt, [1 3])  = sum(tempData1(tt-hh+1:tt, [1 3]), 1);
end

figure('Name', 'PCTL & data, horizon 1 (all)')

for ii = 1

    % pick 1st and last percentile:
    tempM1 = -squeeze(VAR.pct(:, ii, hhp, temppct ));    
    tempM2 = -squeeze(TAR.pct(:, ii, hhp, temppct ));
    tempM3 = -squeeze(TVTP.pct(:,ii, hhp, temppct ));
    
%     subplot(2, 2, ii)
    hold on
    plot(datesn2, tempM3, 'b--', 'LineWidth', 2), axis tight, 
    plot(datesn2, tempM1, 'b', 'LineWidth', 2), axis tight,
    plot(datesn2, tempM2, 'r', 'LineWidth', 2), axis tight, 
%     plot(datesn2, tempData2(T0+hh:T0+T2+hh-1, ii), 'k', 'LineWidth', 2)
    area(datesn2, tempData2(T0+hh:T0+T2+hh-1, ii), 'FaceColor',[.7 .8 .9], 'EdgeColor', [1 1 1])
    hline(0, 'k')
    hold off
    
%     title(vnamestar(ii))
   
    xlim ([min(datesn2(datesind)) max(datesn2(datesind))])
    set(gca, 'XTick', datesn2(datesind));
    set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));

    ylim auto
    
end
% Two percentiles:
legend('VAR^�', 'VAR', 'TAR', 'data', 'Location', 'Best'), legend BOXON


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear temp*
