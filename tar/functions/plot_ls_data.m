%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% log-scores against data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pick horizon (in terms of months):
hh = 12;

% NOTHING TO BE CHANGED FROM HERE ON.

% Get position of the horizon in the results data
hhp = find(horizons==hh);
if isempty(hhp)
    error('There is no such horizon')
end

% Tickmarks for dates:
datesind = 2:48:342;

% Define data and un-do IP and CPI transformation (see arrangedata):
tempData1   = datavar;
tempData1   = tempData1 ./ repmat([12 1 12 1], T,1);

% Get cumulative data
tempData2 = tempData1;
for tt = T0+hh:size(tempData1, 1) 
    tempData2(tt, [1 3])  = sum(tempData1(tt-hh+1:tt, [1 3]), 1);
end


% -------------------------------------------------------------------------
figure('Name', 'LS & data, IP')

x1 = datesn2;
x2 = x1;
y1 = tempData2(T0+hh:T0+T2+hh-1, 1);
y2 = TVTP.ls(:,1,hhp);
y3 = VAR.ls(:,1,hhp);
y4 = TAR.ls(:,1,hhp);

[AX,H1,H2] = plotyy(x1, y1, x2, [y2 y3 y4], @area, @plot);

set(H1, 'FaceColor',[.7 .8 .9], 'EdgeColor', [1 1 1]);
set(H2(1), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
set(H2(2), 'Color', 'b', 'LineWidth', 2);
set(H2(3), 'Color', 'r', 'LineWidth', 2);
hline(0, 'k')

% linkaxes(AX, 'xy') 
set(AX, 'XTick', datesn2(datesind));
set(AX, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
legend('data', 'VAR^�', 'VAR', 'TAR', 'Location', 'SouthWest'), legend BOXOFF


% -------------------------------------------------------------------------
figure('Name', 'LS & data, R')

x1 = datesn2;
x2 = x1;
y1 = tempData2(T0+hh:T0+T2+hh-1, 2);
y2 = TVTP.ls(:,2,hhp);
y3 = VAR.ls(:,2,hhp);
y4 = TAR.ls(:,2,hhp);

[AX,H1,H2] = plotyy(x1, y1, x2, [y2 y3 y4], @area, @plot);

set(H1, 'FaceColor',[.7 .8 .9], 'EdgeColor', [1 1 1]);
set(H2(1), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
set(H2(2), 'Color', 'b', 'LineWidth', 2);
set(H2(3), 'Color', 'r', 'LineWidth', 2);
hline(0, 'k')

% linkaxes(AX, 'xy') 
set(AX, 'XTick', datesn2(datesind));
set(AX, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
legend('data', 'VAR^�', 'VAR', 'TAR', 'Location', 'SouthWest'), legend BOXOFF


% -------------------------------------------------------------------------
figure('Name', 'LS & data, CPI')
y1 = tempData2(T0+hh:T0+T2+hh-1, 3);
y2 = TVTP.ls(:,3,hhp);
y3 = VAR.ls(:,3,hhp);
y4 = TAR.ls(:,3,hhp);

[AX,H1,H2] = plotyy(x1, y1, x2, [y2 y3 y4], @area, @plot);

set(H1, 'FaceColor',[.7 .8 .9], 'EdgeColor', [1 1 1]);
set(H2(1), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
set(H2(2), 'Color', 'b', 'LineWidth', 2);
set(H2(3), 'Color', 'r', 'LineWidth', 2);
hline(0, 'k')

% linkaxes(AX, 'xy') 
set(AX, 'XTick', datesn2(datesind));
set(AX, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
legend('data', 'VAR^�', 'VAR', 'TAR', 'Location', 'SouthWest'), legend BOXOFF


% -------------------------------------------------------------------------
figure('Name', 'LS & data, FIN')
y1 = tempData2(T0+hh:T0+T2+hh-1, 4);
y2 = VAR.ls(:,4,hhp);
y3 = TAR.ls(:,4,hhp);

[AX,H1,H2] = plotyy(x1, y1, x2, [y2 y3], @area, @plot);

set(H1, 'FaceColor',[.7 .8 .9], 'EdgeColor', [1 1 1]);
set(H2(1), 'Color', 'b', 'LineWidth', 2);
set(H2(2), 'Color', 'r', 'LineWidth', 2);
hline(0, 'k')

% linkaxes(AX, 'xy') 
set(AX, 'XTick', datesn2(datesind));
set(AX, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
legend('data', 'VAR', 'TAR', 'Location', 'SouthWest'), legend BOXOFF



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delete temp*
