%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G-W decision rules 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datesind = 2:48:342;

figure('Name', 'Decision rule')

ii = 3;
hh = 4;

tempLS   = -squeeze(GWdr.ls(:,ii,hh)); 
tempRMSE = squeeze(GWdr.rmse(:,ii,hh));
% Signs set up so that for both RMSE and LS "criterion>0 <-> pick TAR" 

hold on
plot(datesn2, tempRMSE, 'b', 'Linewidth', 1.5)
plot(datesn2, tempLS, 'r', 'Linewidth', 1.5)
hline(0, 'k')
hold off

axis tight

ylim([-0.4 0.4])
legend('RMSE', 'LS', 'Location', 'SouthEast'), legend BOXOFF
% title(vnamemsvar(ii))

 
% xlabel('Time')
set(gca, 'XTick', datesn2(datesind));
set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));

% % % % Could also try:
% % % hold on
% % % [AX,H1,H2] = plotyy(datesn2,tempRMSE,datesn2,tempLS,'plot');
% % % hline(0, 'k')
% % % hold off
% % % % But then the 0s on the y axis are not aligned...
% % % 

delete temp*

