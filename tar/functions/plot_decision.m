%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G-W decision rules 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datesind = 2:48:342;

figure('Name', 'Decision, LS, horizon 1')
for ii = 1:4
    subplot(2, 2, ii)
    tempData = squeeze(GWdr.ls(:,ii,1)); 
    plot(datesn2, tempData), axis tight, title(vnamesvar(ii))
    hline(0, 'r')
    hold off
    xlabel('Time')
    set(gca, 'XTick', datesn2(datesind));
    set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
end


figure('Name', 'Decision, LS, horizon 2')
for ii = 1:4
    subplot(2, 2, ii)
    tempData = squeeze(GWdr.ls(:,ii,2)); 
    plot(datesn2, tempData), axis tight, title(vnamesvar(ii))
    hline(0, 'r')
    hold off
    xlabel('Time')
	set(gca, 'XTick', datesn2(datesind));
    set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
end

figure('Name', 'Decision, LS, horizon max')
for ii = 1:4
    subplot(2, 2, ii)
    tempData = squeeze(GWdr.ls(:,ii,end)); 
    plot(datesn2, tempData), axis tight, title(vnamesvar(ii))
    hline(0, 'r')
    hold off
    xlabel('Time')
	set(gca, 'XTick', datesn2(datesind));
    set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
end


figure('Name', 'Decision, WLS L, horizon 1')
for ii = 1:4
    subplot(2, 2, ii)
    tempData = squeeze(GWdr.wlsL(:,ii,1)); 
    plot(tempData), axis tight, title(vnamesvar(ii))
    hline(0, 'r')
    hold off
    xlabel('Time')
end

figure('Name', 'Decision, WLS LR, horizon 1')
for ii = 1:4
    subplot(2, 2, ii)
    tempData = squeeze(GWdr.wlsLR(:,ii,1)); 
    plot(tempData), axis tight, title(vnamesvar(ii))
    hline(0, 'r')
    hold off
    xlabel('Time')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delete temp*
