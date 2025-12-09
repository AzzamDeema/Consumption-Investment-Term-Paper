function ffr_cbj_inflation_plot()

dataFile = 'data.xlsx';

%% data
tbl       = readtable(dataFile);

timeStr   = string(tbl.time);
ffr       = tbl.ffr;
cbj       = tbl.cbj_rate;
inflation = tbl.inflation;

%% build quarterly datetime axis from time string
n = numel(timeStr);
yearVec = zeros(n,1);
qVec    = zeros(n,1);
for t = 1:n
    tok = regexp(timeStr(t),'(\d+):(\d+)','tokens','once');  % YYYY:Q
    yearVec(t) = str2double(tok{1});
    qVec(t)    = str2double(tok{2});
end
time = datetime(yearVec, 3*qVec-2, 1);  % first month of each quarter

%% recession periods
recessions = [
    datetime(2008,1,1)  datetime(2009,6,30)
    datetime(2020,1,1)  datetime(2020,12,31)
];

%% plot
figure;
plot(time, ffr,       'b','LineWidth',1.5); hold on;
plot(time, cbj,       'r','LineWidth',1.5);
plot(time, inflation, 'k','LineWidth',1.5);

xlabel('Time');
ylabel('Percent / annualized rate');
legend({'FFR','CBJ policy rate','Inflation'},'Location','best');
title('FFR, CBJ policy rate, and Inflation with recessions');

%% recession bands
yl = ylim;
for i = 1:size(recessions,1)
    x0 = recessions(i,1); x1 = recessions(i,2);
    patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.4 0.4 0.4], 'EdgeColor','none','FaceAlpha',0.4);
end
ylim(yl);  % restore y-limits
uistack(findobj(gca,'Type','line'),'top');

xlim([time(1) time(end)]);
hold off;

end