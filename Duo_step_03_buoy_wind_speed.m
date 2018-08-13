clear;

dir = '/Volumes/My Passport Pro/SST/Buoy_data_3/';
file_load = [dir,'41002_raw.mat'];
load(file_load);

yr = DATA(:,1); mon = DATA(:,2); dy = DATA(:,3); hr = DATA(:,4);
yr(yr<100) = yr(yr<100) + 1900;
date = (datenum([yr, mon, dy]) - datenum([1970 1 1])) * 24 + round(hr);
date_total = (datenum([2014 12 31]) - datenum([1970 1 1])) * 24 + round(24);

wc = DATA(:,6); wc(wc == 999) = nan;
ws = DATA(:,7); ws(ws == 99)  = nan;

WC = nan(1,date_total); WC(date) = wc;
WS = nan(1,date_total); WS(date) = ws;
date_vec = datevec([1:date_total]'/24 + datenum([1970 1 1]));
logic = date_vec(:,2) == 2 & date_vec(:,3) == 29;

WC(logic) = []; WS(logic) = [];
WC = reshape(WC,24,365,45);
WS = reshape(WS,24,365,45);

U = - sin(WC ./ 180 .* pi) .* WS;
V = - cos(WC ./ 180 .* pi) .* WS; 

lon = -74.835 + 360;



%%
figure(2);clf; hold on;
mon_list = date_vec(1:24:365*24,2);
col = jetCD(12);

for mon = 1:12

    u = nanmean(nanmean(U(:,mon_list == mon,:),3),2);
    v = nanmean(nanmean(V(:,mon_list == mon ,:),3),2);

    % subplot(3,4,mon); hold on;
    % quiver(0,0,nanmean(u),nanmean(v),0,'color',[1 1 1]*.9,'linewi',2);
    patch(u([1:end 1]),v([1:end 1]),[1 1 1]*.7,'facealpha',0.2)
end

for mon = 1:12

    u = nanmean(nanmean(U(:,mon_list == mon,:),3),2);
    v = nanmean(nanmean(V(:,mon_list == mon ,:),3),2);

    plot(u([1:end 1]),v([1:end 1]),'-','linewi',2,'color',[1 1 1]*.7)

    for i = 1:24
        % quiver(0,0,u(i),v(i),0,'color',col(i,:));

        C0_LCL = rem(i + lon./15, 24);
        C0_LCL (C0_LCL < 0.501) = C0_LCL (C0_LCL < 0.501) + 24;
        C0_LCL = fix(C0_LCL - 0.501) + 1;
        C0_LCL (C0_LCL < 1) = C0_LCL (C0_LCL < 1) + 24;
        C0_LCL (C0_LCL > 24) = C0_LCL (C0_LCL > 24) - 24;

        plot(u(i),v(i),'.','markersize',30,'color',col(C0_LCL,:));
    end

    grid on;
    daspect([1 1 1])
    uu(mon) = nanmean(u);
    vv(mon) = nanmean(v);
end

for mon = 1:12
    text(uu(mon),vv(mon),num2str(mon),'color','k',...
        'fontsize',24,'fontweight','bold','horizontalalignment','center');
    text(uu(mon),vv(mon),num2str(mon),'color','w',...
        'fontsize',20,'fontweight','bold','horizontalalignment','center');
end

% plot(uu,vv,'-','linewi',3,'color',[1 1 1]*0.5)
plot(0,0,'rp','markerfacecolor','r','markersize',20)
h = colorbar; set(h,'ytick',[6 12 18 24])
ylabel(h,'local hour')
caxis([.5 24.5])
CDF_panel([-5 5 -5 5],[],{},'U (m/s)','V (m/s)');