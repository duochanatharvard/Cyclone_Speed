clear;
dir_load = '/Volumes/My Passport Pro/GOES16/';

%% *************************************************************************
% Read the reference information
% *************************************************************************
ref = load([dir_load,'IRMA_track.mat']);
temp = ref.iso_time(2:end);
for i = 1:numel(temp)
    a = temp{i};
    aa = [str2num(a(1:4)) str2num(a(6:7)) str2num(a(9:10)),...
          str2num(a(12:13)) str2num(a(15:16)) str2num(a(18:19))];
    ref.time(i) = (datenum(aa) - datenum([2000 01 01 12 0 0])) * 86400;
end
ref.lon = ref.longitude(2:end);
ref.lat = ref.latitude(2:end);
% plot3(ref.longitude,ref.latitude,ref.time,'.')

%% *************************************************************************
% Convert the reference into geostationary projection
% *************************************************************************
[ref.x,ref.y] = Duo_LL2GS(ref.lon,ref.lat,-89.5);

% *************************************************************************
% Check if the information of frame make sense
% *************************************************************************
grid_x = repmat([1:61]',1,61);
grid_y = repmat([1:61],61,1);
POS = [];
ct = 0;
for dy = 248:252
    for hr = 1:24
        
        disp([num2str(dy),'  day   ',num2str(hr),'  hr'])
        
        dir_load_2 = [dir_load,num2str(dy),'/'];
        cdy = CDF_num2str(dy,3);
        chr = CDF_num2str(hr-1,2);
        file_load = [dir_load_2,'OR_ABI-L1b-RadM1-M3C13_G16_s2017',cdy,chr,'.mat'];
        load(file_load,'x','y','t','Rad');
        y = flipud(y);

        ref_x = interp1(ref.time,ref.x,t);
        ref_y = interp1(ref.time,ref.y,t);
        
        for tim = 1:numel(ref_x)
 
            rad     = fliplr(Rad(:,:,tim));  % to be compatible with flipud(y(:,tim));
            pix_x   = discretize(ref_x(tim),x(:,tim));
            pix_y   = discretize(ref_y(tim),y(:,tim));
            
            if ~isnan(pix_x) && ~isnan(pix_y),
                rad_sub = rad(pix_x+[-30:30], pix_y+[-30:30]); 
                sub_x   = x(pix_x+[-30:30],tim);
                sub_y   = y(pix_y+[-30:30],tim);

                temp = rad_sub;
                temp = temp - 20;
                temp(temp<0) = 0;

                w = temp.^10;
                tc_pix_x = nansum(grid_x(:).*w(:)) ./ nansum(w(:));
                tc_pix_y = nansum(grid_y(:).*w(:)) ./ nansum(w(:));

                tc_x = interp1(1:61,sub_x,tc_pix_x);
                tc_y = interp1(1:61,sub_y,tc_pix_y);

                ct = ct + 1;
                POS(ct,:) = [tc_x tc_y t(tim)];
            end
        end
    end
end

save([dir_load,'Geostationary_Projection_Irma_power_10.mat'],'POS','-v7.3')

%% *************************************************************************
% Convert back to lon lat grid
% *************************************************************************
clear;
dir_load = '/Volumes/My Passport Pro/GOES16/';
load([dir_load,'Geostationary_Projection_Irma_power_10.mat']);

file_tab = [dir_load,'ABI-L1b-RadF-M3C13_grids_backup.nc'];
lon = ncread(file_tab,'lon');
lat = ncread(file_tab,'lat');
x = ncread(file_tab,'x');
y = ncread(file_tab,'y');

lon_irma = interp2(x,y,lon,POS(:,1),POS(:,2));
lat_irma = interp2(x,y,lat,POS(:,1),POS(:,2));

ds = distance(lat_irma(1:end-1),lon_irma(1:end-1),lat_irma(2:end),lon_irma(2:end));
dx = lon_irma(2:end) - lon_irma(1:end-1);
dx = dx.*cos(lat_irma(1:end-1) / 180 * pi);
dy = lat_irma(2:end) - lat_irma(1:end-1);
dt = POS(2:end,3) - POS(1:end-1,3);
spd = ds./dt * 3600 * 111;
sx  = dx./dt * 3600 * 111;
sy  = dy./dt * 3600 * 111;

ct = 1:9000;
time_irma = POS(1:9000,3);
spd_irma = [spd(1:9000) sx(1:9000) sy(1:9000)];
time_irma = round((time_irma - time_irma(1)) / 30);
time_fft = 0:8018;
spd_fft = interp1(time_irma,spd_irma,time_fft);

spd_fft(find(isnan(spd_fft(:,1))),:) = repmat(nanmean(spd_fft,1),nnz(isnan(spd_fft(:,1))),1);
a = fft(spd_fft(:,1));
spd_fft = spd_fft - repmat(nanmean(spd_fft,1),size(spd_fft,1),1);

%%
figure(1); clf; hold on;
[P,s,ci] = pmtmPH(spd_fft(:,1),1/2880,3,1);
[P,s,ci] = pmtmPH(spd_fft(:,2),1/2880,3,1);
[P,s,ci] = pmtmPH(spd_fft(:,3),1/2880,3,1);
legend({'Great Circle','Dx','Dy'})
grid on
grid minor
set(gca,'fontsize',16)

fftPH(spd_fft(:,2),1/2880,'r'); 


%% Some analysis srtipt from Peter
clear;
load Irma_tracks.mat;

t = (t_irma - t_irma(1) + 30) / 86400;
pl=find(t<3.5 | isnan(lon_irma) | isnan(lat_irma)); 
t=t(pl);
lon=lon_irma(pl);
lat=lat_irma(pl);

pl=find(~isnan(t) & ~isnan(lon) & ~isnan(lat) & t>0.5);
t=t(pl);
lon=lon(pl);
lat=lat(pl);

ti=t(1):min(diff(t)):t(end);
loni=interp1(t,lon,ti);
lati=interp1(t,lat,ti); 

sm=101;
lons=smoothPH(diff(loni),4*60*60/30)*110*1000/(24*60*60*diff(ti(1:2)));
lats=smoothPH(diff(lati),4*60*60/30)*110*1000/(24*60*60*diff(ti(1:2)));

figure(1); clf; hold on;
plot(lon,lat,'k'); 
plot(loni,lati,'r.'); 
axis tight;
figure(2); clf; hold on;
fftPH(diff(loni),diff(ti(1:2))*24,'k'); 
logPH;
axis tight;

figure(3); clf; hold on;
fftPH(diff(lati),diff(ti(1:2))*24,'r'); 
logPH;
axis tight;

figure(4); clf; 
subplot(211);
plot(ti(2:end),lons,'k');

subplot(212); 
plot(ti(2:end),lats,'k')
axis tight;

figure(5); clf;
plot(lons,lats,'k');

