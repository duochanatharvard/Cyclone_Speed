%% A script that generate the lat lon grid (look up table) 
%  for GOES16 geostationary projection.
%  This version uses iteration algorithum and is more precise than the
%  interpolated grids.
%  D. Chan & P. Chan, Harvard. Jul, 20, 2018

% Reference:
% Combal, Bruno, and Josu¨¦ Noel. Projection of Meteosat Images
% Into World Geodetic System WGS-84 Matching Spot-vegetation Grid.
% EUR-OP, 2009.
% pp.17

clear;

%% Set parameters 
intv = 0.2;        % resolution of the look up grid
intv_2 = 1;       % resolution of the output grids

%% Speficy the geostationary grid for the satellite data.
file = ['ABI-L1b-RadF-2018-008-06-OR_ABI-L1b-RadF-M3C13_G16',...
       '_s20180080615404_e20180080626182_c20180080626239.nc'];

file_save = 'ABI-L1b-RadF-M3C13_grids.nc';

x_sat_0 = ncread(file,'x');
y_sat_0 = ncread(file,'y');
x_sat_0 = x_sat_0(1:intv_2:end);
y_sat_0 = y_sat_0(1:intv_2:end);
[x_sat,y_sat] = meshgrid(x_sat_0,y_sat_0);
rad = ncread(file,'Rad');


%% Set Parameters for the geostationary projection.

disp('Preparing for grid and parameters ...')
r_eq = 6378137;                 % radius of the earth - major axis
r_pol = 6356752.3141;           % radius of the earth - minor axis
h = 35786023 + r_eq;            % distance from satellite to earth center
sub_lon = 0 ./ 180 .* pi;

lon = zeros(size(x_sat)) + sub_lon;
lat = zeros(size(y_sat));

%% Compute a look up table for each longitude and latitude.
disp('Computing a look up table for each longitude and latitude ...')
for i = 1:40

    c_lat = atan((r_pol.^2 ./ r_eq.^2) .* tan(lat));
    rl = r_pol ./ sqrt(1 - (r_eq.^2 - r_pol.^2) ./ r_eq.^2 .* cos(c_lat).^2);
    r1 = h - rl .* cos(c_lat) .* cos(lon - sub_lon);
    r2 = rl .* cos(c_lat) .* sin(lon - sub_lon);
    r3 = rl .* sin(c_lat);
    rn = sqrt(r1.^2 + r2.^2 + r3.^2);
    x = atan(r2./r1);
    y = asin(r3./rn);

    dif_x = x_sat - x;
    dif_y = y_sat - y;

    lambda = 6;
    lon = lon + dif_x * lambda;
    lat = lat + dif_y * lambda;
end

% contourf(dif_x(1:10:end,1:10:end));colorbar;
% contourf(dif_y(1:10:end,1:10:end));colorbar;

lon(isnan(rad)) = nan;
lat(isnan(rad)) = nan;

lon = lon / pi * 180;
lat = lat / pi * 180;

%% Saving data

disp('Saving data ...')
delete(file_save);

nccreate(file_save,'lon','Dimensions',{'x',size(lon,1),'y',size(lon,2)},...
        'Format','netcdf4','datatype','single','deflatelevel',1)
ncwrite(file_save,'lon', single(lon));

nccreate(file_save,'lat','Dimensions',{'x',size(lat,1),'y',size(lat,2)},...
        'Format','netcdf4','datatype','single','deflatelevel',1)
ncwrite(file_save,'lat', single(lat));

nccreate(file_save,'x','Dimensions',{'x',numel(x_sat_0)},...
        'Format','netcdf4','deflatelevel',1)
ncwrite(file_save,'x',x_sat_0);

nccreate(file_save,'y','Dimensions',{'y',numel(y_sat_0)},...
        'Format','netcdf4','deflatelevel',1)
ncwrite(file_save,'y',y_sat_0);

ncwriteatt(file_save,'lon','read_me','The longitude is relative to the longitude_of_projection_origin, which should be read from individual nc files');

