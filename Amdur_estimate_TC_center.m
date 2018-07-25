%% A script that estimates the center of a well-formed tropical cyclone 
% Using Circular Hough Transformation
% Ted Amdur
% July 24, 2018

%Load the lon and lat mesh grid for GOES-16
lonMesh = ncread('ABI-L1b-RadF-M3C13_grids.nc', 'lon');
lonMesh = lonMesh';
latMesh = ncread('ABI-L1b-RadF-M3C13_grids.nc', 'lat');
latMesh = latMesh';
xInd = ncread('ABI-L1b-RadF-M3C13_grids.nc','x');
yInd = ncread('ABI-L1b-RadF-M3C13_grids.nc','y');
%Initialize structure array to hold centers (x,y data), radius, time, longitude of
%projection, and perspecive point height.
center.x=[]; center.y =[]; center.r =[];center.t=[];center.lonc=[];
center.h=[]; center.fracX = []; center.fracY=[];
center.lon=[]; center.lat=[];
center.xPix=[];center.yPix=[];

%----------------------------------------------------------------
% Read in mat files corresponding to one-minute frames of GOES-16 Mesoscale
%----------------------------------------------------------------
irmaStr = 'OR_ABI-L1b-RadM1-M3C13_G16_s2017'; %str not including day/hour
loadflag = 0;
for dayn = 248:253 %corresponding to calendar day of filepaths
    daypath = int2str(dayn);
    for hourn = 0:23 %corresponding to hour of filepath
        hourpath = int2str(hourn);
        if hourn < 10
            totpath = strcat(irmaStr, daypath, '0',hourpath, '.mat');
        else
            totpath = strcat(irmaStr, daypath, hourpath, '.mat');
        end
        try load(totpath);
            loadflag = 1;
        catch
            loadflag = 0;
        end
        
        if loadflag %process data if it can be loaded
            Rad(find(Rad > 1000)) = NaN;
            for ii = 1:length(Rad(1,1,:))
                %process rad to grayscale
                imgLim = max(max(Rad(:,:,ii)));
                img = (Rad(:,:,ii)/imgLim) .* 255;
                [~, circen, cirrad] = CircularHough_Grd(img, [1 50]);
                if ~isempty(circen) %no center may be identifiable or hurricane not in frame
                    xc = x(floor(circen(1,1)),ii); yc = y(floor(circen(1, 2)),ii); rc = cirrad(1);
                    %find fractional distance of x, y
                    fracX = (circen(1,1)-floor(circen(1,1)));
                    fracY = (circen(1,2)-floor(circen(1,2)));
                    tc = t(ii); lonc = longitude_of_projection_origin(ii);
                    hc = perspective_point_height(ii);
                    center.x = [center.x; xc]; center.y = [center.y; yc];
                    center.r = [center.r; rc]; center.t = [center.t; tc];
                    center.lonc = [center.lonc; lonc]; center.h = [center.h; hc];
                    center.fracX=[center.fracX;fracX]; center.fracY=[center.fracY;fracY];
                    center.xPix=[center.xPix; circen(1,1)]; center.yPix=[center.yPix;circen(1,2)];
                end
            end              
        end
    end
end

%----------------------------------------------------------------
% Find geographic centers
%----------------------------------------------------------------
id_x = discretize(center.x, xInd);
id_y = discretize(center.y, flipud(yInd));
lonSet = []; latSet = [];
for n = 1:length(id_x)
    lonSet = [lonSet; lonMesh(id_x(n), id_y(n))];
    latSet = [latSet; latMesh(id_x(n), id_y(n))];
    %add in fractional area
    fracLon = (lonMesh(id_x(n)+1, id_y(n)+1)-lonMesh(id_x(n), id_y(n))).*center.fracX(n);
    fracLat = (latMesh(id_x(n)+1, id_y(n)+1)-latMesh(id_x(n), id_y(n))).*center.fracY(n);
    lonSet(n) = lonSet(n) + fracLon; latSet(n) = latSet(n) + fracLat;
end 
center.lon = lonSet; center.lat = -1.*latSet; %assume lat is symmetrical

save('TC_test1.mat','center')
        
