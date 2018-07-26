% A function that returns geostaionary position based on longitude,
% latitude, and longitude of projection origin
function [x,y] = Duo_LL2GS(lon,lat,sub_lon)

    r_eq = 6378137;                 % radius of the earth - major axis
    r_pol = 6356752.3141;           % radius of the earth - minor axis
    h = 35786023 + r_eq;            % distance from satellite to earth center
    lon = lon ./ 180 .* pi;
    lat = lat ./ 180 .* pi;
    sub_lon = sub_lon ./ 180 .* pi;

    c_lat = atan((r_pol.^2 ./ r_eq.^2) .* tan(lat));
    rl = r_pol ./ sqrt(1 - (r_eq.^2 - r_pol.^2) ./ r_eq.^2 .* cos(c_lat).^2);
    r1 = h - rl .* cos(c_lat) .* cos(lon - sub_lon);
    r2 = rl .* cos(c_lat) .* sin(lon - sub_lon);
    r3 = rl .* sin(c_lat);
    rn = sqrt(r1.^2 + r2.^2 + r3.^2);
    x = atan(r2./r1);
    y = asin(r3./rn);
    
end