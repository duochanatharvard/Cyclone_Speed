function coord = Amdur_track_estimate(time)
%Create a track estimate for each image frame using linear interpolation of NHC
%estimated track location.
%Ted Amdur
% 7/25/2018

%read in NHC center estimates
read_center_file;
J200epoch = datetime('2000-01-01 12:00:00');
%create time comparison, converting to epoch time
tref = IRMAtrack{(2:end),{'iso_time'}};
tref = cell2mat(tref);
tref = datetime(tref(:,1:19));
tref = posixtime(tref) - posixtime(J200epoch);

lonref = IRMAtrack{(2:end),{'longitude'}};
latref = IRMAtrack{(2:end),{'latitude'}};

%For each time value, find time ref below and above for interpolation
coord.lat = zeros(length(time),1);
coord.lon = zeros(length(time),1);
for ii = 1:length(time)
    %Get refernce position before and after query
    indA = find(tref<=time(ii),1,'last');
    indB = find(tref>time(ii), 1, 'first');
    tDiff = tref(indB)-tref(indA);
    
    lonM = (lonref(indB)-lonref(indA))/tDiff;
    latM = (latref(indB)-latref(indA))/tDiff;
    
    x = time(ii)-tref(indA);
    
    %find position as a function of y=mx+b
    coord.lat(ii) = latM*x + latref(indA);
    coord.lon(ii) = lonM*x + lonref(indA);
end
