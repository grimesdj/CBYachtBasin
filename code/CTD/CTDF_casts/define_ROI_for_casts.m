function [perim,Lat0,Lon0,LatInterp,LonInterp,East0,North0,EastInterp,NorthInterp,distInterp] = define_ROI_for_casts(CTDFfiles);
% load('/Users/derekgrimes/git/cbflood/code_v2/FUNCTIONS/regions.mat')
perim  =[ ... 
   34.0586  -77.8927;...
   34.0571  -77.8914;...
   34.0560  -77.8910;...
   34.0534  -77.8909;...
   34.0515  -77.8908;...
   34.0495  -77.8908;...
   34.0476  -77.8908;...
   34.0461  -77.8910;...
   34.0435  -77.8914;...
   34.0404  -77.8921;...
   34.0377  -77.8930;...
   34.0364  -77.8932;...
   34.0362  -77.8923;...
   34.0408  -77.8908;...
   34.0436  -77.8898;...
   34.0452  -77.8893;...
   34.0469  -77.8894;...
   34.0493  -77.8893;...
   34.0514  -77.8891;...
   34.0518  -77.8890;...
   34.0537  -77.8883;...
   34.0547  -77.8882;...
   34.0576  -77.8889;...
   34.0596  -77.8891;...
   34.0598  -77.8909;...
   34.0598  -77.8909;...
   34.0598  -77.8909;...
    34.0586  -77.8927];
%
% Origin for along-channel cooridates
Lat0 = 34.0362;
Lon0 = -77.89275;

Lat = [];
Lon = [];
idx = [];
for ii=1:length(CTDFfiles);
    % 1.0) get CTD file name
    CTDFfile = CTDFfiles{ii};
    % 2.0) load CTDF
    CTD = load(CTDFfile);
    % 3) find all points in ROI:
    inROI = inpolygon(CTD.Latitude,CTD.Longitude,perim(:,1),perim(:,2));
    idx = cat(1,idx,ii*ones(sum(inROI),1));
    Lat = cat(1,Lat,CTD.Latitude(inROI)');
    Lon = cat(1,Lon,CTD.Longitude(inROI)');
end

Lat = cat(1,Lat,perim(:,1));
Lon = cat(1,Lon,perim(:,2));

% 4) use parametric curve to bin-average location
dLat = Lat-Lat0;
dLon = Lon-Lon0;
t = sqrt(dLat.^2 + dLon.^2);

% 4.1) Normalize t to the range [0, 1] (optional but helpful)
tNorm = t / max(t);% or just use t directly

% 4.2) bin average every dt=0.05 for 20 points
dbin = 0.05;
tInterp = [0:dbin:1]';
LatInterp = nan*tInterp;
LonInterp = nan*tInterp;

p = 0.999;
fx = fit(tNorm,Lat,'smoothingspline','SmoothingParam',p);
fy = fit(tNorm,Lon,'smoothingspline','SmoothingParam',p);

LatInterp = fx(tInterp);
LonInterp = fy(tInterp);

figure, geoplot(Lat,Lon,'.')
hold on,geoplot(LatInterp, LonInterp,'-m')
% $$$ hold on,geoplot(CTD.Latitude, CTD.Longitude,'*r')
% $$$ hold on,geoplot(CTD.Latitude(inROI), CTD.Longitude(inROI),'*g')


% $$$ 
% $$$ for ii=1:length(tInterp)
% $$$     idt = tNorm>=(tInterp(ii)-dbin*2) & tNorm<=(tInterp(ii)+dbin*2);
% $$$     LatInterp(ii) = nanmean(Lat(idt));
% $$$     LonInterp(ii) = nanmean(Lon(idt));
% $$$ end
% $$$ ig = ~isnan(LonInterp) & ~isnan(LatInterp);
% $$$ tInterp = tInterp(ig);
% $$$ LatInterp = LatInterp(ig);
% $$$ LonInterp = LonInterp(ig);
%
% 5) convert to meters and get distance along curve:
% define projection:
proj = projcrs(32119);
% Convert to ENU:
[EastInterp,NorthInterp] = projfwd(proj,LatInterp,LonInterp);
[East0,North0] = projfwd(proj,Lat0,Lon0);

distInterp = sqrt( (EastInterp-East0).^2 + (NorthInterp-North0).^2 );









% $$$ % 1. Create a cumulative chord-length parameter t
% $$$ dLat = perim(:,1)-Lat0;
% $$$ dLon = perim(:,2)-Lon0;
% $$$ t = sqrt(dLat.^2 + dLon.^2);
% $$$ 
% $$$ % 2. Normalize t to the range [0, 1] (optional but helpful)
% $$$ tNorm = t / max(t(end));% or just use t directly
% $$$ 
% $$$ % 2.1 bin average every dt=0.1;
% $$$ dbin = 0.1;
% $$$ tInterp = [0:dbin:1]';
% $$$ LatInterp = nan*tInterp;
% $$$ LonInterp = nan*tInterp;
% $$$ for ii=1:length(tInterp)
% $$$     idt = tNorm>=(tInterp(ii)-dbin) & tNorm<=(tInterp(ii)+dbin);
% $$$     LatInterp(ii) = mean(perim(idt,1));
% $$$     LonInterp(ii) = mean(perim(idt,2));
% $$$ end
% $$$ % now up-sample:
% $$$ tmp = linspace(0, 1,20);
% $$$ X   = interp1(tInterp,[LatInterp,LonInterp],tmp);
% $$$ tInterp = tmp;
% $$$ LatInterp = X(:,1);
% $$$ LonInterp = X(:,2);
% $$$ % define projection:
% $$$ proj = projcrs(32119);
% $$$ % Convert to ENU:
% $$$ [EastInterp,NorthInterp] = projfwd(proj,LatInterp,LonInterp);
% $$$ [East0,North0] = projfwd(proj,Lat0,Lon0);
% $$$ distInterp = sqrt( (EastInterp-East0).^2 + (NorthInterp-North0).^2 );
% $$$ %

end

