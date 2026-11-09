function inEclipse = CheckEclipse(r, rSatSun)
%inEclipse = CheckEclipse(r, rSun)
%   Determines if current orbital position is in eclipse using a simplified
%   check. Does not consider penumbra.
%
%   Inputs
%   r           (3, 1)  Satellite orbital position (earth -> sat) (ECI) (km)
%   rSatSun     (3, 1)  Sun position (sat -> sun) (ECI) (km)
%
%   Outputs
%   inEclipse   (bool)  If satellite is in eclipse

inEclipse = false;
REarth = 6378; % earth max radius (km)

% earth -> sat
uEarthSat = r ./ norm(r);

% sun to earth
rSunEarth = -rSatSun-r;
uSunEarth = rSunEarth ./ norm(rSunEarth);


% is in sun side
if dot(uEarthSat, -uSunEarth) > 0
    return
end


% perpendicular vector
v = cross(uSunEarth, [1; 0; 0]);
if norm(v) < eps
    v = cross(uSunEarth, [0; 1; 0]);
end

v = v ./ norm(v);
v = v .* REarth;

% sun->earth outer radius
rSunEarthEdge = rSunEarth+v;
uSunEarthEdge = rSunEarthEdge ./ norm(rSunEarthEdge);

% angle between sun->earth and sun->earth outer radius
earthBounds = acos(dot(uSunEarth, uSunEarthEdge));

% sun -> sat
uSunSat = -rSatSun ./ norm(rSatSun);

% angle between sun->satellite and sun->earth
satAng = acos(dot(uSunSat, uSunEarth));

% if satAng is less than earth bounds then in eclipse
if satAng < earthBounds
    inEclipse = true;
end