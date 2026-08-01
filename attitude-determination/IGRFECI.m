function bECI = IGRFECI(rECI, jD)
%IGRFECI Computes the IGRF magnetic field in ECI coordinates
%   bECI = IGRFECI(rECI, jD)
%
%   Inputs
%   rECI    (3, 1)  Position in the ECI frame (m)
%   jD      (1, 1)  Julian date (UTC)
%
%   Outputs
%   bECI    (3, 1)  Magnetic field vector in the ECI frame (T)

reduction = 'IAU-2000/2006';

% convert JD to UTC datetime object
utc = datetime(jD, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
% remove TimeZone field since decyear can't handle it
utc.TimeZone = "";

% convert ECI to LLA (latitude, longitude, altitude)
lla = eci2lla(rECI', datevec(utc), reduction);
lat = lla(1); % (deg)
lon = lla(2); % (deg)
alt = lla(3); % (m)

% get magnetic field in NED (North-East-Down)
bNED = igrfmagm(alt, lat, lon, decyear(utc));
bNED = bNED' * 1e-9; % convert to column vector and T

% convert NED to ECEF (Earth-Centred-Earth-Fixed)
ECEF2NED = dcmecef2ned(lat, lon);
bECEF = ECEF2NED' * bNED;

% convert ECEF to ECI
ECI2ECEF = dcmeci2ecef(reduction, utc);
bECI = ECI2ECEF' * bECEF;

end