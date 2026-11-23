function rho = AtmosphericDensity(rECI, jD)
%AtmosphericDensity Computes the local air density
%   rho = AtmosphericDensity(rECI, jD)
%
%   Input
%   rECI    (3, 1)  Satellite's ECI position (km)
%   jD      (1, 1)  Julian date
%
%   Output
%   rho     (1, 1)  Atmospheric density (kg/m^3)

% convert eci to lla to get altitude
utc = datetime(jD, 'ConvertFrom', 'juliandate', 'TimeZone', 'UTC');
lla = eci2lla(rECI'.*1e3, datevec(utc), 'IAU-2000/2006');
h = lla(3); % m

h = h*1e-3; % km

% markley table 11.1
if h >= 300 && h < 350
    h0 = 300;
    rho0 = 2.418e-11;
    H = 52.5;
elseif h >= 350 && h < 400
    h0 = 350;
    rho0 = 9.158e-12;
    H = 56.4;
elseif h >= 400 && h < 450
    h0 = 400;
    rho0 = 3.725e-12;
    H = 59.4;
elseif h >= 450 && h < 500
    h0 = 450;
    rho0 = 1.585e-12;
    H = 62.2;
elseif h >= 500 && h < 600
    h0 = 500;
    rho0 = 6.967e-13;
    H = 65.8;
end

rho = rho0*exp(-(h-h0)/H);

end