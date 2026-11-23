function rho = AtmosphericDensity(rECI, jD)
%AtmosphericDensity Computes the local air density
%   rho = AtmosphericDensity(rECI, jD)
%
%   Currently uses a very simple exponential decay model. See JB2008 for a
%   high fidelity model.
%
%   Inputs
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
persistent tbl
if isempty(tbl)
    %      h0 (km)   rho0 (kg/m^3)   H (km)
    tbl = [
        0      1.225       8.44;
        25     3.899e-2    6.49;
        30     1.774e-2    6.75;
        35     8.279e-3    7.07;
        40     3.972e-3    7.47;
        45     1.995e-3    7.83;
        50     1.057e-3    7.95;
        55     5.821e-4    7.73;
        60     3.206e-4    7.29;
        65     1.718e-4    6.81;
        70     8.770e-5    6.33;
        75     4.178e-5    6.00;
        80     1.905e-5    5.70;
        85     8.337e-6    5.41;
        90     3.396e-6    5.38;
        95     1.343e-6    5.74;
        100    5.297e-7    6.15;
        110    9.661e-8    8.06;
        120    2.438e-8    11.6;
        130    8.484e-9    16.1;
        140    3.845e-9    20.6;
        150    2.070e-9    24.6;
        160    1.224e-9    26.3;
        180    5.464e-10   33.2;
        200    2.789e-10   38.5;
        250    7.248e-11   46.9;
        300    2.418e-11   52.5;
        350    9.158e-12   56.4;
        400    3.725e-12   59.4;
        450    1.585e-12   62.2;
        500    6.967e-13   65.8;
        600    1.454e-13   79.0;
        700    3.614e-14  109.0;
        800    1.170e-14  164.0;
        900    5.245e-15  225.0;
        1000   3.019e-15  268.0;
    ];
end

if h < min(tbl(:, 1)) || h > max(tbl(:, 1))
    fprintf(['Altitude %.1f km is outside the supported range of 0-1000 km ' ...
        'for the atmospheric model. Consider another atmospheric model.\n'], h);
    rho = 0;
    return
end

index = discretize(h, tbl(:, 1));

h0 = tbl(index, 1);
rho0 = tbl(index, 2);
H = tbl(index, 3);

rho = rho0*exp(-(h-h0)/H);

end