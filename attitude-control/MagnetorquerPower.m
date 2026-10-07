function P_MTQ = MagnetorquerPower(T)

%% Program Summary

% INPUTS
% T                     (3, 1) control torque in 3 dimensions, should be
% lower than 5.36 uNm for each value (x, y, z)

% OUTPUT
% P_MTQ                 (3, 1) power vector: amount of power to be supplied
% to each orthogonal magnetorquer (x, y, z)

%% Finding power using linear relationship between torque and power

% T = [-6.72e-7, -1.6977e-7, 2.69e-8]; % use for testing


% torque = magnetic dipole moment x earths magnetic field
% t = M x b
% Nm = Am^2 x T

maxDipoleStrength = 0.39; % Am^2
minDipoleStrength = 0;

magFieldStrength = 2.5e-5*(6371/(6371 + 500))^3; % orbit height = 500 km
% at the minimum magnetic field, e.g. at the equator

maxTorque = maxDipoleStrength * magFieldStrength;

maxTorque3U = maxTorque*[1, 1, 1];

torquePerWatt = maxTorque3U/1.5; % Nm/W


% Calculate the power required for each torque component
P = (T ./ torquePerWatt);
negative = [0,0,0]; %flag if power values are negative

for i = 1:3
    if P(i) < 0
        negative(i) = 1;
        P(i) = -P(i);
    end
end

P = P + 0.25; % correct for idle power offset of 0.25W

P_MTQ = P;

end