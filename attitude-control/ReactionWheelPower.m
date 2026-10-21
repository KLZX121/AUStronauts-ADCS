function P_RWA = ReactionWheelPower(T)

%% Program Summary

% INPUTS
% T                     (3, 1) control torque in 3 dimensions, should be
% lower than 200 uNm for each value (x, y, z)

% OUTPUT
% P_RWA                 (3, 1) power vector: amount of power to be supplied
% to each orthogonal reaction wheel (x, y, z)

%% Finding power using linear relationship between torque and power

% T = [-6.72e-6, -1.6977e-6, 2.69e-7]; % use for testing (larger torques
% than MTQ testing)

maxTorque1U = 200e-6; % 200 mNm (Nm e-3)

maxTorque3U = maxTorque1U*[1, 1, 1];

torquePerWatt = maxTorque3U/0.3; % Nm/W

% Calculate the power required for each torque component
P = (T ./ torquePerWatt);
negative = [0,0,0]; %flag if power values are negative

for i = 1:3
    if P(i) < 0
        negative(i) = 1;
        P(i) = -P(i);
    end
end

P = P + 0.1; % correct for idle power offset of 0.1W

P_RWA = P;

end