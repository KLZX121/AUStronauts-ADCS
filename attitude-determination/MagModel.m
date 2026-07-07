function uB = MagModel(x, jD0)
%MagModel Magnetometer Model
%   uB = MagModel(x, jD0)
%   
%   Inputs
%   x   (:, 1)  state vector
%   jD0         Julian Date of epoch
%   
%   Outputs
%   uB  (3, 1)  unit vector magnetic field measurement (body frame)
%   
%   TODO: improve with IGRF
%   TODO: work out scale, bias, noise, quantization
%   TODO: rewrite function with only orbital position

% set data structure
dMag = MeasMagnetometerEarth;
dMag.jD = jD0;
% dMag.quantization = 1e-5;
% default values:
% d.kR - indices of ECI position vector in state vector (1:3)
% d.kQ - indices of attitude quaternion in state vector (7:10)
% d.qBToS - rotation from body frame to sensor (none)
% d.scale, d.bias, d.noise - sensor characteristics (ideal)
% d.quantization - LSB (1e-8)

% get measurements in integer counts, then convert integer counts to tesla
bCount = MeasMagnetometerEarth( x, dMag );
bMeas = bCount*dMag.quantization; % [T]
% convert to unit vector
uB = bMeas./norm(bMeas);


end
