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
d_mag = MeasMagnetometerEarth;
d_mag.jD = jD0;
% d_mag.quantization = 1e-20;
% default values:
% d.kR - indices of ECI position vector in state vector (1:3)
% d.kQ - indices of attitude quaternion in state vector (7:10)
% d.qBToS - rotation from body frame to sensor (none)
% d.scale, d.bias, d.noise - sensor characteristics (ideal)
% d.quantization - LSB (1e-8)

% get measurements in integer counts, then convert integer counts to tesla
b_count = MeasMagnetometerEarth( x, d_mag );
b_meas = b_count*d_mag.quantization; % [T]
% convert to unit vector
uB = b_meas./norm(b_meas);


end
