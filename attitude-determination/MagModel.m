function vB = MagModel(x, jD0)
%MagModel Magnetometer Model
%   vB = MagModel(x, jD0)
%   
%   Inputs
%   x   (:, 1)  state vector
%   jD0         Julian Date of epoch
%   
%   Outputs
%   vB  (3, 1)  vector magnetic field measurement (T) (body frame)
%   
%   TODO: improve with IGRF
%   TODO: work out scale, bias, noise, quantization
%   TODO: rewrite function with only orbital position

% set data structure
dMag = MeasMagnetometerEarth;
dMag.jD = jD0;
dMag.quantization = 1e-20;
% default values:
% d.kR - indices of ECI position vector in state vector (1:3)
% d.kQ - indices of attitude quaternion in state vector (7:10)
% d.qBToS - rotation from body frame to sensor (none)
% d.scale, d.bias, d.noise - sensor characteristics (ideal)
% d.quantization - LSB (1e-8)

% get measurements in integer counts, then convert integer counts to tesla

% changed line 66-67 from
% bBody   = QForm( x(d.kQ), bECI );
% bSensor = QForm( d.qBToS, bBody );
% to
% bBody   = QToDCM(x(d.kQ))*bECI;
% bSensor = QToDCM(d.qBToS)*bBody;
bCount = MeasMagnetometerEarth( x, dMag );

vB = bCount*dMag.quantization; % [T]

end
