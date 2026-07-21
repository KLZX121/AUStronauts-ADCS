function vB = MagModel(x, d)
%MagModel Magnetometer Model
%   vB = MagModel(x, d)
%   
%   Inputs
%   x               (:, 1)      state vector
%   d                           data struct
%   .jD             (1, 1)      Julian Date of Epoch
%   .kR             (1, 3)      indices of ECI pos in state vector
%   .kQ             (1, 4)      indices of attitude quaternion
%   .qBToS          (4, 1)      quaternion for body to sensor frame
%   .bias           (3, 1)      biases
%   .scale          (3, 1)      scale factors
%   .noise          (3, 1)      1-sigma noise
%   .quantization   (1, 1)      Least Significant Bit (LSB)
%   
%   Outputs
%   vB  (3, 1)  vector magnetic field measurement (T) (body frame)
%   
%   TODO: improve with IGRF
%   TODO: work out scale, bias, noise, quantization
%   TODO: rewrite function with only orbital position

% get measurements in integer counts, then convert integer counts to tesla

% changed line 66-67 from
% bBody   = QForm( x(d.kQ), bECI );
% bSensor = QForm( d.qBToS, bBody );
% to
% bBody   = QToDCM(x(d.kQ))*bECI;
% bSensor = QToDCM(d.qBToS)*bBody;
bCount = MeasMagnetometerEarth(x, d);

vB = bCount*d.quantization; % [T]

end
