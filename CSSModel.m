function uS = CSSModel(x, uS_ref)
%CSSModel Coarse Sun Sensor Model
%   uS = CSSModel(x, uS_ref)
%   
%   Inputs
%   x       (:, 1)  state vector
%   uS_ref  (3, 1)  reference unit sun vector (reference frame)
%   
%   Outputs
%   uS      (3, 1)  measured unit sun vector (body frame)
%   
%   TODO: analyse sun sensor script
%   TODO: add line of sight for sun sensor (earth blocks sun)
%   TODO: add fov of sensor
%   TODO: rewrite function for orbital parameters

% set data structure
d_sun = MeasSunSensorAnalog;
% sensor unit vectors (6 sensors)
d_sun.u = [1  -1  0  0  0  0;
           0   0  1 -1  0  0;
           0   0  0  0  1 -1];
% sensor coefficients
% m(1,1) + m(1,2)*cos(theta) + ... + m(1,n)*cos(theta)^(n-1)
% set to pure cosine response
d_sun.m = repmat([0 1], [6, 1]);
% sun reference vector
d_sun.uSunECI = uS_ref;
% gaussian noise level (none)
d_sun.noise = zeros(1, 6);
% default values:
% d.kQ - indices of quaternion in state vector (7:10)

% SCT FIXES in MeasSunSensorAnalog.m: 
% changed line 65 from
% p = length(d.m);
% to
% p = size(d.m, 2);
% changed line 68 from
% y = d.m(k,1);
% to
% y(k) = d.m(k,1);
% changed line 67 from
% c = cos(uSunBody'*d.u(:,k));
% to
% c = uSunBody'*d.u(:,k);
y_sun = MeasSunSensorAnalog(x, d_sun);
% trim negative values
y_sun = max(0, y_sun);
% compute sun vector
uS = d_sun.u*y_sun;

end