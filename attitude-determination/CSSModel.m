function uS = CSSModel(x, uSRef)
%CSSModel Coarse Sun Sensor Model
%   uS = CSSModel(x, uS_ref)
%   
%   Inputs
%   x       (:, 1)  state vector
%   uSRef  (3, 1)  reference unit sun vector (reference frame)
%   
%   Outputs
%   uS      (3, 1)  measured unit sun vector (body frame)
%   
%   TODO: analyse sun sensor script
%   TODO: add line of sight for sun sensor (earth blocks sun)
%   TODO: add fov of sensor
%   TODO: rewrite function for orbital parameters

% set data structure
dSun = MeasSunSensorAnalog;
% sensor unit vectors (6 sensors)
dSun.u = [1  -1  0  0  0  0;
           0   0  1 -1  0  0;
           0   0  0  0  1 -1];
% sensor coefficients
% m(1,1) + m(1,2)*cos(theta) + ... + m(1,n)*cos(theta)^(n-1)
% set to pure cosine response
dSun.m = repmat([0 1], [6, 1]);
% sun reference vector
dSun.uSunECI = uSRef;
% gaussian noise level (none)
dSun.noise = zeros(1, 6);
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
ySun = MeasSunSensorAnalog(x, dSun);
% trim negative values
ySun = max(0, ySun);
% compute sun vector
uS = dSun.u*ySun;

end