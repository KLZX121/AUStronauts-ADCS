function T_c = DetumbleTorque(x, d, jD)

%% Program Summary

%   INPUTS
%   x               (:, 1)  state vector
%   d                         Input parameters (struct)
%       .iWSat      (1, 3)  Indices of satellite ang vel in state vector
%   jD              (1, 1)  Julian date for detumbling
%
%   OUTPUT
%   T_c               (3, 1) control torque vector (x,y,z)

if isempty(which('IGRFECI'))
    addpath(genpath('../attitude-determination/'));
end

w = x(d.iWSat); % current satellite angular velocities

% Local magnetic field
B = IGRFECI(x(d.iR),jD);

% Gain for controller
k = 4*pi/(90*60)*(1+sin(51.6*pi/180))*0.0067;
m = cross(k/norm(B)*w,B/norm(B));
T_c = cross(m,B);

end