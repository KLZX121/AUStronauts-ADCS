function [T, m] = DetumbleTorque(x, d, env)

%% Program Summary

%   INPUTS
%   x               (:, 1)  State vector
%   d               struct  Input parameters
%   .iWSat          (1, 3)  Indices of satellite ang vel in state vector
%   env             struct  Environmental Data
%   .B              (3, 1)  Local magnetic field (ECI) (T)
%
%   OUTPUT
%   T               (3, 1)  Control torque vector (x,y,z) (N m)
%   m               (3, 1)  Magnetorquer dipole moment (A m^2)

w = x(d.iWSat); % current satellite angular velocities

% Gain for controller
k = 4*pi/(90*60)*(1+sin(51.6*pi/180))*0.0067;
m = cross(k/norm(env.B)*w,env.B/norm(env.B));
T = cross(m,env.B);

end