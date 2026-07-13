function wSatDot = EulerDynamics(x, d)
%EulerDynamics Euler Rigid Body Dynamics Equation
%   wSatDot = EulerDynamics(wSat, wRW, d)
%   
%   INPUT
%   x               (:, 1)  State vector
%   d                       Input parameters (struct)
%           .TExt   (3, 1)  Total external applied torque (N m)
%           .TRW    (3, 1)  Total reaction wheel torque (N m)
%           .ISat   (3, 3)  Moment of Inertia tensor of satellite (kg m^2)
%           .IRW    (1, 1)  Moment of Inertia of Reaction Wheels (kg m^2)
%           .iWSat  (1, 3)  Indices of satellite ang vel in state vector
%           .iWRW   (1, 3)  Indices of rw ang vel in state vector
%   
%   OUTPUT
%   wSatDot         (3, 1)  Angular velocity ODE (rad/s)
%   
%   TODO: account for time varying I (subsystem deployment)

% state variables
wSat = x(d.iWSat);  % (3, 1) angular velocity of satellite (rad/s)
wRW = x(d.iWRW);     % (3, 1) angular velocity of reaction wheels (rad/s)

% ode
hT      = d.ISat*wSat + d.IRW*wRW; % total angular momentum (kg m^2/s^2)
wSatDot = d.ISat \ (d.TExt - d.TRW - cross(wSat, hT));

end