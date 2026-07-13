function wRWDot = RWDynamics(wSatDot, d)

% RWDynamics Reaction Wheel Dynamics Equation
%   wRWDot = RWDynamics(x, d)
%   
%   INPUT
%   .wSatDot        (3, 1)  satellite angular velocity ODE (rad/s)/s
%   d               Input parameters (struct)
%           .TRW    (3, 1)  Torque of reaction wheels (N m)
%           .IRW    (1, 1)  Moment of Inertia of Reaction Wheel (kg m^2)
%
%   
%   OUTPUT
%   wRWDot          (3, 1)  Reaction wheel speeds rel. to sat ODE (rad/s)/s

wRWDot = d.TRW./d.IRW - wSatDot;

end